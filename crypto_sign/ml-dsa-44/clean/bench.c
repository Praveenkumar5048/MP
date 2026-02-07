/*
 * Simple benchmark that measures per-operation timing and memory usage for
 * PQClean ml-dsa-44 (clean) implementation. For reliable peak RSS reporting
 * we run each measured operation in a child process (fork) and collect
 * getrusage(RUSAGE_SELF).ru_maxrss from the child process.
 *
 * Outputs CSV-like lines for keygen/sign/verify with:
 *  - iterations
 *  - avg/min/max time in microseconds
 *  - heap_delta_bytes (mallinfo.uordblks difference)
 *  - peak_rss_kb (ru_maxrss)
 *
 * Usage: ./bench [iterations]
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <sys/wait.h>
#include <malloc.h>
#include "api.h"
#include <pthread.h>

static inline uint64_t now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
}

static size_t heap_used_bytes(void) {
#if defined(__GLIBC__)
    struct mallinfo mi = mallinfo();
    return (size_t)mi.uordblks; /* total allocated bytes */
#else
    void *brk = sbrk(0);
    return (size_t)brk;
#endif
}

/* Read a few memory fields from /proc/self/status if available. Values are
 * returned in kilobytes. If a field is not found, its value is -1.
 */
struct proc_status_fields {
    long VmPeak_kb;
    long VmHWM_kb;
    long VmRSS_kb;
    long VmData_kb;
    long VmStk_kb;
};

static struct proc_status_fields read_proc_status(void) {
    struct proc_status_fields ret = { -1, -1, -1, -1, -1 };
    FILE *f = fopen("/proc/self/status", "r");
    if (!f) return ret;
    char line[256];
    while (fgets(line, sizeof(line), f)) {
        long v;
        if (sscanf(line, "VmPeak: %ld kB", &v) == 1) ret.VmPeak_kb = v;
        else if (sscanf(line, "VmHWM: %ld kB", &v) == 1) ret.VmHWM_kb = v;
        else if (sscanf(line, "VmRSS: %ld kB", &v) == 1) ret.VmRSS_kb = v;
        else if (sscanf(line, "VmData: %ld kB", &v) == 1) ret.VmData_kb = v;
        else if (sscanf(line, "VmStk: %ld kB", &v) == 1) ret.VmStk_kb = v;
    }
    fclose(f);
    return ret;
}

/* Fill stack region with a pattern (touch pages) so we can detect overwritten
 * areas after calling an operation. Returns 0 on success.
 */
static int fill_stack_pattern(unsigned char pattern) {
    pthread_attr_t attr;
    size_t stack_size = 0;
    void *stack_addr = NULL;
    if (pthread_getattr_np(pthread_self(), &attr) != 0) {
        return -1;
    }
    pthread_attr_getstack(&attr, &stack_addr, &stack_size);
    pthread_attr_destroy(&attr);
    if (stack_addr == NULL || stack_size == 0) return -1;

    /* stack_addr is the lowest address of the stack region; stack grows down
     * on most Linux systems, so the usable stack top is stack_addr + stack_size
     */
    unsigned char *stack_low = (unsigned char *)stack_addr;
    unsigned char *stack_high = stack_low + stack_size;
    /* Current stack pointer */
    unsigned char marker;
    unsigned char *curr_sp = (unsigned char *)&marker;

    /* Leave a safety margin near the current frame to avoid overwriting
     * local variables, canaries, or saved registers. Touch every 64 bytes
     * from (curr_sp + safety_margin) up to stack_high with the pattern. */
    const size_t STACK_SAFETY_MARGIN = 4096; /* bytes */
    unsigned char *p = curr_sp + STACK_SAFETY_MARGIN;
    for (; p < stack_high; p += 64) {
        /* volatile write to avoid optimization */
        volatile unsigned char *vp = p;
        *vp = pattern;
    }
    /* Ensure last page touched */
    if (p >= stack_high && (stack_high - 1) >= curr_sp) {
        volatile unsigned char *vp = stack_high - 1;
        *vp = pattern;
    }
    return 0;
}

/* After calling an operation, scan from current SP upwards to find the first
 * byte that differs from the pattern; that indicates stack growth into the
 * filled region. Returns bytes of stack used (approx) or -1 on error.
 */
static long measure_stack_highwater_bytes(unsigned char pattern) {
    pthread_attr_t attr;
    size_t stack_size = 0;
    void *stack_addr = NULL;
    if (pthread_getattr_np(pthread_self(), &attr) != 0) {
        return -1;
    }
    pthread_attr_getstack(&attr, &stack_addr, &stack_size);
    pthread_attr_destroy(&attr);
    if (stack_addr == NULL || stack_size == 0) return -1;

    unsigned char *stack_low = (unsigned char *)stack_addr;
    unsigned char *stack_high = stack_low + stack_size;
    unsigned char marker;
    unsigned char *curr_sp = (unsigned char *)&marker;

    const size_t STACK_SAFETY_MARGIN = 4096; /* bytes */
    unsigned char *p = curr_sp + STACK_SAFETY_MARGIN;
    if (p < stack_low) p = stack_low;
    for (; p < stack_high; p++) {
        if (*p != pattern) {
            /* overwritten region from p .. stack_high-1 */
            return (long)(stack_high - p);
        }
    }
    return 0;
}

/* Operation wrappers */
static void op_keygen_once(uint8_t *pk, uint8_t *sk) {
    PQCLEAN_MLDSA44_CLEAN_crypto_sign_keypair(pk, sk);
}
static void op_sign_once(const uint8_t *m, size_t mlen, uint8_t *sig, size_t *siglen, const uint8_t *sk) {
    PQCLEAN_MLDSA44_CLEAN_crypto_sign_signature(sig, siglen, m, mlen, sk);
}
static void op_verify_once(const uint8_t *sig, size_t siglen, const uint8_t *m, size_t mlen, const uint8_t *pk) {
    PQCLEAN_MLDSA44_CLEAN_crypto_sign_verify(sig, siglen, m, mlen, pk);
}

/* Run a loop of "iterations" of the provided operation in the current process
 * and print a single CSV-like result line to stdout. This is intended to be
 * called from a child process after fork() so we can collect per-operation
 * peak RSS via getrusage(RUSAGE_SELF).
 */
static void run_keygen_loop(int iterations) {
    uint8_t *pk = malloc(PQCLEAN_MLDSA44_CLEAN_CRYPTO_PUBLICKEYBYTES);
    uint8_t *sk = malloc(PQCLEAN_MLDSA44_CLEAN_CRYPTO_SECRETKEYBYTES);
    if (!pk || !sk) { perror("alloc"); exit(1); }

    uint64_t min_ns = (uint64_t)-1, max_ns = 0, total_ns = 0;
    struct proc_status_fields before = read_proc_status();
    size_t heap_before = heap_used_bytes();
    /* Fill stack with pattern so we can measure high-water after operations */
    fill_stack_pattern(0xAA);
    long max_stack_used_bytes = 0;
    for (int i = 0; i < iterations; i++) {
        uint64_t t0 = now_ns();
        op_keygen_once(pk, sk);
        uint64_t t1 = now_ns();
        uint64_t dt = t1 - t0;
        if (dt < min_ns) min_ns = dt;
        if (dt > max_ns) max_ns = dt;
        total_ns += dt;
        /* measure stack high-water after this iteration */
        long stack_used = measure_stack_highwater_bytes(0xAA);
        if (stack_used > max_stack_used_bytes) max_stack_used_bytes = stack_used;
    }
    size_t heap_after = heap_used_bytes();
    struct rusage ru; getrusage(RUSAGE_SELF, &ru);
    long peak_rss_kb = ru.ru_maxrss;
    struct proc_status_fields after = read_proc_status();

    long VmData_delta_bytes = -1;
    if (before.VmData_kb != -1 && after.VmData_kb != -1) {
        VmData_delta_bytes = (after.VmData_kb - before.VmData_kb) * 1024L;
        if (VmData_delta_bytes < 0) VmData_delta_bytes = 0;
    }

    printf("keygen,iterations=%d,avg_us=%.2f,min_us=%.2f,max_us=%.2f,heap_delta_bytes=%zu,VmData_delta_bytes=%ld,stack_peak_bytes=%ld,VmStk_kb=%ld,VmHWM_kb=%ld,VmPeak_kb=%ld,peak_rss_kb=%ld\n",
           iterations,
           (double)total_ns / iterations / 1000.0,
           (double)min_ns / 1000.0,
           (double)max_ns / 1000.0,
           (heap_after >= heap_before) ? (heap_after - heap_before) : 0,
           VmData_delta_bytes,
        max_stack_used_bytes,
        after.VmStk_kb,
           after.VmHWM_kb,
           after.VmPeak_kb,
           peak_rss_kb);
    fflush(stdout);
    free(pk); free(sk);
}

static void run_sign_loop(int iterations) {
    size_t mlen = 1024;
    uint8_t *m = malloc(mlen);
    uint8_t *sig = malloc(PQCLEAN_MLDSA44_CLEAN_CRYPTO_BYTES + 16);
    size_t siglen = 0;
    uint8_t *pk = malloc(PQCLEAN_MLDSA44_CLEAN_CRYPTO_PUBLICKEYBYTES);
    uint8_t *sk = malloc(PQCLEAN_MLDSA44_CLEAN_CRYPTO_SECRETKEYBYTES);
    if (!m || !sig || !pk || !sk) { perror("alloc"); exit(1); }
    /* seed */
    for (size_t i = 0; i < mlen; i++) m[i] = (uint8_t)(i & 0xFF);
    /* ensure keys exist */
    PQCLEAN_MLDSA44_CLEAN_crypto_sign_keypair(pk, sk);

    uint64_t min_ns = (uint64_t)-1, max_ns = 0, total_ns = 0;
    struct proc_status_fields before = read_proc_status();
    size_t heap_before = heap_used_bytes();
    fill_stack_pattern(0xAA);
    long max_stack_used_bytes = 0;
    for (int i = 0; i < iterations; i++) {
        uint64_t t0 = now_ns();
        op_sign_once(m, mlen, sig, &siglen, sk);
        uint64_t t1 = now_ns();
        uint64_t dt = t1 - t0;
        if (dt < min_ns) min_ns = dt;
        if (dt > max_ns) max_ns = dt;
        total_ns += dt;
        long stack_used = measure_stack_highwater_bytes(0xAA);
        if (stack_used > max_stack_used_bytes) max_stack_used_bytes = stack_used;
    }
    size_t heap_after = heap_used_bytes();
    struct rusage ru; getrusage(RUSAGE_SELF, &ru);
    long peak_rss_kb = ru.ru_maxrss;
    struct proc_status_fields after = read_proc_status();

    long VmData_delta_bytes = -1;
    if (before.VmData_kb != -1 && after.VmData_kb != -1) {
        VmData_delta_bytes = (after.VmData_kb - before.VmData_kb) * 1024L;
        if (VmData_delta_bytes < 0) VmData_delta_bytes = 0;
    }

    printf("sign,iterations=%d,avg_us=%.2f,min_us=%.2f,max_us=%.2f,heap_delta_bytes=%zu,VmData_delta_bytes=%ld,stack_peak_bytes=%ld,VmStk_kb=%ld,VmHWM_kb=%ld,VmPeak_kb=%ld,peak_rss_kb=%ld\n",
           iterations,
           (double)total_ns / iterations / 1000.0,
           (double)min_ns / 1000.0,
           (double)max_ns / 1000.0,
           (heap_after >= heap_before) ? (heap_after - heap_before) : 0,
           VmData_delta_bytes,
        max_stack_used_bytes,
        after.VmStk_kb,
           after.VmHWM_kb,
           after.VmPeak_kb,
           peak_rss_kb);
    fflush(stdout);
    free(m); free(sig); free(pk); free(sk);
}

static void run_verify_loop(int iterations) {
    size_t mlen = 1024;
    uint8_t *m = malloc(mlen);
    uint8_t *sig = malloc(PQCLEAN_MLDSA44_CLEAN_CRYPTO_BYTES + 16);
    size_t siglen = 0;
    uint8_t *pk = malloc(PQCLEAN_MLDSA44_CLEAN_CRYPTO_PUBLICKEYBYTES);
    uint8_t *sk = malloc(PQCLEAN_MLDSA44_CLEAN_CRYPTO_SECRETKEYBYTES);
    if (!m || !sig || !pk || !sk) { perror("alloc"); exit(1); }
    for (size_t i = 0; i < mlen; i++) m[i] = (uint8_t)(i & 0xFF);
    PQCLEAN_MLDSA44_CLEAN_crypto_sign_keypair(pk, sk);
    PQCLEAN_MLDSA44_CLEAN_crypto_sign_signature(sig, &siglen, m, mlen, sk);

    uint64_t min_ns = (uint64_t)-1, max_ns = 0, total_ns = 0;
    struct proc_status_fields before = read_proc_status();
    size_t heap_before = heap_used_bytes();
    fill_stack_pattern(0xAA);
    long max_stack_used_bytes = 0;
    for (int i = 0; i < iterations; i++) {
        uint64_t t0 = now_ns();
        op_verify_once(sig, siglen, m, mlen, pk);
        uint64_t t1 = now_ns();
        uint64_t dt = t1 - t0;
        if (dt < min_ns) min_ns = dt;
        if (dt > max_ns) max_ns = dt;
        total_ns += dt;
        long stack_used = measure_stack_highwater_bytes(0xAA);
        if (stack_used > max_stack_used_bytes) max_stack_used_bytes = stack_used;
    }
    size_t heap_after = heap_used_bytes();
    struct rusage ru; getrusage(RUSAGE_SELF, &ru);
    long peak_rss_kb = ru.ru_maxrss;
    struct proc_status_fields after = read_proc_status();

    long VmData_delta_bytes = -1;
    if (before.VmData_kb != -1 && after.VmData_kb != -1) {
        VmData_delta_bytes = (after.VmData_kb - before.VmData_kb) * 1024L;
        if (VmData_delta_bytes < 0) VmData_delta_bytes = 0;
    }

    printf("verify,iterations=%d,avg_us=%.2f,min_us=%.2f,max_us=%.2f,heap_delta_bytes=%zu,VmData_delta_bytes=%ld,stack_peak_bytes=%ld,VmStk_kb=%ld,VmHWM_kb=%ld,VmPeak_kb=%ld,peak_rss_kb=%ld\n",
           iterations,
           (double)total_ns / iterations / 1000.0,
           (double)min_ns / 1000.0,
           (double)max_ns / 1000.0,
           (heap_after >= heap_before) ? (heap_after - heap_before) : 0,
           VmData_delta_bytes,
        max_stack_used_bytes,
        after.VmStk_kb,
           after.VmHWM_kb,
           after.VmPeak_kb,
           peak_rss_kb);
    fflush(stdout);
    free(m); free(sig); free(pk); free(sk);
}

static int spawn_and_wait(void (*child_fn)(int), int iterations) {
    pid_t pid = fork();
    if (pid < 0) { perror("fork"); return -1; }
    if (pid == 0) {
        /* child */
        child_fn(iterations);
        _exit(0);
    }
    /* parent */
    int status = 0;
    if (waitpid(pid, &status, 0) < 0) { perror("waitpid"); return -1; }
    if (WIFEXITED(status)) return WEXITSTATUS(status);
    return -1;
}

int main(int argc, char **argv) {
    int iterations = 10;
    if (argc > 1) iterations = atoi(argv[1]);

    printf("# PQClean ml-dsa-44 benchmark (fork per operation)\n");
    printf("# iterations=%d\n", iterations);

    /* Run keygen in a child process to collect peak RSS separately */
    if (spawn_and_wait(run_keygen_loop, iterations) != 0) {
        fprintf(stderr, "keygen child failed\n");
    }
    if (spawn_and_wait(run_sign_loop, iterations) != 0) {
        fprintf(stderr, "sign child failed\n");
    }
    if (spawn_and_wait(run_verify_loop, iterations) != 0) {
        fprintf(stderr, "verify child failed\n");
    }
    return 0;
}
