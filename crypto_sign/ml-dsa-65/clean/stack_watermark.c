#define _GNU_SOURCE
#include <errno.h>
#include <pthread.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "api.h"

#define STACK_SIZE (256*1024)
#define WATERMARK 0xAA
#define WATERMARK_GUARD_GAP 4096

typedef enum {
    JOB_EMPTY = 0,
    JOB_KEYPAIR = 1,
    JOB_SIGN = 2,
    JOB_VERIFY = 3
} job_t;

typedef struct {
    job_t job;
    int rc;
    int grows_down;
    const uint8_t *stack_base;
    size_t stack_size;
    size_t peak;
} measure_ctx_t;

static uint8_t g_pk[PQCLEAN_MLDSA65_CLEAN_CRYPTO_PUBLICKEYBYTES];
static uint8_t g_sk[PQCLEAN_MLDSA65_CLEAN_CRYPTO_SECRETKEYBYTES];
static uint8_t g_sig[PQCLEAN_MLDSA65_CLEAN_CRYPTO_BYTES];
static size_t g_siglen;
static const uint8_t g_msg[] = "Runtime stack watermark benchmark for ML-DSA-65";
static const size_t g_msglen = sizeof(g_msg) - 1;

static int run_algorithm(job_t job) {
    switch (job) {
        case JOB_EMPTY:
            return 0;
        case JOB_KEYPAIR:
            return PQCLEAN_MLDSA65_CLEAN_crypto_sign_keypair(g_pk, g_sk);
        case JOB_SIGN:
            return PQCLEAN_MLDSA65_CLEAN_crypto_sign_signature(g_sig, &g_siglen, g_msg, g_msglen, g_sk);
        case JOB_VERIFY:
            return PQCLEAN_MLDSA65_CLEAN_crypto_sign_verify(g_sig, g_siglen, g_msg, g_msglen, g_pk);
        default:
            return -1;
    }
}

static void *thread_func(void *arg) {
    measure_ctx_t *ctx = (measure_ctx_t *)arg;
    char marker;
    uintptr_t sp = (uintptr_t)&marker;
    uintptr_t low = (uintptr_t)ctx->stack_base;
    uintptr_t high = low + ctx->stack_size;
    uintptr_t fill_start;
    size_t fill_len;

    ctx->grows_down = (sp - low) > (high - sp);

    if (ctx->grows_down) {
        if (sp <= low + WATERMARK_GUARD_GAP) {
            ctx->rc = EFAULT;
            return NULL;
        }
        fill_start = low;
        fill_len = (size_t)(sp - low - WATERMARK_GUARD_GAP);
    } else {
        if (sp + WATERMARK_GUARD_GAP >= high) {
            ctx->rc = EFAULT;
            return NULL;
        }
        fill_start = sp + WATERMARK_GUARD_GAP;
        fill_len = (size_t)(high - fill_start);
    }

    memset((void *)fill_start, WATERMARK, fill_len);

    ctx->rc = run_algorithm(ctx->job);

    if (ctx->rc != 0) {
        return NULL;
    }

    if (ctx->grows_down) {
        size_t first_non = fill_len;
        const uint8_t *buf = (const uint8_t *)fill_start;
        for (size_t i = 0; i < fill_len; i++) {
            if (buf[i] != WATERMARK) {
                first_non = i;
                break;
            }
        }
        ctx->peak = first_non == fill_len ? 0 : (fill_len - first_non);
    } else {
        size_t suffix_watermark = 0;
        const uint8_t *buf = (const uint8_t *)fill_start;
        while (suffix_watermark < fill_len && buf[fill_len - 1 - suffix_watermark] == WATERMARK) {
            suffix_watermark++;
        }
        ctx->peak = fill_len - suffix_watermark;
    }

    return NULL;
}

static int measure_stack(job_t job, size_t *peak) {
    pthread_attr_t attr;
    pthread_t thread;
    measure_ctx_t ctx;
    long page_sz = sysconf(_SC_PAGESIZE);
    if (page_sz <= 0) {
        page_sz = 4096;
    }

    uint8_t *stack = aligned_alloc((size_t)page_sz, STACK_SIZE);
    if (stack == NULL) {
        return ENOMEM;
    }

    memset(&ctx, 0, sizeof(ctx));
    ctx.job = job;
    ctx.stack_base = stack;
    ctx.stack_size = STACK_SIZE;

    int rc = pthread_attr_init(&attr);
    if (rc != 0) {
        free(stack);
        return rc;
    }

    rc = pthread_attr_setstack(&attr, stack, STACK_SIZE);
    if (rc != 0) {
        pthread_attr_destroy(&attr);
        free(stack);
        return rc;
    }

    rc = pthread_create(&thread, &attr, thread_func, &ctx);
    pthread_attr_destroy(&attr);
    if (rc != 0) {
        free(stack);
        return rc;
    }

    rc = pthread_join(thread, NULL);
    if (rc != 0) {
        free(stack);
        return rc;
    }

    if (ctx.rc != 0) {
        free(stack);
        return EFAULT;
    }

    *peak = ctx.peak;
    free(stack);

    return 0;
}

int main() {
    size_t peak_baseline = 0;
    size_t peak_keypair = 0;
    size_t peak_sign = 0;
    size_t peak_verify = 0;
    size_t crypto_keypair = 0;
    size_t crypto_sign = 0;
    size_t crypto_verify = 0;
    int rc;

    rc = measure_stack(JOB_EMPTY, &peak_baseline);
    if (rc != 0) {
        fprintf(stderr, "baseline measurement failed: %s\n", strerror(rc));
        return 1;
    }

    rc = measure_stack(JOB_KEYPAIR, &peak_keypair);
    if (rc != 0) {
        fprintf(stderr, "keypair measurement failed: %s\n", strerror(rc));
        return 1;
    }

    if (PQCLEAN_MLDSA65_CLEAN_crypto_sign_keypair(g_pk, g_sk) != 0) {
        fprintf(stderr, "setup keypair failed\n");
        return 1;
    }

    rc = measure_stack(JOB_SIGN, &peak_sign);
    if (rc != 0) {
        fprintf(stderr, "sign measurement failed: %s\n", strerror(rc));
        return 1;
    }

    if (PQCLEAN_MLDSA65_CLEAN_crypto_sign_signature(g_sig, &g_siglen, g_msg, g_msglen, g_sk) != 0) {
        fprintf(stderr, "setup signature failed\n");
        return 1;
    }

    rc = measure_stack(JOB_VERIFY, &peak_verify);
    if (rc != 0) {
        fprintf(stderr, "verify measurement failed: %s\n", strerror(rc));
        return 1;
    }

    crypto_keypair = peak_keypair > peak_baseline ? peak_keypair - peak_baseline : 0;
    crypto_sign = peak_sign > peak_baseline ? peak_sign - peak_baseline : 0;
    crypto_verify = peak_verify > peak_baseline ? peak_verify - peak_baseline : 0;

    printf("Crypto stack usage keypair: %zu bytes\n", crypto_keypair);
    printf("Crypto stack usage sign   : %zu bytes\n", crypto_sign);
    printf("Crypto stack usage verify : %zu bytes\n", crypto_verify);

    return 0;
}