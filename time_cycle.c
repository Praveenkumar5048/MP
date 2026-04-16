#define _POSIX_C_SOURCE 200809L

#include <x86intrin.h>
#include <cpuid.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>

#include "crypto_sign/ml-dsa-65/clean/api.h"

#define ITERS_KEYGEN 300
#define ITERS_SIGN 1000
#define ITERS_VERIFY 1000
#define MSG_LEN 64

static double now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1e9 + (double)ts.tv_nsec;
}

static uint64_t rdtsc_start(void) {
    unsigned int a, b, c, d;
    __cpuid(0, a, b, c, d);
    return __rdtsc();
}

static uint64_t rdtsc_end(void) {
    unsigned int aux;
    uint64_t t = __rdtscp(&aux);
    unsigned int a, b, c, d;
    __cpuid(0, a, b, c, d);
    return t;
}

int main(void) {
    uint8_t pk[PQCLEAN_MLDSA65_CLEAN_CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[PQCLEAN_MLDSA65_CLEAN_CRYPTO_SECRETKEYBYTES];
    uint8_t msg[MSG_LEN];
    uint8_t sig[PQCLEAN_MLDSA65_CLEAN_CRYPTO_BYTES];
    size_t siglen = 0;

    for (size_t i = 0; i < MSG_LEN; i++) {
        msg[i] = (uint8_t)(i * 3 + 1);
    }

    if (PQCLEAN_MLDSA65_CLEAN_crypto_sign_keypair(pk, sk) != 0) {
        fprintf(stderr, "initial keypair failed\n");
        return 1;
    }

    for (int i = 0; i < 20; i++) {
        PQCLEAN_MLDSA65_CLEAN_crypto_sign_keypair(pk, sk);
        PQCLEAN_MLDSA65_CLEAN_crypto_sign_signature(sig, &siglen, msg, MSG_LEN, sk);
        PQCLEAN_MLDSA65_CLEAN_crypto_sign_verify(sig, siglen, msg, MSG_LEN, pk);
    }

    double t0, t1;
    uint64_t c0, c1;

    t0 = now_ns();
    c0 = rdtsc_start();
    for (int i = 0; i < ITERS_KEYGEN; i++) {
        if (PQCLEAN_MLDSA65_CLEAN_crypto_sign_keypair(pk, sk) != 0) {
            fprintf(stderr, "keypair failed at iter %d\n", i);
            return 1;
        }
    }
    c1 = rdtsc_end();
    t1 = now_ns();
    double keygen_ns = (t1 - t0) / ITERS_KEYGEN;
    double keygen_cycles = (double)(c1 - c0) / ITERS_KEYGEN;

    t0 = now_ns();
    c0 = rdtsc_start();
    for (int i = 0; i < ITERS_SIGN; i++) {
        if (PQCLEAN_MLDSA65_CLEAN_crypto_sign_signature(sig, &siglen, msg, MSG_LEN, sk) != 0) {
            fprintf(stderr, "sign failed at iter %d\n", i);
            return 1;
        }
    }
    c1 = rdtsc_end();
    t1 = now_ns();
    double sign_ns = (t1 - t0) / ITERS_SIGN;
    double sign_cycles = (double)(c1 - c0) / ITERS_SIGN;

    if (PQCLEAN_MLDSA65_CLEAN_crypto_sign_signature(sig, &siglen, msg, MSG_LEN, sk) != 0) {
        fprintf(stderr, "prep sign failed\n");
        return 1;
    }

    t0 = now_ns();
    c0 = rdtsc_start();
    for (int i = 0; i < ITERS_VERIFY; i++) {
        if (PQCLEAN_MLDSA65_CLEAN_crypto_sign_verify(sig, siglen, msg, MSG_LEN, pk) != 0) {
            fprintf(stderr, "verify failed at iter %d\n", i);
            return 1;
        }
    }
    c1 = rdtsc_end();
    t1 = now_ns();
    double verify_ns = (t1 - t0) / ITERS_VERIFY;
    double verify_cycles = (double)(c1 - c0) / ITERS_VERIFY;

    printf("ML-DSA-65 timing and cycles on current branch\n");
    printf("Iterations: keygen=%d sign=%d verify=%d\n", ITERS_KEYGEN, ITERS_SIGN, ITERS_VERIFY);
    printf("KeyGen : %.3f ms/op | %.0f cycles/op\n", keygen_ns / 1e6, keygen_cycles);
    printf("Sign   : %.3f ms/op | %.0f cycles/op\n", sign_ns / 1e6, sign_cycles);
    printf("Verify : %.3f ms/op | %.0f cycles/op\n", verify_ns / 1e6, verify_cycles);

    return 0;
}