#define _POSIX_C_SOURCE 200809L
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

    t0 = now_ns();
    for (int i = 0; i < ITERS_KEYGEN; i++) {
        if (PQCLEAN_MLDSA65_CLEAN_crypto_sign_keypair(pk, sk) != 0) {
            fprintf(stderr, "keypair failed at iter %d\n", i);
            return 1;
        }
    }
    t1 = now_ns();
    double keygen_ns = (t1 - t0) / ITERS_KEYGEN;

    t0 = now_ns();
    for (int i = 0; i < ITERS_SIGN; i++) {
        if (PQCLEAN_MLDSA65_CLEAN_crypto_sign_signature(sig, &siglen, msg, MSG_LEN, sk) != 0) {
            fprintf(stderr, "sign failed at iter %d\n", i);
            return 1;
        }
    }
    t1 = now_ns();
    double sign_ns = (t1 - t0) / ITERS_SIGN;

    if (PQCLEAN_MLDSA65_CLEAN_crypto_sign_signature(sig, &siglen, msg, MSG_LEN, sk) != 0) {
        fprintf(stderr, "prep sign failed\n");
        return 1;
    }

    t0 = now_ns();
    for (int i = 0; i < ITERS_VERIFY; i++) {
        if (PQCLEAN_MLDSA65_CLEAN_crypto_sign_verify(sig, siglen, msg, MSG_LEN, pk) != 0) {
            fprintf(stderr, "verify failed at iter %d\n", i);
            return 1;
        }
    }
    t1 = now_ns();
    double verify_ns = (t1 - t0) / ITERS_VERIFY;

    printf("ML-DSA-65 timing on current branch\n");
    printf("Iterations: keygen=%d sign=%d verify=%d\n", ITERS_KEYGEN, ITERS_SIGN, ITERS_VERIFY);
    printf("KeyGen : %.3f ms/op\n", keygen_ns / 1e6);
    printf("Sign   : %.3f ms/op\n", sign_ns / 1e6);
    printf("Verify : %.3f ms/op\n", verify_ns / 1e6);

    return 0;
}
