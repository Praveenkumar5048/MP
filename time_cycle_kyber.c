#define _POSIX_C_SOURCE 200809L

#include <x86intrin.h>
#include <cpuid.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>

#include "crypto_kem/ml-kem-768/clean/api.h"

#define ITERS_KEYGEN 1000
#define ITERS_ENCAP 3000
#define ITERS_DECAP 3000

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
    uint8_t pk[PQCLEAN_MLKEM768_CLEAN_CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[PQCLEAN_MLKEM768_CLEAN_CRYPTO_SECRETKEYBYTES];
    uint8_t ct[PQCLEAN_MLKEM768_CLEAN_CRYPTO_CIPHERTEXTBYTES];
    uint8_t ss_enc[PQCLEAN_MLKEM768_CLEAN_CRYPTO_BYTES];
    uint8_t ss_dec[PQCLEAN_MLKEM768_CLEAN_CRYPTO_BYTES];

    double t0, t1;
    uint64_t c0, c1;

    if (PQCLEAN_MLKEM768_CLEAN_crypto_kem_keypair(pk, sk) != 0) {
        fprintf(stderr, "initial keypair failed\n");
        return 1;
    }
    if (PQCLEAN_MLKEM768_CLEAN_crypto_kem_enc(ct, ss_enc, pk) != 0) {
        fprintf(stderr, "initial encaps failed\n");
        return 1;
    }
    if (PQCLEAN_MLKEM768_CLEAN_crypto_kem_dec(ss_dec, ct, sk) != 0) {
        fprintf(stderr, "initial decaps failed\n");
        return 1;
    }

    for (int i = 0; i < 50; i++) {
        PQCLEAN_MLKEM768_CLEAN_crypto_kem_keypair(pk, sk);
        PQCLEAN_MLKEM768_CLEAN_crypto_kem_enc(ct, ss_enc, pk);
        PQCLEAN_MLKEM768_CLEAN_crypto_kem_dec(ss_dec, ct, sk);
    }

    t0 = now_ns();
    c0 = rdtsc_start();
    for (int i = 0; i < ITERS_KEYGEN; i++) {
        if (PQCLEAN_MLKEM768_CLEAN_crypto_kem_keypair(pk, sk) != 0) {
            fprintf(stderr, "keypair failed at iter %d\n", i);
            return 1;
        }
    }
    c1 = rdtsc_end();
    t1 = now_ns();
    double keygen_ns = (t1 - t0) / ITERS_KEYGEN;
    double keygen_cycles = (double)(c1 - c0) / ITERS_KEYGEN;

    if (PQCLEAN_MLKEM768_CLEAN_crypto_kem_keypair(pk, sk) != 0) {
        fprintf(stderr, "prep keypair failed\n");
        return 1;
    }

    t0 = now_ns();
    c0 = rdtsc_start();
    for (int i = 0; i < ITERS_ENCAP; i++) {
        if (PQCLEAN_MLKEM768_CLEAN_crypto_kem_enc(ct, ss_enc, pk) != 0) {
            fprintf(stderr, "encaps failed at iter %d\n", i);
            return 1;
        }
    }
    c1 = rdtsc_end();
    t1 = now_ns();
    double enc_ns = (t1 - t0) / ITERS_ENCAP;
    double enc_cycles = (double)(c1 - c0) / ITERS_ENCAP;

    if (PQCLEAN_MLKEM768_CLEAN_crypto_kem_enc(ct, ss_enc, pk) != 0) {
        fprintf(stderr, "prep encaps failed\n");
        return 1;
    }

    t0 = now_ns();
    c0 = rdtsc_start();
    for (int i = 0; i < ITERS_DECAP; i++) {
        if (PQCLEAN_MLKEM768_CLEAN_crypto_kem_dec(ss_dec, ct, sk) != 0) {
            fprintf(stderr, "decaps failed at iter %d\n", i);
            return 1;
        }
    }
    c1 = rdtsc_end();
    t1 = now_ns();
    double dec_ns = (t1 - t0) / ITERS_DECAP;
    double dec_cycles = (double)(c1 - c0) / ITERS_DECAP;

    printf("ML-KEM-768 timing and cycles on current branch\n");
    printf("Iterations: keygen=%d encaps=%d decaps=%d\n", ITERS_KEYGEN, ITERS_ENCAP, ITERS_DECAP);
    printf("KeyGen : %.3f ms/op | %.0f cycles/op\n", keygen_ns / 1e6, keygen_cycles);
    printf("Encaps : %.3f ms/op | %.0f cycles/op\n", enc_ns / 1e6, enc_cycles);
    printf("Decaps : %.3f ms/op | %.0f cycles/op\n", dec_ns / 1e6, dec_cycles);

    return 0;
}