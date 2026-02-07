#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "../crypto_sign/ml-dsa-65/clean/api.h"
#include "../common/randombytes.h"

#define NUM_TESTS 100

int main(void) {
    uint8_t pk[PQCLEAN_MLDSA65_CLEAN_CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[PQCLEAN_MLDSA65_CLEAN_CRYPTO_SECRETKEYBYTES];
    uint8_t sig[PQCLEAN_MLDSA65_CLEAN_CRYPTO_BYTES];
    uint8_t msg[128];
    size_t siglen;
    int ret;
    int failures = 0;

    printf("=== ML-DSA-65 Correctness Test (%d iterations) ===\n\n", NUM_TESTS);

    for (int i = 0; i < NUM_TESTS; i++) {
        // Generate random message
        randombytes(msg, sizeof(msg));
        
        // Generate keypair
        ret = PQCLEAN_MLDSA65_CLEAN_crypto_sign_keypair(pk, sk);
        if (ret != 0) {
            printf("Test %d: Keypair generation FAILED\n", i);
            failures++;
            continue;
        }

        // Sign
        ret = PQCLEAN_MLDSA65_CLEAN_crypto_sign_signature(sig, &siglen, msg, sizeof(msg), sk);
        if (ret != 0) {
            printf("Test %d: Signing FAILED\n", i);
            failures++;
            continue;
        }

        // Verify
        ret = PQCLEAN_MLDSA65_CLEAN_crypto_sign_verify(sig, siglen, msg, sizeof(msg), pk);
        if (ret != 0) {
            printf("Test %d: Verification FAILED\n", i);
            failures++;
            continue;
        }

        // Verify with wrong message should fail
        msg[0] ^= 0xFF;
        ret = PQCLEAN_MLDSA65_CLEAN_crypto_sign_verify(sig, siglen, msg, sizeof(msg), pk);
        if (ret == 0) {
            printf("Test %d: Verification with wrong message should have FAILED\n", i);
            failures++;
        }
        msg[0] ^= 0xFF;  // Restore

        if ((i + 1) % 10 == 0) {
            printf("  Completed %d/%d tests...\n", i + 1, NUM_TESTS);
        }
    }

    printf("\n========================================\n");
    if (failures == 0) {
        printf("All %d tests PASSED!\n", NUM_TESTS);
    } else {
        printf("%d out of %d tests FAILED!\n", failures, NUM_TESTS);
    }
    printf("========================================\n");

    return failures > 0 ? 1 : 0;
}
