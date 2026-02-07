#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "../crypto_sign/ml-dsa-65/clean/api.h"
#include "../common/randombytes.h"

// Message to sign
#define MSG_LEN 32
static uint8_t msg[MSG_LEN] = "Test message for ML-DSA-65!!!!";

int main(void) {
    uint8_t pk[PQCLEAN_MLDSA65_CLEAN_CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[PQCLEAN_MLDSA65_CLEAN_CRYPTO_SECRETKEYBYTES];
    uint8_t sig[PQCLEAN_MLDSA65_CLEAN_CRYPTO_BYTES];
    size_t siglen;
    int ret;

    printf("=== ML-DSA-65 (Dilithium 3) Memory Test ===\n\n");

    printf("Key and Signature Sizes:\n");
    printf("  Public Key:  %d bytes\n", PQCLEAN_MLDSA65_CLEAN_CRYPTO_PUBLICKEYBYTES);
    printf("  Secret Key:  %d bytes\n", PQCLEAN_MLDSA65_CLEAN_CRYPTO_SECRETKEYBYTES);
    printf("  Signature:   %d bytes\n", PQCLEAN_MLDSA65_CLEAN_CRYPTO_BYTES);
    printf("\n");

    // Generate keypair
    printf("Generating keypair...\n");
    ret = PQCLEAN_MLDSA65_CLEAN_crypto_sign_keypair(pk, sk);
    if (ret != 0) {
        printf("Keypair generation failed!\n");
        return 1;
    }
    printf("Keypair generated successfully.\n\n");

    // Sign the message
    printf("Signing message...\n");
    ret = PQCLEAN_MLDSA65_CLEAN_crypto_sign_signature(sig, &siglen, msg, MSG_LEN, sk);
    if (ret != 0) {
        printf("Signing failed!\n");
        return 1;
    }
    printf("Signature created successfully. Length: %zu bytes\n\n", siglen);

    // Verify the signature
    printf("Verifying signature...\n");
    ret = PQCLEAN_MLDSA65_CLEAN_crypto_sign_verify(sig, siglen, msg, MSG_LEN, pk);
    if (ret != 0) {
        printf("Verification failed!\n");
        return 1;
    }
    printf("Signature verified successfully.\n\n");

    printf("=== Test completed successfully ===\n");
    return 0;
}
