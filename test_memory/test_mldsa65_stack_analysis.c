#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "../crypto_sign/ml-dsa-65/clean/api.h"
#include "../crypto_sign/ml-dsa-65/clean/params.h"
#include "../crypto_sign/ml-dsa-65/clean/poly.h"
#include "../crypto_sign/ml-dsa-65/clean/polyvec.h"
#include "../common/randombytes.h"
#include "../common/fips202.h"

// Message to sign
#define MSG_LEN 32
static uint8_t msg[MSG_LEN] = "Test message for ML-DSA-65!!!!";

int main(void) {
    printf("=== ML-DSA-65 (Dilithium 3) Stack Memory Analysis ===\n\n");

    printf("Algorithm Parameters:\n");
    printf("  N (polynomial degree):    %d\n", N);
    printf("  Q (modulus):              %d\n", Q);
    printf("  K (rows):                 %d\n", K);
    printf("  L (columns):              %d\n", L);
    printf("\n");

    // Calculate sizes
    size_t poly_size = sizeof(poly);  // N * sizeof(int32_t)
    size_t polyvecl_size = sizeof(polyvecl);  // L * poly_size
    size_t polyveck_size = sizeof(polyveck);  // K * poly_size

    printf("Data Structure Sizes:\n");
    printf("  poly (1 polynomial):      %zu bytes (%d coeffs × %zu bytes)\n", 
           poly_size, N, sizeof(int32_t));
    printf("  polyvecl (L=%d polys):    %zu bytes\n", L, polyvecl_size);
    printf("  polyveck (K=%d polys):    %zu bytes\n", K, polyveck_size);
    printf("\n");

    printf("=== SIGNING FUNCTION STACK VARIABLES ===\n\n");

    // Variables in crypto_sign_signature_ctx
    size_t seedbuf_size = 2 * SEEDBYTES + TRBYTES + RNDBYTES + 2 * CRHBYTES;
    size_t mat_size = K * polyvecl_size;  // polyvecl mat[K]
    size_t s1_size = polyvecl_size;
    size_t y_size = polyvecl_size;
    size_t z_size = polyvecl_size;
    size_t t0_size = polyveck_size;
    size_t s2_size = polyveck_size;
    size_t w1_size = polyveck_size;
    size_t w0_size = polyveck_size;
    size_t h_size = polyveck_size;
    size_t cp_size = poly_size;
    size_t shake_ctx_size = sizeof(shake256incctx);

    printf("Byte arrays:\n");
    printf("  seedbuf:                  %zu bytes\n", seedbuf_size);
    printf("\n");

    printf("Matrix A (mat[K]):\n");
    printf("  polyvecl mat[%d]:         %zu bytes  <-- LARGEST!\n", K, mat_size);
    printf("\n");

    printf("polyvecl vectors (L=%d polynomials each):\n", L);
    printf("  s1:                       %zu bytes\n", s1_size);
    printf("  y:                        %zu bytes\n", y_size);
    printf("  z:                        %zu bytes\n", z_size);
    printf("  Subtotal polyvecl:        %zu bytes\n", s1_size + y_size + z_size);
    printf("\n");

    printf("polyveck vectors (K=%d polynomials each):\n", K);
    printf("  t0:                       %zu bytes\n", t0_size);
    printf("  s2:                       %zu bytes\n", s2_size);
    printf("  w1:                       %zu bytes  <-- TARGET for recompute optimization\n", w1_size);
    printf("  w0:                       %zu bytes  <-- TARGET for recompute optimization\n", w0_size);
    printf("  h:                        %zu bytes\n", h_size);
    printf("  Subtotal polyveck:        %zu bytes\n", t0_size + s2_size + w1_size + w0_size + h_size);
    printf("\n");

    printf("Other:\n");
    printf("  cp (challenge poly):      %zu bytes\n", cp_size);
    printf("  shake256incctx:           %zu bytes\n", shake_ctx_size);
    printf("\n");

    size_t total_signing = seedbuf_size + mat_size + 
                           s1_size + y_size + z_size +
                           t0_size + s2_size + w1_size + w0_size + h_size +
                           cp_size + shake_ctx_size;

    printf("========================================\n");
    printf("TOTAL SIGNING STACK:        %zu bytes (%.2f KB)\n", total_signing, total_signing / 1024.0);
    printf("========================================\n\n");

    printf("=== MEMORY SAVINGS WITH RECOMPUTE OPTIMIZATION ===\n\n");
    printf("If we recompute w = A·y instead of storing w0:\n");
    printf("  Savings from w0:          %zu bytes (%.2f KB)\n", w0_size, w0_size / 1024.0);
    printf("  New total:                %zu bytes (%.2f KB)\n", 
           total_signing - w0_size, (total_signing - w0_size) / 1024.0);
    printf("  Reduction:                %.1f%%\n", 100.0 * w0_size / total_signing);
    printf("\n");

    // Actually run signing to verify correctness
    printf("=== FUNCTIONAL TEST ===\n\n");
    
    uint8_t pk[PQCLEAN_MLDSA65_CLEAN_CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[PQCLEAN_MLDSA65_CLEAN_CRYPTO_SECRETKEYBYTES];
    uint8_t sig[PQCLEAN_MLDSA65_CLEAN_CRYPTO_BYTES];
    size_t siglen;
    int ret;

    printf("Generating keypair... ");
    ret = PQCLEAN_MLDSA65_CLEAN_crypto_sign_keypair(pk, sk);
    printf("%s\n", ret == 0 ? "OK" : "FAILED");

    printf("Signing message...    ");
    ret = PQCLEAN_MLDSA65_CLEAN_crypto_sign_signature(sig, &siglen, msg, MSG_LEN, sk);
    printf("%s\n", ret == 0 ? "OK" : "FAILED");

    printf("Verifying signature.. ");
    ret = PQCLEAN_MLDSA65_CLEAN_crypto_sign_verify(sig, siglen, msg, MSG_LEN, pk);
    printf("%s\n", ret == 0 ? "OK" : "FAILED");

    printf("\n=== DONE ===\n");
    return 0;
}
