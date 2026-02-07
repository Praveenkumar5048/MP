#ifndef PQCLEAN_MLDSA65_CLEAN_SMALLNTT_H
#define PQCLEAN_MLDSA65_CLEAN_SMALLNTT_H

#include "params.h"
#include <stdint.h>

/**
 * Small-modulus NTT for memory-efficient c*s1 and c*s2 computation.
 * 
 * Key insight: The challenge polynomial c has TAU=49 coefficients in {-1,0,1}
 * and secret polynomials s1,s2 have coefficients in [-ETA, ETA] = [-4, 4].
 * The product c*s has coefficients bounded by TAU * ETA = 49 * 4 = 196.
 * 
 * We use q_small = 7681 = 15*512 + 1, which:
 *   - Is prime and supports 512th roots of unity (needed for negacyclic NTT)
 *   - Is much larger than 2 * 196 = 392, so no wraparound occurs
 *   - Fits in 13 bits, allowing 16-bit coefficient storage with headroom
 *
 * This reduces NTT buffer memory by 50% compared to 32-bit coefficients.
 *
 * Memory savings per polynomial: 256 * (4 - 2) = 512 bytes
 * For c*s1 (5 polys) + c*s2 (6 polys): potential savings of ~5.5 KB
 */

#define SMALL_Q 7681

/**
 * 16-bit polynomial for small-modulus operations
 */
typedef struct {
    int16_t coeffs[N];
} smallpoly;

/**
 * Convert standard poly (32-bit coeffs) to smallpoly (16-bit coeffs).
 * Input coefficients must be in range [-SMALL_Q/2, SMALL_Q/2].
 * For ML-DSA-65: s coefficients in [-4, 4], c coefficients in {-1, 0, 1}.
 *
 * @param b Output smallpoly
 * @param a Input array of 32-bit coefficients
 */
void poly_to_smallpoly(smallpoly *b, const int32_t *a);

/**
 * Convert smallpoly back to standard poly (32-bit coeffs).
 * Output is centered in [-(SMALL_Q-1)/2, (SMALL_Q-1)/2].
 *
 * @param b Output array of 32-bit coefficients
 * @param a Input smallpoly
 */
void smallpoly_to_poly(int32_t *b, const smallpoly *a);

/**
 * Forward NTT on smallpoly.
 * Uses twist-based negacyclic NTT with omega = 2028 (256th root of unity).
 * After this, polynomial is in NTT domain.
 *
 * @param a Polynomial to transform (modified in place)
 */
void smallpoly_ntt(smallpoly *a);

/**
 * Inverse NTT on smallpoly.
 * Includes division by N and untwist step.
 * After this, polynomial is back in normal domain.
 *
 * @param a Polynomial to transform (modified in place)
 */
void smallpoly_invntt(smallpoly *a);

/**
 * Pointwise multiplication in NTT domain.
 * c[i] = a[i] * b[i] mod SMALL_Q
 *
 * @param c Output polynomial (in NTT domain)
 * @param a First input polynomial (in NTT domain)
 * @param b Second input polynomial (in NTT domain)
 */
void smallpoly_pointwise(smallpoly *c, const smallpoly *a, const smallpoly *b);

/**
 * Full polynomial multiplication using small NTT.
 * Computes c = a * b mod (X^N + 1) mod SMALL_Q
 *
 * Input a and b should be in normal domain (not NTT).
 * Output c is in normal domain.
 *
 * @param c Output polynomial
 * @param a First input polynomial
 * @param b Second input polynomial
 */
void smallntt_poly_mul(smallpoly *c, const smallpoly *a, const smallpoly *b);

/*************************************************
 * Integration functions for ML-DSA-65 signing
 *************************************************/

/**
 * Convert challenge polynomial to small NTT form.
 * Can be reused for multiple c*s multiplications.
 *
 * @param c_ntt Output challenge in small NTT domain
 * @param c Input challenge coefficients (TAU non-zero ±1 entries)
 */
void smallntt_challenge_to_ntt(smallpoly *c_ntt, const int32_t *c);

/**
 * Compute c*s for a single polynomial using small NTT.
 *
 * @param result Output coefficients in [-TAU*ETA, TAU*ETA]
 * @param c_ntt Challenge already in small NTT domain
 * @param s Secret polynomial coefficients in [-ETA, ETA]
 */
void smallntt_mul_poly(int32_t *result, const smallpoly *c_ntt, const int32_t *s);

/**
 * Compute c*s1 for polyvecl (L=5 polynomials) using small NTT.
 * 
 * Replaces the standard:
 *   poly_ntt(&cp);
 *   polyvecl_ntt(&s1);
 *   polyvecl_pointwise_poly_montgomery(&z, &cp, &s1);
 *   polyvecl_invntt_tomont(&z);
 *
 * @param z Output: L polynomials, each coefficient in [-TAU*ETA, TAU*ETA]
 * @param c Challenge polynomial coefficients
 * @param s1 L secret polynomials, each coefficient in [-ETA, ETA]
 */
void smallntt_polyvecl_mul(int32_t z[L][N], const int32_t c[N], const int32_t s1[L][N]);

/**
 * Compute c*s2 for polyveck (K=6 polynomials) using small NTT.
 *
 * @param h Output: K polynomials, each coefficient in [-TAU*ETA, TAU*ETA]
 * @param c Challenge polynomial coefficients
 * @param s2 K secret polynomials, each coefficient in [-ETA, ETA]
 */
void smallntt_polyveck_mul(int32_t h[K][N], const int32_t c[N], const int32_t s2[K][N]);

#endif /* PQCLEAN_MLDSA65_CLEAN_SMALLNTT_H */
