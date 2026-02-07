/**
 * Small-modulus NTT for memory-efficient c*s1 and c*s2 computation in ML-DSA-65.
 * 
 * Uses twist-based negacyclic NTT with q=7681 = 15*512 + 1.
 * This allows 16-bit coefficient storage instead of 32-bit.
 * 
 * Algorithm:
 *   Forward NTT: twist by psi^i, then cyclic NTT with omega
 *   Inverse NTT: inverse cyclic NTT with omega, then untwist by psi^(-i)
 * 
 * Where psi is a primitive 512th root of unity and omega = psi^2.
 */

#include "smallntt.h"
#include "params.h"
#include <stdint.h>

/* Twist factors: psi^i for i = 0..255 where psi = 7146 */
static const int16_t psi_powers[N] = {
        1,  7146,  2028,  5722,  3449,  5906,  4862,  2689,
     5413,  7463,  1415,  3394,  4607,   856,  2900,    62,
     5235,  2840,  1438,  6451,  5165,  1885,  5417,  5323,
     1846,  3239,  3041,  1437,  6986,  3137,  3844,  1968,
     7098,  4665,   550,  5309,  1655,  5571,  7424,  6918,
     1112,  4198,  4603,  2996,  2469,   217,  6801,  2259,
     5033,  3376,  6556,  2757,  7438,  7109,  6461,  7496,
     6803,  1189,  1408,  7139,  5773,  6888,  1800,  4806,
     1925,  7060,  1952,   296,  2941,  1170,  3892,  7012,
     4589,  2805,  4801,  4600,  4601,  4066,  6094,  4135,
     7584,  5809,  2990,  5679,  3411,  3193,  4608,   321,
     4928,  5784,  1003,  1065,  6300,  1459,  2897,  1667,
     6832,  1036,  6453,  4095,  5941,  1499,  4540,  5977,
     5282,   738,  4582,  6550,  5967,  2951,  3501,  1129,
     2784,   674,   417,  7335,   766,  4964,  1886,  4882,
     7351,  7568,  6688,  1266,  6299,  1994,   869,  3626,
     3383,  2811,  1591,  1406,   528,  1717,  3125,  2583,
      675,  7563,  1682,  6488,   732,   111,  2063,  2359,
     5300,  6470,  2681,  2012,  6601,  1725,  6526,  3445,
      365,  4431,  2844,  6979,  6882,  5010,   319,  5998,
     1728,  4921,  1848,  2169,  7097,  5200,  6203,  7268,
     5887,  7346,  2562,  4229,  3380,  4416,  3188,  7283,
     5543,  7042,  3901,  2197,  7479,   536,  5118,  3987,
     2273,  5224,  1044,  2173,  4957,  5631,  6048,  5702,
     6468,  3751,  5637,  2838,  2508,  2395,  1402,  2668,
     1286,  3280,  4149,    94,  3477,  6288,   198,  1604,
     2132,  3849,  6974,  1876,  2551,  2433,  4115,  2922,
     3654,  3765,  5828,   506,  5806,  4595,  7276,  1607,
      527,  2252,  1097,  4542,  4907,  1657,  4501,  3799,
     3000,   329,   648,  6646,   693,  5614,  7462,  1950,
     1366,  6566,  5088,  4675,  2881,  2546,  5108,  1656,
     5036,  1771,  4959,  4561,  2423,  1784,  5685,   201
};

/* Untwist factors: psi^(-i) for i = 0..255 */
static const int16_t psi_inv_powers[N] = {
        1,  7480,  1996,  5897,  5258,  3120,  2722,  5910,
     2645,  6025,  2573,  5135,  4800,  3006,  2593,  1115,
     6315,  5731,   219,  2067,  6988,  1035,  7033,  7352,
     4681,  3882,  3180,  6024,  2774,  3139,  6584,  5429,
     7154,  6074,   405,  3086,  1875,  7175,  1853,  3916,
     4027,  4759,  3566,  5248,  5130,  5805,   707,  3832,
     5549,  6077,  7483,  1393,  4204,  7587,  3532,  4401,
     6395,  5013,  6279,  5286,  5173,  4843,  2044,  3930,
     1213,  1979,  1633,  2050,  2724,  5508,  6637,  2457,
     5408,  3694,  2563,  7145,   202,  5484,  3780,   639,
     2138,   398,  4493,  3265,  4301,  3452,  5119,   335,
     1794,   413,  1478,  2481,   584,  5512,  5833,  2760,
     5953,  1683,  7362,  2671,   799,   702,  4837,  3250,
     7316,  4236,  1155,  5956,  1080,  5669,  5000,  1211,
     2381,  5322,  5618,  7570,  6949,  1193,  5999,   118,
     7006,  5098,  4556,  5964,  7153,  6275,  6090,  4870,
     4298,  4055,  6812,  5687,  1382,  6415,   993,   113,
      330,  2799,  5795,  2717,  6915,   346,  7264,  7007,
     4897,  6552,  4180,  4730,  1714,  1131,  3099,  6943,
     2399,  1704,  3141,  6182,  1740,  3586,  1228,  6645,
      849,  6014,  4784,  6222,  1381,  6616,  6678,  1897,
     2753,  7360,  3073,  4488,  4270,  2002,  4691,  1872,
       97,  3546,  1587,  3615,  3080,  3081,  2880,  4876,
     3092,   669,  3789,  6511,  4740,  7385,  5729,   621,
     5756,  2875,  5881,   793,  1908,   542,  6273,  6492,
      878,   185,  1220,   572,   243,  4924,  1125,  4305,
     2648,  5422,   880,  7464,  5212,  4685,  3078,  3483,
     6569,   763,   257,  2110,  6026,  2372,  7131,  3016,
      583,  5713,  3837,  4544,   695,  6244,  4640,  4442,
     5835,  2358,  2264,  5796,  2516,  1230,  6243,  4841,
     2446,  7619,  4781,  6825,  3074,  4287,  6266,   218,
     2268,  4992,  2819,  1775,  4232,  1959,  5653,   535
};

/* NTT twiddles: omega^k for k = 0..255 where omega = 2028 */
static const int16_t omega_powers[N] = {
        1,  2028,  3449,  4862,  5413,  1415,  4607,  2900,
     5235,  1438,  5165,  5417,  1846,  3041,  6986,  3844,
     7098,   550,  1655,  7424,  1112,  4603,  2469,  6801,
     5033,  6556,  7438,  6461,  6803,  1408,  5773,  1800,
     1925,  1952,  2941,  3892,  4589,  4801,  4601,  6094,
     7584,  2990,  3411,  4608,  4928,  1003,  6300,  2897,
     6832,  6453,  5941,  4540,  5282,  4582,  5967,  3501,
     2784,   417,   766,  1886,  7351,  6688,  6299,   869,
     3383,  1591,   528,  3125,   675,  1682,   732,  2063,
     5300,  2681,  6601,  6526,   365,  2844,  6882,   319,
     1728,  1848,  7097,  6203,  5887,  2562,  3380,  3188,
     5543,  3901,  7479,  5118,  2273,  1044,  4957,  6048,
     6468,  5637,  2508,  1402,  1286,  4149,  3477,   198,
     2132,  6974,  2551,  4115,  3654,  5828,  5806,  7276,
      527,  1097,  4907,  4501,  3000,   648,   693,  7462,
     1366,  5088,  2881,  5108,  5036,  4959,  2423,  5685,
     7680,  5653,  4232,  2819,  2268,  6266,  3074,  4781,
     2446,  6243,  2516,  2264,  5835,  4640,   695,  3837,
      583,  7131,  6026,   257,  6569,  3078,  5212,   880,
     2648,  1125,   243,  1220,   878,  6273,  1908,  5881,
     5756,  5729,  4740,  3789,  3092,  2880,  3080,  1587,
       97,  4691,  4270,  3073,  2753,  6678,  1381,  4784,
      849,  1228,  1740,  3141,  2399,  3099,  1714,  4180,
     4897,  7264,  6915,  5795,   330,   993,  1382,  6812,
     4298,  6090,  7153,  4556,  7006,  5999,  6949,  5618,
     2381,  5000,  1080,  1155,  7316,  4837,   799,  7362,
     5953,  5833,   584,  1478,  1794,  5119,  4301,  4493,
     2138,  3780,   202,  2563,  5408,  6637,  2724,  1633,
     1213,  2044,  5173,  6279,  6395,  3532,  4204,  7483,
     5549,   707,  5130,  3566,  4027,  1853,  1875,   405,
     7154,  6584,  2774,  3180,  4681,  7033,  6988,   219,
     6315,  2593,  4800,  2573,  2645,  2722,  5258,  1996
};

/* N^(-1) mod q = 7651 */
#define SMALL_N_INV 7651

/*
 * Bit-reverse an 8-bit number
 */
static unsigned int bitrev8(unsigned int x) {
    x = ((x & 0xaa) >> 1) | ((x & 0x55) << 1);
    x = ((x & 0xcc) >> 2) | ((x & 0x33) << 2);
    x = ((x & 0xf0) >> 4) | ((x & 0x0f) << 4);
    return x;
}

/*
 * Reduce a 32-bit value mod SMALL_Q to 16-bit
 */
static int16_t reduce_small(int32_t a) {
    int32_t t = a % SMALL_Q;
    if (t < 0) {
        t += SMALL_Q;
    }
    return (int16_t)t;
}

/*
 * Convert polynomial from standard (32-bit) to small (16-bit) format
 * Input coefficients must be in range [-eta, eta] = [-4, 4]
 */
void poly_to_smallpoly(smallpoly *b, const int32_t *a) {
    unsigned int i;
    for (i = 0; i < N; i++) {
        int32_t coeff = a[i];
        /* Reduce mod SMALL_Q, handling negative values */
        if (coeff < 0) {
            coeff += SMALL_Q;
        }
        b->coeffs[i] = (int16_t)coeff;
    }
}

/*
 * Convert polynomial from small (16-bit) to standard (32-bit) format
 * Produces centered output in range [-(q-1)/2, (q-1)/2]
 */
void smallpoly_to_poly(int32_t *b, const smallpoly *a) {
    unsigned int i;
    for (i = 0; i < N; i++) {
        int16_t coeff = a->coeffs[i];
        /* Reduce to [0, q-1] first */
        while (coeff < 0) {
            coeff += SMALL_Q;
        }
        while (coeff >= SMALL_Q) {
            coeff -= SMALL_Q;
        }
        /* Center the result */
        if (coeff > SMALL_Q / 2) {
            b[i] = coeff - SMALL_Q;
        } else {
            b[i] = coeff;
        }
    }
}

/*
 * Forward NTT using twist-based negacyclic approach.
 * Step 1: Multiply by psi^i (twist)
 * Step 2: Bit-reversal permutation
 * Step 3: Cyclic NTT with omega
 */
void smallpoly_ntt(smallpoly *a) {
    unsigned int i, j, k;
    unsigned int len, start;
    int32_t u, v;
    int16_t w;
    
    /* Step 1: Twist - multiply each coefficient by psi^i */
    for (i = 0; i < N; i++) {
        int32_t tmp = (int32_t)a->coeffs[i] * (int32_t)psi_powers[i];
        a->coeffs[i] = reduce_small(tmp);
    }
    
    /* Step 2: Bit-reversal permutation */
    for (i = 0; i < N; i++) {
        j = bitrev8(i);
        if (i < j) {
            int16_t tmp = a->coeffs[i];
            a->coeffs[i] = a->coeffs[j];
            a->coeffs[j] = tmp;
        }
    }
    
    /* Step 3: Cooley-Tukey NTT butterflies with omega */
    for (len = 2; len <= N; len <<= 1) {
        unsigned int half = len >> 1;
        unsigned int step = N / len;
        
        for (start = 0; start < N; start += len) {
            for (k = 0; k < half; k++) {
                w = omega_powers[k * step];
                u = a->coeffs[start + k];
                v = ((int32_t)a->coeffs[start + k + half] * (int32_t)w) % SMALL_Q;
                a->coeffs[start + k] = reduce_small(u + v);
                a->coeffs[start + k + half] = reduce_small(u - v);
            }
        }
    }
}

/*
 * Inverse NTT using twist-based negacyclic approach.
 * Step 1: Gentleman-Sande INTT butterflies with omega
 * Step 2: Bit-reversal permutation
 * Step 3: Divide by N
 * Step 4: Multiply by psi^(-i) (untwist)
 */
void smallpoly_invntt(smallpoly *a) {
    unsigned int i, j, k;
    unsigned int len, start;
    int32_t t, u, v;
    int16_t w;
    
    /* Step 1: Gentleman-Sande INTT butterflies */
    for (len = N; len >= 2; len >>= 1) {
        unsigned int half = len >> 1;
        unsigned int step = N / len;
        
        for (start = 0; start < N; start += len) {
            for (k = 0; k < half; k++) {
                /* omega^(-k*step) = omega^(N - k*step) */
                unsigned int widx = (N - k * step) % N;
                w = omega_powers[widx];
                u = a->coeffs[start + k];
                v = a->coeffs[start + k + half];
                a->coeffs[start + k] = reduce_small(u + v);
                t = ((int32_t)(u - v + SMALL_Q) * (int32_t)w) % SMALL_Q;
                a->coeffs[start + k + half] = reduce_small(t);
            }
        }
    }
    
    /* Step 2: Bit-reversal permutation */
    for (i = 0; i < N; i++) {
        j = bitrev8(i);
        if (i < j) {
            int16_t tmp = a->coeffs[i];
            a->coeffs[i] = a->coeffs[j];
            a->coeffs[j] = tmp;
        }
    }
    
    /* Step 3: Divide by N and Step 4: Untwist */
    for (i = 0; i < N; i++) {
        /* Divide by N */
        int32_t tmp = ((int32_t)a->coeffs[i] * SMALL_N_INV) % SMALL_Q;
        /* Untwist: multiply by psi^(-i) */
        tmp = (tmp * (int32_t)psi_inv_powers[i]) % SMALL_Q;
        a->coeffs[i] = reduce_small(tmp);
    }
}

/*
 * Pointwise multiplication in NTT domain
 */
void smallpoly_pointwise(smallpoly *c, const smallpoly *a, const smallpoly *b) {
    unsigned int i;
    for (i = 0; i < N; i++) {
        int32_t tmp = (int32_t)a->coeffs[i] * (int32_t)b->coeffs[i];
        c->coeffs[i] = reduce_small(tmp);
    }
}

/*
 * Full polynomial multiplication using small NTT:
 *   c = a * b mod (X^N + 1) mod SMALL_Q
 * 
 * Input a and b should be in coefficient domain.
 * Output c is in coefficient domain.
 */
void smallntt_poly_mul(smallpoly *c, const smallpoly *a, const smallpoly *b) {
    smallpoly a_ntt, b_ntt;
    unsigned int i;
    
    /* Copy inputs */
    for (i = 0; i < N; i++) {
        a_ntt.coeffs[i] = a->coeffs[i];
        b_ntt.coeffs[i] = b->coeffs[i];
    }
    
    /* Transform to NTT domain */
    smallpoly_ntt(&a_ntt);
    smallpoly_ntt(&b_ntt);
    
    /* Pointwise multiplication */
    smallpoly_pointwise(c, &a_ntt, &b_ntt);
    
    /* Transform back */
    smallpoly_invntt(c);
}

/*************************************************
 * Integration functions for ML-DSA-65
 *************************************************/

/*
 * Convert challenge polynomial (from normal domain) to small NTT form.
 * Challenge c has TAU=49 non-zero entries in {-1, 1}.
 */
void smallntt_challenge_to_ntt(smallpoly *c_ntt, const int32_t *c) {
    unsigned int i;
    
    /* Convert to smallpoly */
    for (i = 0; i < N; i++) {
        int32_t coeff = c[i];
        if (coeff < 0) {
            coeff += SMALL_Q;
        }
        c_ntt->coeffs[i] = (int16_t)coeff;
    }
    
    /* Transform to NTT domain */
    smallpoly_ntt(c_ntt);
}

/*
 * Compute c*s for a single polynomial using small NTT.
 * c_ntt: challenge already in small NTT domain
 * s: secret polynomial in normal domain (coefficients in [-ETA, ETA])
 * result: output in normal domain (coefficients in [-TAU*ETA, TAU*ETA])
 */
void smallntt_mul_poly(int32_t *result, const smallpoly *c_ntt, const int32_t *s) {
    smallpoly s_small, prod;
    unsigned int i;
    
    /* Convert secret to smallpoly and transform to NTT */
    for (i = 0; i < N; i++) {
        int32_t coeff = s[i];
        if (coeff < 0) {
            coeff += SMALL_Q;
        }
        s_small.coeffs[i] = (int16_t)coeff;
    }
    smallpoly_ntt(&s_small);
    
    /* Pointwise multiplication in NTT domain */
    smallpoly_pointwise(&prod, c_ntt, &s_small);
    
    /* Transform back to normal domain */
    smallpoly_invntt(&prod);
    
    /* Convert to 32-bit centered representation */
    for (i = 0; i < N; i++) {
        int16_t coeff = prod.coeffs[i];
        /* Reduce to [0, q-1] */
        while (coeff < 0) {
            coeff += SMALL_Q;
        }
        while (coeff >= SMALL_Q) {
            coeff -= SMALL_Q;
        }
        /* Center to [-(q-1)/2, (q-1)/2] */
        if (coeff > SMALL_Q / 2) {
            result[i] = coeff - SMALL_Q;
        } else {
            result[i] = coeff;
        }
    }
}

/*
 * Compute c*s1 for polyvecl (L=5 polynomials) using small NTT.
 * Challenge c should be in normal domain.
 * s1 should be in normal domain (coefficients in [-ETA, ETA]).
 * Result z is in normal domain.
 * 
 * This replaces:
 *   poly_ntt(&cp);
 *   polyvecl_ntt(&s1);
 *   polyvecl_pointwise_poly_montgomery(&z, &cp, &s1);
 *   polyvecl_invntt_tomont(&z);
 */
void smallntt_polyvecl_mul(int32_t z[L][N], const int32_t c[N], const int32_t s1[L][N]) {
    smallpoly c_ntt;
    unsigned int i;
    
    /* Convert challenge to small NTT form (done once) */
    smallntt_challenge_to_ntt(&c_ntt, c);
    
    /* Multiply challenge with each polynomial in s1 */
    for (i = 0; i < L; i++) {
        smallntt_mul_poly(z[i], &c_ntt, s1[i]);
    }
}

/*
 * Compute c*s2 for polyveck (K=6 polynomials) using small NTT.
 * Same as above but for K polynomials.
 */
void smallntt_polyveck_mul(int32_t h[K][N], const int32_t c[N], const int32_t s2[K][N]) {
    smallpoly c_ntt;
    unsigned int i;
    
    /* Convert challenge to small NTT form (done once) */
    smallntt_challenge_to_ntt(&c_ntt, c);
    
    /* Multiply challenge with each polynomial in s2 */
    for (i = 0; i < K; i++) {
        smallntt_mul_poly(h[i], &c_ntt, s2[i]);
    }
}

