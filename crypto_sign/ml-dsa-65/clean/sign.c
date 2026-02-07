#include "fips202.h"
#include "packing.h"
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "randombytes.h"
#include "sign.h"
#include "smallntt.h"
#include "symmetric.h"
#include <stdint.h>

/*************************************************
* Name:        PQCLEAN_MLDSA65_CLEAN_crypto_sign_keypair
*
* Description: Generates public and private key.
*
* Arguments:   - uint8_t *pk: pointer to output public key (allocated
*                             array of PQCLEAN_MLDSA65_CLEAN_CRYPTO_PUBLICKEYBYTES bytes)
*              - uint8_t *sk: pointer to output private key (allocated
*                             array of PQCLEAN_MLDSA65_CLEAN_CRYPTO_SECRETKEYBYTES bytes)
*
* Returns 0 (success)
**************************************************/
int PQCLEAN_MLDSA65_CLEAN_crypto_sign_keypair(uint8_t *pk, uint8_t *sk) {
    uint8_t seedbuf[2 * SEEDBYTES + CRHBYTES];
    uint8_t tr[TRBYTES];
    const uint8_t *rho, *rhoprime, *key;
    polyvecl mat[K];
    polyvecl s1, s1hat;
    polyveck s2, t1, t0;

    /* Get randomness for rho, rhoprime and key */
    randombytes(seedbuf, SEEDBYTES);
    seedbuf[SEEDBYTES + 0] = K;
    seedbuf[SEEDBYTES + 1] = L;
    shake256(seedbuf, 2 * SEEDBYTES + CRHBYTES, seedbuf, SEEDBYTES + 2);
    rho = seedbuf;
    rhoprime = rho + SEEDBYTES;
    key = rhoprime + CRHBYTES;

    /* Expand matrix */
    PQCLEAN_MLDSA65_CLEAN_polyvec_matrix_expand(mat, rho);

    /* Sample short vectors s1 and s2 */
    PQCLEAN_MLDSA65_CLEAN_polyvecl_uniform_eta(&s1, rhoprime, 0);
    PQCLEAN_MLDSA65_CLEAN_polyveck_uniform_eta(&s2, rhoprime, L);

    /* Matrix-vector multiplication */
    s1hat = s1;
    PQCLEAN_MLDSA65_CLEAN_polyvecl_ntt(&s1hat);
    PQCLEAN_MLDSA65_CLEAN_polyvec_matrix_pointwise_montgomery(&t1, mat, &s1hat);
    PQCLEAN_MLDSA65_CLEAN_polyveck_reduce(&t1);
    PQCLEAN_MLDSA65_CLEAN_polyveck_invntt_tomont(&t1);

    /* Add error vector s2 */
    PQCLEAN_MLDSA65_CLEAN_polyveck_add(&t1, &t1, &s2);

    /* Extract t1 and write public key */
    PQCLEAN_MLDSA65_CLEAN_polyveck_caddq(&t1);
    PQCLEAN_MLDSA65_CLEAN_polyveck_power2round(&t1, &t0, &t1);
    PQCLEAN_MLDSA65_CLEAN_pack_pk(pk, rho, &t1);

    /* Compute H(rho, t1) and write secret key */
    shake256(tr, TRBYTES, pk, PQCLEAN_MLDSA65_CLEAN_CRYPTO_PUBLICKEYBYTES);
    PQCLEAN_MLDSA65_CLEAN_pack_sk(sk, rho, tr, key, &t0, &s1, &s2);

    return 0;
}

/*************************************************
* Name:        crypto_sign_signature
*
* Description: Computes signature.
*
* Arguments:   - uint8_t *sig:   pointer to output signature (of length PQCLEAN_MLDSA65_CLEAN_CRYPTO_BYTES)
*              - size_t *siglen: pointer to output length of signature
*              - uint8_t *m:     pointer to message to be signed
*              - size_t mlen:    length of message
*              - uint8_t *ctx:   pointer to context string
*              - size_t ctxlen:  length of context string
*              - uint8_t *sk:    pointer to bit-packed secret key
*
* Returns 0 (success) or -1 (context string too long)
**************************************************/
int PQCLEAN_MLDSA65_CLEAN_crypto_sign_signature_ctx(uint8_t *sig,
        size_t *siglen,
        const uint8_t *m,
        size_t mlen,
        const uint8_t *ctx,
        size_t ctxlen,
        const uint8_t *sk) {
    unsigned int i, n;
    uint8_t seedbuf[2 * SEEDBYTES + TRBYTES + RNDBYTES + 2 * CRHBYTES];
    uint8_t *rho, *tr, *key, *mu, *rhoprime, *rnd;
    uint16_t nonce = 0;
    polyvecl mat[K], s1, y;
    /* OPTIMIZATION: Removed polyvecl z - using streaming computation instead */
    polyveck t0, s2, w1, w0, h;
    poly cp;
    poly tmp;  /* Single poly buffer for streaming z computation */
    shake256incctx state;

    if (ctxlen > 255) {
        return -1;
    }

    rho = seedbuf;
    tr = rho + SEEDBYTES;
    key = tr + TRBYTES;
    rnd = key + SEEDBYTES;
    mu = rnd + RNDBYTES;
    rhoprime = mu + CRHBYTES;
    PQCLEAN_MLDSA65_CLEAN_unpack_sk(rho, tr, key, &t0, &s1, &s2, sk);

    /* Compute mu = CRH(tr, 0, ctxlen, ctx, msg) */
    mu[0] = 0;
    mu[1] = (uint8_t)ctxlen;
    shake256_inc_init(&state);
    shake256_inc_absorb(&state, tr, TRBYTES);
    shake256_inc_absorb(&state, mu, 2);
    shake256_inc_absorb(&state, ctx, ctxlen);
    shake256_inc_absorb(&state, m, mlen);
    shake256_inc_finalize(&state);
    shake256_inc_squeeze(mu, CRHBYTES, &state);
    shake256_inc_ctx_release(&state);

    randombytes(rnd, RNDBYTES);
    shake256(rhoprime, CRHBYTES, key, SEEDBYTES + RNDBYTES + CRHBYTES);

    /* Expand matrix and transform vectors */
    PQCLEAN_MLDSA65_CLEAN_polyvec_matrix_expand(mat, rho);
    /* Note: s1 and s2 are kept in normal domain for small NTT optimization */
    /* Only t0 is transformed to NTT domain for c*t0 computation */
    PQCLEAN_MLDSA65_CLEAN_polyveck_ntt(&t0);

rej:
    /* Sample intermediate vector y */
    PQCLEAN_MLDSA65_CLEAN_polyvecl_uniform_gamma1(&y, rhoprime, nonce++);

    /* Matrix-vector multiplication: compute w = A*NTT(y) */
    /* Use tmp poly as scratch for each polynomial NTT */
    for (i = 0; i < L; ++i) {
        tmp = y.vec[i];
        PQCLEAN_MLDSA65_CLEAN_poly_ntt(&tmp);
        /* Accumulate into w1 row by row */
        for (unsigned int j = 0; j < K; ++j) {
            poly acc;
            PQCLEAN_MLDSA65_CLEAN_poly_pointwise_montgomery(&acc, &mat[j].vec[i], &tmp);
            if (i == 0) {
                w1.vec[j] = acc;
            } else {
                PQCLEAN_MLDSA65_CLEAN_poly_add(&w1.vec[j], &w1.vec[j], &acc);
            }
        }
    }
    PQCLEAN_MLDSA65_CLEAN_polyveck_reduce(&w1);
    PQCLEAN_MLDSA65_CLEAN_polyveck_invntt_tomont(&w1);

    /* Decompose w and call the random oracle */
    PQCLEAN_MLDSA65_CLEAN_polyveck_caddq(&w1);
    PQCLEAN_MLDSA65_CLEAN_polyveck_decompose(&w1, &w0, &w1);
    PQCLEAN_MLDSA65_CLEAN_polyveck_pack_w1(sig, &w1);

    shake256_inc_init(&state);
    shake256_inc_absorb(&state, mu, CRHBYTES);
    shake256_inc_absorb(&state, sig, K * POLYW1_PACKEDBYTES);
    shake256_inc_finalize(&state);
    shake256_inc_squeeze(sig, CTILDEBYTES, &state);
    shake256_inc_ctx_release(&state);
    PQCLEAN_MLDSA65_CLEAN_poly_challenge(&cp, sig);

    /*
     * COMBINED OPTIMIZATIONS:
     * 1. Streaming z computation with early rejection (eliminates polyvecl z ~5KB)
     * 2. Small NTT for c*s1 and c*s2 (16-bit coefficients, q=7681)
     *
     * Two-pass approach for z:
     * Pass 1: Compute each z[i] using small NTT, check norm immediately
     *         Abort as soon as any coefficient exceeds bound (early rejection)
     * Pass 2: If all pass, recompute z[i] and pack directly into signature
     */
    {
        smallpoly c_ntt_small;
        
        /* Convert challenge to small NTT domain once - reused for all c*s multiplications */
        smallntt_challenge_to_ntt(&c_ntt_small, cp.coeffs);
        
        /* Pass 1: Stream compute z = y + c*s1 and check norm with early rejection */
        for (i = 0; i < L; ++i) {
            /* Compute tmp = c * s1[i] using small NTT */
            smallntt_mul_poly(tmp.coeffs, &c_ntt_small, s1.vec[i].coeffs);
            /* Add y[i]: z[i] = c*s1[i] + y[i] */
            for (unsigned int j = 0; j < N; ++j) {
                tmp.coeffs[j] += y.vec[i].coeffs[j];
            }
            PQCLEAN_MLDSA65_CLEAN_poly_reduce(&tmp);
            /* Check norm - early rejection if any polynomial fails */
            if (PQCLEAN_MLDSA65_CLEAN_poly_chknorm(&tmp, GAMMA1 - BETA)) {
                goto rej;
            }
        }
        
        /* Check that subtracting cs2 does not change high bits of w and low bits
         * do not reveal secret information - using small NTT for c*s2 */
        for (i = 0; i < K; ++i) {
            smallntt_mul_poly(h.vec[i].coeffs, &c_ntt_small, s2.vec[i].coeffs);
        }
    }
    PQCLEAN_MLDSA65_CLEAN_polyveck_sub(&w0, &w0, &h);
    PQCLEAN_MLDSA65_CLEAN_polyveck_reduce(&w0);
    if (PQCLEAN_MLDSA65_CLEAN_polyveck_chknorm(&w0, GAMMA2 - BETA)) {
        goto rej;
    }

    /* Compute hints for w1 - use standard NTT for c*t0 (t0 has larger coefficients) */
    PQCLEAN_MLDSA65_CLEAN_poly_ntt(&cp);
    PQCLEAN_MLDSA65_CLEAN_polyveck_pointwise_poly_montgomery(&h, &cp, &t0);
    PQCLEAN_MLDSA65_CLEAN_polyveck_invntt_tomont(&h);
    PQCLEAN_MLDSA65_CLEAN_polyveck_reduce(&h);
    if (PQCLEAN_MLDSA65_CLEAN_polyveck_chknorm(&h, GAMMA2)) {
        goto rej;
    }

    PQCLEAN_MLDSA65_CLEAN_polyveck_add(&w0, &w0, &h);
    n = PQCLEAN_MLDSA65_CLEAN_polyveck_make_hint(&h, &w0, &w1);
    if (n > OMEGA) {
        goto rej;
    }

    /*
     * Write signature with streaming z packing
     * 
     * Pass 2: Recompute z using small NTT and pack directly into signature buffer
     * sig layout: c (CTILDEBYTES) | z (L * POLYZ_PACKEDBYTES) | h (OMEGA + K)
     * 
     * Note: sig already contains c at the beginning from earlier squeeze
     * Note: cp was modified by poly_ntt() for c*t0, so we reconstruct challenge
     */
    {
        uint8_t *sig_z = sig + CTILDEBYTES;  /* Point to z portion of signature */
        smallpoly c_ntt_small;
        
        /* Reconstruct challenge from sig and convert to small NTT domain */
        PQCLEAN_MLDSA65_CLEAN_poly_challenge(&cp, sig);
        smallntt_challenge_to_ntt(&c_ntt_small, cp.coeffs);
        
        /* Recompute each z[i] using small NTT and pack directly */
        for (i = 0; i < L; ++i) {
            /* Recompute tmp = c * s1[i] using small NTT */
            smallntt_mul_poly(tmp.coeffs, &c_ntt_small, s1.vec[i].coeffs);
            /* Add y[i] */
            for (unsigned int j = 0; j < N; ++j) {
                tmp.coeffs[j] += y.vec[i].coeffs[j];
            }
            PQCLEAN_MLDSA65_CLEAN_poly_reduce(&tmp);
            /* Pack directly into signature buffer */
            PQCLEAN_MLDSA65_CLEAN_polyz_pack(sig_z + i * POLYZ_PACKEDBYTES, &tmp);
        }
        
        /* Pack h (hint vector) */
        {
            uint8_t *sig_h = sig + CTILDEBYTES + L * POLYZ_PACKEDBYTES;
            unsigned int j, k_idx;
            
            /* Zero out h portion */
            for (i = 0; i < OMEGA + K; ++i) {
                sig_h[i] = 0;
            }
            
            k_idx = 0;
            for (i = 0; i < K; ++i) {
                for (j = 0; j < N; ++j) {
                    if (h.vec[i].coeffs[j] != 0) {
                        sig_h[k_idx++] = (uint8_t) j;
                    }
                }
                sig_h[OMEGA + i] = (uint8_t) k_idx;
            }
        }
    }
    
    *siglen = PQCLEAN_MLDSA65_CLEAN_CRYPTO_BYTES;
    return 0;
}

/*************************************************
* Name:        crypto_sign
*
* Description: Compute signed message.
*
* Arguments:   - uint8_t *sm: pointer to output signed message (allocated
*                             array with PQCLEAN_MLDSA65_CLEAN_CRYPTO_BYTES + mlen bytes),
*                             can be equal to m
*              - size_t *smlen: pointer to output length of signed
*                               message
*              - const uint8_t *m: pointer to message to be signed
*              - size_t mlen: length of message
*              - const uint8_t *ctx: pointer to context string
*              - size_t ctxlen: length of context string
*              - const uint8_t *sk: pointer to bit-packed secret key
*
* Returns 0 (success) or -1 (context string too long)
**************************************************/
int PQCLEAN_MLDSA65_CLEAN_crypto_sign_ctx(uint8_t *sm,
        size_t *smlen,
        const uint8_t *m,
        size_t mlen,
        const uint8_t *ctx,
        size_t ctxlen,
        const uint8_t *sk) {
    int ret;
    size_t i;

    for (i = 0; i < mlen; ++i) {
        sm[PQCLEAN_MLDSA65_CLEAN_CRYPTO_BYTES + mlen - 1 - i] = m[mlen - 1 - i];
    }
    ret = PQCLEAN_MLDSA65_CLEAN_crypto_sign_signature_ctx(sm, smlen, sm + PQCLEAN_MLDSA65_CLEAN_CRYPTO_BYTES, mlen, ctx, ctxlen, sk);
    *smlen += mlen;
    return ret;
}

/*************************************************
* Name:        crypto_sign_verify
*
* Description: Verifies signature.
*
* Arguments:   - uint8_t *m: pointer to input signature
*              - size_t siglen: length of signature
*              - const uint8_t *m: pointer to message
*              - size_t mlen: length of message
*              - const uint8_t *ctx: pointer to context string
*              - size_t ctxlen: length of context string
*              - const uint8_t *pk: pointer to bit-packed public key
*
* Returns 0 if signature could be verified correctly and -1 otherwise
**************************************************/
int PQCLEAN_MLDSA65_CLEAN_crypto_sign_verify_ctx(const uint8_t *sig,
        size_t siglen,
        const uint8_t *m,
        size_t mlen,
        const uint8_t *ctx,
        size_t ctxlen,
        const uint8_t *pk) {
    unsigned int i;
    uint8_t buf[K * POLYW1_PACKEDBYTES];
    uint8_t rho[SEEDBYTES];
    uint8_t mu[CRHBYTES];
    uint8_t c[CTILDEBYTES];
    uint8_t c2[CTILDEBYTES];
    poly cp;
    polyvecl mat[K], z;
    polyveck t1, w1, h;
    shake256incctx state;

    if (ctxlen > 255 || siglen != PQCLEAN_MLDSA65_CLEAN_CRYPTO_BYTES) {
        return -1;
    }

    PQCLEAN_MLDSA65_CLEAN_unpack_pk(rho, &t1, pk);
    if (PQCLEAN_MLDSA65_CLEAN_unpack_sig(c, &z, &h, sig)) {
        return -1;
    }
    if (PQCLEAN_MLDSA65_CLEAN_polyvecl_chknorm(&z, GAMMA1 - BETA)) {
        return -1;
    }

    /* Compute CRH(H(rho, t1), msg) */
    shake256(mu, TRBYTES, pk, PQCLEAN_MLDSA65_CLEAN_CRYPTO_PUBLICKEYBYTES);
    shake256_inc_init(&state);
    shake256_inc_absorb(&state, mu, TRBYTES);
    mu[0] = 0;
    mu[1] = (uint8_t)ctxlen;
    shake256_inc_absorb(&state, mu, 2);
    shake256_inc_absorb(&state, ctx, ctxlen);
    shake256_inc_absorb(&state, m, mlen);
    shake256_inc_finalize(&state);
    shake256_inc_squeeze(mu, CRHBYTES, &state);
    shake256_inc_ctx_release(&state);

    /* Matrix-vector multiplication; compute Az - c2^dt1 */
    PQCLEAN_MLDSA65_CLEAN_poly_challenge(&cp, c);
    PQCLEAN_MLDSA65_CLEAN_polyvec_matrix_expand(mat, rho);

    PQCLEAN_MLDSA65_CLEAN_polyvecl_ntt(&z);
    PQCLEAN_MLDSA65_CLEAN_polyvec_matrix_pointwise_montgomery(&w1, mat, &z);

    PQCLEAN_MLDSA65_CLEAN_poly_ntt(&cp);
    PQCLEAN_MLDSA65_CLEAN_polyveck_shiftl(&t1);
    PQCLEAN_MLDSA65_CLEAN_polyveck_ntt(&t1);
    PQCLEAN_MLDSA65_CLEAN_polyveck_pointwise_poly_montgomery(&t1, &cp, &t1);

    PQCLEAN_MLDSA65_CLEAN_polyveck_sub(&w1, &w1, &t1);
    PQCLEAN_MLDSA65_CLEAN_polyveck_reduce(&w1);
    PQCLEAN_MLDSA65_CLEAN_polyveck_invntt_tomont(&w1);

    /* Reconstruct w1 */
    PQCLEAN_MLDSA65_CLEAN_polyveck_caddq(&w1);
    PQCLEAN_MLDSA65_CLEAN_polyveck_use_hint(&w1, &w1, &h);
    PQCLEAN_MLDSA65_CLEAN_polyveck_pack_w1(buf, &w1);

    /* Call random oracle and verify challenge */
    shake256_inc_init(&state);
    shake256_inc_absorb(&state, mu, CRHBYTES);
    shake256_inc_absorb(&state, buf, K * POLYW1_PACKEDBYTES);
    shake256_inc_finalize(&state);
    shake256_inc_squeeze(c2, CTILDEBYTES, &state);
    shake256_inc_ctx_release(&state);
    for (i = 0; i < CTILDEBYTES; ++i) {
        if (c[i] != c2[i]) {
            return -1;
        }
    }

    return 0;
}

/*************************************************
* Name:        crypto_sign_open
*
* Description: Verify signed message.
*
* Arguments:   - uint8_t *m: pointer to output message (allocated
*                            array with smlen bytes), can be equal to sm
*              - size_t *mlen: pointer to output length of message
*              - const uint8_t *sm: pointer to signed message
*              - size_t smlen: length of signed message
*              - const uint8_t *ctx: pointer to context tring
*              - size_t ctxlen: length of context string
*              - const uint8_t *pk: pointer to bit-packed public key
*
* Returns 0 if signed message could be verified correctly and -1 otherwise
**************************************************/
int PQCLEAN_MLDSA65_CLEAN_crypto_sign_open_ctx(uint8_t *m,
        size_t *mlen,
        const uint8_t *sm,
        size_t smlen,
        const uint8_t *ctx,
        size_t ctxlen,
        const uint8_t *pk) {
    size_t i;

    if (smlen < PQCLEAN_MLDSA65_CLEAN_CRYPTO_BYTES) {
        goto badsig;
    }

    *mlen = smlen - PQCLEAN_MLDSA65_CLEAN_CRYPTO_BYTES;
    if (PQCLEAN_MLDSA65_CLEAN_crypto_sign_verify_ctx(sm, PQCLEAN_MLDSA65_CLEAN_CRYPTO_BYTES, sm + PQCLEAN_MLDSA65_CLEAN_CRYPTO_BYTES, *mlen, ctx, ctxlen, pk)) {
        goto badsig;
    } else {
        /* All good, copy msg, return 0 */
        for (i = 0; i < *mlen; ++i) {
            m[i] = sm[PQCLEAN_MLDSA65_CLEAN_CRYPTO_BYTES + i];
        }
        return 0;
    }

badsig:
    /* Signature verification failed */
    *mlen = 0;
    for (i = 0; i < smlen; ++i) {
        m[i] = 0;
    }

    return -1;
}

int PQCLEAN_MLDSA65_CLEAN_crypto_sign_signature(uint8_t *sig,
        size_t *siglen,
        const uint8_t *m,
        size_t mlen,
        const uint8_t *sk) {
    return PQCLEAN_MLDSA65_CLEAN_crypto_sign_signature_ctx(sig, siglen, m, mlen, NULL, 0, sk);
}

int PQCLEAN_MLDSA65_CLEAN_crypto_sign(uint8_t *sm,
                                      size_t *smlen,
                                      const uint8_t *m,
                                      size_t mlen,
                                      const uint8_t *sk) {
    return PQCLEAN_MLDSA65_CLEAN_crypto_sign_ctx(sm, smlen, m, mlen, NULL, 0, sk);
}

int PQCLEAN_MLDSA65_CLEAN_crypto_sign_verify(const uint8_t *sig,
        size_t siglen,
        const uint8_t *m,
        size_t mlen,
        const uint8_t *pk) {
    return PQCLEAN_MLDSA65_CLEAN_crypto_sign_verify_ctx(sig, siglen, m, mlen, NULL, 0, pk);
}

int PQCLEAN_MLDSA65_CLEAN_crypto_sign_open(uint8_t *m,
        size_t *mlen,
        const uint8_t *sm, size_t smlen,
        const uint8_t *pk) {
    return PQCLEAN_MLDSA65_CLEAN_crypto_sign_open_ctx(m, mlen, sm, smlen, NULL, 0, pk);
}
