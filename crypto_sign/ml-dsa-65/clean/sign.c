#include "fips202.h"
#include "packing.h"
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "randombytes.h"
#include "sign.h"
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
    unsigned int i, j, n;
    uint8_t seedbuf[2 * SEEDBYTES + TRBYTES + RNDBYTES + 2 * CRHBYTES];
    uint8_t *rho, *tr, *key, *mu, *rhoprime, *rnd;
    uint16_t nonce = 0;
    uint16_t nonce_y;  /* Save nonce for y regeneration */
    /* ALL OPTIMIZATIONS COMBINED:
     * - Removed polyvecl mat[K] (~30 KiB): A streamed on-the-fly from rho
     * - Removed polyvecl y (~5 KiB): y generated one poly at a time
     * - Removed polyvecl z (~5 KiB): z streamed with early rejection
     * - Removed polyveck w0 (~6 KiB): w recomputed from saved nonce */
    polyvecl s1;
    polyveck t0, s2, w1, h;
    poly cp;
    poly tmp;  /* Single poly buffer for streaming A/y/z computation */
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

    /* Transform vectors to NTT domain (no matrix expansion!) */
    PQCLEAN_MLDSA65_CLEAN_polyvecl_ntt(&s1);
    PQCLEAN_MLDSA65_CLEAN_polyveck_ntt(&s2);
    PQCLEAN_MLDSA65_CLEAN_polyveck_ntt(&t0);

rej:
    /* Save nonce for y regeneration later */
    nonce_y = nonce++;

    /*
     * STREAMING A AND y OPTIMIZATION:
     * Instead of storing full matrix A (K*L*1KiB = 30 KiB) and full y (L*1KiB = 5 KiB),
     * generate each A[i][j] and y[j] on-the-fly.
     *
     * For each column j: generate y[j], NTT it, then for each row i:
     *   generate A[i][j], multiply with NTT(y[j]), accumulate into w1[i]
     * This way only one poly of y is live at a time.
     */
    for (j = 0; j < L; ++j) {
        /* Generate y[j] and NTT it */
        PQCLEAN_MLDSA65_CLEAN_poly_uniform_gamma1(&tmp, rhoprime, (uint16_t)(L * nonce_y + j));
        PQCLEAN_MLDSA65_CLEAN_poly_ntt(&tmp);
        for (i = 0; i < K; ++i) {
            poly a_ij;
            PQCLEAN_MLDSA65_CLEAN_poly_uniform(&a_ij, rho, (uint16_t)((i << 8) + j));
            PQCLEAN_MLDSA65_CLEAN_poly_pointwise_montgomery(&a_ij, &a_ij, &tmp);
            if (j == 0) {
                w1.vec[i] = a_ij;
            } else {
                PQCLEAN_MLDSA65_CLEAN_poly_add(&w1.vec[i], &w1.vec[i], &a_ij);
            }
        }
    }
    PQCLEAN_MLDSA65_CLEAN_polyveck_reduce(&w1);
    PQCLEAN_MLDSA65_CLEAN_polyveck_invntt_tomont(&w1);

    /* Decompose w - discard low bits (w0), we'll recompute w later when needed */
    PQCLEAN_MLDSA65_CLEAN_polyveck_caddq(&w1);
    PQCLEAN_MLDSA65_CLEAN_polyveck_decompose(&w1, &h, &w1);  /* h = temp low bits, discarded */
    PQCLEAN_MLDSA65_CLEAN_polyveck_pack_w1(sig, &w1);

    shake256_inc_init(&state);
    shake256_inc_absorb(&state, mu, CRHBYTES);
    shake256_inc_absorb(&state, sig, K * POLYW1_PACKEDBYTES);
    shake256_inc_finalize(&state);
    shake256_inc_squeeze(sig, CTILDEBYTES, &state);
    shake256_inc_ctx_release(&state);
    PQCLEAN_MLDSA65_CLEAN_poly_challenge(&cp, sig);
    PQCLEAN_MLDSA65_CLEAN_poly_ntt(&cp);

    /*
     * ALL OPTIMIZATIONS COMBINED:
     * 1. Streaming z: compute z[i] one at a time, check norm with early rejection
     * 2. Streaming y: regenerate y[i] per-polynomial (no polyvecl y stored)
     * 3. Streaming A: regenerate A[i][j] on-the-fly (no mat[K] stored)
     * 4. Recompute w: regenerate w=A*y from saved nonce instead of storing w0
     */

    /* Pass 1: Stream compute z = y + c*s1 and check norm with early rejection */
    for (i = 0; i < L; ++i) {
        PQCLEAN_MLDSA65_CLEAN_poly_pointwise_montgomery(&tmp, &cp, &s1.vec[i]);
        PQCLEAN_MLDSA65_CLEAN_poly_invntt_tomont(&tmp);
        {
            poly yi;
            PQCLEAN_MLDSA65_CLEAN_poly_uniform_gamma1(&yi, rhoprime, (uint16_t)(L * nonce_y + i));
            PQCLEAN_MLDSA65_CLEAN_poly_add(&tmp, &tmp, &yi);
        }
        PQCLEAN_MLDSA65_CLEAN_poly_reduce(&tmp);
        if (PQCLEAN_MLDSA65_CLEAN_poly_chknorm(&tmp, GAMMA1 - BETA)) {
            goto rej;
        }
    }

    /* RECOMPUTE w = A*y to get w0 (low bits) - streaming A and y */
    for (j = 0; j < L; ++j) {
        PQCLEAN_MLDSA65_CLEAN_poly_uniform_gamma1(&tmp, rhoprime, (uint16_t)(L * nonce_y + j));
        PQCLEAN_MLDSA65_CLEAN_poly_ntt(&tmp);
        for (i = 0; i < K; ++i) {
            poly a_ij;
            PQCLEAN_MLDSA65_CLEAN_poly_uniform(&a_ij, rho, (uint16_t)((i << 8) + j));
            PQCLEAN_MLDSA65_CLEAN_poly_pointwise_montgomery(&a_ij, &a_ij, &tmp);
            if (j == 0) {
                h.vec[i] = a_ij;
            } else {
                PQCLEAN_MLDSA65_CLEAN_poly_add(&h.vec[i], &h.vec[i], &a_ij);
            }
        }
    }
    PQCLEAN_MLDSA65_CLEAN_polyveck_reduce(&h);
    PQCLEAN_MLDSA65_CLEAN_polyveck_invntt_tomont(&h);
    PQCLEAN_MLDSA65_CLEAN_polyveck_caddq(&h);
    PQCLEAN_MLDSA65_CLEAN_polyveck_decompose(&w1, &h, &h);  /* w1 = HighBits, h = LowBits = w0 */

    /* Compute c*s2 and subtract from h (which holds w0) */
    PQCLEAN_MLDSA65_CLEAN_polyveck_pointwise_poly_montgomery(&w1, &cp, &s2);
    PQCLEAN_MLDSA65_CLEAN_polyveck_invntt_tomont(&w1);
    PQCLEAN_MLDSA65_CLEAN_polyveck_sub(&h, &h, &w1);
    PQCLEAN_MLDSA65_CLEAN_polyveck_reduce(&h);
    if (PQCLEAN_MLDSA65_CLEAN_polyveck_chknorm(&h, GAMMA2 - BETA)) {
        goto rej;
    }

    /* Compute c*t0 and check its norm */
    PQCLEAN_MLDSA65_CLEAN_polyveck_pointwise_poly_montgomery(&w1, &cp, &t0);
    PQCLEAN_MLDSA65_CLEAN_polyveck_invntt_tomont(&w1);
    PQCLEAN_MLDSA65_CLEAN_polyveck_reduce(&w1);
    if (PQCLEAN_MLDSA65_CLEAN_polyveck_chknorm(&w1, GAMMA2)) {
        goto rej;
    }
    PQCLEAN_MLDSA65_CLEAN_polyveck_add(&h, &h, &w1);  /* h = w0 - c*s2 + c*t0 */

    /* Recompute w for hints - streaming A and y */
    for (j = 0; j < L; ++j) {
        PQCLEAN_MLDSA65_CLEAN_poly_uniform_gamma1(&tmp, rhoprime, (uint16_t)(L * nonce_y + j));
        PQCLEAN_MLDSA65_CLEAN_poly_ntt(&tmp);
        for (i = 0; i < K; ++i) {
            poly a_ij;
            PQCLEAN_MLDSA65_CLEAN_poly_uniform(&a_ij, rho, (uint16_t)((i << 8) + j));
            PQCLEAN_MLDSA65_CLEAN_poly_pointwise_montgomery(&a_ij, &a_ij, &tmp);
            if (j == 0) {
                w1.vec[i] = a_ij;
            } else {
                PQCLEAN_MLDSA65_CLEAN_poly_add(&w1.vec[i], &w1.vec[i], &a_ij);
            }
        }
    }
    PQCLEAN_MLDSA65_CLEAN_polyveck_reduce(&w1);
    PQCLEAN_MLDSA65_CLEAN_polyveck_invntt_tomont(&w1);
    PQCLEAN_MLDSA65_CLEAN_polyveck_caddq(&w1);

    /* Compute hints one polynomial at a time using single poly temporary */
    {
        poly w_high_i;
        n = 0;
        for (i = 0; i < K; ++i) {
            PQCLEAN_MLDSA65_CLEAN_poly_decompose(&w_high_i, &w1.vec[i], &w1.vec[i]);
            n += PQCLEAN_MLDSA65_CLEAN_poly_make_hint(&w1.vec[i], &h.vec[i], &w_high_i);
        }
        if (n > OMEGA) {
            goto rej;
        }
        h = w1;
    }

    /* Pass 2: Recompute z and pack directly into signature buffer */
    {
        uint8_t *sig_z = sig + CTILDEBYTES;

        for (i = 0; i < L; ++i) {
            poly yi;
            PQCLEAN_MLDSA65_CLEAN_poly_pointwise_montgomery(&tmp, &cp, &s1.vec[i]);
            PQCLEAN_MLDSA65_CLEAN_poly_invntt_tomont(&tmp);
            PQCLEAN_MLDSA65_CLEAN_poly_uniform_gamma1(&yi, rhoprime, (uint16_t)(L * nonce_y + i));
            PQCLEAN_MLDSA65_CLEAN_poly_add(&tmp, &tmp, &yi);
            PQCLEAN_MLDSA65_CLEAN_poly_reduce(&tmp);
            PQCLEAN_MLDSA65_CLEAN_polyz_pack(sig_z + i * POLYZ_PACKEDBYTES, &tmp);
        }

        /* Pack h (hint vector) */
        {
            uint8_t *sig_h = sig + CTILDEBYTES + L * POLYZ_PACKEDBYTES;
            unsigned int k_idx;

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
