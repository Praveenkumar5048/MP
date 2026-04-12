#include "indcpa.h"
#include "ntt.h"
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "randombytes.h"
#include "symmetric.h"
#include <stddef.h>
#include <stdint.h>
#include <string.h>

/*************************************************
* Name:        pack_pk
*
* Description: Serialize the public key as concatenation of the
*              serialized vector of polynomials pk
*              and the public seed used to generate the matrix A.
*
* Arguments:   uint8_t *r:           pointer to the output serialized public key
*              polyvec *pk:          pointer to the input public-key polyvec
*              const uint8_t *seed:  pointer to the input public seed
**************************************************/
static void pack_pk(uint8_t r[KYBER_INDCPA_PUBLICKEYBYTES],
                    polyvec *pk,
                    const uint8_t seed[KYBER_SYMBYTES]) {
    PQCLEAN_MLKEM768_CLEAN_polyvec_tobytes(r, pk);
    memcpy(r + KYBER_POLYVECBYTES, seed, KYBER_SYMBYTES);
}

/*************************************************
* Name:        unpack_pk
*
* Description: De-serialize public key from a byte array;
*              approximate inverse of pack_pk
*
* Arguments:   - polyvec *pk:              pointer to output public-key polynomial vector
*              - uint8_t *seed:            pointer to output seed to generate matrix A
*              - const uint8_t *packedpk:  pointer to input serialized public key
**************************************************/
static void unpack_pk(polyvec *pk,
                      uint8_t seed[KYBER_SYMBYTES],
                      const uint8_t packedpk[KYBER_INDCPA_PUBLICKEYBYTES]) {
    PQCLEAN_MLKEM768_CLEAN_polyvec_frombytes(pk, packedpk);
    memcpy(seed, packedpk + KYBER_POLYVECBYTES, KYBER_SYMBYTES);
}

/*************************************************
* Name:        pack_sk
*
* Description: Serialize the secret key
*
* Arguments:   - uint8_t *r:    pointer to output serialized secret key
*              - polyvec *sk:   pointer to input vector of polynomials (secret key)
**************************************************/
static void pack_sk(uint8_t r[KYBER_INDCPA_SECRETKEYBYTES], polyvec *sk) {
    PQCLEAN_MLKEM768_CLEAN_polyvec_tobytes(r, sk);
}

/*************************************************
* Name:        unpack_sk
*
* Description: De-serialize the secret key; inverse of pack_sk
*
* Arguments:   - polyvec *sk:          pointer to output vector of polynomials (secret key)
*              - const uint8_t *packedsk: pointer to input serialized secret key
**************************************************/
static void unpack_sk(polyvec *sk, const uint8_t packedsk[KYBER_INDCPA_SECRETKEYBYTES]) {
    PQCLEAN_MLKEM768_CLEAN_polyvec_frombytes(sk, packedsk);
}

/*************************************************
* Name:        unpack_ciphertext
*
* Description: De-serialize and decompress ciphertext from a byte array;
*              approximate inverse of pack_ciphertext
*
* Arguments:   - polyvec *b:      pointer to the output vector of polynomials b
*              - poly *v:         pointer to the output polynomial v
*              - const uint8_t *c: pointer to the input serialized ciphertext
**************************************************/
static void unpack_ciphertext(polyvec *b, poly *v,
                               const uint8_t c[KYBER_INDCPA_BYTES]) {
    PQCLEAN_MLKEM768_CLEAN_polyvec_decompress(b, c);
    PQCLEAN_MLKEM768_CLEAN_poly_decompress(v, c + KYBER_POLYVECCOMPRESSEDBYTES);
}

/*************************************************
* Name:        poly_compress_u_slot
*
* Description: Compress and serialize one polynomial into its slot of the
*              u-component of the ciphertext.  The pointer r must already
*              be advanced to the correct slot offset by the caller.
*
*              This is a standalone extraction of the 10-bit branch from
*              polyvec_compress (KYBER_POLYVECCOMPRESSEDBYTES/KYBER_K == 320
*              for k=2,3,4), so it handles exactly the 10-bits-per-coefficient
*              case used by all three parameter sets.
*
* Arguments:   - uint8_t *r:    pointer to output byte stream (slot start)
*              - const poly *a: pointer to input polynomial
**************************************************/
static void poly_compress_u_slot(uint8_t *r, const poly *a) {
    unsigned int j, k;
    uint64_t d0;
    uint16_t t[4];

    for (j = 0; j < KYBER_N / 4; j++) {
        for (k = 0; k < 4; k++) {
            t[k] = a->coeffs[4 * j + k];
            t[k] += ((int16_t)t[k] >> 15) & KYBER_Q;
            d0    = t[k];
            d0  <<= 10;
            d0   += 1665;
            d0   *= 1290167;
            d0  >>= 32;
            t[k]  = d0 & 0x3ff;
        }

        r[0] = (uint8_t)(t[0] >> 0);
        r[1] = (uint8_t)((t[0] >> 8) | (t[1] << 2));
        r[2] = (uint8_t)((t[1] >> 6) | (t[2] << 4));
        r[3] = (uint8_t)((t[2] >> 4) | (t[3] << 6));
        r[4] = (uint8_t)(t[3] >> 2);
        r   += 5;
    }
}

/*************************************************
* Name:        rej_uniform
*
* Description: Run rejection sampling on uniform random bytes to generate
*              uniform random integers mod q
*
* Arguments:   - int16_t *r:        pointer to output buffer
*              - unsigned int len:   requested number of 16-bit integers (uniform mod q)
*              - const uint8_t *buf: pointer to input buffer (assumed uniform random)
*              - unsigned int buflen: length of input buffer in bytes
*
* Returns number of sampled 16-bit integers (at most len)
**************************************************/
static unsigned int rej_uniform(int16_t *r,
                                unsigned int len,
                                const uint8_t *buf,
                                unsigned int buflen) {
    unsigned int ctr, pos;
    uint16_t val0, val1;

    ctr = pos = 0;
    while (ctr < len && pos + 3 <= buflen) {
        val0 = ((buf[pos + 0] >> 0) | ((uint16_t)buf[pos + 1] << 8)) & 0xFFF;
        val1 = ((buf[pos + 1] >> 4) | ((uint16_t)buf[pos + 2] << 4)) & 0xFFF;
        pos += 3;

        if (val0 < KYBER_Q) {
            r[ctr++] = val0;
        }
        if (ctr < len && val1 < KYBER_Q) {
            r[ctr++] = val1;
        }
    }

    return ctr;
}

#define gen_a(A,B)  PQCLEAN_MLKEM768_CLEAN_gen_matrix(A,B,0)
#define gen_at(A,B) PQCLEAN_MLKEM768_CLEAN_gen_matrix(A,B,1)

#define GEN_MATRIX_NBLOCKS \
    ((12*KYBER_N/8*(1 << 12)/KYBER_Q + XOF_BLOCKBYTES)/XOF_BLOCKBYTES)

/*************************************************
* Name:        gen_matrix_entry
*
* Description: Generate a single polynomial A[i][j] (or A^T[i][j]) on the
*              fly from the public seed.  Avoids allocating the full k×k
*              matrix on the stack — the key enabler for row-at-a-time enc.
*
* Arguments:   - poly *aij:             pointer to output polynomial
*              - const uint8_t *seed:   public seed
*              - int transposed:        0 → A[i][j],  1 → A^T[i][j]
*              - unsigned int i, j:     row / column indices
**************************************************/
static void gen_matrix_entry(poly *aij,
                             const uint8_t seed[KYBER_SYMBYTES],
                             int transposed,
                             unsigned int i,
                             unsigned int j) {
    unsigned int ctr, buflen;
    uint8_t buf[GEN_MATRIX_NBLOCKS * XOF_BLOCKBYTES];
    xof_state state;

    /* For A^T: A^T[i][j] = A[j][i], so swap the domain-separation bytes. */
    if (transposed) {
        xof_absorb(&state, seed, (uint8_t)i, (uint8_t)j);
    } else {
        xof_absorb(&state, seed, (uint8_t)j, (uint8_t)i);
    }

    xof_squeezeblocks(buf, GEN_MATRIX_NBLOCKS, &state);
    buflen = GEN_MATRIX_NBLOCKS * XOF_BLOCKBYTES;
    ctr = rej_uniform(aij->coeffs, KYBER_N, buf, buflen);

    while (ctr < KYBER_N) {
        xof_squeezeblocks(buf, 1, &state);
        buflen = XOF_BLOCKBYTES;
        ctr += rej_uniform(aij->coeffs + ctr, KYBER_N - ctr, buf, buflen);
    }

    xof_ctx_release(&state);
}

/*************************************************
* Name:        PQCLEAN_MLKEM768_CLEAN_gen_matrix
*
* Description: Deterministically generate matrix A (or its transpose)
*              from a seed. Entries look uniformly random. Uses rejection
*              sampling on the output of a XOF.
*
* Arguments:   - polyvec *a:           pointer to output matrix A
*              - const uint8_t *seed:  pointer to input seed
*              - int transposed:       boolean — generate A^T instead of A
*
* Note: kept for keypair generation which still needs the full matrix.
*       Encryption uses gen_matrix_entry instead.
**************************************************/
void PQCLEAN_MLKEM768_CLEAN_gen_matrix(polyvec *a,
                                       const uint8_t seed[KYBER_SYMBYTES],
                                       int transposed) {
    unsigned int ctr, i, j, buflen;
    uint8_t buf[GEN_MATRIX_NBLOCKS * XOF_BLOCKBYTES];
    xof_state state;

    for (i = 0; i < KYBER_K; i++) {
        for (j = 0; j < KYBER_K; j++) {
            if (transposed) {
                xof_absorb(&state, seed, (uint8_t)i, (uint8_t)j);
            } else {
                xof_absorb(&state, seed, (uint8_t)j, (uint8_t)i);
            }

            xof_squeezeblocks(buf, GEN_MATRIX_NBLOCKS, &state);
            buflen = GEN_MATRIX_NBLOCKS * XOF_BLOCKBYTES;
            ctr = rej_uniform(a[i].vec[j].coeffs, KYBER_N, buf, buflen);

            while (ctr < KYBER_N) {
                xof_squeezeblocks(buf, 1, &state);
                buflen = XOF_BLOCKBYTES;
                ctr += rej_uniform(a[i].vec[j].coeffs + ctr,
                                   KYBER_N - ctr, buf, buflen);
            }
            xof_ctx_release(&state);
        }
    }
}

/*************************************************
* Name:        PQCLEAN_MLKEM768_CLEAN_indcpa_keypair_derand
*
* Description: Generates public and private key for the CPA-secure
*              public-key encryption scheme underlying Kyber
*
* Arguments:   - uint8_t *pk:          pointer to output public key
*              - uint8_t *sk:          pointer to output private key
*              - const uint8_t *coins: pointer to input randomness
*                                      (of length KYBER_SYMBYTES bytes)
**************************************************/
void PQCLEAN_MLKEM768_CLEAN_indcpa_keypair_derand(
        uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
        uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES],
        const uint8_t coins[KYBER_SYMBYTES]) {
    unsigned int i;
    uint8_t buf[2 * KYBER_SYMBYTES];
    const uint8_t *publicseed = buf;
    const uint8_t *noiseseed  = buf + KYBER_SYMBYTES;
    uint8_t nonce = 0;
    polyvec a[KYBER_K], e, pkpv, skpv;

    memcpy(buf, coins, KYBER_SYMBYTES);
    buf[KYBER_SYMBYTES] = KYBER_K;
    hash_g(buf, buf, KYBER_SYMBYTES + 1);

    gen_a(a, publicseed);

    for (i = 0; i < KYBER_K; i++) {
        PQCLEAN_MLKEM768_CLEAN_poly_getnoise_eta1(&skpv.vec[i], noiseseed, nonce++);
    }
    for (i = 0; i < KYBER_K; i++) {
        PQCLEAN_MLKEM768_CLEAN_poly_getnoise_eta1(&e.vec[i], noiseseed, nonce++);
    }

    PQCLEAN_MLKEM768_CLEAN_polyvec_ntt(&skpv);
    PQCLEAN_MLKEM768_CLEAN_polyvec_ntt(&e);

    for (i = 0; i < KYBER_K; i++) {
        PQCLEAN_MLKEM768_CLEAN_polyvec_basemul_acc_montgomery(&pkpv.vec[i],
                                                              &a[i], &skpv);
        PQCLEAN_MLKEM768_CLEAN_poly_tomont(&pkpv.vec[i]);
    }

    PQCLEAN_MLKEM768_CLEAN_polyvec_add(&pkpv, &pkpv, &e);
    PQCLEAN_MLKEM768_CLEAN_polyvec_reduce(&pkpv);

    pack_sk(sk, &skpv);
    pack_pk(pk, &pkpv, publicseed);
}

/*************************************************
* Name:        PQCLEAN_MLKEM768_CLEAN_indcpa_enc
*
* Description: Encryption function of the CPA-secure public-key encryption
*              scheme underlying Kyber.
*
* Stack optimisation (this version)
* ----------------------------------
* Original code held four full polyvecs simultaneously:
*   sp (r̂), pkpv (t̂), ep (e₁), b (accumulator)  →  4 × k × 512 B
*
* This version eliminates ep and b entirely by processing u one output
* polynomial at a time while preserving the original arithmetic order:
*
*   u[i] = iNTT( A^T[i] ∘ r̂ ) + e₁[i]
*
* For each row i, we compute iNTT(A^T[i]∘r̂) into one scratch polynomial,
* add the freshly sampled e₁[i] in normal domain, reduce, and compress
* directly into the output ciphertext. Scratch is reused for row i+1.
*
* Locals inside the per-row block:
*   u_i  — row accumulator in normal domain after iNTT     (512 B, 1 poly)
*   e1_i — one e₁ polynomial                               (512 B, 1 poly)
*   aij  — one matrix entry generated on the fly           (512 B, 1 poly)
*   t    — one basemul result before accumulation          (512 B, 1 poly)
*
* Nonce assignment (unchanged from spec):
*   0 … k-1        : r̂  components (eta1)
*   k … 2k-1       : e₁ components (eta2)   ← sampled here, per row
*   2k              : e₂ scalar    (eta2)
*
* Stack savings vs. original (both ep and b gone):
*   k=2  −2 048 B  |  k=3  −3 072 B  |  k=4  −4 096 B
*
* Arguments:   - uint8_t *c:           pointer to output ciphertext
*              - const uint8_t *m:     pointer to input message
*              - const uint8_t *pk:    pointer to input public key
*              - const uint8_t *coins: pointer to input random coins
**************************************************/
void PQCLEAN_MLKEM768_CLEAN_indcpa_enc(uint8_t c[KYBER_INDCPA_BYTES],
                                       const uint8_t m[KYBER_INDCPA_MSGBYTES],
                                       const uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                                       const uint8_t coins[KYBER_SYMBYTES]) {
    unsigned int i, j;
    uint8_t seed[KYBER_SYMBYTES];
    uint8_t nonce = 0;
    polyvec sp, pkpv;           /* r̂ and t̂ — still needed in full */
    poly v, k, epp;             /* scalar v, message poly, e₂ */
    poly u_i, e1_i, aij, t;     /* per-row scratch */
    const size_t u_slot_bytes = KYBER_POLYVECCOMPRESSEDBYTES / KYBER_K;

    /* ------------------------------------------------------------------ */
    /* Unpack public key and encode message                                */
    /* ------------------------------------------------------------------ */
    unpack_pk(&pkpv, seed, pk);
    PQCLEAN_MLKEM768_CLEAN_poly_frommsg(&k, m);

    /* ------------------------------------------------------------------ */
    /* Sample r, compute r̂ = NTT(r)  — nonces 0 … k-1                   */
    /* ------------------------------------------------------------------ */
    for (i = 0; i < KYBER_K; i++) {
        PQCLEAN_MLKEM768_CLEAN_poly_getnoise_eta1(sp.vec + i, coins, nonce++);
    }
    PQCLEAN_MLKEM768_CLEAN_polyvec_ntt(&sp);

    /* ------------------------------------------------------------------ */
    /* Compute u one row at a time  — nonces k … 2k-1                     */
    /*                                                                     */
    /* For row i:                                                          */
    /*   u_i ← NTT(e₁[i])  (Montgomery-scaled)                           */
    /*       + Σ_j  A^T[i][j] ∘ r̂[j]  (also Montgomery-scaled)          */
    /*   then apply one iNTT_tomont and compress directly into c.         */
    /* ------------------------------------------------------------------ */
    for (i = 0; i < KYBER_K; i++) {

        /* (a) Sample e₁[i] in normal domain. */
        PQCLEAN_MLKEM768_CLEAN_poly_getnoise_eta2(&e1_i, coins, nonce++);

        /*
         * (b) Accumulate A^T[i][j] ∘ r̂[j] into u_i for all j.
         *     gen_matrix_entry generates exactly one polynomial on the fly —
         *     no full matrix is ever allocated.
         */
        gen_matrix_entry(&aij, seed, /*transposed=*/1, i, 0);
        PQCLEAN_MLKEM768_CLEAN_poly_basemul_montgomery(&u_i, &aij, &sp.vec[0]);
        for (j = 1; j < KYBER_K; j++) {
            gen_matrix_entry(&aij, seed, /*transposed=*/1, i, j);
            PQCLEAN_MLKEM768_CLEAN_poly_basemul_montgomery(&t, &aij, &sp.vec[j]);
            PQCLEAN_MLKEM768_CLEAN_poly_add(&u_i, &u_i, &t);
        }
        PQCLEAN_MLKEM768_CLEAN_poly_reduce(&u_i);

        /*
         * (c) Bring row product back to normal domain and add e₁[i],
         *     preserving baseline semantics exactly.
         */
        PQCLEAN_MLKEM768_CLEAN_poly_invntt_tomont(&u_i);
        PQCLEAN_MLKEM768_CLEAN_poly_add(&u_i, &u_i, &e1_i);
        PQCLEAN_MLKEM768_CLEAN_poly_reduce(&u_i);

        /*
         * (d) Write compressed u[i] directly into the output ciphertext.
         *     No full polyvec u is ever allocated.
         */
        poly_compress_u_slot(c + i * u_slot_bytes, &u_i);
    }

    /* ------------------------------------------------------------------ */
    /* Compute v = iNTT(t̂ᵀ ∘ r̂) + e₂ + msg  — nonce 2k                 */
    /* ------------------------------------------------------------------ */
    PQCLEAN_MLKEM768_CLEAN_poly_getnoise_eta2(&epp, coins, nonce++);

    PQCLEAN_MLKEM768_CLEAN_polyvec_basemul_acc_montgomery(&v, &pkpv, &sp);
    PQCLEAN_MLKEM768_CLEAN_poly_invntt_tomont(&v);

    PQCLEAN_MLKEM768_CLEAN_poly_add(&v, &v, &epp);
    PQCLEAN_MLKEM768_CLEAN_poly_add(&v, &v, &k);
    PQCLEAN_MLKEM768_CLEAN_poly_reduce(&v);

    PQCLEAN_MLKEM768_CLEAN_poly_compress(c + KYBER_POLYVECCOMPRESSEDBYTES, &v);
}

/*************************************************
* Name:        PQCLEAN_MLKEM768_CLEAN_indcpa_dec
*
* Description: Decryption function of the CPA-secure public-key encryption
*              scheme underlying Kyber.
*
* Arguments:   - uint8_t *m:       pointer to output decrypted message
*              - const uint8_t *c: pointer to input ciphertext
*              - const uint8_t *sk: pointer to input secret key
**************************************************/
void PQCLEAN_MLKEM768_CLEAN_indcpa_dec(uint8_t m[KYBER_INDCPA_MSGBYTES],
                                       const uint8_t c[KYBER_INDCPA_BYTES],
                                       const uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES]) {
    polyvec b, skpv;
    poly v, mp;

    unpack_ciphertext(&b, &v, c);
    unpack_sk(&skpv, sk);

    PQCLEAN_MLKEM768_CLEAN_polyvec_ntt(&b);
    PQCLEAN_MLKEM768_CLEAN_polyvec_basemul_acc_montgomery(&mp, &skpv, &b);
    PQCLEAN_MLKEM768_CLEAN_poly_invntt_tomont(&mp);

    PQCLEAN_MLKEM768_CLEAN_poly_sub(&mp, &v, &mp);
    PQCLEAN_MLKEM768_CLEAN_poly_reduce(&mp);

    PQCLEAN_MLKEM768_CLEAN_poly_tomsg(m, &mp);
}