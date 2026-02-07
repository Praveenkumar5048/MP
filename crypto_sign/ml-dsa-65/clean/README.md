# ML-DSA-65 Clean Implementation - Small NTT Optimization

## Overview
This implementation uses a smaller prime for challenge-secret multiplications to reduce stack usage.

## Optimization: Small NTT (q=7681)

**Problem:** Standard NTT uses q=8380417 (23-bit) requiring 32-bit coefficients.

**Solution:** For c*s1 and c*s2 computations, use q=7681 (13-bit) with 16-bit coefficients.

**Why it works:**
- Challenge c has τ=49 non-zero ±1 coefficients
- Secret s1/s2 coefficients are in {-4, ..., 4}
- Product bound: |c*s|∞ ≤ τ × 4 = 196
- This fits in q=7681 (13-bit prime with 256th root of unity)

**Benefit:** `smallpoly` uses 512 bytes vs `poly`'s 1024 bytes

## Stack Usage

| Function | Bytes | vs Baseline |
|----------|-------|-------------|
| KeyGen | 59,776 | 0 |
| **Sign** | **78,864** | **-440 (-0.55%)** |
| Verify | 56,680 | 0 |

Baseline Sign was 79,304 bytes.

## Files
- `sign.c` - Uses small NTT for c*s1 and c*s2
- `smallntt.c/h` - Small modulus NTT implementation

## Verification
- All functional tests pass
- Testvectors match reference implementation
