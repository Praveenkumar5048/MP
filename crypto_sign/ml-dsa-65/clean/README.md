# ML-DSA-65 Clean Implementation - Streaming Z Optimization

## Overview
This implementation includes a stack optimization that reduces signing memory usage by ~4KB.

## Optimization: Streaming Z Computation

**Problem:** The standard implementation stores the full `polyvecl z` vector (5,120 bytes for L=5).

**Solution:** Two-pass streaming approach that eliminates z storage:
- **Pass 1:** Compute each z[i] = y[i] + c*s1[i], check norm immediately (early rejection)
- **Pass 2:** Recompute z[i] and pack directly into signature buffer

This trades ~2× computation of c*s1 for significant stack savings.

## Stack Usage

| Function | Bytes | vs Baseline |
|----------|-------|-------------|
| KeyGen | 59,696 | -80 |
| **Sign** | **75,280** | **-4,024 (-5.1%)** |
| Verify | 56,400 | -280 |

Baseline Sign was 79,304 bytes.

## Files Modified
- `sign.c` - Streaming z computation in signature generation

## Verification
- All functional tests pass
- Testvectors match reference implementation
