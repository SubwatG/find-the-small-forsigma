# Sigma(70n+1) > Sigma(70n) — Smallest n

## Result

```
n = 48,049,097  (8 digits)
z = 70n + 1 = 3,363,436,791  (10 digits)
```

## Factorization

**z = 70n+1:**
```
3⁴ × 11 × 13 × 17 × 19 × 29 × 31
```
- Distinct primes: 7
- Gap prime used: **3** (blocked primes are 2,5,7)

**70n:**
```
2 × 5 × 7 × 48,049,097
```
where 48,049,097 is prime!

## Sigma Values

| Quantity | Value |
|----------|-------|
| σ(z) | 7,034,221,440 |
| σ(70n) | 6,980,835,360 |
| **σ(z)/σ(70n)** | **1.0076665540** |
| A(z) | 2.0887439614 |
| A(70n) | 2.0728515967 |

## Comparison

| Modulus | Blocked | Gap | n digits | Ratio |
|---------|---------|-----|----------|-------|
| 30 | {2,3,5} | — | 65 | 1.0017 |
| 42 | {2,3,7} | 5 | 42 | 1.0111 |
| **70** | **{2,5,7}** | **3** | **8** | **1.0077** |

Having 3 available makes this case almost trivial — only 7 primes needed,
and the smallest n is just 8 digits.
