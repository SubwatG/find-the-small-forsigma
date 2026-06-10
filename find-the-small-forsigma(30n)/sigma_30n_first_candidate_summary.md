# Smallest Known Candidate: σ(30n+1) > σ(30n)

## Status

**Current best candidate found** — minimality NOT yet rigorously proven.
Believed to be the smallest possible.

## Result — Candidate

```
n = 13032232299953752123721014703967652574762133589869955544193341473
z = 30n + 1 = 390966968998612563711630441119029577242864007696098666325800244191
```

- n digits: **65**
- z digits: **66**

## Factorization

**z = 30n + 1 (66 digits):**
```
7³ × 11² × 13² × 17² × 19² × 23² × 29 × 31 × 37 × 41 × 43 × 47 
× 53 × 59 × 61 × 67 × 71 × 73 × 79 × 83 × 89 × 97 × 101 × 103 
× 107 × 109 × 113 × 127 × 131 × 137 × 139 × 149 × 151
```
- Distinct primes: **33**
- Prime powers: 7³, 11², 13², 17², 19², 23² (rest exponent 1)
- Smallest allowed prime: 7
- Largest prime used: 151

**z - 1 = 30n (66 digits):**
```
2 × 3 × 5 × 13654233307799957 × 954446288281093901941669500832420518262112112989
```

**n = (z-1)/30 (65 digits):**
```
13654233307799957 × 954446288281093901941669500832420518262112112989
```

n has only TWO large prime factors (the best case for us).

## Sigma Values

| Quantity | Value |
|----------|-------|
| σ(z) | 939,891,231,856,361,256,945,954,360,893,619,510,395,124,991,191,794,319,360,000,000,000 |
| σ(z-1) | 938,320,725,596,670,221,628,045,814,924,431,925,183,077,678,405,897,218,852,154,318,240 |
| **σ(z)/σ(z-1)** | **1.0016737414** |
| A(z) = σ(z)/z | 2.4040169794 |
| A(z-1) = σ(z-1)/(z-1) | 2.4000000000 |

## Why This Candidate Works

1. **A(z) = 2.404** — built from 33 primes from 7 to 151 with powers on the smallest ones
2. **A(z-1) = 2.400 exactly** — n has only 2 prime factors, so A(30n) = σ(30)/30 × A(n) = 2.4 × (1+1/p)(1+1/q)/(1×1) ≈ 2.4
3. The gap of 0.004 in abundancy is just enough to make σ(z) > σ(z-1)

## Search History

| ZMAX | n digits | A(z) | A(z-1) | Ratio | Status |
|------|----------|------|--------|-------|--------|
| 10⁸⁰ | 79 | 2.490 | 2.401 | 1.037 | ✅ Larger |
| 10⁷⁰ | 69 | 2.421 | 2.400 | 1.008 | ✅ Larger |
| **10⁶⁶** | **65** | **2.404** | **2.400** | **1.002** | ✅ **Smallest found** |
| 10⁶⁵ | — | 2.391 | 2.400 | 0.996 | ❌ No success |

## Search Parameters (for 10^66 candidate)

- Method: Beam search diagnostic
- ZMAX: 10^66
- pmax: 300 (primes up to 293)
- K: 2000 (beam width per residue)
- λ (lam): 0.0002 (size penalty)
- η (eta): 0 (no z-1 penalty)
- Date: 2026-06-07

## Key Observation for Minimality

At ZMAX = 10^65, the beam search found candidates with A(z) up to 2.391
but ALL had σ(z)/σ(z-1) < 1 (best ratio 0.996).

This suggests that any z < 10^66 with z ≡ 1 (mod 30) cannot have A(z) ≥ 2.40
(because there isn't enough "budget" to fit the required 32+ primes).

To rigorously prove minimality, we would adapt Greg Martin's Claim:

**Claim (adapted for σ):** For m not divisible by p₁,…,p_r with at most s distinct
prime factors, the abundancy is bounded above by the product of σ(p)/p
for the s smallest allowed primes.

With r=3 (blocks 2,3,5) and appropriate s, we can show that any m < z 
with m ≡ 1 (mod 30) must have A(m) < A(z) ≈ 2.404, hence cannot be
a smaller counterexample.

## Comparison with Other Moduli

| Modulus | Smallest n | n digits | z structure |
|---------|-----------|----------|-------------|
| 6n | 4,802,693,729 | 10 | 5³×7×11×13×17×19×23×31 (8 primes) |
| 30n | ~1.3×10⁶⁴ | **65** | 33 primes from 7 to 151 |
