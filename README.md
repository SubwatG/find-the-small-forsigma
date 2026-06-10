# find-the-small-forsigma

Beam search for the smallest positive integer \(n\) such that

\[
\sigma(an+1) > \sigma(an)
\]

for squarefree moduli \(a = 30, 42, 70\).

This repository accompanies the paper  
**"The Smallest Solution of \(\sigma(an+1) > \sigma(an)\) for \(a = 30, 42,\) and \(70\)".**

## Results

| \(a\) | Blocked primes | Gap | Smallest \(n\) | Digits | \(\frac{\sigma(an+1)}{\sigma(an)}\) |
|-------|---------------|-----|----------------|--------|-------------------------------------|
| 30    | {2,3,5}       | —   | 1.3032×10⁶⁴    | 65     | 1.00167                             |
| 42    | {2,3,7}       | 5   | 6.6874×10⁴¹    | 42     | 1.01111                             |
| 70    | {2,5,7}       | 3   | 48,049,097     | 8      | 1.00767                             |

Having a *gap prime* dramatically reduces the size of the smallest \(n\).

## Repository Structure

```
find-the-small-forsigma/
├── README.md
├── LICENSE
├── .gitignore
├── paper/
│   └── paper.tex                      # LaTeX source of the paper
├── src/
│   └── beam_search.py                 # Beam search implementation
├── results/
│   ├── 30/
│   │   ├── candidate_summary.md
│   │   └── minimality_proof.md
│   ├── 42/
│   │   └── candidate_summary.md
│   └── 70/
│       └── candidate_summary.md
```

## Reproducing the Results

The main script is `src/beam_search.py`.

### Prerequisites

```bash
pip install sympy
```

### Beam Search (Example for modulus 30)

```python
from src.beam_search import beam_search, test_candidate

candidates = beam_search(
    ZMAX=10**80,
    K=2000,
    pmax=300,
    lam=0.0002,
    TARGET_M=30,
    verbose=True
)

for z, fac in candidates:
    ok, ratio = test_candidate(z)
    if ok:
        print(f"Found: n = {(z-1)//30}, ratio = {ratio}")
```

### Parameters used in the paper

| \(a\) | \(Z_{\max}\) | \(K\) | \(p_{\max}\) | \(\lambda\) |
|-------|-------------|-------|-------------|------------|
| 30    | \(10^{80}\) | 2000  | 300         | 0.0002     |
| 42    | \(10^{44}\) | 2000  | 300         | 0.0001     |
| 70    | \(10^{11}\) | 2000  | 200         | 0.0002     |

## Citation

If you use this code or the results, please cite:

```
K. Subwattanachai,
"The Smallest Solution of σ(an+1) > σ(an) for a = 30, 42, and 70",
preprint, 2026.
```

## License

MIT License — see [LICENSE](LICENSE).

## Contact

Kittipong Subwattanachai  
Department of Mathematics, Kasetsart University  
subwattanachai.k@gmail.com
