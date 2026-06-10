# find-the-small-forsigma

Beam search for the smallest positive integer \(n\) such that
\[
\sigma(an+1) > \sigma(an)
\]
for squarefree moduli \(a = 30, 42, 70\).

This repository contains the code and results accompanying the paper
**"The Smallest Solution of \(\sigma(an+1) > \sigma(an)\) for \(a = 30, 42,\) and \(70\)"**.

## Main Results

| \(a\) | Blocked primes | Gap prime | Smallest \(n\) | Digits | \(\sigma(an+1)/\sigma(an)\) |
|-------|----------------|-----------|----------------|--------|-----------------------------|
| 30    | \{2,3,5\}      | —         | 1.303×10⁶⁴     | 65     | 1.00167                     |
| 42    | \{2,3,7\}      | 5         | 6.687×10⁴¹     | 42     | 1.01111                     |
| 70    | \{2,5,7\}      | 3         | 48 049 097     | 8      | 1.00767                     |

The dramatic difference in size is explained by the presence or absence of
**gap primes** — small primes that do not divide \(a\) and can therefore be
used to construct a highly abundant \(an+1\).

## Repository Structure

```
find-the-small-forsigma/
├── find-the-small-forsigma(30n)/
│   ├── hermes_sigma_30n_experiment.py   # Beam search implementation
│   ├── sigma_30n_first_candidate_summary.md
│   └── sigma_30n_minimality_proof.md    # Full mathematical proof
├── sigma_42n_first_candidate_summary.md
├── sigma_70n_first_candidate_summary.md
├── paper.tex
├── paper.pdf
└── README.md
```

## Reproducing the Results

The main script is
`find-the-small-forsigma(30n)/hermes_sigma_30n_experiment.py`.

Typical usage for modulus 30:
```bash
python3 hermes_sigma_30n_experiment.py beam-search \
  --a 30 --ZMAX 1e80 --K 2000 --pmax 300 --lam 0.0002
```

For other moduli, change the `--a` parameter and adjust `ZMAX` accordingly
(see Table 1 in the paper).

## License

MIT License

## Citation

If you use this code or the results, please cite the accompanying paper:

```
K. Subwattanachai,
"The Smallest Solution of σ(an+1) > σ(an) for a = 30, 42, and 70",
preprint, 2026.
```

## Contact

Kittipong Subwattanachai  
Department of Mathematics, Kasetsart University  
subwattanachai.k@gmail.com
