#!/usr/bin/env python3
"""Batch exact searches for sigma(Mn) versus sigma(Mn + 1).

Here M ranges over squarefree products of three distinct primes:

    M = p*q*r, p < q < r, M <= M_max.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


DEFAULT_M_MAX = 10_000
DEFAULT_N = 100_000
DEFAULT_SIEVE_LIMIT = 50_000_000
M30_EXTERNAL_NOTE = (
    "Existing 30n segmented work reports no original-direction hit for "
    "1 <= n <= 10,000,000,000."
)
M30_EXTERNAL_SOURCE = "find-the-small-forsigma(30n)/Experimental Plan for (30n).md"


def parse_int(value: str) -> int:
    """Parse CLI integers, accepting forms like 10^6 and 1_000_000."""
    cleaned = value.strip().replace("_", "")
    if "^" in cleaned:
        base, exponent = cleaned.split("^", 1)
        return int(base) ** int(exponent)
    return int(cleaned)


def prime_sieve(limit: int) -> list[int]:
    if limit < 2:
        return []
    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0:2] = b"\x00\x00"
    for p in range(2, int(limit**0.5) + 1):
        if is_prime[p]:
            start = p * p
            is_prime[start : limit + 1 : p] = b"\x00" * (((limit - start) // p) + 1)
    return [p for p in range(2, limit + 1) if is_prime[p]]


def parse_int_list(value: str) -> list[int]:
    items = []
    for part in value.split(","):
        part = part.strip()
        if part:
            items.append(parse_int(part))
    return items


def enumerate_squarefree_triples(mmax: int) -> list[dict[str, Any]]:
    """Enumerate p*q*r <= mmax using early breaks, without all combinations."""
    primes = prime_sieve(mmax)
    rows: list[dict[str, Any]] = []
    count = len(primes)
    for i, p in enumerate(primes):
        if i + 2 >= count or p * primes[i + 1] * primes[i + 2] > mmax:
            break
        for j in range(i + 1, count):
            q = primes[j]
            if j + 1 >= count or p * q * primes[j + 1] > mmax:
                break
            for k in range(j + 1, count):
                r = primes[k]
                M = p * q * r
                if M > mmax:
                    break
                rows.append({"M": M, "factors": [p, q, r]})
    rows.sort(key=lambda row: (row["M"], row["factors"]))
    return rows


def smallest_prime_factors(limit: int) -> list[int]:
    spf = list(range(limit + 1))
    if limit >= 0:
        spf[0] = 0
    if limit >= 1:
        spf[1] = 1
    for p in range(2, int(limit**0.5) + 1):
        if spf[p] == p:
            for multiple in range(p * p, limit + 1, p):
                if spf[multiple] == multiple:
                    spf[multiple] = p
    return spf


def factor_from_spf(n: int, spf: list[int]) -> dict[int, int]:
    fac: dict[int, int] = {}
    while n > 1:
        p = spf[n]
        fac[p] = fac.get(p, 0) + 1
        n //= p
    return fac


def is_prime(n: int) -> bool:
    """Deterministic Miller-Rabin for the 32-bit sized values used here."""
    if n < 2:
        return False
    small_primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for p in small_primes:
        if n == p:
            return True
        if n % p == 0:
            return False

    d = n - 1
    s = 0
    while d % 2 == 0:
        s += 1
        d //= 2

    for a in (2, 3, 5, 7):
        if a >= n:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                break
        else:
            return False
    return True


def factor_by_trial_division(n: int, primes: list[int]) -> dict[int, int]:
    if is_prime(n):
        return {n: 1}

    fac: dict[int, int] = {}
    residual = n
    for p in primes:
        if p * p > residual:
            break
        while residual % p == 0:
            fac[p] = fac.get(p, 0) + 1
            residual //= p
        if residual > 1 and is_prime(residual):
            fac[residual] = fac.get(residual, 0) + 1
            residual = 1
            break
    if residual > 1:
        fac[residual] = fac.get(residual, 0) + 1
    return fac


def sigma_from_factorization(fac: dict[int, int]) -> int:
    result = 1
    for p, exponent in fac.items():
        result *= (p ** (exponent + 1) - 1) // (p - 1)
    return result


def sigma_Mn(n: int, factors: list[int], spf: list[int]) -> int:
    fac = factor_from_spf(n, spf)
    for p in factors:
        fac[p] = fac.get(p, 0) + 1
    return sigma_from_factorization(fac)


def sigma_values_for_progression(
    M: int,
    nmax: int,
    factor_primes: list[int],
) -> list[int]:
    """Compute sigma(Mn + 1) for 1 <= n <= nmax by sieving a progression."""
    residuals = [M * n + 1 for n in range(1, nmax + 1)]
    sigma_values = [1] * nmax
    max_value = M * nmax + 1

    for p in factor_primes:
        if p * p > max_value:
            break
        if M % p == 0:
            continue
        residue = (-pow(M % p, -1, p)) % p
        first_n = residue if residue != 0 else p
        if first_n > nmax:
            continue
        for n in range(first_n, nmax + 1, p):
            idx = n - 1
            if residuals[idx] % p != 0:
                continue
            term = 1
            power = 1
            while residuals[idx] % p == 0:
                residuals[idx] //= p
                power *= p
                term += power
            sigma_values[idx] *= term

    for idx, residual in enumerate(residuals):
        if residual > 1:
            sigma_values[idx] *= 1 + residual
    return sigma_values


def make_hit(n: int, M: int, sigma_mn: int, sigma_z: int) -> dict[str, Any]:
    return {
        "n": n,
        "z": M * n + 1,
        "sigma_Mn": sigma_mn,
        "sigma_z": sigma_z,
        "ratio": sigma_z / sigma_mn,
        "found_under_batch_bound": True,
    }


def sigma_sieve(m: int) -> list[int]:
    sigma = [0] * (m + 1)
    for divisor in range(1, m + 1):
        for multiple in range(divisor, m + 1, divisor):
            sigma[multiple] += divisor
    return sigma


def first_hits_sieve(
    sigma: list[int],
    M: int,
    nmax: int,
    include_reverse: bool,
) -> tuple[dict[str, Any] | None, dict[str, Any] | None]:
    original = None
    reverse = None
    for n in range(1, nmax + 1):
        Mn = M * n
        z = Mn + 1
        sigma_mn = sigma[Mn]
        sigma_z = sigma[z]
        if original is None and sigma_mn < sigma_z:
            original = make_hit(n, M, sigma_mn, sigma_z)
        if include_reverse and reverse is None and sigma_z < sigma_mn:
            reverse = make_hit(n, M, sigma_mn, sigma_z)
        if original is not None and (not include_reverse or reverse is not None):
            break
    return original, reverse


def first_hits_factored(
    M: int,
    factors: list[int],
    nmax: int,
    spf: list[int],
    factor_primes: list[int],
    include_reverse: bool,
) -> tuple[dict[str, Any] | None, dict[str, Any] | None]:
    original = None
    reverse = None
    sigma_z_values = sigma_values_for_progression(M, nmax, factor_primes)
    for n in range(1, nmax + 1):
        sigma_mn = sigma_Mn(n, factors, spf)
        sigma_z = sigma_z_values[n - 1]

        if original is None and sigma_mn < sigma_z:
            original = make_hit(n, M, sigma_mn, sigma_z)
        if include_reverse and reverse is None and sigma_z < sigma_mn:
            reverse = make_hit(n, M, sigma_mn, sigma_z)
        if original is not None and (not include_reverse or reverse is not None):
            break
    return original, reverse


def add_minimality_metadata(row: dict[str, Any], nmax: int) -> None:
    row["checked_n_min"] = 1
    row["checked_n_max"] = nmax
    row["minimality_note"] = (
        f"Exact brute-force search checked 1 <= n <= {nmax}. "
        "A reported hit is first under this checked bound."
    )
    if row["M"] == 30 and row["original"] is None:
        row["external_original_note"] = M30_EXTERNAL_NOTE
        row["external_original_source"] = M30_EXTERNAL_SOURCE


def run_batch_sieve(
    cases: list[dict[str, Any]],
    requested_mmax: int,
    nmax: int,
    include_reverse: bool,
    excluded_m: list[int] | None = None,
) -> dict[str, Any]:
    max_m = max(row["M"] for row in cases) * nmax + 1
    sigma = sigma_sieve(max_m)
    rows = []
    for case in cases:
        original, reverse = first_hits_sieve(
            sigma, case["M"], nmax, include_reverse=include_reverse
        )
        row = {
            "M": case["M"],
            "factors": case["factors"],
            "N": nmax,
            "original": original,
            "reverse": reverse,
        }
        add_minimality_metadata(row, nmax)
        rows.append(row)
    return {
        "generated": datetime.now(timezone.utc).isoformat(),
        "method": "sieve",
        "M_max": requested_mmax,
        "max_tested_M": max(row["M"] for row in cases),
        "N": nmax,
        "include_reverse": include_reverse,
        "excluded_M": excluded_m or [],
        "case_count": len(cases),
        "max_sieve_value": max_m,
        "results": rows,
    }


def run_batch_factored(
    cases: list[dict[str, Any]],
    requested_mmax: int,
    nmax: int,
    include_reverse: bool,
    excluded_m: list[int] | None = None,
) -> dict[str, Any]:
    max_m = max(row["M"] for row in cases) * nmax + 1
    factor_primes = prime_sieve(int(max_m**0.5) + 1)
    spf = smallest_prime_factors(nmax)
    rows = []
    for case in cases:
        original, reverse = first_hits_factored(
            case["M"],
            case["factors"],
            nmax,
            spf,
            factor_primes,
            include_reverse=include_reverse,
        )
        row = {
            "M": case["M"],
            "factors": case["factors"],
            "N": nmax,
            "original": original,
            "reverse": reverse,
        }
        add_minimality_metadata(row, nmax)
        rows.append(row)
    return {
        "generated": datetime.now(timezone.utc).isoformat(),
        "method": "factor",
        "M_max": requested_mmax,
        "max_tested_M": max(row["M"] for row in cases),
        "N": nmax,
        "include_reverse": include_reverse,
        "excluded_M": excluded_m or [],
        "case_count": len(cases),
        "max_sieve_value": None,
        "max_tested_value": max_m,
        "factor_prime_limit": int(max_m**0.5) + 1,
        "results": rows,
    }


def run_batch(
    mmax: int,
    nmax: int,
    include_reverse: bool = True,
    exclude_m: set[int] | None = None,
    method: str = "auto",
    sieve_limit: int = DEFAULT_SIEVE_LIMIT,
) -> dict[str, Any]:
    cases = enumerate_squarefree_triples(mmax)
    excluded = sorted(exclude_m or [])
    if exclude_m:
        cases = [case for case in cases if case["M"] not in exclude_m]
    if not cases:
        raise ValueError(f"no squarefree triple-prime M values found for M <= {mmax}")
    max_m = max(row["M"] for row in cases) * nmax + 1
    if method == "sieve":
        return run_batch_sieve(
            cases,
            mmax,
            nmax,
            include_reverse=include_reverse,
            excluded_m=excluded,
        )
    if method == "factor":
        return run_batch_factored(
            cases,
            mmax,
            nmax,
            include_reverse=include_reverse,
            excluded_m=excluded,
        )
    if method != "auto":
        raise ValueError(f"unknown method: {method}")
    if max_m <= sieve_limit:
        return run_batch_sieve(
            cases,
            mmax,
            nmax,
            include_reverse=include_reverse,
            excluded_m=excluded,
        )
    return run_batch_factored(
        cases,
        mmax,
        nmax,
        include_reverse=include_reverse,
        excluded_m=excluded,
    )


def stringify_large_ints(value: Any) -> Any:
    if isinstance(value, bool) or value is None:
        return value
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return value
    if isinstance(value, dict):
        return {str(k): stringify_large_ints(v) for k, v in value.items()}
    if isinstance(value, list):
        return [stringify_large_ints(item) for item in value]
    return value


def write_json(path: Path, result: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(stringify_large_ints(result), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def hit_text(hit: dict[str, Any] | None, key: str) -> str:
    if hit is None:
        return "not found"
    return str(hit[key])


def factor_text(factors: list[int]) -> str:
    return "*".join(str(p) for p in factors)


def write_report(path: Path, result: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    nmax = result["N"]
    max_value_label = (
        f"- Maximum sieved value: `{result['max_sieve_value']}`"
        if result.get("max_sieve_value") is not None
        else f"- Maximum tested value: `{result['max_tested_value']}`"
    )
    missing_original = [row for row in result["results"] if row["original"] is None]
    missing_reverse = [
        row
        for row in result["results"]
        if result["include_reverse"] and row["reverse"] is None
    ]
    row30 = next((row for row in result["results"] if row["M"] == 30), None)
    nearby30 = [
        row for row in result["results"] if row["M"] in {30, 42, 66, 70, 78, 102}
    ]

    lines = [
        "# Batch Sigma Experiment for Squarefree `M = p*q*r`",
        "",
        f"- Generated: {result['generated']}",
        f"- Method: `{result['method']}`",
        f"- Exact brute-force bound: `1 <= n <= {nmax}`",
        f"- Maximum M bound requested: `M <= {result['M_max']}`",
        f"- Maximum M actually tested: `{result['max_tested_M']}`",
        f"- Total tested M values: `{result['case_count']}`",
        f"- Excluded M values: `{', '.join(str(M) for M in result['excluded_M']) or 'none'}`",
        max_value_label,
        "",
        "## Problems",
        "",
        "```text",
        "Original: sigma(Mn) < sigma(Mn + 1)",
        "Reverse:  sigma(Mn + 1) < sigma(Mn)",
        "M = p*q*r, where p < q < r are distinct primes",
        "```",
        "",
        "All results below are exact first hits found within the stated brute-force",
        "bound. A missing hit means only that no hit was found for `n <= N`;",
        "it is not a proof of nonexistence beyond the checked bound.",
        "",
        "## Results",
        "",
        "| M | factors | original n | original z | sigma(Mn) | sigma(Mn + 1) | reverse n | reverse z | reverse sigma(Mn) | reverse sigma(Mn + 1) |",
        "|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]

    for row in result["results"]:
        original = row["original"]
        reverse = row["reverse"]
        lines.append(
            f"| {row['M']} | {factor_text(row['factors'])} | "
            f"{hit_text(original, 'n')} | {hit_text(original, 'z')} | "
            f"{hit_text(original, 'sigma_Mn')} | {hit_text(original, 'sigma_z')} | "
            f"{hit_text(reverse, 'n')} | {hit_text(reverse, 'z')} | "
            f"{hit_text(reverse, 'sigma_Mn')} | {hit_text(reverse, 'sigma_z')} |"
        )

    lines.extend(["", "## Not Found Under Bound", ""])
    if not missing_original and not missing_reverse:
        lines.append("- No missing hits under the checked bound.")
    else:
        lines.append("| direction | M | factors | note |")
        lines.append("|:---|---:|:---|:---|")
        for row in missing_original:
            note = f"No original hit found for `1 <= n <= {nmax}`."
            if row["M"] == 30:
                note += f" {M30_EXTERNAL_NOTE}"
            lines.append(
                f"| original | {row['M']} | {factor_text(row['factors'])} | {note} |"
            )
        for row in missing_reverse:
            lines.append(
                f"| reverse | {row['M']} | {factor_text(row['factors'])} | "
                f"No reverse hit found for `1 <= n <= {nmax}`. |"
            )

    lines.extend(["", "## Around `M = 30`", ""])
    if nearby30:
        lines.append(
            "| M | factors | original n | original z | reverse n | reverse z |"
        )
        lines.append("|---:|:---|---:|---:|---:|---:|")
        for row in nearby30:
            lines.append(
                f"| {row['M']} | {factor_text(row['factors'])} | "
                f"{hit_text(row['original'], 'n')} | {hit_text(row['original'], 'z')} | "
                f"{hit_text(row['reverse'], 'n')} | {hit_text(row['reverse'], 'z')} |"
            )
    if row30 and row30.get("external_original_note"):
        lines.extend(
            [
                "",
                f"- For `M = 30`, this report keeps the batch result separate from the existing segmented evidence: {M30_EXTERNAL_NOTE}",
                f"- Source: `{M30_EXTERNAL_SOURCE}`",
            ]
        )

    original_found = result["case_count"] - len(missing_original)
    reverse_found = (
        result["case_count"] - len(missing_reverse)
        if result["include_reverse"]
        else "not tested"
    )
    lines.extend(
        [
            "",
            "## Observations",
            "",
            f"- Original-direction hits found under the bound: `{original_found}` of `{result['case_count']}`.",
            f"- Reverse-direction hits found under the bound: `{reverse_found}`.",
            "- Cases unresolved under this first bound are candidates for a later segmented or candidate-search pass.",
        ]
    )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Batch exact searches for squarefree M=p*q*r."
    )
    parser.add_argument(
        "--M-max",
        dest="mmax",
        type=parse_int,
        default=DEFAULT_M_MAX,
        help="Maximum squarefree triple-prime M to include.",
    )
    parser.add_argument(
        "--N",
        type=parse_int,
        default=DEFAULT_N,
        help="Maximum n to test. Default: 100000.",
    )
    parser.add_argument(
        "--exclude-M",
        type=parse_int_list,
        default=[],
        help="Comma-separated M values to omit from the batch, e.g. 30.",
    )
    reverse_group = parser.add_mutually_exclusive_group()
    reverse_group.add_argument(
        "--include-reverse",
        dest="include_reverse",
        action="store_true",
        default=True,
        help="Test both directions. This is the default.",
    )
    reverse_group.add_argument(
        "--original-only",
        dest="include_reverse",
        action="store_false",
        help="Only test sigma(Mn) < sigma(Mn + 1).",
    )
    parser.add_argument(
        "--method",
        choices=("auto", "sieve", "factor"),
        default="auto",
        help="Search method. auto uses sieve for small max values and factorization otherwise.",
    )
    parser.add_argument(
        "--sieve-limit",
        type=parse_int,
        default=DEFAULT_SIEVE_LIMIT,
        help="Maximum max(MN+1) for auto to use a single sigma sieve.",
    )
    parser.add_argument("--json", type=Path, help="Optional JSON summary path.")
    parser.add_argument("--report", type=Path, help="Optional Markdown report path.")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    result = run_batch(
        args.mmax,
        args.N,
        include_reverse=args.include_reverse,
        exclude_m=set(args.exclude_M),
        method=args.method,
        sieve_limit=args.sieve_limit,
    )
    print(json.dumps(stringify_large_ints(result), sort_keys=True))

    if args.json:
        write_json(args.json, result)
    if args.report:
        write_report(args.report, result)
    return 0


if __name__ == "__main__":
    sys.exit(main())
