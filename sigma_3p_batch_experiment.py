#!/usr/bin/env python3
"""Batch exact searches for sigma(3pn) versus sigma(3pn + 1)."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


DEFAULT_PRIMES = [5, 7, 11, 13, 17, 19, 23, 29, 31]
DEFAULT_SIEVE_LIMIT = 50_000_000
MULTIPLIER = 3


def parse_int(value: str) -> int:
    """Parse CLI integers, accepting forms like 10^6 and 1_000_000."""
    cleaned = value.strip().replace("_", "")
    if "^" in cleaned:
        base, exponent = cleaned.split("^", 1)
        return int(base) ** int(exponent)
    return int(cleaned)


def parse_prime_list(value: str) -> list[int]:
    primes = []
    for part in value.split(","):
        part = part.strip()
        if not part:
            continue
        primes.append(parse_int(part))
    if not primes:
        raise argparse.ArgumentTypeError("at least one prime must be provided")
    return primes


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


def factor_by_trial_division(n: int, primes: list[int]) -> dict[int, int]:
    fac: dict[int, int] = {}
    residual = n
    for p in primes:
        if p * p > residual:
            break
        while residual % p == 0:
            fac[p] = fac.get(p, 0) + 1
            residual //= p
    if residual > 1:
        fac[residual] = fac.get(residual, 0) + 1
    return fac


def sigma_from_factorization(fac: dict[int, int]) -> int:
    result = 1
    for p, exponent in fac.items():
        result *= (p ** (exponent + 1) - 1) // (p - 1)
    return result


def sigma_3pn(n: int, p: int, spf: list[int]) -> int:
    fac = factor_from_spf(n, spf)
    fac[MULTIPLIER] = fac.get(MULTIPLIER, 0) + 1
    fac[p] = fac.get(p, 0) + 1
    return sigma_from_factorization(fac)


def make_hit(n: int, M: int, sigma_Mn: int, sigma_z: int) -> dict[str, Any]:
    return {
        "n": n,
        "z": M * n + 1,
        "sigma_Mn": sigma_Mn,
        "sigma_z": sigma_z,
        "ratio": sigma_z / sigma_Mn,
    }


def sigma_sieve(m: int) -> list[int]:
    sigma = [0] * (m + 1)
    for divisor in range(1, m + 1):
        for multiple in range(divisor, m + 1, divisor):
            sigma[multiple] += divisor
    return sigma


def first_hit_for_direction(
    sigma: list[int],
    M: int,
    nmax: int,
    reverse: bool,
) -> dict[str, Any] | None:
    for n in range(1, nmax + 1):
        Mn = M * n
        z = Mn + 1
        sigma_Mn = sigma[Mn]
        sigma_z = sigma[z]
        success = sigma_z < sigma_Mn if reverse else sigma_Mn < sigma_z
        if success:
            return make_hit(n, M, sigma_Mn, sigma_z)
    return None


def first_hits_factored(
    p: int,
    nmax: int,
    spf: list[int],
    factor_primes: list[int],
) -> tuple[dict[str, Any] | None, dict[str, Any] | None]:
    M = MULTIPLIER * p
    original = None
    reverse = None

    for n in range(1, nmax + 1):
        sigma_Mn = sigma_3pn(n, p, spf)
        z = M * n + 1
        sigma_z = sigma_from_factorization(factor_by_trial_division(z, factor_primes))

        if original is None and sigma_Mn < sigma_z:
            original = make_hit(n, M, sigma_Mn, sigma_z)
            original["found_under_batch_bound"] = True
        if reverse is None and sigma_z < sigma_Mn:
            reverse = make_hit(n, M, sigma_Mn, sigma_z)
            reverse["found_under_batch_bound"] = True
        if original is not None and reverse is not None:
            break

    return original, reverse


def run_batch_sieve(primes: list[int], nmax: int) -> dict[str, Any]:
    max_m = MULTIPLIER * max(primes) * nmax + 1
    sigma = sigma_sieve(max_m)
    rows = []

    for p in primes:
        M = MULTIPLIER * p
        original = first_hit_for_direction(sigma, M, nmax, reverse=False)
        reverse = first_hit_for_direction(sigma, M, nmax, reverse=True)
        if original is not None:
            original["found_under_batch_bound"] = True

        if reverse is not None:
            reverse["found_under_batch_bound"] = True

        rows.append(
            {
                "p": p,
                "M": M,
                "N": nmax,
                "original": original,
                "reverse": reverse,
                "minimality": (
                    f"Exact brute-force search checked 1 <= n <= {nmax}. "
                    "A reported hit is first under this bound."
                ),
            }
        )

    return {
        "generated": datetime.now(timezone.utc).isoformat(),
        "method": "sieve",
        "N": nmax,
        "primes": primes,
        "max_sieve_value": max_m,
        "results": rows,
    }


def run_batch_factored(primes: list[int], nmax: int) -> dict[str, Any]:
    max_m = MULTIPLIER * max(primes) * nmax + 1
    factor_primes = prime_sieve(int(max_m**0.5) + 1)
    spf = smallest_prime_factors(nmax)
    rows = []

    for p in primes:
        original, reverse = first_hits_factored(p, nmax, spf, factor_primes)
        rows.append(
            {
                "p": p,
                "M": MULTIPLIER * p,
                "N": nmax,
                "original": original,
                "reverse": reverse,
                "minimality": (
                    f"Exact factorization search checked 1 <= n <= {nmax}. "
                    "A reported hit is first under this bound unless an external "
                    "minimality source is listed."
                ),
            }
        )

    return {
        "generated": datetime.now(timezone.utc).isoformat(),
        "method": "factor",
        "N": nmax,
        "primes": primes,
        "max_sieve_value": None,
        "max_tested_value": max_m,
        "factor_prime_limit": int(max_m**0.5) + 1,
        "results": rows,
    }


def run_batch(
    primes: list[int],
    nmax: int,
    method: str = "auto",
    sieve_limit: int = DEFAULT_SIEVE_LIMIT,
) -> dict[str, Any]:
    max_m = MULTIPLIER * max(primes) * nmax + 1
    if method == "sieve":
        return run_batch_sieve(primes, nmax)
    if method == "factor":
        return run_batch_factored(primes, nmax)
    if method != "auto":
        raise ValueError(f"unknown method: {method}")
    if max_m <= sieve_limit:
        return run_batch_sieve(primes, nmax)
    return run_batch_factored(primes, nmax)


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


def write_report(path: Path, result: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    nmax = result["N"]
    prime_cases = (
        f"{len(result['primes'])} primes, min `p = {min(result['primes'])}`, "
        f"max `p = {max(result['primes'])}`"
        if len(result["primes"]) > 40
        else f"`{', '.join(str(p) for p in result['primes'])}`"
    )
    max_value_label = (
        f"- Maximum sieved value: `{result['max_sieve_value']}`"
        if result.get("max_sieve_value") is not None
        else f"- Maximum tested value: `{result['max_tested_value']}`"
    )
    missing_original = [
        str(row["p"]) for row in result["results"] if row["original"] is None
    ]
    lines = [
        "# Batch Sigma Experiment for `M = 3p`",
        "",
        f"- Generated: {result['generated']}",
        f"- Method: `{result['method']}`",
        f"- Exact brute-force bound: `1 <= n <= {nmax}`",
        f"- Prime cases: {prime_cases}",
        max_value_label,
        "",
        "## Problems",
        "",
        "```text",
        "Original: sigma(Mn) < sigma(Mn + 1)",
        "Reverse:  sigma(Mn + 1) < sigma(Mn)",
        "M = 3p",
        "```",
        "",
        "All results below are exact first hits found within the stated brute-force",
        "bound. A missing hit would mean only that no hit was found for `n <= N`;",
        "it would not be a proof of nonexistence.",
        "",
        "## Results",
        "",
        "| p | M = 3p | original n | original z | sigma(Mn) | sigma(Mn + 1) | reverse n | reverse z | reverse sigma(Mn) | reverse sigma(Mn + 1) |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]

    for row in result["results"]:
        original = row["original"]
        reverse = row["reverse"]
        lines.append(
            f"| {row['p']} | {row['M']} | "
            f"{hit_text(original, 'n')} | {hit_text(original, 'z')} | "
            f"{hit_text(original, 'sigma_Mn')} | {hit_text(original, 'sigma_z')} | "
            f"{hit_text(reverse, 'n')} | {hit_text(reverse, 'z')} | "
            f"{hit_text(reverse, 'sigma_Mn')} | {hit_text(reverse, 'sigma_z')} |"
        )

    lines.extend(
        [
            "",
            "## Observations",
            "",
            "- The original direction is immediate in this batch: every tested prime has first hit `n = 1`.",
            "- The reverse direction is also uniform in this batch: every tested prime has first hit `n = 2`.",
        ]
    )
    if missing_original:
        lines.append(
            "- Original-direction hits were not found under the bound for: "
            + ", ".join(missing_original)
            + "."
        )
    else:
        lines.append("- Every original-direction case in this report has a known first hit.")

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Batch exact searches for sigma(3pn) versus sigma(3pn + 1)."
    )
    parser.add_argument(
        "--primes",
        type=parse_prime_list,
        default=DEFAULT_PRIMES,
        help="Comma-separated prime list. Default: 5,7,11,13,17,19,23,29,31.",
    )
    parser.add_argument(
        "--prime-max",
        type=parse_int,
        help="Use all primes p with 3 < p <= PRIME_MAX. Overrides --primes.",
    )
    parser.add_argument("--N", type=parse_int, required=True, help="Maximum n to test.")
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
        help="Maximum max(3pN+1) for auto to use a single sigma sieve.",
    )
    parser.add_argument("--json", type=Path, help="Optional JSON summary path.")
    parser.add_argument("--report", type=Path, help="Optional Markdown report path.")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    primes = [p for p in prime_sieve(args.prime_max) if p > 3] if args.prime_max else args.primes
    result = run_batch(primes, args.N, method=args.method, sieve_limit=args.sieve_limit)
    print(json.dumps(stringify_large_ints(result), sort_keys=True))

    if args.json:
        write_json(args.json, result)
    if args.report:
        write_report(args.report, result)
    return 0


if __name__ == "__main__":
    sys.exit(main())
