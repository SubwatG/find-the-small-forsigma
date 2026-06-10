#!/usr/bin/env python3
"""CLI tools for the Hermes sigma 30n experiment.

Searches for candidates with sigma(30n) < sigma(30n + 1), written in terms of
z = 30n + 1 so that the exact test is sigma(z) > sigma(z - 1).
"""


from __future__ import annotations

import argparse
import glob
import heapq
import itertools
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from math import gcd, log
from pathlib import Path
from typing import Any

from sympy import factorint


PRIME_REPO_URL = "https://github.com/ByteJoseph/primes"
TARGET_M = 30
_SEGMENTED_WORKER_PRIMES: list[int] | None = None


@dataclass(frozen=True)
class Candidate:
    z: int
    log_a: float
    factorization: dict[int, int]


def parse_int(value: str) -> int:
    """Parse CLI integers, accepting forms like 10^6 and 1_000_000."""
    cleaned = value.strip().replace("_", "")
    if "^" in cleaned:
        base, exponent = cleaned.split("^", 1)
        return int(base) ** int(exponent)
    return int(cleaned)


def parse_int_grid(value: str) -> list[int]:
    return [parse_int(part) for part in value.split(",") if part.strip()]


def parse_float_grid(value: str) -> list[float]:
    return [float(part.strip()) for part in value.split(",") if part.strip()]


def parse_factorization_spec(value: str | None) -> dict[int, int]:
    if value is None or not value.strip():
        return {}

    fac: dict[int, int] = {}
    for item in value.split(","):
        item = item.strip()
        if not item:
            continue
        if ":" in item:
            p_text, a_text = item.split(":", 1)
            p, a = int(p_text), int(a_text)
        else:
            p, a = int(item), 1
        if p <= 1 or a < 1:
            raise ValueError("--seed-fac and --require-primes entries must use p[:a] with p > 1 and a >= 1")
        fac[p] = max(fac.get(p, 0), a)
    return fac


def validate_seed_factorization(fac: dict[int, int]) -> None:
    if any(gcd(p, TARGET_M) != 1 for p in fac):
        raise ValueError("--seed-fac/--require-primes cannot include primes dividing 30")


def load_primes(path: str | Path = "primes.json", pmax: int | None = None) -> list[int]:
    prime_path = Path(path)
    if not prime_path.exists():
        raise FileNotFoundError(
            f"Prime source not found: {prime_path}\n"
            f"Download primes.json from {PRIME_REPO_URL} and pass it with "
            "--primes-json PATH."
        )

    with prime_path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    primes = sorted(int(p) for p in data.keys())
    primes = [p for p in primes if gcd(p, TARGET_M) == 1]
    if pmax is not None:
        primes = [p for p in primes if p <= pmax]
    return primes


def load_all_primes(path: str | Path = "primes.json", pmax: int | None = None) -> list[int]:
    prime_path = Path(path)
    if not prime_path.exists():
        raise FileNotFoundError(
            f"Prime source not found: {prime_path}\n"
            f"Download primes.json from {PRIME_REPO_URL} and pass it with "
            "--primes-json PATH."
        )

    with prime_path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    primes = sorted(int(p) for p in data.keys())
    if pmax is not None:
        primes = [p for p in primes if p <= pmax]
    return primes


def sigma_from_factorization(fac: dict[int, int]) -> int:
    result = 1
    for p, a in fac.items():
        result *= (p ** (a + 1) - 1) // (p - 1)
    return result


def log_a_prime_power(p: int, a: int) -> float:
    sigma_pa = (p ** (a + 1) - 1) // (p - 1)
    return log(sigma_pa) - a * log(p)


def log_a_from_factorization(fac: dict[int, int]) -> float:
    return sum(log_a_prime_power(p, a) for p, a in fac.items())


def small_factor_penalty(y: int, small_primes: list[int]) -> float:
    total = 0.0
    for q in small_primes:
        if TARGET_M % q == 0:
            continue
        if y % q == 0:
            total += log(q / (q - 1))
    return total


def score_candidate(
    z: int,
    log_a: float,
    lam: float,
    eta: float,
    small_primes: list[int],
) -> float:
    return log_a - lam * log(z) - eta * small_factor_penalty(z - 1, small_primes)


def beam_search(
    primes: list[int],
    zmax: int,
    k: int,
    lam: float,
    eta: float,
    small_primes: list[int],
    k_smallest: int = 0,
    seed_fac: dict[int, int] | None = None,
    min_log_a: float | None = None,
) -> list[Candidate]:
    seed_fac = seed_fac or {}
    seed_z = 1
    for p, a in seed_fac.items():
        seed_z *= p**a
    if seed_z > zmax:
        return []

    seed_log_a = log_a_from_factorization(seed_fac) if seed_fac else 0.0
    beam = [Candidate(seed_z, seed_log_a, seed_fac.copy())]
    primes = [p for p in primes if p not in seed_fac]

    for p in primes:
        new_candidates: list[Candidate] = []

        for cand in beam:
            new_candidates.append(cand)

            pp = p
            exponent = 1
            while cand.z * pp <= zmax:
                new_fac = cand.factorization.copy()
                new_fac[p] = exponent
                new_z = cand.z * pp
                new_log_a = cand.log_a + log_a_prime_power(p, exponent)
                new_candidates.append(Candidate(new_z, new_log_a, new_fac))

                exponent += 1
                pp *= p

        buckets: dict[int, list[Candidate]] = {
            residue: [] for residue in range(TARGET_M)
        }
        for cand in new_candidates:
            buckets[cand.z % TARGET_M].append(cand)

        beam = []
        for residue in range(TARGET_M):
            by_score = sorted(
                buckets[residue],
                key=lambda c: score_candidate(c.z, c.log_a, lam, eta, small_primes),
                reverse=True,
            )
            kept = {cand.z: cand for cand in by_score[:k]}
            if k_smallest > 0:
                for cand in sorted(buckets[residue], key=lambda c: c.z)[:k_smallest]:
                    kept[cand.z] = cand
            beam.extend(kept.values())

    final = [cand for cand in beam if cand.z % TARGET_M == 1]
    if min_log_a is not None:
        final = [cand for cand in final if cand.log_a >= min_log_a]
    return final


def beam_search_diagnostic(
    primes: list[int],
    zmax: int,
    k: int,
    lam: float,
    eta: float,
    small_primes: list[int],
    k_smallest: int = 0,
    seed_fac: dict[int, int] | None = None,
    min_log_a: float | None = None,
    verbose: bool = False,
    top_n_diag: int = 10,
) -> tuple[list[Candidate], list[Candidate]]:
    """Beam search returning (final, top_all) where top_all is top-N by
    abundancy regardless of residue constraint.  Also prints diagnostic
    info when verbose=True.
    """
    import math as _math

    seed_fac = seed_fac or {}
    seed_z = 1
    for p, a in seed_fac.items():
        seed_z *= p**a
    if seed_z > zmax:
        return [], []

    seed_log_a = log_a_from_factorization(seed_fac) if seed_fac else 0.0
    beam = [Candidate(seed_z, seed_log_a, seed_fac.copy())]
    primes = [p for p in primes if p not in seed_fac]

    for idx, p in enumerate(primes):
        new_candidates: list[Candidate] = []
        for cand in beam:
            new_candidates.append(cand)
            pp = p
            exponent = 1
            while cand.z * pp <= zmax:
                new_fac = cand.factorization.copy()
                new_fac[p] = exponent
                new_z = cand.z * pp
                new_log_a = cand.log_a + log_a_prime_power(p, exponent)
                new_candidates.append(Candidate(new_z, new_log_a, new_fac))
                exponent += 1
                pp *= p

        buckets: dict[int, list[Candidate]] = {
            residue: [] for residue in range(TARGET_M)}
        for cand in new_candidates:
            buckets[cand.z % TARGET_M].append(cand)

        beam = []
        for residue in range(TARGET_M):
            by_score = sorted(
                buckets[residue],
                key=lambda c: score_candidate(c.z, c.log_a, lam, eta, small_primes),
                reverse=True,
            )
            kept = {cand.z: cand for cand in by_score[:k]}
            if k_smallest > 0:
                for cand in sorted(
                    buckets[residue], key=lambda c: c.z
                )[:k_smallest]:
                    kept[cand.z] = cand
            beam.extend(kept.values())

        if verbose and (idx + 1) % 5 == 0:
            top_ab = sorted(beam, key=lambda c: c.log_a, reverse=True)[:3]
            print(
                f"  prime {idx+1}/{len(primes)} (p={p}): "
                f"beam_size={len(beam)}, "
                f"top_logA={top_ab[0].log_a:.5f} "
                f"(A={_math.exp(top_ab[0].log_a):.4f}), "
                f"z_digits={len(str(top_ab[0].z))}"
            )

    final = [cand for cand in beam if cand.z % TARGET_M == 1]
    if min_log_a is not None:
        final = [cand for cand in final if cand.log_a >= min_log_a]

    top_all = sorted(beam, key=lambda c: c.log_a, reverse=True)[:top_n_diag]
    return final, top_all


def frontier_search(
    primes: list[int],
    zmax: int,
    frontier_size: int,
    lam: float,
    eta: float,
    small_primes: list[int],
    k_smallest: int = 0,
    seed_fac: dict[int, int] | None = None,
    min_log_a: float | None = None,
) -> list[Candidate]:
    seed_fac = seed_fac or {}
    seed_z = 1
    for p, a in seed_fac.items():
        seed_z *= p**a
    if seed_z > zmax:
        return []

    seed_log_a = log_a_from_factorization(seed_fac) if seed_fac else 0.0
    frontier = [Candidate(seed_z, seed_log_a, seed_fac.copy())]
    primes = [p for p in primes if p not in seed_fac]
    per_residue = max(1, frontier_size // TARGET_M)

    for p in primes:
        expanded: list[Candidate] = []
        for cand in frontier:
            expanded.append(cand)
            pp = p
            exponent = 1
            while cand.z * pp <= zmax:
                new_fac = cand.factorization.copy()
                new_fac[p] = exponent
                expanded.append(
                    Candidate(
                        cand.z * pp,
                        cand.log_a + log_a_prime_power(p, exponent),
                        new_fac,
                    )
                )
                exponent += 1
                pp *= p

        buckets: dict[int, list[Candidate]] = {
            residue: [] for residue in range(TARGET_M)
        }
        for cand in expanded:
            buckets[cand.z % TARGET_M].append(cand)

        next_frontier_by_z: dict[int, Candidate] = {}
        for residue, bucket in buckets.items():
            if not bucket:
                continue
            scored = heapq.nlargest(
                per_residue,
                bucket,
                key=lambda c: score_candidate(c.z, c.log_a, lam, eta, small_primes),
            )
            for cand in scored:
                next_frontier_by_z[cand.z] = cand
            if k_smallest > 0:
                for cand in heapq.nsmallest(k_smallest, bucket, key=lambda c: c.z):
                    next_frontier_by_z[cand.z] = cand

        frontier = list(next_frontier_by_z.values())

    final = [cand for cand in frontier if cand.z % TARGET_M == 1]
    if min_log_a is not None:
        final = [cand for cand in final if cand.log_a >= min_log_a]
    return final


def sigma_sieve(m: int) -> list[int]:
    sigma = [0] * (m + 1)
    for divisor in range(1, m + 1):
        for multiple in range(divisor, m + 1, divisor):
            sigma[multiple] += divisor
    return sigma


def brute_force_first(nmax: int) -> dict[str, Any] | None:
    m = TARGET_M * nmax + 1
    sigma = sigma_sieve(m)
    for n in range(1, nmax + 1):
        if sigma[TARGET_M * n] < sigma[TARGET_M * n + 1]:
            return {
                "M": TARGET_M,
                "n": n,
                "z": TARGET_M * n + 1,
                "sigma_Mn": sigma[TARGET_M * n],
                "sigma_z": sigma[TARGET_M * n + 1],
                "ratio": sigma[TARGET_M * n + 1] / sigma[TARGET_M * n],
            }
    return None


def initialize_segmented_worker(primes: list[int]) -> None:
    global _SEGMENTED_WORKER_PRIMES
    _SEGMENTED_WORKER_PRIMES = primes


def sigma_values_for_linear_block(
    n_start: int,
    block_len: int,
    multiplier: int,
    offset: int,
    primes: list[int],
) -> list[int]:
    values = [multiplier * (n_start + i) + offset for i in range(block_len)]
    residuals = values.copy()
    sigmas = [1] * block_len
    max_value = values[-1]

    for p in primes:
        if p * p > max_value:
            break

        if multiplier % p == 0:
            if offset % p != 0:
                continue
            first_index = 0
            step = 1
        else:
            residue = (-offset * pow(multiplier, -1, p)) % p
            first_n = n_start + ((residue - n_start) % p)
            first_index = first_n - n_start
            step = p

        for idx in range(first_index, block_len, step):
            exponent = 0
            power = 1
            while residuals[idx] % p == 0:
                residuals[idx] //= p
                exponent += 1
                power *= p
            if exponent:
                sigmas[idx] *= (power * p - 1) // (p - 1)

    for idx, residual in enumerate(residuals):
        if residual > 1:
            sigmas[idx] *= residual + 1

    return sigmas


def check_segmented_block(block: tuple[int, int], primes: list[int]) -> dict[str, Any]:
    block_start, block_end = block
    block_len = block_end - block_start + 1
    sigma_Mn = sigma_values_for_linear_block(block_start, block_len, TARGET_M, 0, primes)
    sigma_Mn_plus_1 = sigma_values_for_linear_block(
        block_start,
        block_len,
        TARGET_M,
        1,
        primes,
    )

    for idx in range(block_len):
        if sigma_Mn[idx] < sigma_Mn_plus_1[idx]:
            n = block_start + idx
            return {
                "M": TARGET_M,
                "n": n,
                "z": TARGET_M * n + 1,
                "sigma_Mn": sigma_Mn[idx],
                "sigma_z": sigma_Mn_plus_1[idx],
                "block_start": block_start,
                "block_end": block_end,
                "found": True,
            }

    return {
        "M": TARGET_M,
        "block_start": block_start,
        "block_end": block_end,
        "checked": block_len,
        "found": False,
    }


def check_segmented_block_worker(block: tuple[int, int]) -> dict[str, Any]:
    if _SEGMENTED_WORKER_PRIMES is None:
        raise RuntimeError("segmented worker was not initialized with primes")
    return check_segmented_block(block, _SEGMENTED_WORKER_PRIMES)


def iter_blocks(n_start: int, n_end: int, block_size: int) -> list[tuple[int, int]]:
    blocks = []
    current = n_start
    while current <= n_end:
        block_end = min(current + block_size - 1, n_end)
        blocks.append((current, block_end))
        current = block_end + 1
    return blocks


def write_progress_row(handle: Any, row: dict[str, Any]) -> None:
    progress_row = row.copy()
    progress_row["timestamp"] = datetime.now(timezone.utc).isoformat()
    handle.write(json.dumps(stringify_large_ints(progress_row), sort_keys=True) + "\n")
    handle.flush()


def segmented_first(
    n_start: int,
    n_end: int,
    block_size: int,
    primes: list[int],
    progress_jsonl: Path | None = None,
    workers: int = 1,
) -> dict[str, Any] | None:
    if n_start < 1:
        raise ValueError("n_start must be >= 1")
    if n_end < n_start:
        raise ValueError("n_end must be >= n_start")
    if block_size < 1:
        raise ValueError("block_size must be >= 1")
    if workers < 1:
        raise ValueError("workers must be >= 1")

    progress_handle = None
    if progress_jsonl is not None:
        progress_jsonl.parent.mkdir(parents=True, exist_ok=True)
        progress_handle = progress_jsonl.open("a", encoding="utf-8")

    try:
        blocks = iter_blocks(n_start, n_end, block_size)
        if workers == 1:
            results = (check_segmented_block(block, primes) for block in blocks)
            for row in results:
                if progress_handle is not None:
                    write_progress_row(progress_handle, row)
                if row.get("found"):
                    return row
        else:
            with ProcessPoolExecutor(
                max_workers=workers,
                initializer=initialize_segmented_worker,
                initargs=(primes,),
            ) as executor:
                for row in executor.map(check_segmented_block_worker, blocks, chunksize=1):
                    if progress_handle is not None:
                        write_progress_row(progress_handle, row)
                    if row.get("found"):
                        return row
    finally:
        if progress_handle is not None:
            progress_handle.close()

    return None


def read_progress_ranges(path: Path) -> list[tuple[int, int]]:
    if not path.exists():
        return []

    ranges = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            ranges.append((int(row["block_start"]), int(row["block_end"])))
    return ranges


def merge_ranges(ranges: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not ranges:
        return []

    merged: list[list[int]] = []
    for start, end in sorted(ranges):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        elif end > merged[-1][1]:
            merged[-1][1] = end
    return [(start, end) for start, end in merged]


def find_unchecked_ranges(
    progress_jsonl: Path,
    n_start: int,
    n_end: int,
) -> list[tuple[int, int]]:
    checked = []
    for start, end in read_progress_ranges(progress_jsonl):
        overlap_start = max(start, n_start)
        overlap_end = min(end, n_end)
        if overlap_start <= overlap_end:
            checked.append((overlap_start, overlap_end))

    merged = merge_ranges(checked)
    gaps = []
    cursor = n_start
    for start, end in merged:
        if cursor < start:
            gaps.append((cursor, start - 1))
        cursor = max(cursor, end + 1)
    if cursor <= n_end:
        gaps.append((cursor, n_end))
    return gaps


def read_progress_rows(paths: list[Path]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in paths:
        if not path.exists():
            continue
        with path.open("r", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip():
                    continue
                try:
                    row = json.loads(line)
                except json.JSONDecodeError as exc:
                    print(
                        f"warning: skipping malformed JSONL row in {path}:{line_number}: {exc}",
                        file=sys.stderr,
                    )
                    continue
                row["_source"] = str(path)
                rows.append(row)
    return rows


def summarize_progress_coverage(
    paths: list[Path],
    n_start: int,
    n_end: int,
) -> dict[str, Any]:
    rows = read_progress_rows(paths)
    checked_ranges: list[tuple[int, int]] = []
    success_rows = []

    for row in rows:
        if "block_start" not in row or "block_end" not in row:
            continue
        block_start = int(row["block_start"])
        block_end = int(row["block_end"])
        checked_ranges.append((block_start, block_end))
        if row.get("found") is True or row.get("success") is True:
            success_rows.append(row)

    clipped = []
    for start, end in checked_ranges:
        overlap_start = max(start, n_start)
        overlap_end = min(end, n_end)
        if overlap_start <= overlap_end:
            clipped.append((overlap_start, overlap_end))

    merged = merge_ranges(clipped)
    gaps = []
    cursor = n_start
    highest_contiguous = n_start - 1
    first_gap_seen = False
    for start, end in merged:
        if start > cursor:
            gaps.append((cursor, start - 1))
            if not first_gap_seen:
                highest_contiguous = cursor - 1
                first_gap_seen = True
        if start <= cursor:
            if not first_gap_seen:
                highest_contiguous = max(highest_contiguous, end)
        cursor = max(cursor, end + 1)
    if cursor <= n_end:
        gaps.append((cursor, n_end))
        if not first_gap_seen:
            highest_contiguous = cursor - 1
    elif not first_gap_seen:
        highest_contiguous = n_end

    smallest_success = None
    if success_rows:
        smallest_success = min(success_rows, key=lambda row: int(row.get("n", row["block_start"])))

    return {
        "M": TARGET_M,
        "n_start": n_start,
        "n_end": n_end,
        "progress_files": [str(path) for path in paths],
        "progress_file_count": len(paths),
        "progress_rows": len(rows),
        "checked_ranges": merged,
        "gaps": gaps,
        "complete": not gaps,
        "highest_contiguous_from_start": highest_contiguous,
        "success_count": len(success_rows),
        "smallest_success": smallest_success,
    }


def read_jsonl_rows(paths: list[Path]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in paths:
        if not path.exists():
            continue
        with path.open("r", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip():
                    continue
                try:
                    row = json.loads(line)
                except json.JSONDecodeError as exc:
                    print(
                        f"warning: skipping malformed JSONL row in {path}:{line_number}: {exc}",
                        file=sys.stderr,
                    )
                    continue
                row["_source"] = str(path)
                rows.append(row)
    return rows


def numeric_value(row: dict[str, Any], key: str, default: float = 0.0) -> float:
    value = row.get(key, default)
    if value is None:
        return default
    return float(value)


def int_value(row: dict[str, Any], key: str, default: int = 0) -> int:
    value = row.get(key, default)
    if value is None:
        return default
    return int(value)


def stringify_factorization(fac: dict[int, int]) -> dict[str, str]:
    return {str(p): str(a) for p, a in sorted(fac.items())}


def stringify_large_ints(value: Any) -> Any:
    if isinstance(value, bool) or value is None:
        return value
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return value
    if isinstance(value, dict):
        return {str(k): stringify_large_ints(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [stringify_large_ints(item) for item in value]
    return value


def test_candidate(z: int, fac_z: dict[int, int]) -> dict[str, Any]:
    sigma_z = sigma_from_factorization(fac_z)
    fac_z_minus_1 = {int(p): int(a) for p, a in factorint(z - 1).items()}
    sigma_z_minus_1 = sigma_from_factorization(fac_z_minus_1)
    ratio = sigma_z / sigma_z_minus_1

    return {
        "M": TARGET_M,
        "z": z,
        "n": (z - 1) // TARGET_M if z % TARGET_M == 1 else None,
        "residue_mod_30": z % TARGET_M,
        "fac_z": fac_z,
        "fac_z_minus_1": fac_z_minus_1,
        "sigma_z": sigma_z,
        "sigma_z_minus_1": sigma_z_minus_1,
        "success": sigma_z > sigma_z_minus_1,
        "ratio": ratio,
        "A_z": sigma_z / z,
        "A_z_minus_1": sigma_z_minus_1 / (z - 1),
    }


def write_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(stringify_large_ints(row), sort_keys=True) + "\n")


def append_jsonl_row(path: Path, row: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(stringify_large_ints(row), sort_keys=True) + "\n")
        handle.flush()


def append_progress_log(path: Path | None, message: str) -> None:
    if path is None:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(message.rstrip() + "\n")
        handle.flush()


def checkpoint_completed_runs(path: Path | None) -> set[int]:
    if path is None or not path.exists():
        return set()
    completed: set[int] = set()
    for row in read_jsonl_rows([path]):
        if row.get("status") == "complete" and "grid_run" in row:
            completed.add(int_value(row, "grid_run"))
    return completed


def write_report(
    path: Path,
    command: str,
    params: dict[str, Any],
    rows: list[dict[str, Any]],
    minimality_note: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    successes = [row for row in rows if row.get("success")]
    best = min(successes, key=lambda row: int(row["z"])) if successes else None

    lines = [
        "# Hermes Sigma Experiment Report",
        "",
        f"- Generated: {datetime.now(timezone.utc).isoformat()}",
        f"- Command: `{command}`",
        f"- Results written: {len(rows)}",
        f"- Successful candidates: {len(successes)}",
        f"- Minimality: {minimality_note}",
        "",
        "## Parameters",
        "",
        "```json",
        json.dumps(stringify_large_ints(params), indent=2, sort_keys=True),
        "```",
        "",
        "## Smallest Successful Candidate",
        "",
    ]

    if best is None:
        lines.append("No successful candidate was found in this run.")
    elif "sigma_z_minus_1" not in best:
        lines.extend(
            [
                f"- n: `{best.get('n')}`",
                f"- z: `{best.get('z')}`",
                f"- M: `{best.get('M')}`",
                f"- sigma(Mn): `{best.get('sigma_Mn')}`",
                f"- sigma(z): `{best.get('sigma_z')}`",
                "",
                "Full row:",
                "",
                "```json",
                json.dumps(stringify_large_ints(best), indent=2, sort_keys=True),
                "```",
            ]
        )
    else:
        lines.extend(
            [
                f"- n: `{best['n']}`",
                f"- M: `{best.get('M', TARGET_M)}`",
                f"- z: `{best['z']}`",
                f"- sigma(z): `{best['sigma_z']}`",
                f"- sigma(z - 1): `{best['sigma_z_minus_1']}`",
                f"- sigma ratio: `{best['ratio']}`",
                f"- A(z): `{best['A_z']}`",
                f"- A(z - 1): `{best['A_z_minus_1']}`",
                f"- factorization z: `{best['fac_z']}`",
                f"- factorization z - 1: `{best['fac_z_minus_1']}`",
            ]
        )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def enrich_result(
    result: dict[str, Any],
    params: dict[str, Any],
    log_a: float | None = None,
) -> dict[str, Any]:
    enriched = result.copy()
    if log_a is not None:
        enriched["logA_z"] = log_a
    enriched["params"] = params
    enriched["fac_z"] = stringify_factorization(enriched["fac_z"])
    enriched["fac_z_minus_1"] = stringify_factorization(enriched["fac_z_minus_1"])
    return enriched


def command_brute(args: argparse.Namespace) -> int:
    result = brute_force_first(args.N)
    params = {"M": TARGET_M, "N": args.N}
    rows: list[dict[str, Any]] = []

    if result is not None:
        rows.append({**result, "success": True, "params": params})
        minimality_note = f"Brute force proves this is the first hit for n <= {args.N}."
        print(json.dumps(stringify_large_ints(result), sort_keys=True))
    else:
        minimality_note = f"Brute force found no hit for n <= {args.N}."
        print("None")

    if args.jsonl:
        write_jsonl(args.jsonl, rows)
    if args.report:
        write_report(args.report, "brute", params, rows, minimality_note)
    return 0


def command_test(args: argparse.Namespace) -> int:
    z = args.z
    fac_z = {int(p): int(a) for p, a in factorint(z).items()}
    result = test_candidate(z, fac_z)
    params = {"M": TARGET_M, "z": z}
    row = enrich_result(result, params, log_a_from_factorization(fac_z))
    minimality_note = "Explicit z test only; minimality not proven."

    print(json.dumps(stringify_large_ints(row), sort_keys=True))
    if z % TARGET_M != 1:
        print("Warning: z is not congruent to 1 modulo 30; n is not integral.")

    if args.jsonl:
        write_jsonl(args.jsonl, [row])
    if args.report:
        write_report(args.report, "test", params, [row], minimality_note)
    return 0


def command_beam(args: argparse.Namespace) -> int:
    primes = load_primes(args.primes_json, args.pmax)
    seed_fac = parse_factorization_spec(args.seed_fac)
    for p, a in parse_factorization_spec(args.require_primes).items():
        seed_fac[p] = max(seed_fac.get(p, 0), a)
    if any(gcd(p, TARGET_M) != 1 for p in seed_fac):
        raise ValueError("--seed-fac/--require-primes cannot include primes dividing 30")

    small_primes = [p for p in primes if p <= args.penalty_qmax]
    min_log_a = log(args.min_A) if args.min_A is not None else args.min_logA
    params = {
        "M": TARGET_M,
        "ZMAX": args.ZMAX,
        "K": args.K,
        "k_smallest": args.k_smallest,
        "pmax": args.pmax,
        "lambda": args.lam,
        "eta": args.eta,
        "penalty_qmax": args.penalty_qmax,
        "primes_json": str(args.primes_json),
        "max_tests": args.max_tests,
        "near_ratio": args.near_ratio,
        "min_logA": min_log_a,
        "min_A": args.min_A,
        "seed_fac": stringify_factorization(seed_fac),
    }

    candidates = beam_search(
        primes,
        args.ZMAX,
        args.K,
        args.lam,
        args.eta,
        small_primes,
        args.k_smallest,
        seed_fac,
        min_log_a,
    )
    candidates.sort(
        key=lambda c: score_candidate(c.z, c.log_a, args.lam, args.eta, small_primes),
        reverse=True,
    )

    rows: list[dict[str, Any]] = []
    for cand in candidates[: args.max_tests]:
        if cand.z <= 1:
            continue
        result = test_candidate(cand.z, cand.factorization)
        if result["success"] or result["ratio"] >= args.near_ratio:
            rows.append(enrich_result(result, params, cand.log_a))

    rows.sort(key=lambda row: (not row["success"], -float(row["ratio"]), int(row["z"])))
    successes = [row for row in rows if row["success"]]
    minimality_note = "Candidate search only; minimality not proven."

    print(
        json.dumps(
            stringify_large_ints(
                {
                    "generated_candidates": len(candidates),
                    "tested_candidates": min(len(candidates), args.max_tests),
                    "reported_rows": len(rows),
                    "successful_candidates": len(successes),
                    "best": rows[0] if rows else None,
                }
            ),
            sort_keys=True,
        )
    )

    if args.jsonl:
        write_jsonl(args.jsonl, rows)
    if args.report:
        write_report(args.report, "beam", params, rows, minimality_note)
    return 0


def command_coverage_summary(args: argparse.Namespace) -> int:
    paths: list[Path] = []
    progress_globs = args.progress_glob or ["segmented_30n_*_progress.jsonl"]
    for pattern in progress_globs:
        paths.extend(Path(path) for path in glob.glob(pattern))
    paths = sorted(set(paths))
    result = summarize_progress_coverage(paths, args.n_start, args.n_end)

    print(json.dumps(stringify_large_ints(result), indent=2, sort_keys=True))
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(
            json.dumps(stringify_large_ints(result), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    return 0


def run_beam_grid_job(job: dict[str, Any]) -> dict[str, Any]:
    primes = load_primes(job["primes_json"], job["pmax"])
    small_primes = [p for p in primes if p <= job["penalty_qmax"]]
    candidates = beam_search(
        primes,
        job["zmax"],
        job["k"],
        job["lam"],
        job["eta"],
        small_primes,
        job["k_smallest"],
        job["seed_fac"],
        job["min_log_a"],
    )
    candidates.sort(
        key=lambda c: score_candidate(c.z, c.log_a, job["lam"], job["eta"], small_primes),
        reverse=True,
    )
    params = {
        "M": TARGET_M,
        "ZMAX": job["zmax"],
        "K": job["k"],
        "k_smallest": job["k_smallest"],
        "pmax": job["pmax"],
        "lambda": job["lam"],
        "eta": job["eta"],
        "penalty_qmax": job["penalty_qmax"],
        "primes_json": str(job["primes_json"]),
        "near_ratio": job["near_ratio"],
        "min_logA": job["min_log_a"],
        "min_A": job["min_A"],
        "seed_fac": stringify_factorization(job["seed_fac"]),
        "grid_run": job["grid_run"],
    }

    rows = []
    tested = 0
    for cand in candidates[: job["max_tests_per_run"]]:
        if cand.z <= 1:
            continue
        tested += 1
        result = test_candidate(cand.z, cand.factorization)
        if result["success"] or result["ratio"] >= job["near_ratio"]:
            rows.append(enrich_result(result, params, cand.log_a))

    return {
        "grid_run": job["grid_run"],
        "generated_candidates": len(candidates),
        "tested_candidates": tested,
        "rows": rows,
    }


def command_beam_grid(args: argparse.Namespace) -> int:
    if args.workers < 1:
        raise ValueError("--workers must be >= 1")
    if args.resume_checkpoint and args.checkpoint_jsonl is None:
        raise ValueError("--resume-checkpoint requires --checkpoint-jsonl")
    zmax_values = parse_int_grid(args.ZMAX)
    k_values = parse_int_grid(args.K)
    pmax_values = parse_int_grid(args.pmax)
    lam_values = parse_float_grid(args.lam)
    eta_values = parse_float_grid(args.eta)
    penalty_values = parse_int_grid(args.penalty_qmax)
    seed_fac = parse_factorization_spec(args.seed_fac)
    for p, a in parse_factorization_spec(args.require_primes).items():
        seed_fac[p] = max(seed_fac.get(p, 0), a)
    validate_seed_factorization(seed_fac)
    min_log_a = log(args.min_A) if args.min_A is not None else args.min_logA

    rows_by_z: dict[int, dict[str, Any]] = {}
    total_generated = 0
    total_tested = 0
    jobs = []
    for run_count, (zmax, k, pmax, lam, eta, penalty_qmax) in enumerate(
        itertools.product(
            zmax_values,
            k_values,
            pmax_values,
            lam_values,
            eta_values,
            penalty_values,
        ),
        start=1,
    ):
        jobs.append(
            {
                "grid_run": run_count,
                "primes_json": args.primes_json,
                "zmax": zmax,
                "k": k,
                "k_smallest": args.k_smallest,
                "pmax": pmax,
                "lam": lam,
                "eta": eta,
                "penalty_qmax": penalty_qmax,
                "max_tests_per_run": args.max_tests_per_run,
                "near_ratio": args.near_ratio,
                "min_log_a": min_log_a,
                "min_A": args.min_A,
                "seed_fac": seed_fac,
            }
        )

    current_run_ids = {int(job["grid_run"]) for job in jobs}
    completed_runs = checkpoint_completed_runs(args.checkpoint_jsonl) if args.resume_checkpoint else set()
    completed_runs &= current_run_ids
    pending_jobs = [job for job in jobs if int(job["grid_run"]) not in completed_runs]
    total_runs = len(jobs)

    def record_job_result(job_result: dict[str, Any]) -> None:
        nonlocal total_generated, total_tested
        total_generated += int(job_result["generated_candidates"])
        total_tested += int(job_result["tested_candidates"])
        for row in job_result["rows"]:
            z = int(row["z"])
            old = rows_by_z.get(z)
            if old is None or numeric_value(row, "ratio") > numeric_value(old, "ratio"):
                rows_by_z[z] = row
        if args.checkpoint_jsonl:
            append_jsonl_row(
                args.checkpoint_jsonl,
                {
                    "status": "complete",
                    "grid_run": job_result["grid_run"],
                    "generated_candidates": job_result["generated_candidates"],
                    "tested_candidates": job_result["tested_candidates"],
                    "reported_rows": len(job_result["rows"]),
                    "rows": job_result["rows"],
                    "completed_at": datetime.now(timezone.utc).isoformat(),
                },
            )
        append_progress_log(
            args.progress_log,
            (
                f"{datetime.now(timezone.utc).isoformat()} complete "
                f"grid_run={job_result['grid_run']}/{total_runs} "
                f"generated={job_result['generated_candidates']} "
                f"tested={job_result['tested_candidates']} "
                f"reported={len(job_result['rows'])}"
            ),
        )

    append_progress_log(
        args.progress_log,
        (
            f"{datetime.now(timezone.utc).isoformat()} start "
            f"runs={total_runs} pending={len(pending_jobs)} "
            f"skipped_from_checkpoint={len(completed_runs)} workers={args.workers}"
        ),
    )

    if args.resume_checkpoint and args.checkpoint_jsonl:
        for checkpoint_row in read_jsonl_rows([args.checkpoint_jsonl]):
            if checkpoint_row.get("status") != "complete":
                continue
            if int_value(checkpoint_row, "grid_run") not in completed_runs:
                continue
            total_generated += int_value(checkpoint_row, "generated_candidates")
            total_tested += int_value(checkpoint_row, "tested_candidates")
            for row in checkpoint_row.get("rows", []):
                if "z" not in row:
                    continue
                z = int_value(row, "z")
                old = rows_by_z.get(z)
                if old is None or numeric_value(row, "ratio") > numeric_value(old, "ratio"):
                    rows_by_z[z] = row

    if args.workers == 1:
        for job in pending_jobs:
            append_progress_log(
                args.progress_log,
                f"{datetime.now(timezone.utc).isoformat()} launch grid_run={job['grid_run']}/{total_runs}",
            )
            record_job_result(run_beam_grid_job(job))
    elif pending_jobs:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            future_to_job = {
                executor.submit(run_beam_grid_job, job): job
                for job in pending_jobs
            }
            for future in as_completed(future_to_job):
                job = future_to_job[future]
                try:
                    record_job_result(future.result())
                except Exception as exc:
                    if args.checkpoint_jsonl:
                        append_jsonl_row(
                            args.checkpoint_jsonl,
                            {
                                "status": "error",
                                "grid_run": job["grid_run"],
                                "error": repr(exc),
                                "completed_at": datetime.now(timezone.utc).isoformat(),
                            },
                        )
                    append_progress_log(
                        args.progress_log,
                        (
                            f"{datetime.now(timezone.utc).isoformat()} error "
                            f"grid_run={job['grid_run']}/{total_runs} error={exc!r}"
                        ),
                    )
                    raise

    rows = list(rows_by_z.values())
    rows.sort(key=lambda row: (not row["success"], -numeric_value(row, "ratio"), int_value(row, "z")))
    successes = [row for row in rows if row.get("success")]

    if args.jsonl:
        write_jsonl(args.jsonl, rows)
    if args.report:
        write_candidate_summary_report(
            args.report,
            rows,
            {
                "runs": len(jobs),
                "workers": args.workers,
                "generated_candidates": total_generated,
                "tested_candidates": total_tested,
                "deduped_reported_rows": len(rows),
                "successful_candidates": len(successes),
            },
            args.top,
        )

    print(
        json.dumps(
            stringify_large_ints(
                {
                    "runs": len(jobs),
                    "pending_runs": len(pending_jobs),
                    "skipped_from_checkpoint": len(completed_runs),
                    "workers": args.workers,
                    "generated_candidates": total_generated,
                    "tested_candidates": total_tested,
                    "reported_rows": len(rows),
                    "successful_candidates": len(successes),
                    "best": rows[0] if rows else None,
                    "jsonl": str(args.jsonl) if args.jsonl else None,
                    "report": str(args.report) if args.report else None,
                    "checkpoint_jsonl": str(args.checkpoint_jsonl) if args.checkpoint_jsonl else None,
                    "progress_log": str(args.progress_log) if args.progress_log else None,
                }
            ),
            sort_keys=True,
        )
    )
    return 0


def write_candidate_summary_report(
    path: Path,
    rows: list[dict[str, Any]],
    params: dict[str, Any],
    top: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    successes = [row for row in rows if row.get("success")]
    shown = rows[:top]
    lines = [
        "# Hermes Sigma Candidate Summary",
        "",
        f"- Generated: {datetime.now(timezone.utc).isoformat()}",
        f"- Rows: {len(rows)}",
        f"- Successful candidates: {len(successes)}",
        "",
        "## Parameters",
        "",
        "```json",
        json.dumps(stringify_large_ints(params), indent=2, sort_keys=True),
        "```",
        "",
        "## Top Candidates",
        "",
    ]
    if not shown:
        lines.append("No candidates met the reporting threshold.")
    else:
        lines.append("| rank | success | n | z | ratio | A(z) | A(z-1) | factorization z |")
        lines.append("|---:|:---:|---:|---:|---:|---:|---:|---|")
        for idx, row in enumerate(shown, start=1):
            lines.append(
                "| "
                f"{idx} | {row.get('success')} | {row.get('n')} | {row.get('z')} | "
                f"{row.get('ratio')} | {row.get('A_z')} | {row.get('A_z_minus_1')} | "
                f"`{row.get('fac_z')}` |"
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def command_candidate_summary(args: argparse.Namespace) -> int:
    paths: list[Path] = []
    for pattern in args.jsonl:
        paths.extend(Path(path) for path in glob.glob(pattern))
    rows_by_z: dict[int, dict[str, Any]] = {}
    for row in read_jsonl_rows(sorted(set(paths))):
        if "z" not in row or "ratio" not in row:
            continue
        if (not row.get("success")) and numeric_value(row, "ratio") < args.near_ratio:
            continue
        z = int_value(row, "z")
        old = rows_by_z.get(z)
        if old is None or numeric_value(row, "ratio") > numeric_value(old, "ratio"):
            rows_by_z[z] = row

    rows = list(rows_by_z.values())
    rows.sort(key=lambda row: (not row.get("success"), -numeric_value(row, "ratio"), int_value(row, "z")))
    summary = {
        "input_files": [str(path) for path in sorted(set(paths))],
        "reported_rows": len(rows),
        "successful_candidates": sum(1 for row in rows if row.get("success")),
        "best": rows[0] if rows else None,
    }
    print(json.dumps(stringify_large_ints(summary), indent=2, sort_keys=True))
    if args.report:
        write_candidate_summary_report(args.report, rows, summary, args.top)
    if args.jsonl_out:
        write_jsonl(args.jsonl_out, rows)
    return 0


def command_frontier_search(args: argparse.Namespace) -> int:
    if args.frontier_size < 1:
        raise ValueError("--frontier-size must be >= 1")
    if args.max_tests < 1:
        raise ValueError("--max-tests must be >= 1")
    seed_fac = parse_factorization_spec(args.seed_fac)
    for p, a in parse_factorization_spec(args.require_primes).items():
        seed_fac[p] = max(seed_fac.get(p, 0), a)
    validate_seed_factorization(seed_fac)
    min_log_a = log(args.min_A) if args.min_A is not None else args.min_logA

    primes = load_primes(args.primes_json, args.pmax)
    small_primes = [p for p in primes if p <= args.penalty_qmax]
    candidates = frontier_search(
        primes,
        args.ZMAX,
        args.frontier_size,
        args.lam,
        args.eta,
        small_primes,
        args.k_smallest,
        seed_fac,
        min_log_a,
    )
    candidates.sort(
        key=lambda c: score_candidate(c.z, c.log_a, args.lam, args.eta, small_primes),
        reverse=True,
    )

    params = {
        "M": TARGET_M,
        "ZMAX": args.ZMAX,
        "frontier_size": args.frontier_size,
        "k_smallest": args.k_smallest,
        "pmax": args.pmax,
        "lambda": args.lam,
        "eta": args.eta,
        "penalty_qmax": args.penalty_qmax,
        "primes_json": str(args.primes_json),
        "max_tests": args.max_tests,
        "near_ratio": args.near_ratio,
        "min_logA": min_log_a,
        "min_A": args.min_A,
        "seed_fac": stringify_factorization(seed_fac),
    }

    rows = []
    for cand in candidates[: args.max_tests]:
        if cand.z <= 1:
            continue
        result = test_candidate(cand.z, cand.factorization)
        if result["success"] or result["ratio"] >= args.near_ratio:
            rows.append(enrich_result(result, params, cand.log_a))
    rows.sort(key=lambda row: (not row["success"], -numeric_value(row, "ratio"), int_value(row, "z")))
    successes = [row for row in rows if row.get("success")]

    if args.jsonl:
        write_jsonl(args.jsonl, rows)
    if args.report:
        write_candidate_summary_report(
            args.report,
            rows,
            {
                "generated_candidates": len(candidates),
                "tested_candidates": min(len(candidates), args.max_tests),
                "reported_rows": len(rows),
                "successful_candidates": len(successes),
                **params,
            },
            args.top,
        )

    print(
        json.dumps(
            stringify_large_ints(
                {
                    "generated_candidates": len(candidates),
                    "tested_candidates": min(len(candidates), args.max_tests),
                    "reported_rows": len(rows),
                    "successful_candidates": len(successes),
                    "best": rows[0] if rows else None,
                    "jsonl": str(args.jsonl) if args.jsonl else None,
                    "report": str(args.report) if args.report else None,
                }
            ),
            sort_keys=True,
        )
    )
    return 0


def command_segmented(args: argparse.Namespace) -> int:
    effective_n_start = args.n_start
    effective_n_end = args.n_end

    if args.resume:
        if args.progress_jsonl is None:
            raise ValueError("--resume requires --progress-jsonl")
        gaps = find_unchecked_ranges(args.progress_jsonl, args.n_start, args.n_end)
        if not gaps:
            params = {
                "M": TARGET_M,
                "n_start": args.n_start,
                "n_end": args.n_end,
                "block_size": args.block_size,
                "workers": args.workers,
                "progress_jsonl": str(args.progress_jsonl),
                "resume": True,
            }
            print(
                json.dumps(
                    stringify_large_ints(
                        {
                            "complete": True,
                            "message": "No unchecked ranges remain.",
                            "params": params,
                        }
                    ),
                    sort_keys=True,
                )
            )
            return 0
        effective_n_start, effective_n_end = gaps[0]

    max_value = TARGET_M * effective_n_end + 1
    prime_limit = int(max_value**0.5) + 1
    primes = load_all_primes(args.primes_json, prime_limit)
    params = {
        "n_start": effective_n_start,
        "n_end": effective_n_end,
        "requested_n_start": args.n_start,
        "requested_n_end": args.n_end,
        "M": TARGET_M,
        "block_size": args.block_size,
        "workers": args.workers,
        "prime_limit": prime_limit,
        "primes_json": str(args.primes_json),
        "progress_jsonl": str(args.progress_jsonl) if args.progress_jsonl else None,
        "resume": args.resume,
    }

    result = segmented_first(
        effective_n_start,
        effective_n_end,
        args.block_size,
        primes,
        args.progress_jsonl,
        args.workers,
    )
    rows: list[dict[str, Any]] = []

    if result is not None:
        row = {**result, "success": True, "params": params}
        rows.append(row)
        minimality_note = (
            f"Segmented exact check found the first hit in checked range "
            f"{effective_n_start} <= n <= {effective_n_end}."
        )
        print(json.dumps(stringify_large_ints(row), sort_keys=True))
    else:
        minimality_note = (
            f"Segmented exact check found no hit for "
            f"{effective_n_start} <= n <= {effective_n_end}."
        )
        print("None")

    if args.jsonl:
        write_jsonl(args.jsonl, rows)
    if args.report:
        write_report(args.report, "segmented", params, rows, minimality_note)
    return 0


def command_progress_check(args: argparse.Namespace) -> int:
    gaps = find_unchecked_ranges(args.progress_jsonl, args.n_start, args.n_end)
    checked_ranges = merge_ranges(
        [
            (max(start, args.n_start), min(end, args.n_end))
            for start, end in read_progress_ranges(args.progress_jsonl)
            if max(start, args.n_start) <= min(end, args.n_end)
        ]
    )
    result = {
        "complete": not gaps,
        "M": TARGET_M,
        "n_start": args.n_start,
        "n_end": args.n_end,
        "checked_ranges": checked_ranges,
        "gaps": gaps,
        "next_n_start": gaps[0][0] if gaps else None,
        "next_n_end": gaps[0][1] if gaps else None,
        "progress_jsonl": str(args.progress_jsonl),
    }
    print(json.dumps(stringify_large_ints(result), indent=2, sort_keys=True))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Hermes sigma experiment CLI for sigma(30n) < sigma(30n + 1)."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    brute = subparsers.add_parser("brute", help="Run a small brute-force sigma sieve.")
    brute.add_argument("--N", type=parse_int, required=True, help="Maximum n to test.")
    brute.add_argument("--jsonl", type=Path, help="Optional JSONL output path.")
    brute.add_argument("--report", type=Path, help="Optional Markdown report path.")
    brute.set_defaults(func=command_brute)

    beam = subparsers.add_parser("beam", help="Run beam-search candidate generation.")
    beam.add_argument("--primes-json", type=Path, default=Path("primes.json"))
    beam.add_argument("--ZMAX", type=parse_int, required=True)
    beam.add_argument("--K", type=int, required=True)
    beam.add_argument("--k-smallest", type=int, default=0)
    beam.add_argument("--pmax", type=int, required=True)
    beam.add_argument("--lam", type=float, default=0.02)
    beam.add_argument("--eta", type=float, default=0.5)
    beam.add_argument("--penalty-qmax", type=int, default=100)
    beam.add_argument("--max-tests", type=int, default=100)
    beam.add_argument("--near-ratio", type=float, default=0.98)
    beam.add_argument("--min-logA", type=float)
    beam.add_argument("--min-A", type=float)
    beam.add_argument("--seed-fac", help="Comma-separated p[:a] factorization to seed every candidate.")
    beam.add_argument("--require-primes", help="Comma-separated p[:a] entries required in every candidate.")
    beam.add_argument("--jsonl", type=Path, default=Path("hermes_sigma_30n_results.jsonl"))
    beam.add_argument("--report", type=Path, default=Path("hermes_sigma_30n_report.md"))
    beam.set_defaults(func=command_beam)

    coverage_summary = subparsers.add_parser(
        "coverage-summary",
        help="Summarize segmented progress coverage across progress JSONL files.",
    )
    coverage_summary.add_argument(
        "--progress-glob",
        action="append",
    )
    coverage_summary.add_argument("--n-start", type=parse_int, default=1)
    coverage_summary.add_argument("--n-end", type=parse_int, required=True)
    coverage_summary.add_argument("--json", type=Path, help="Optional JSON summary output path.")
    coverage_summary.set_defaults(func=command_coverage_summary)

    beam_grid = subparsers.add_parser(
        "beam-grid",
        help="Run a grid of beam-search candidate discovery experiments.",
    )
    beam_grid.add_argument("--primes-json", type=Path, default=Path("primes.json"))
    beam_grid.add_argument("--ZMAX", required=True, help="Comma-separated integer grid.")
    beam_grid.add_argument("--K", required=True, help="Comma-separated integer grid.")
    beam_grid.add_argument("--k-smallest", type=int, default=0)
    beam_grid.add_argument("--pmax", required=True, help="Comma-separated integer grid.")
    beam_grid.add_argument("--lam", default="0.02", help="Comma-separated float grid.")
    beam_grid.add_argument("--eta", default="0.5", help="Comma-separated float grid.")
    beam_grid.add_argument("--penalty-qmax", default="100", help="Comma-separated integer grid.")
    beam_grid.add_argument("--workers", type=int, default=1)
    beam_grid.add_argument("--max-tests-per-run", type=int, default=100)
    beam_grid.add_argument("--near-ratio", type=float, default=0.95)
    beam_grid.add_argument("--min-logA", type=float)
    beam_grid.add_argument("--min-A", type=float)
    beam_grid.add_argument("--seed-fac")
    beam_grid.add_argument("--require-primes")
    beam_grid.add_argument("--top", type=int, default=25)
    beam_grid.add_argument("--jsonl", type=Path, default=Path("beam_grid_30n_candidates.jsonl"))
    beam_grid.add_argument("--report", type=Path, default=Path("beam_grid_30n_candidates_report.md"))
    beam_grid.add_argument(
        "--checkpoint-jsonl",
        type=Path,
        help="Append one checkpoint row after each completed grid run.",
    )
    beam_grid.add_argument(
        "--resume-checkpoint",
        action="store_true",
        help="Skip grid runs already marked complete in --checkpoint-jsonl.",
    )
    beam_grid.add_argument(
        "--progress-log",
        type=Path,
        help="Append human-readable progress messages as grid runs finish.",
    )
    beam_grid.set_defaults(func=command_beam_grid)

    candidate_summary = subparsers.add_parser(
        "candidate-summary",
        help="Summarize successful and near-miss candidate JSONL rows.",
    )
    candidate_summary.add_argument("--jsonl", action="append", required=True)
    candidate_summary.add_argument("--near-ratio", type=float, default=0.95)
    candidate_summary.add_argument("--top", type=int, default=25)
    candidate_summary.add_argument("--report", type=Path)
    candidate_summary.add_argument("--jsonl-out", type=Path)
    candidate_summary.set_defaults(func=command_candidate_summary)

    frontier = subparsers.add_parser(
        "frontier-search",
        help="Run residue-aware frontier candidate discovery.",
    )
    frontier.add_argument("--primes-json", type=Path, default=Path("primes.json"))
    frontier.add_argument("--ZMAX", type=parse_int, required=True)
    frontier.add_argument("--pmax", type=int, required=True)
    frontier.add_argument("--frontier-size", type=int, default=200000)
    frontier.add_argument("--k-smallest", type=int, default=200)
    frontier.add_argument("--lam", type=float, default=0.005)
    frontier.add_argument("--eta", type=float, default=0.25)
    frontier.add_argument("--penalty-qmax", type=int, default=1000)
    frontier.add_argument("--max-tests", type=int, default=5000)
    frontier.add_argument("--near-ratio", type=float, default=0.95)
    frontier.add_argument("--min-logA", type=float)
    frontier.add_argument("--min-A", type=float)
    frontier.add_argument("--seed-fac")
    frontier.add_argument("--require-primes")
    frontier.add_argument("--top", type=int, default=25)
    frontier.add_argument("--jsonl", type=Path, default=Path("frontier_30n_candidates.jsonl"))
    frontier.add_argument("--report", type=Path, default=Path("frontier_30n_candidates_report.md"))
    frontier.set_defaults(func=command_frontier_search)

    explicit = subparsers.add_parser("test", help="Exactly test one explicit z.")
    explicit.add_argument("--z", type=parse_int, required=True)
    explicit.add_argument("--jsonl", type=Path, help="Optional JSONL output path.")
    explicit.add_argument("--report", type=Path, help="Optional Markdown report path.")
    explicit.set_defaults(func=command_test)

    segmented = subparsers.add_parser(
        "segmented",
        help="Run an exact segmented check over a range of n.",
    )
    segmented.add_argument("--primes-json", type=Path, default=Path("primes.json"))
    segmented.add_argument("--n-start", type=parse_int, default=1)
    segmented.add_argument("--n-end", type=parse_int, required=True)
    segmented.add_argument("--block-size", type=parse_int, default=100000)
    segmented.add_argument("--workers", type=int, default=1)
    segmented.add_argument(
        "--resume",
        action="store_true",
        help="Find the first unchecked range from --progress-jsonl and run that range.",
    )
    segmented.add_argument("--progress-jsonl", type=Path)
    segmented.add_argument("--jsonl", type=Path, help="Optional JSONL output path.")
    segmented.add_argument("--report", type=Path, help="Optional Markdown report path.")
    segmented.set_defaults(func=command_segmented)

    progress_check = subparsers.add_parser(
        "progress-check",
        help="Inspect segmented progress coverage and list unchecked ranges.",
    )
    progress_check.add_argument("--progress-jsonl", type=Path, required=True)
    progress_check.add_argument("--n-start", type=parse_int, default=1)
    progress_check.add_argument("--n-end", type=parse_int, required=True)
    progress_check.set_defaults(func=command_progress_check)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except FileNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
