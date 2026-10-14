"""Build trait-wise GWAS P-value background cutoffs for normalized evidence scoring."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.project_paths import PROCESSED_DIR
from src.utils.gwas_index import chrom_parquet_path, default_index_dir


TRAITS = [
    "Heading_date",
    "Plant_height",
    "Num_panicles",
    "Num_effective_panicles",
    "Yield",
    "Grain_weight",
    "Spikelet_length",
    "Grain_length",
    "Grain_width",
    "Grain_thickness",
]

TOP_FRACTIONS = [0.001, 0.01, 0.05]


def collect_trait_pvalues(index_dir: Path, trait: str) -> np.ndarray:
    arrays = []
    for chrom in range(1, 13):
        path = chrom_parquet_path(index_dir, trait, chrom)
        if not path.exists():
            continue
        df = pd.read_parquet(path, columns=["P"])
        values = pd.to_numeric(df["P"], errors="coerce").dropna().to_numpy(dtype="float64")
        if values.size:
            arrays.append(values)
    if not arrays:
        return np.array([], dtype="float64")
    return np.concatenate(arrays)


def build_background(index_dir: Path, traits: list[str]) -> dict:
    payload = {
        "method": "trait-wise empirical P-value distribution",
        "evidence_levels": {
            "Strong": "candidate SNP P value is within the trait-specific top 0.1% genome-wide signals",
            "Moderate": "candidate SNP P value is within the trait-specific top 1% genome-wide signals",
            "Weak": "candidate SNP P value is within the trait-specific top 5% genome-wide signals",
        },
        "traits": {},
    }
    for trait in traits:
        pvalues = collect_trait_pvalues(index_dir, trait)
        if pvalues.size == 0:
            continue
        cutoffs = {f"top_{fraction:g}": float(np.quantile(pvalues, fraction)) for fraction in TOP_FRACTIONS}
        payload["traits"][trait] = {
            "total_snps": int(pvalues.size),
            "min_p": float(np.min(pvalues)),
            "median_p": float(np.median(pvalues)),
            "p_cutoffs": cutoffs,
        }
        print(
            trait,
            f"n={pvalues.size}",
            " ".join([f"top {fraction * 100:g}% <= {cutoffs[f'top_{fraction:g}']:.3e}" for fraction in TOP_FRACTIONS]),
        )
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description="Build trait-wise GWAS background P-value cutoffs.")
    parser.add_argument("--index-dir", default=str(default_index_dir(PROCESSED_DIR)))
    parser.add_argument("--out", default=str(PROCESSED_DIR / "gwas_trait_background.json"))
    parser.add_argument("--traits", nargs="*", default=TRAITS)
    args = parser.parse_args()

    index_dir = Path(args.index_dir)
    out_path = Path(args.out)
    payload = build_background(index_dir, args.traits)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(out_path)


if __name__ == "__main__":
    main()
