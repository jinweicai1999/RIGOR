#!/usr/bin/env python3
"""Build a fast per-trait/per-chromosome Parquet index for local GWAS scans."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.agents.gwas_agent import GwasAgent
from src.project_paths import PROCESSED_DIR
from src.utils.gwas_index import (
    default_index_dir,
    trait_source_path,
    write_index_metadata,
    write_trait_parquet_index,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build GWAS Parquet index.")
    parser.add_argument("--gwas-dir", default=str(PROCESSED_DIR), help="Directory containing gwas_<Trait>.assoc.linear files.")
    parser.add_argument("--index-dir", default=None, help="Output index directory. Default: <gwas-dir>/gwas_parquet_index")
    parser.add_argument("--traits", nargs="*", default=None, help="Traits to index. Default: all RIGOR GWAS traits.")
    parser.add_argument("--chunksize", type=int, default=500_000, help="Rows per text chunk while converting.")
    parser.add_argument("--compression", default="zstd", choices=["zstd", "snappy", "gzip", "brotli", "none"], help="Parquet compression.")
    parser.add_argument("--force", action="store_true", help="Rebuild existing trait index files.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    gwas_dir = Path(args.gwas_dir)
    index_dir = Path(args.index_dir) if args.index_dir else default_index_dir(gwas_dir)
    traits = args.traits or GwasAgent().traits
    compression = None if args.compression == "none" else args.compression

    index_dir.mkdir(parents=True, exist_ok=True)
    summaries = []
    missing = []
    started = datetime.now()

    print(f"[GWAS index] source dir: {gwas_dir}")
    print(f"[GWAS index] output dir: {index_dir}")
    print(f"[GWAS index] traits: {', '.join(traits)}")

    for idx, trait in enumerate(traits, start=1):
        source = trait_source_path(gwas_dir, trait)
        trait_dir = index_dir / trait
        existing = trait_dir.exists() and any(trait_dir.glob("chr*.parquet"))
        if existing and not args.force:
            print(f"[{idx}/{len(traits)}] {trait}: existing index found, skip. Use --force to rebuild.")
            summaries.append({"trait": trait, "skipped_existing": True, "index_dir": str(trait_dir)})
            continue
        if not source.exists():
            print(f"[{idx}/{len(traits)}] {trait}: missing source file {source}")
            missing.append(str(source))
            continue
        if trait_dir.exists() and args.force:
            for path in trait_dir.glob("chr*.parquet"):
                path.unlink()
        print(f"[{idx}/{len(traits)}] {trait}: indexing {source.name} ...", flush=True)
        summary = write_trait_parquet_index(
            source,
            trait=trait,
            index_dir=index_dir,
            chunksize=args.chunksize,
            compression=compression,
        )
        summaries.append(summary)
        print(f"[{idx}/{len(traits)}] {trait}: indexed {summary['indexed_rows']:,} rows.")

    metadata = {
        "created_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "elapsed_seconds": round((datetime.now() - started).total_seconds(), 3),
        "gwas_dir": str(gwas_dir),
        "index_dir": str(index_dir),
        "traits": traits,
        "summaries": summaries,
        "missing_sources": missing,
        "layout": "trait/chrXX.parquet",
        "query_columns": ["CHR", "SNP", "BP", "BETA", "P"],
    }
    write_index_metadata(index_dir, metadata)
    print("[GWAS index] done.")
    print(json.dumps(metadata, ensure_ascii=False, indent=2)[:4000])
    return 0 if not missing else 1


if __name__ == "__main__":
    raise SystemExit(main())
