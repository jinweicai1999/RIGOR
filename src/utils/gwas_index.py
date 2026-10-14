"""Fast local GWAS window queries backed by per-trait/per-chromosome Parquet files."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Iterable

import pandas as pd

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
except Exception:  # pragma: no cover - handled by caller at runtime.
    pa = None
    pq = None


GWAS_COLUMNS = ["CHR", "SNP", "BP", "BETA", "P"]
INDEX_VERSION = 1


def default_index_dir(gwas_dir: str | os.PathLike[str]) -> Path:
    return Path(gwas_dir) / "gwas_parquet_index"


def trait_source_path(gwas_dir: str | os.PathLike[str], trait: str) -> Path:
    return Path(gwas_dir) / f"gwas_{trait}.assoc.linear"


def trait_index_dir(index_dir: str | os.PathLike[str], trait: str) -> Path:
    return Path(index_dir) / trait


def chrom_parquet_path(index_dir: str | os.PathLike[str], trait: str, chrom: int | str) -> Path:
    chrom_int = int(chrom)
    return trait_index_dir(index_dir, trait) / f"chr{chrom_int:02d}.parquet"


def metadata_path(index_dir: str | os.PathLike[str]) -> Path:
    return Path(index_dir) / "_metadata.json"


def read_index_metadata(index_dir: str | os.PathLike[str]) -> dict:
    path = metadata_path(index_dir)
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def write_index_metadata(index_dir: str | os.PathLike[str], metadata: dict) -> None:
    path = metadata_path(index_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {"index_version": INDEX_VERSION, **metadata}
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def available_indexed_traits(index_dir: str | os.PathLike[str], traits: Iterable[str]) -> list[str]:
    root = Path(index_dir)
    found = []
    for trait in traits:
        trait_dir = root / trait
        if trait_dir.exists() and any(trait_dir.glob("chr*.parquet")):
            found.append(trait)
    return found


def has_trait_index(index_dir: str | os.PathLike[str], trait: str) -> bool:
    trait_dir = trait_index_dir(index_dir, trait)
    return trait_dir.exists() and any(trait_dir.glob("chr*.parquet"))


def normalize_gwas_chunk(chunk: pd.DataFrame) -> pd.DataFrame:
    """Normalize PLINK GWAS columns to compact, typed columns used by the index."""
    missing = [col for col in GWAS_COLUMNS if col not in chunk.columns]
    if missing:
        raise ValueError(f"GWAS input missing required columns: {missing}")
    out = chunk.loc[:, GWAS_COLUMNS].copy()
    out["CHR"] = pd.to_numeric(out["CHR"], errors="coerce").astype("Int16")
    out["BP"] = pd.to_numeric(out["BP"], errors="coerce").astype("Int64")
    out["BETA"] = pd.to_numeric(out["BETA"], errors="coerce").astype("float32")
    out["P"] = pd.to_numeric(out["P"], errors="coerce").astype("float64")
    out["SNP"] = out["SNP"].astype(str)
    out = out.dropna(subset=["CHR", "BP", "BETA", "P"])
    return out


def write_trait_parquet_index(
    source_path: str | os.PathLike[str],
    trait: str,
    index_dir: str | os.PathLike[str],
    chunksize: int = 500_000,
    compression: str = "zstd",
) -> dict:
    """Stream one large PLINK association file into chr-specific Parquet files."""
    if pa is None or pq is None:
        raise RuntimeError("pyarrow is required to build the GWAS Parquet index.")

    source = Path(source_path)
    if not source.exists():
        raise FileNotFoundError(source)

    out_dir = trait_index_dir(index_dir, trait)
    out_dir.mkdir(parents=True, exist_ok=True)

    writers: dict[int, pq.ParquetWriter] = {}
    row_counts: dict[str, int] = {}
    try:
        reader = pd.read_csv(
            source,
            sep=r"\s+",
            chunksize=chunksize,
            usecols=GWAS_COLUMNS,
            engine="c",
        )
        total_rows = 0
        for chunk in reader:
            chunk = normalize_gwas_chunk(chunk)
            if chunk.empty:
                continue
            for chrom, sub in chunk.groupby("CHR", sort=False):
                chrom_int = int(chrom)
                table = pa.Table.from_pandas(sub.reset_index(drop=True), preserve_index=False)
                if chrom_int not in writers:
                    path = chrom_parquet_path(index_dir, trait, chrom_int)
                    writers[chrom_int] = pq.ParquetWriter(
                        path,
                        table.schema,
                        compression=compression,
                        use_dictionary=["SNP"],
                        write_statistics=True,
                    )
                writers[chrom_int].write_table(table)
                row_counts[f"chr{chrom_int:02d}"] = row_counts.get(f"chr{chrom_int:02d}", 0) + len(sub)
                total_rows += len(sub)
    finally:
        for writer in writers.values():
            writer.close()

    return {
        "trait": trait,
        "source_file": str(source),
        "indexed_rows": int(total_rows),
        "chromosomes": row_counts,
    }


def query_trait_window(
    index_dir: str | os.PathLike[str],
    trait: str,
    chrom: int,
    start: int,
    end: int,
    p_threshold: float = 0.05,
) -> pd.DataFrame:
    """Query one trait/chromosome Parquet shard for SNPs in a coordinate window."""
    if pq is None:
        raise RuntimeError("pyarrow is required for fast GWAS index queries.")
    path = chrom_parquet_path(index_dir, trait, chrom)
    if not path.exists():
        return pd.DataFrame(columns=GWAS_COLUMNS)
    filters = [
        ("BP", ">=", int(start)),
        ("BP", "<=", int(end)),
        ("P", "<", float(p_threshold)),
    ]
    table = pq.read_table(path, columns=GWAS_COLUMNS, filters=filters)
    if table.num_rows == 0:
        return pd.DataFrame(columns=GWAS_COLUMNS)
    return table.to_pandas()


def query_trait_chrom_significant(
    index_dir: str | os.PathLike[str],
    trait: str,
    chrom: int,
    p_threshold: float = 0.05,
) -> pd.DataFrame:
    """Load significant SNPs for one trait/chromosome shard once for batch window matching."""
    if pq is None:
        raise RuntimeError("pyarrow is required for fast GWAS index queries.")
    path = chrom_parquet_path(index_dir, trait, chrom)
    if not path.exists():
        return pd.DataFrame(columns=GWAS_COLUMNS)
    table = pq.read_table(path, columns=GWAS_COLUMNS, filters=[("P", "<", float(p_threshold))])
    if table.num_rows == 0:
        return pd.DataFrame(columns=GWAS_COLUMNS)
    return table.to_pandas()
