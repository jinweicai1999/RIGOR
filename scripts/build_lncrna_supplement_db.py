#!/usr/bin/env python3
from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SUPPLEMENT_XLSX = PROJECT_ROOT / "Supplement Tables.xlsx"
NOVEL_GTF = PROJECT_ROOT / "data" / "ncrna" / "novel_lncRNA.gtf"
OUT_DIR = PROJECT_ROOT / "data" / "ncrna" / "lncrna_supplement"


def normalize_chr(text: str) -> str:
    text = str(text).strip().lower().replace("chr.", "").replace("chr", "")
    return text.lstrip("0") or "0"


def clean_text(value):
    if pd.isna(value):
        return ""
    return str(value).strip()


def load_project_intervals():
    intervals = defaultdict(list)
    with open(NOVEL_GTF, "r") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "transcript":
                continue
            match = re.search(r'gene_id "([^"]+)"', parts[8])
            if not match:
                continue
            gene_id = match.group(1)
            intervals[normalize_chr(parts[0])].append((int(parts[3]), int(parts[4]), gene_id))
    for chrom in intervals:
        intervals[chrom].sort()
    return intervals


def best_map(intervals, chrom: str, start: int, end: int):
    qlen = end - start + 1
    best = None
    for proj_start, proj_end, project_id in intervals.get(chrom, []):
        if proj_end < start:
            continue
        if proj_start > end:
            break
        overlap = min(end, proj_end) - max(start, proj_start) + 1
        if overlap <= 0:
            continue
        score = overlap / min(qlen, proj_end - proj_start + 1)
        if best is None or score > best[0]:
            best = (score, project_id)
    return best


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    intervals = load_project_intervals()

    s2 = pd.read_excel(
        SUPPLEMENT_XLSX,
        sheet_name="Sheet2",
        header=1,
        usecols=["lncRNA ID1", "Chr", "Start", "End", "Tissue2", "lncRNA class3", "Gene ID4", "Gene name5"],
    ).dropna(subset=["lncRNA ID1", "Chr", "Start", "End"])

    mappings = []
    for row in s2.itertuples(index=False):
        supplement_id = clean_text(row[0])
        chrom = normalize_chr(row[1])
        start = int(row[2])
        end = int(row[3])
        match = best_map(intervals, chrom, start, end)
        if not match or match[0] < 0.8:
            continue
        mappings.append(
            {
                "project_lnc_id": match[1],
                "supplement_lnc_id": supplement_id,
                "overlap_score": round(float(match[0]), 6),
                "chromosome": clean_text(row[1]),
                "start": start,
                "end": end,
                "tissue": clean_text(row[4]),
                "lncrna_class": clean_text(row[5]),
                "linked_gene_id": clean_text(row[6]),
                "linked_gene_name": clean_text(row[7]),
            }
        )
    mapping_df = pd.DataFrame(mappings).sort_values(
        ["project_lnc_id", "overlap_score", "supplement_lnc_id"], ascending=[True, False, True]
    )
    mapping_df.to_csv(OUT_DIR / "mapping.tsv", sep="\t", index=False)

    mapped_ids = set(mapping_df["supplement_lnc_id"].tolist())

    s5 = pd.read_excel(
        SUPPLEMENT_XLSX,
        sheet_name="Sheet5",
        header=1,
        usecols=[
            "lncRNA",
            "Chr",
            "Start",
            "End",
            "Function1",
            "Partner mRNA2",
            "Pearson correlation coefficient3",
            "P_value4",
            "State class6",
            "Tissue7",
        ],
    ).dropna(subset=["lncRNA", "Partner mRNA2"])
    s5["supplement_lnc_id"] = s5["lncRNA"].astype(str).str.strip()
    s5 = s5[s5["supplement_lnc_id"].isin(mapped_ids)].copy()
    s5["correlation"] = pd.to_numeric(s5["Pearson correlation coefficient3"], errors="coerce")
    s5["p_value"] = pd.to_numeric(s5["P_value4"], errors="coerce")
    s5["regulation"] = s5["Function1"].map(clean_text)
    s5["partner_gene_id"] = s5["Partner mRNA2"].map(clean_text)
    s5["state_class"] = s5["State class6"].map(clean_text)
    s5["tissue"] = s5["Tissue7"].map(clean_text)
    coexp_df = s5[
        ["supplement_lnc_id", "regulation", "partner_gene_id", "correlation", "p_value", "state_class", "tissue"]
    ].drop_duplicates()
    coexp_df.to_csv(OUT_DIR / "coexpression_pairs.tsv", sep="\t", index=False)

    s9 = pd.read_excel(
        SUPPLEMENT_XLSX,
        sheet_name="Sheet9",
        header=1,
        usecols=["lncRNA", "SNP", "Trait", "P-value1"],
    )
    s9["supplement_lnc_id"] = s9["lncRNA"].astype(str).str.strip()
    s9 = s9[s9["supplement_lnc_id"].isin(mapped_ids)].copy()
    s9["trait"] = s9["Trait"].map(clean_text)
    s9["snp"] = s9["SNP"].map(clean_text)
    s9["p_value"] = pd.to_numeric(s9["P-value1"], errors="coerce")
    gwas_hits_df = s9[["supplement_lnc_id", "snp", "trait", "p_value"]].drop_duplicates()
    gwas_hits_df.to_csv(OUT_DIR / "gwas_hits.tsv", sep="\t", index=False)

    s10 = pd.read_excel(
        SUPPLEMENT_XLSX,
        sheet_name="Sheet10",
        header=1,
        usecols=["lncRNA", "SNP", "PVE (%)", "Trait"],
    )
    s10["supplement_lnc_id"] = s10["lncRNA"].astype(str).str.strip()
    s10 = s10[s10["supplement_lnc_id"].isin(mapped_ids)].copy()
    s10["trait"] = s10["Trait"].map(clean_text)
    s10["snp"] = s10["SNP"].map(clean_text)
    s10["pve"] = pd.to_numeric(s10["PVE (%)"], errors="coerce")
    gwas_pve_df = s10[["supplement_lnc_id", "snp", "trait", "pve"]].drop_duplicates()
    gwas_pve_df.to_csv(OUT_DIR / "gwas_pve.tsv", sep="\t", index=False)

    summary = pd.DataFrame(
        [
            {"metric": "mapped_project_lncRNAs", "value": int(mapping_df["project_lnc_id"].nunique())},
            {"metric": "mapped_supplement_lncRNAs", "value": int(mapping_df["supplement_lnc_id"].nunique())},
            {"metric": "mapping_rows", "value": int(len(mapping_df))},
            {"metric": "coexpression_rows", "value": int(len(coexp_df))},
            {"metric": "gwas_hit_rows", "value": int(len(gwas_hits_df))},
            {"metric": "gwas_pve_rows", "value": int(len(gwas_pve_df))},
        ]
    )
    summary.to_csv(OUT_DIR / "summary.tsv", sep="\t", index=False)

    print(f"Wrote internal lncRNA supplement tables to {OUT_DIR}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
