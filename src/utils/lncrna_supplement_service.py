from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import pandas as pd

from src.project_paths import NCRNA_DIR


class LncRnaSupplementService:
    """Load and query internally standardized lncRNA evidence tables."""

    def __init__(
        self,
        data_dir: str | None = None,
    ):
        self.data_dir = Path(data_dir) if data_dir else NCRNA_DIR / "lncrna_supplement"

        self._loaded = False
        self._available = False
        self._supplement_to_project: dict[str, dict[str, Any]] = {}
        self._project_to_supplement: dict[str, list[dict[str, Any]]] = defaultdict(list)
        self.sheet5_df: pd.DataFrame | None = None
        self.sheet9_df: pd.DataFrame | None = None
        self.sheet10_df: pd.DataFrame | None = None

    @property
    def available(self) -> bool:
        if not self._loaded:
            self.load()
        return self._available

    def load(self):
        if self._loaded:
            return
        self._loaded = True

        required = [
            self.data_dir / "mapping.tsv",
            self.data_dir / "coexpression_pairs.tsv",
            self.data_dir / "gwas_hits.tsv",
            self.data_dir / "gwas_pve.tsv",
        ]
        if not all(path.exists() for path in required):
            self._available = False
            return

        self._load_mapping_table()
        self._load_support_tables()
        self._available = True

    def _load_mapping_table(self):
        mapping_df = pd.read_csv(self.data_dir / "mapping.tsv", sep="\t")
        for row in mapping_df.itertuples(index=False):
            payload = {
                "supplement_id": row.supplement_lnc_id,
                "project_id": row.project_lnc_id,
                "overlap_score": float(row.overlap_score),
                "chromosome": self._clean_text(row.chromosome),
                "start": int(row.start),
                "end": int(row.end),
                "tissue": self._clean_text(row.tissue),
                "class": self._clean_text(row.lncrna_class),
                "linked_gene_id": self._clean_text(row.linked_gene_id),
                "linked_gene_name": self._clean_text(row.linked_gene_name),
            }
            self._supplement_to_project[payload["supplement_id"]] = payload
            self._project_to_supplement[payload["project_id"]].append(payload)

        for project_id in self._project_to_supplement:
            self._project_to_supplement[project_id].sort(
                key=lambda item: (-float(item["overlap_score"]), item["supplement_id"])
            )

    def _load_support_tables(self):
        self.sheet5_df = pd.read_csv(self.data_dir / "coexpression_pairs.tsv", sep="\t")
        self.sheet5_df["supplement_id"] = self.sheet5_df["supplement_lnc_id"].astype(str).str.strip()
        self.sheet5_df["regulation"] = self.sheet5_df["regulation"].map(self._clean_text)
        self.sheet5_df["partner_gene"] = self.sheet5_df["partner_gene_id"].astype(str).str.strip()
        self.sheet5_df["correlation"] = pd.to_numeric(self.sheet5_df["correlation"], errors="coerce")
        self.sheet5_df["p_value"] = pd.to_numeric(self.sheet5_df["p_value"], errors="coerce")
        self.sheet5_df["state_class"] = self.sheet5_df["state_class"].map(self._clean_text)
        self.sheet5_df["tissue"] = self.sheet5_df["tissue"].map(self._clean_text)

        self.sheet9_df = pd.read_csv(self.data_dir / "gwas_hits.tsv", sep="\t")
        self.sheet9_df["supplement_id"] = self.sheet9_df["supplement_lnc_id"].astype(str).str.strip()
        self.sheet9_df["trait"] = self.sheet9_df["trait"].map(self._clean_text)
        self.sheet9_df["snp"] = self.sheet9_df["snp"].astype(str).str.strip()
        self.sheet9_df["p_value"] = pd.to_numeric(self.sheet9_df["p_value"], errors="coerce")

        self.sheet10_df = pd.read_csv(self.data_dir / "gwas_pve.tsv", sep="\t")
        self.sheet10_df["supplement_id"] = self.sheet10_df["supplement_lnc_id"].astype(str).str.strip()
        self.sheet10_df["trait"] = self.sheet10_df["trait"].map(self._clean_text)
        self.sheet10_df["snp"] = self.sheet10_df["snp"].astype(str).str.strip()
        self.sheet10_df["pve"] = pd.to_numeric(self.sheet10_df["pve"], errors="coerce")

    def get_mapping(self, project_lnc_id: str) -> list[dict[str, Any]]:
        if not self.available:
            return []
        return list(self._project_to_supplement.get(project_lnc_id, []))

    def get_coexpression_evidence(self, project_lnc_id: str, top_n: int = 10) -> dict[str, Any]:
        if not self.available or self.sheet5_df is None:
            return {"available": False, "pairs": [], "summary": {}, "partner_gene_ids": []}

        supplement_ids = [item["supplement_id"] for item in self.get_mapping(project_lnc_id)]
        if not supplement_ids:
            return {"available": False, "pairs": [], "summary": {}, "partner_gene_ids": []}

        df = self.sheet5_df[self.sheet5_df["supplement_id"].isin(supplement_ids)].copy()
        if df.empty:
            return {"available": False, "pairs": [], "summary": {}, "partner_gene_ids": []}

        df = df.sort_values(["supplement_id", "partner_gene", "tissue", "correlation"], ascending=[True, True, True, False])
        df = df.drop_duplicates(subset=["supplement_id", "partner_gene", "tissue", "regulation", "correlation"])
        df["abs_corr"] = df["correlation"].abs()
        ranked = df.sort_values(["abs_corr", "p_value"], ascending=[False, True], na_position="last")

        top_pairs = []
        for row in ranked.head(top_n).itertuples(index=False):
            top_pairs.append(
                {
                    "supplement_id": row.supplement_id,
                    "partner_gene": row.partner_gene,
                    "regulation": row.regulation,
                    "correlation": self._safe_round(row.correlation, 4),
                    "p_value": self._to_native_number(row.p_value),
                    "state_class": row.state_class,
                    "tissue": row.tissue,
                }
            )

        summary = {
            "mapped_supplement_ids": supplement_ids,
            "pair_count": int(len(df)),
            "unique_partner_count": int(df["partner_gene"].nunique()),
            "cis_count": int((df["regulation"].str.contains("cis", case=False, na=False)).sum()),
            "trans_count": int((df["regulation"].str.contains("trans", case=False, na=False)).sum()),
            "positive_pair_count": int((df["correlation"] > 0).sum()),
            "negative_pair_count": int((df["correlation"] < 0).sum()),
            "tissues": [x for x, _ in Counter(df["tissue"].dropna()).most_common(6)],
            "state_classes": [x for x, _ in Counter(df["state_class"].dropna()).most_common(6)],
        }
        partner_gene_ids = ranked["partner_gene"].dropna().astype(str).drop_duplicates().tolist()
        return {
            "available": True,
            "pairs": top_pairs,
            "summary": summary,
            "partner_gene_ids": partner_gene_ids,
        }

    def get_gwas_evidence(self, project_lnc_id: str, top_n: int = 10) -> dict[str, Any]:
        if not self.available or self.sheet9_df is None or self.sheet10_df is None:
            return {"available": False, "hits": [], "trait_summary": [], "pve_hits": [], "summary": {}}

        supplement_ids = [item["supplement_id"] for item in self.get_mapping(project_lnc_id)]
        if not supplement_ids:
            return {"available": False, "hits": [], "trait_summary": [], "pve_hits": [], "summary": {}}

        hit_df = self.sheet9_df[self.sheet9_df["supplement_id"].isin(supplement_ids)].copy()
        pve_df = self.sheet10_df[self.sheet10_df["supplement_id"].isin(supplement_ids)].copy()
        if hit_df.empty and pve_df.empty:
            return {"available": False, "hits": [], "trait_summary": [], "pve_hits": [], "summary": {}}

        trait_summary = []
        if not hit_df.empty:
            grouped = (
                hit_df.groupby("trait", dropna=False)
                .agg(min_p_value=("p_value", "min"), snp_count=("snp", "nunique"))
                .reset_index()
                .sort_values("min_p_value", ascending=True)
            )
            for row in grouped.head(top_n).itertuples(index=False):
                trait_summary.append(
                    {
                        "trait": row.trait,
                        "min_p_value": self._to_native_number(row.min_p_value),
                        "snp_count": int(row.snp_count),
                    }
                )

        top_hits = []
        for row in hit_df.sort_values("p_value", ascending=True).head(top_n).itertuples(index=False):
            top_hits.append(
                {
                    "supplement_id": row.supplement_id,
                    "trait": row.trait,
                    "snp": row.snp,
                    "p_value": self._to_native_number(row.p_value),
                }
            )

        top_pve_hits = []
        if not pve_df.empty:
            for row in pve_df.sort_values("pve", ascending=False).head(top_n).itertuples(index=False):
                top_pve_hits.append(
                    {
                        "supplement_id": row.supplement_id,
                        "trait": row.trait,
                        "snp": row.snp,
                        "pve": self._safe_round(row.pve, 3),
                    }
                )

        summary = {
            "mapped_supplement_ids": supplement_ids,
            "direct_snp_hit_count": int(len(hit_df)),
            "direct_trait_count": int(hit_df["trait"].nunique()) if not hit_df.empty else 0,
            "pve_record_count": int(len(pve_df)),
            "best_p_value": self._to_native_number(hit_df["p_value"].min()) if not hit_df.empty else None,
            "best_pve": self._to_native_number(pve_df["pve"].max()) if not pve_df.empty else None,
        }
        return {
            "available": True,
            "hits": top_hits,
            "trait_summary": trait_summary,
            "pve_hits": top_pve_hits,
            "summary": summary,
        }

    def _clean_text(self, value: Any) -> str:
        if pd.isna(value):
            return ""
        return str(value).strip()

    def _safe_round(self, value: Any, digits: int) -> float | None:
        try:
            return round(float(value), digits)
        except Exception:
            return None

    def _to_native_number(self, value: Any) -> float | int | None:
        try:
            if pd.isna(value):
                return None
            value = float(value)
            if value.is_integer():
                return int(value)
            return value
        except Exception:
            return None
