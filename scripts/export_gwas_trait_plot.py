import argparse
import math
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42})

os.environ["APP_LANG"] = "en"

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.agents.gwas_agent import GwasAgent
from src.project_paths import ANNOTATION_DIR
from src.utils.gene_annotator import RiceGeneAnnotator


TRAIT_ORDER = [
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


def _save_single_row_heatmap(
    df: pd.DataFrame,
    value_col: str,
    value_label: str,
    cmap: str,
    title: str,
    stem: Path,
) -> tuple[Path, Path]:
    labels = [str(x).replace("_", " ") for x in df["Trait"]]
    values = df[value_col].astype(float).to_numpy().reshape(1, -1)
    text_col = "display_text_p" if value_col == "display_neg_log10_p" else "display_text_beta"
    text_values = [df[text_col].astype(str).tolist()] if text_col in df.columns else [[f"{value:.2f}" if value_col == "neg_log10_p" else f"{value:.3g}" for value in values[0]]]
    vmax = max(1.0, float(values.max()))

    fig, ax = plt.subplots(figsize=(11.8, 2.6))
    image = ax.imshow(values, aspect="auto", cmap=cmap, vmin=0, vmax=vmax)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=9.5)
    ax.set_yticks([0])
    ax.set_yticklabels([str(df["Gene"].iloc[0])], fontsize=10.5)

    for idx, text in enumerate(text_values[0]):
        ax.text(idx, 0, text, ha="center", va="center", fontsize=11.5, fontweight="bold", color="#0f172a")

    ax.set_xlabel("Trait", fontsize=11)
    ax.set_ylabel("Gene", fontsize=11)
    ax.set_title(title, fontsize=13, fontweight="bold", pad=12)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    cbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.025)
    cbar.set_label(value_label, fontsize=10)
    cbar.ax.tick_params(labelsize=8.5)

    fig.tight_layout()
    svg_path = stem.with_suffix(".svg")
    pdf_path = stem.with_suffix(".pdf")
    fig.savefig(svg_path, format="svg", bbox_inches="tight")
    fig.savefig(pdf_path, format="pdf", bbox_inches="tight")
    plt.close(fig)
    return svg_path, pdf_path


def complete_trait_grid(df: pd.DataFrame) -> pd.DataFrame:
    gene = str(df["Gene"].iloc[0]) if not df.empty and "Gene" in df.columns else "Target"
    if "GWAS_Query_Mode" in df.columns:
        df = df.drop(columns=["GWAS_Query_Mode"])
    lookup = df.drop_duplicates(subset=["Trait"], keep="first").set_index("Trait", drop=False)
    rows = []
    for trait in TRAIT_ORDER:
        if trait in lookup.index:
            item = lookup.loc[trait].to_dict()
        else:
            item = {
                "Gene": gene,
                "Trait": trait,
                "P_Value": float("nan"),
                "Effect_Size(BETA)": 0.0,
                "neg_log10_p": 0.0,
                "abs_beta": 0.0,
                "Trait_Normalized_Evidence": "Not significant",
                "Trait_Top_Percent_Max": "",
            }
        is_core = str(item.get("Trait_Normalized_Evidence", "")) in {"Strong", "Moderate"}
        item["display_abs_beta"] = float(item.get("abs_beta") or 0.0) if is_core else 0.0
        item["display_neg_log10_p"] = float(item.get("neg_log10_p") or 0.0) if is_core else 0.0
        item["display_text_beta"] = f"{item['display_abs_beta']:.3g}" if is_core else "-"
        item["display_text_p"] = f"{item['display_neg_log10_p']:.2f}" if is_core else "-"
        rows.append(item)
    return pd.DataFrame(rows)


def export_trait_plot(gene: str, out_dir: Path) -> tuple[Path, Path, Path, Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)

    annotator = RiceGeneAnnotator(
        str(ANNOTATION_DIR / "IRGSP-1.0_representative" / "transcripts.gff"),
        str(ANNOTATION_DIR / "RAP-MSU.txt"),
    )
    annotator.build_database()
    info = annotator.query(gene)
    if not info:
        raise SystemExit(f"Gene not found: {gene}")

    agent = GwasAgent()
    df = agent.scan_network_gwas([gene], annotator)
    if df.empty:
        raise SystemExit(f"No significant GWAS records found for {gene}")

    df = df.copy()
    df["neg_log10_p"] = df["P_Value"].apply(lambda x: -math.log10(max(float(x), 1e-300)))
    df["abs_beta"] = df["Effect_Size(BETA)"].astype(float).abs()
    if "GWAS_Query_Mode" in df.columns:
        df = df.drop(columns=["GWAS_Query_Mode"])
    df["Trait"] = pd.Categorical(df["Trait"], categories=TRAIT_ORDER, ordered=True)
    df = df.sort_values("Trait").reset_index(drop=True)

    data_path = out_dir / f"{gene}_GWAS_trait_effect_data.csv"
    df.to_csv(data_path, index=False)
    suggestive_path = out_dir / f"{gene}_GWAS_suggestive_weak_data.csv"
    weak_df = df[df["Trait_Normalized_Evidence"] == "Weak"] if "Trait_Normalized_Evidence" in df.columns else df.iloc[0:0]
    weak_df.to_csv(suggestive_path, index=False)
    plot_df = complete_trait_grid(df)

    beta_svg, beta_pdf = _save_single_row_heatmap(
        plot_df,
        "display_abs_beta",
        "|BETA|",
        "Blues",
        f"{gene} GWAS effect-size magnitude (|BETA|)",
        out_dir / f"{gene}_GWAS_absBETA_heatmap",
    )
    p_svg, p_pdf = _save_single_row_heatmap(
        plot_df,
        "display_neg_log10_p",
        "-log10(P)",
        "YlOrRd",
        f"{gene} GWAS association significance (-log10(P))",
        out_dir / f"{gene}_GWAS_neglog10P_heatmap",
    )

    return beta_svg, beta_pdf, p_svg, p_pdf, data_path, suggestive_path


def main():
    parser = argparse.ArgumentParser(description="Export a RIGOR-style GWAS trait effect plot.")
    parser.add_argument("--gene", required=True, help="RAP gene ID, e.g. Os07g0622200")
    parser.add_argument(
        "--out-dir",
        default=str(PROJECT_ROOT / "deliverables" / "gwas_trait_effect"),
        help="Output directory.",
    )
    args = parser.parse_args()

    for path in export_trait_plot(args.gene, Path(args.out_dir)):
        print(path)


if __name__ == "__main__":
    main()
