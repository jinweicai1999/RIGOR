from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
IN_DIR = PROJECT_ROOT / "deliverables" / "gwas_trait_effect"
OUT_DIR = IN_DIR / "validation_vertical"

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

TRAIT_LABELS = {
    "Heading_date": "Heading date",
    "Plant_height": "Plant height",
    "Num_panicles": "Panicle number",
    "Num_effective_panicles": "Effective panicle number",
    "Yield": "Yield",
    "Grain_weight": "Grain weight",
    "Spikelet_length": "Spikelet length",
    "Grain_length": "Grain length",
    "Grain_width": "Grain width",
    "Grain_thickness": "Grain thickness",
}

GENES = [
    ("Os07g0622200", "DS1 / DAHPS2"),
    ("Os04g0396500", "ORL4 / LAX2"),
]

CORE_CLASSES = {"Strong", "Moderate"}


plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.7,
    }
)


def load_gene_grid(gene: str) -> pd.DataFrame:
    path = IN_DIR / f"{gene}_GWAS_trait_effect_data.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    source = pd.read_csv(path)
    lookup = source.drop_duplicates(subset=["Trait"], keep="first").set_index("Trait", drop=False)

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
            }

        is_core = str(item.get("Trait_Normalized_Evidence", "")) in CORE_CLASSES
        abs_beta = abs(float(item.get("Effect_Size(BETA)") or 0.0))
        neg_log10_p = float(item.get("neg_log10_p") or 0.0)
        item["Trait_label"] = TRAIT_LABELS.get(trait, trait.replace("_", " "))
        item["plot_abs_beta"] = abs_beta if is_core else 0.0
        item["plot_neg_log10_p"] = neg_log10_p if is_core else 0.0
        item["text_abs_beta"] = f"{abs_beta:.2g}" if is_core else "-"
        item["text_neg_log10_p"] = f"{neg_log10_p:.1f}" if is_core else "-"
        rows.append(item)

    return pd.DataFrame(rows)


def draw_vertical_panel(ax, df: pd.DataFrame, value_col: str, text_col: str, title: str, cmap: str, vmax: float):
    values = df[value_col].astype(float).to_numpy().reshape(-1, 1)
    image = ax.imshow(values, cmap=cmap, vmin=0, vmax=max(vmax, 1.0), aspect="auto")

    ax.set_xticks([0])
    ax.set_xticklabels([title], fontsize=8.5)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["Trait_label"], fontsize=7.8)
    ax.tick_params(length=0)

    for row_idx, row in df.iterrows():
        ax.text(
            0,
            row_idx,
            str(row[text_col]),
            ha="center",
            va="center",
            fontsize=7.8,
            color="#111827",
        )

    for spine in ax.spines.values():
        spine.set_visible(False)
    return image


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    grids = [(gene, label, load_gene_grid(gene)) for gene, label in GENES]
    beta_vmax = max(float(df["plot_abs_beta"].max()) for _, _, df in grids)
    p_vmax = max(float(df["plot_neg_log10_p"].max()) for _, _, df in grids)

    fig, axes = plt.subplots(
        nrows=2,
        ncols=2,
        figsize=(5.8, 8.2),
        gridspec_kw={"width_ratios": [1, 1], "wspace": 0.46, "hspace": 0.38},
    )

    beta_images = []
    p_images = []
    for row_idx, (_, display_name, df) in enumerate(grids):
        beta_images.append(
            draw_vertical_panel(
                axes[row_idx, 0],
                df,
                "plot_abs_beta",
                "text_abs_beta",
                "|BETA|",
                "Blues",
                beta_vmax,
            )
        )
        p_images.append(
            draw_vertical_panel(
                axes[row_idx, 1],
                df,
                "plot_neg_log10_p",
                "text_neg_log10_p",
                "-log10(P)",
                "YlOrRd",
                p_vmax,
            )
        )
        axes[row_idx, 0].text(
            -0.62,
            1.07,
            display_name,
            transform=axes[row_idx, 0].transAxes,
            ha="left",
            va="bottom",
            fontsize=9.5,
            fontweight="bold",
        )

    axes[0, 0].set_title("Effect size magnitude", fontsize=10, pad=10, fontweight="bold")
    axes[0, 1].set_title("Association significance", fontsize=10, pad=10, fontweight="bold")

    cbar_beta = fig.colorbar(beta_images[0], ax=axes[:, 0], fraction=0.045, pad=0.04)
    cbar_beta.ax.set_title("|BETA|", fontsize=8.5, pad=6)
    cbar_beta.ax.tick_params(labelsize=7.5, length=2)

    cbar_p = fig.colorbar(p_images[0], ax=axes[:, 1], fraction=0.045, pad=0.04)
    cbar_p.ax.set_title("-log10(P)", fontsize=8.5, pad=6)
    cbar_p.ax.tick_params(labelsize=7.5, length=2)

    stem = OUT_DIR / "DS1_ORL4_combined_vertical_GWAS_heatmaps"
    fig.savefig(stem.with_suffix(".svg"), format="svg", bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), format="pdf", bbox_inches="tight")
    plt.close(fig)

    print(stem.with_suffix(".svg"))
    print(stem.with_suffix(".pdf"))


if __name__ == "__main__":
    main()
