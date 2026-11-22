from pathlib import Path
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

sys.path.insert(0, str(Path(__file__).resolve().parent))
from get_mcra_and_pmoa_tpm import (
    FIGDIR,
    TABDIR,
    NIDA,
    VIDM,
    SEASON_COLOR,
    CURATED_BLOCKS,
    load,
    sample_order_from_meta,
)


def style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    })


def column_guides(ax, so):
    n_samp = len(so)
    n_sed = int((so["matrix"] == "sediment").sum())
    ax.axvline(n_sed - 0.5, color="black", lw=1.6, zorder=4)
    for i in range(1, n_samp):
        prev, cur = so.iloc[i - 1], so.iloc[i]
        if cur["matrix"] != prev["matrix"]:
            continue
        if cur["region"] != prev["region"]:
            ax.axvline(i - 0.5, color="#4a5568", lw=0.8, ls="--", zorder=3)
        elif cur["depth"] != prev["depth"]:
            ax.axvline(i - 0.5, color="#718096", lw=0.9, zorder=3)
        elif cur["season"] != prev["season"]:
            ax.axvline(i - 0.5, color="#e2e8f0", lw=0.7, zorder=2)


def main():
    data = load()
    ko = data["ko"]
    meta = data["meta"]
    pmoA_tpm = data["pmoA_tpm"]
    sample_cols = data["sample_cols"]

    style()
    so = sample_order_from_meta(meta, sample_cols)
    sample_order = list(so.index)

    mat = ko.drop(columns=["process"]).copy()
    mat.loc["pmoA (phylogenetic)"] = pmoA_tpm.reindex(mat.columns)
    if "pmoA-amoA (K10944)" in mat.index:
        mat = mat.drop(index=["pmoA-amoA (K10944)"])

    genes = [g for _, block in CURATED_BLOCKS for g in block]
    missing = [g for g in genes if g not in mat.index]
    if missing:
        raise KeyError(f"curated genes missing from KO table: {missing}")

    M = mat.loc[genes, sample_order].astype(float)
    L = np.log10(M + 0.1)
    M.to_csv(TABDIR / "curated_heatmap_TPM.csv")
    L.to_csv(TABDIR / "curated_heatmap_log10TPM.csv")

    cmap = LinearSegmentedColormap.from_list(
        "wg", ["#f7f7f7", "#fee08b", "#fc8d59", "#d73027", "#7a0177"]
    )
    vmin, vmax = np.log10(0.1), np.log10(max(80.0, float(M.to_numpy().max()) * 1.05))

    heights = [max(len(block), 2.5) for _, block in CURATED_BLOCKS]
    n_samp = len(sample_order)
    n_sed = int((so["matrix"] == "sediment").sum())
    nida_sed = int(((so["matrix"] == "sediment") & (so["region"] == "Nida")).sum())
    nida_wat = int(((so["matrix"] == "water") & (so["region"] == "Nida")).sum())
    n_surf = int(((so["matrix"] == "water") & (so["depth"] == "surface")).sum())
    n_bot = int(((so["matrix"] == "water") & (so["depth"] == "bottom")).sum())

    fig = plt.figure(figsize=(15.4, 9.6))
    outer = GridSpec(
        2, 2, figure=fig,
        width_ratios=[1, 0.035],
        height_ratios=[0.85, 10],
        hspace=0.02, wspace=0.05,
        left=0.18, right=0.92, top=0.93, bottom=0.10,
    )
    gs_heat = GridSpecFromSubplotSpec(
        len(CURATED_BLOCKS), 1,
        subplot_spec=outer[1, 0],
        height_ratios=heights,
        hspace=0.18,
    )

    ax_head = fig.add_subplot(outer[0, 0])
    ax_head.set_xlim(-0.5, n_samp - 0.5)
    ax_head.set_ylim(0, 1)
    ax_head.axis("off")

    def span(x0, x1, y, text, **kw):
        ax_head.text((x0 + x1) / 2 - 0.5, y, text, ha="center", va="center",
                     clip_on=False, **kw)

    span(0, n_sed, 0.82, "SEDIMENT", fontsize=16, fontweight="bold", color="black")
    span(n_sed, n_samp, 0.82, "WATER", fontsize=16, fontweight="bold", color="black")
    span(0, nida_sed, 0.48, "Muddy", fontsize=12, fontweight="bold", color=NIDA)
    span(nida_sed, n_sed, 0.48, "Sandy", fontsize=12, fontweight="bold", color=VIDM)
    span(n_sed, n_sed + nida_wat, 0.48, "Muddy", fontsize=12, fontweight="bold", color=NIDA)
    span(n_sed + nida_wat, n_samp, 0.48, "Sandy", fontsize=12, fontweight="bold", color=VIDM)
    span(n_sed, n_sed + n_surf, 0.12, "Surface", fontsize=9, color=NIDA)
    span(n_sed + n_surf, n_sed + n_surf + n_bot, 0.12, "Bottom", fontsize=9, color=NIDA)
    span(n_sed + nida_wat, n_samp, 0.12, "Surface", fontsize=9, color=VIDM)
    ax_head.axvline(n_sed - 0.5, color="black", lw=1.6, ymin=0, ymax=1)
    ax_head.axvline(nida_sed - 0.5, color="#4a5568", lw=0.8, ls="--", ymin=0, ymax=0.70)
    ax_head.axvline(n_sed + nida_wat - 0.5, color="#4a5568", lw=0.8, ls="--", ymin=0, ymax=0.70)
    ax_head.axvline(n_sed + n_surf - 0.5, color="#718096", lw=0.8, ymin=0, ymax=0.28)

    im = None
    heat_axes = []
    for i, (title, block) in enumerate(CURATED_BLOCKS):
        ax = fig.add_subplot(gs_heat[i], sharex=ax_head)
        heat_axes.append(ax)
        sub = L.loc[block, sample_order]
        im = ax.imshow(
            sub.values, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
        )
        ax.set_yticks(range(len(block)))
        ax.set_yticklabels(
            [g.replace(" (phylogenetic)", "") for g in block], fontsize=9.5,
        )
        column_guides(ax, so)
        ax.tick_params(length=0)
        ax.tick_params(axis="y", pad=2)
        if i == 0:
            ax.set_ylim(len(block) - 0.5, -0.78)
            ax.scatter(
                np.arange(n_samp), np.full(n_samp, -0.68),
                c=[SEASON_COLOR[s] for s in so["season"]],
                s=20, marker="s", clip_on=False, linewidths=0, zorder=5,
            )
        if i < len(CURATED_BLOCKS) - 1:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xticks(range(n_samp))
            ax.set_xticklabels(sample_order, rotation=90, fontsize=7)

    for i, (ax, (title, _)) in enumerate(zip(heat_axes, CURATED_BLOCKS)):
        bb = ax.get_position()
        y = bb.y0 + (0.42 if i == 0 else 0.50) * (bb.y1 - bb.y0)
        fig.text(
            0.072, y, title,
            rotation=90, ha="center", va="center",
            fontsize=9, fontweight="bold", color="#1a202c",
            linespacing=1.15, transform=fig.transFigure,
        )

    cax = fig.add_subplot(outer[1, 1])
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("log10(TPM + 0.1)", fontsize=9)

    fig.suptitle("Methane-cycle marker genes", fontsize=13, y=0.98)

    fig.legend(
        handles=[
            Patch(facecolor=SEASON_COLOR["spring"], label="Spring"),
            Patch(facecolor=SEASON_COLOR["summer"], label="Summer"),
            Patch(facecolor=SEASON_COLOR["autumn"], label="Autumn"),
        ],
        loc="upper right", bbox_to_anchor=(0.91, 0.98),
        frameon=False, fontsize=8, ncol=3, title=None,
    )
    out = FIGDIR / "fig2_process_heatmaps.png"
    fig.savefig(out, dpi=180, bbox_inches="tight", pad_inches=0.15)
    fig.savefig(FIGDIR / "fig2_process_heatmaps.pdf", bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    print("Saved", out)


if __name__ == "__main__":
    main()
