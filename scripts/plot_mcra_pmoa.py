from pathlib import Path
from shutil import copy2
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from get_mcra_and_pmoa_tpm import (
    FIGDIR,
    TABDIR,
    ROOT,
    NIDA,
    VIDM,
    NIDA_SURF,
    NIDA_BOT,
    SEASON_ORDER,
    SEASON_POS,
    load,
)

PAPER = ROOT / "final_results" / "PAPER"
HUE_OFF = {1: [0.0], 2: [-0.18, 0.18], 3: [-0.24, 0.0, 0.24]}
HUE_W = {1: 0.38, 2: 0.28, 3: 0.20}
WATER_HUE = ["Muddy surface", "Muddy bottom", "Sandy"]
WATER_X = {
    "Muddy surface": HUE_OFF[3][0],
    "Muddy bottom": HUE_OFF[3][1],
    "Sandy": HUE_OFF[3][2],
}
WATER_PAIRS = [
    ("Muddy surface", "Sandy"),
    ("Muddy surface", "Muddy bottom"),
]


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


def boxes(ax, sub, hue_col, hue_order, hue_colors, rng, hatches=None):
    n_h = len(hue_order)
    offsets = HUE_OFF[n_h]
    width = HUE_W[n_h]
    for season in SEASON_ORDER:
        x0 = SEASON_POS[season]
        for i, hue in enumerate(hue_order):
            vals = sub[(sub["season"] == season) & (sub[hue_col] == hue)]["TPM"].to_numpy()
            if len(vals) == 0:
                continue
            pos = x0 + offsets[i]
            color = hue_colors[hue]
            hatch = None if hatches is None else hatches.get(hue)
            if len(vals) >= 2:
                bp = ax.boxplot(
                    [vals], positions=[pos], widths=width,
                    patch_artist=True, showfliers=False, zorder=2,
                    manage_ticks=False,
                    medianprops={"color": "#1a202c", "lw": 1.2},
                    whiskerprops={"color": "#4a5568", "lw": 0.9},
                    capprops={"color": "#4a5568", "lw": 0.9},
                    boxprops={"lw": 0.8, "color": "#2d3748"},
                )
                bp["boxes"][0].set_facecolor(color)
                bp["boxes"][0].set_alpha(0.45)
                if hatch:
                    bp["boxes"][0].set_hatch(hatch)
            jitter = rng.uniform(-0.05, 0.05, size=len(vals))
            ax.scatter(
                np.full(len(vals), pos) + jitter, vals,
                s=28, c=color, edgecolors="white", linewidths=0.5,
                zorder=3, alpha=0.95,
            )


def water_series(row):
    if row["region"] == "Nida":
        return "Muddy surface" if row["depth"] == "surface" else "Muddy bottom"
    return "Sandy"


def stars(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "n=1"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def holm(p_by_key):
    items = [(k, p) for k, p in p_by_key.items() if p is not None]
    out = {k: None for k in p_by_key}
    m = len(items)
    if m == 0:
        return out
    running = 0.0
    for i, (k, p) in enumerate(sorted(items, key=lambda t: t[1])):
        running = max(running, min(1.0, p * (m - i)))
        out[k] = running
    return out


def welch(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 2 or len(y) < 2:
        return None
    t, p = stats.ttest_ind(x, y, equal_var=False)
    return float(t), float(p), int(len(x)), int(len(y)), float(x.mean()), float(y.mean())


def sediment_stats(sub, gene):
    rows = []
    raw = {}
    for season in SEASON_ORDER:
        mud = sub.loc[
            (sub["season"] == season) & (sub["region_label"] == "Muddy"), "TPM"
        ].to_numpy()
        san = sub.loc[
            (sub["season"] == season) & (sub["region_label"] == "Sandy"), "TPM"
        ].to_numpy()
        got = welch(mud, san)
        if got is None:
            raw[season] = None
            continue
        t, p, n1, n2, m1, m2 = got
        raw[season] = p
        rows.append({
            "compartment": "sediment",
            "gene": gene,
            "test": "Welch t-test",
            "contrast": f"Muddy vs Sandy ({season})",
            "n_muddy": n1,
            "n_sandy": n2,
            "mean_muddy": m1,
            "mean_sandy": m2,
            "statistic": t,
            "p": p,
        })
    adj = holm(raw)
    for row in rows:
        season = row["contrast"].split("(")[1].rstrip(")")
        row["p_holm"] = adj[season]
        row["stars"] = stars(adj[season])
    return pd.DataFrame(rows), adj


def water_stats(wat, gene="pmoA"):
    rows = []
    raw_by_pair = {pair: {} for pair in WATER_PAIRS}
    for a, b in WATER_PAIRS:
        for season in SEASON_ORDER:
            got = welch(
                wat.loc[(wat["season"] == season) & (wat["series"] == a), "TPM"].to_numpy(),
                wat.loc[(wat["season"] == season) & (wat["series"] == b), "TPM"].to_numpy(),
            )
            if got is None:
                raw_by_pair[(a, b)][season] = None
                rows.append({
                    "compartment": "water",
                    "gene": gene,
                    "test": "Welch t-test",
                    "contrast": f"{a} vs {b} ({season})",
                    "n_muddy": 1,
                    "n_sandy": 0,
                    "mean_muddy": np.nan,
                    "mean_sandy": np.nan,
                    "statistic": np.nan,
                    "p": np.nan,
                    "p_holm": np.nan,
                    "stars": "n=1",
                })
                continue
            t, p, n1, n2, m1, m2 = got
            raw_by_pair[(a, b)][season] = p
            rows.append({
                "compartment": "water",
                "gene": gene,
                "test": "Welch t-test",
                "contrast": f"{a} vs {b} ({season})",
                "n_muddy": n1,
                "n_sandy": n2,
                "mean_muddy": m1,
                "mean_sandy": m2,
                "statistic": t,
                "p": p,
            })
    adj_by_pair = {pair: holm(raw) for pair, raw in raw_by_pair.items()}
    within = {s: [] for s in SEASON_ORDER}
    for row in rows:
        if row.get("stars") == "n=1":
            continue
        a, rest = row["contrast"].split(" vs ")
        b, season = rest.split(" (")
        season = season.rstrip(")")
        padj = adj_by_pair[(a, b)][season]
        row["p_holm"] = padj
        row["stars"] = stars(padj)
        if padj is not None and padj < 0.05:
            within[season].append((a, b, padj))
    return pd.DataFrame(rows), within


def sig_bracket(ax, x1, x2, y, text, h):
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color="#2d3748", lw=0.8, clip_on=False)
    ax.text((x1 + x2) / 2, y + h, text, ha="center", va="bottom", fontsize=9, color="#1a202c")


def mark_seasons(ax, sub, x_left, x_right, adj_map, ymax):
    pad = ymax * 0.045
    h = ymax * 0.025
    used = 0.0
    for season in SEASON_ORDER:
        padj = adj_map.get(season)
        if padj is None or padj >= 0.05:
            continue
        x0 = SEASON_POS[season]
        y = sub.loc[sub["season"] == season, "TPM"].max() + pad
        sig_bracket(ax, x0 + x_left, x0 + x_right, y, stars(padj), h)
        used = max(used, y + h + ymax * 0.06)
    return used


def mark_water_pairs(ax, sub, pairs_by_season, ymax):
    pad = ymax * 0.04
    h = ymax * 0.022
    step = ymax * 0.08
    used = 0.0
    for season in SEASON_ORDER:
        pairs = sorted(
            pairs_by_season.get(season, []),
            key=lambda t: abs(WATER_X[t[1]] - WATER_X[t[0]]),
        )
        if not pairs:
            continue
        y = float(sub.loc[sub["season"] == season, "TPM"].max()) + pad
        x0 = SEASON_POS[season]
        for a, b, p in pairs:
            sig_bracket(ax, x0 + WATER_X[a], x0 + WATER_X[b], y, stars(p), h)
            y += step
        used = max(used, y + ymax * 0.04)
    return used


def draw_sediment_gene(ax, sub, gene, rng, p_map, letter=None):
    ymax = max(sub["TPM"].max() * 1.22, 0.05)
    decorate(ax, f"{gene} TPM", letter)
    boxes(ax, sub, "region", ["Nida", "Vidmares"], {"Nida": NIDA, "Vidmares": VIDM}, rng)
    needed = mark_seasons(ax, sub, HUE_OFF[2][0], HUE_OFF[2][1], p_map, ymax)
    ax.set_ylim(0, max(ymax, needed) * 1.08)


def decorate(ax, ylabel, letter=None):
    ax.set_xlim(-0.55, 2.55)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["Spring", "Summer", "Autumn"], fontsize=9)
    ax.yaxis.grid(True, ls=":", lw=0.5, color="#e2e8f0", zorder=0)
    ax.set_axisbelow(True)
    ax.set_ylabel(ylabel, fontsize=9)
    if letter:
        ax.text(-0.08, 1.01, letter, transform=ax.transAxes, fontsize=16,
                fontweight="bold", ha="left", va="bottom")


def save_named(fig, stem):
    PAPER.mkdir(parents=True, exist_ok=True)
    out_png = FIGDIR / f"{stem}.png"
    out_pdf = FIGDIR / f"{stem}.pdf"
    fig.savefig(out_png, dpi=300, bbox_inches="tight", pad_inches=0.08)
    fig.savefig(out_pdf, bbox_inches="tight", pad_inches=0.08)
    for dest in (PAPER, ROOT / "final_results"):
        dest.mkdir(parents=True, exist_ok=True)
        try:
            copy2(out_png, dest / f"{stem}.png")
            copy2(out_pdf, dest / f"{stem}.pdf")
        except PermissionError:
            copy2(out_png, dest / f"{stem}.png")
    print("Saved", out_png)


def main():
    data = load()
    long = data["long"]
    long.to_csv(TABDIR / "indicator_TPM_long.csv", index=False)
    (
        long.groupby(["gene", "panel", "region_label", "depth", "season"], observed=True)["TPM"]
        .agg(mean="mean", std="std", n="count")
        .reset_index()
        .to_csv(TABDIR / "indicator_TPM_means.csv", index=False)
    )

    style()
    rng = np.random.default_rng(1)
    site_leg = [
        Patch(facecolor=NIDA, alpha=0.55, edgecolor="#2d3748", label="Muddy"),
        Patch(facecolor=VIDM, alpha=0.55, edgecolor="#2d3748", label="Sandy"),
    ]

    sed = long[(long["panel"] == "sediment") & (long["gene"].isin(["mcrA", "pmoA"]))]
    wat = long[(long["gene"] == "pmoA") & (long["panel"] != "sediment")].copy()
    wat["series"] = wat.apply(water_series, axis=1)

    stat_frames = []
    adj = {}
    for gene, sub, key in (
        ("mcrA", sed[sed["gene"] == "mcrA"], "sediment_mcrA"),
        ("pmoA", sed[sed["gene"] == "pmoA"], "sediment_pmoA"),
    ):
        tab, adj_map = sediment_stats(sub, gene)
        stat_frames.append(tab)
        adj[key] = adj_map
    wtab, w_within = water_stats(wat)
    stat_frames.append(wtab)

    stats_tab = pd.concat(stat_frames, ignore_index=True)
    stats_path = TABDIR / "indicator_site_stats.csv"
    stats_tab.to_csv(stats_path, index=False)
    try:
        copy2(stats_path, ROOT / "final_results" / "indicator_site_stats.csv")
    except PermissionError:
        pass

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 5.0), sharex=True)
    fig.subplots_adjust(left=0.10, right=0.98, top=0.90, bottom=0.10, wspace=0.28)
    for ax, gene, letter, key in (
        (axes[0], "mcrA", "A", "sediment_mcrA"),
        (axes[1], "pmoA", "B", "sediment_pmoA"),
    ):
        draw_sediment_gene(ax, sed[sed["gene"] == gene], gene, rng, adj[key], letter)
    axes[1].legend(handles=site_leg, loc="upper right", frameon=False, fontsize=8,
                   handlelength=1.2, handletextpad=0.4, borderaxespad=0.3)
    save_named(fig, "fig1_indicator_genes")
    plt.close(fig)

    fig_w, ax_w = plt.subplots(figsize=(6.6, 5.1))
    fig_w.subplots_adjust(left=0.14, right=0.97, top=0.92, bottom=0.12)
    ymax_w = max(wat["TPM"].max() * 1.22, 0.05)
    decorate(ax_w, "pmoA TPM")
    boxes(
        ax_w, wat, "series",
        ["Muddy surface", "Muddy bottom", "Sandy"],
        {"Muddy surface": NIDA_SURF, "Muddy bottom": NIDA_BOT, "Sandy": VIDM},
        rng,
    )
    needed = mark_water_pairs(ax_w, wat, w_within, ymax_w)
    ax_w.set_ylim(0, max(ymax_w, needed) * 1.08)
    ax_w.legend(
        handles=[
            Patch(facecolor=NIDA_SURF, alpha=0.55, edgecolor="#2d3748", label="Muddy surface"),
            Patch(facecolor=NIDA_BOT, alpha=0.55, edgecolor="#2d3748", label="Muddy bottom"),
            Patch(facecolor=VIDM, alpha=0.55, edgecolor="#2d3748", label="Sandy"),
        ],
        loc="upper left", bbox_to_anchor=(0.06, 0.98),
        frameon=False, fontsize=8,
        handlelength=1.2, handletextpad=0.4, borderaxespad=0.0,
    )
    save_named(fig_w, "fig1_indicator_genes_water")
    plt.close(fig_w)
    print(stats_tab.to_string(index=False))


if __name__ == "__main__":
    main()
