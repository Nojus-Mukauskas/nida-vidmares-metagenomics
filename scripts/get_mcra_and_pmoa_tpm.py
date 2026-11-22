from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(r"C:\Users\nojus\OneDrive\Desktop\Bubble project")
PACK = ROOT / "DATA" / "squeezemeta" / "FINAL_PCoA_analysis_pack"
OUT = ROOT / "seasonal_analysis"
FIGDIR = OUT / "results" / "figures"
TABDIR = OUT / "results" / "tables"

SEASON_ORDER = ["spring", "summer", "autumn"]
SEASON_POS = {s: i for i, s in enumerate(SEASON_ORDER)}
REGION_LABEL = {"Nida": "Muddy", "Vidmares": "Sandy"}
NIDA = "#2c5282"
VIDM = "#dd6b20"
NIDA_SURF = "#63b3ed"
NIDA_BOT = "#2c5282"
SEASON_COLOR = {"spring": "#43a047", "summer": "#fb8c00", "autumn": "#1e88e5"}

PROJECT_SAMPLES = {
    "sed_n": ["PN1", "PN2", "PN3", "VN1", "VN2", "VN3", "RN1", "RN2", "RN3"],
    "sed_v": ["PV1", "PV2", "PV3", "VV1", "VV2", "VV3", "RV1", "RV2", "RV3"],
    "water_autumn_N_V": [
        "RNS1", "RNS2", "RNS3", "RNB1", "RNB2", "RNB3",
        "RVF1", "RVF2", "RVF3", "PNB1", "PNB2", "PNB3", "PNS2",
    ],
    "water_spring_vidmares": ["PVF1", "PVF2", "PVF3"],
    "water_summer": [
        "VVF1", "VVF2", "VVF3", "VNS1", "VNS2", "VNS3", "VNB1", "VNB2", "VNB3",
    ],
}

CURATED_BLOCKS = [
    ("Methyl-CoM reductase\n(methanogenesis)",
     ["mcrA (K00399)", "mcrB (K00401)", "mcrG (K00402)"]),
    ("Acetoclastic\nmethanogenesis",
     ["cdhC (K00193)", "cdhD (K00194)", "cdhE (K00197)"]),
    ("Hydrogenotrophic\nmethanogenesis",
     ["mtd (K00319)", "mtrH (K00584)"]),
    ("Methylotrophic\nmethanogenesis",
     ["mtaA (K14080)", "mtbB (K16178)"]),
    ("Aerobic\nmethanotrophy",
     ["pmoA (phylogenetic)"]),
]


def load_pmoa_guilds(path: Path) -> dict[str, str]:
    guilds = {}
    text = path.read_text(encoding="utf-8")
    for chunk in text.split("pmoA_amoA|")[1:]:
        fields = []
        for line in chunk.strip().splitlines():
            fields.extend(line.split("\t"))
        if "|" in fields[0]:
            orf = fields[0].split("|")[-1]
        else:
            orf = fields[0]
        if len(fields) >= 5:
            guilds[orf] = fields[4]
    return guilds


def parse_orftable_abund(path: Path, samples: list[str]) -> pd.DataFrame:
    n = len(samples)
    rows, kos = {}, {}
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            orf = parts[0]
            kos[orf] = parts[9] if len(parts) > 9 else ""
            numeric = parts[-(4 * n + 1):-1]
            if len(numeric) != 4 * n:
                raise ValueError(
                    f"{path.name} {orf}: expected {4 * n} numeric fields, got {len(numeric)}"
                )
            rows[orf] = [float(x) for x in numeric[:n]]
    df = pd.DataFrame.from_dict(rows, orient="index", columns=samples)
    df.attrs["ko"] = kos
    return df


def habitat_panel(row: pd.Series) -> str:
    if row["matrix"] == "sediment":
        return "sediment"
    if row["region"] == "Nida":
        return "nida_water"
    return "vidmares_water"


def sample_order_from_meta(meta: pd.DataFrame, sample_cols: list[str]) -> pd.DataFrame:
    depth_rank = {"sediment": 0, "surface": 1, "bottom": 2, "filtered": 3}
    so = meta.loc[sample_cols].copy()
    so["matrix"] = pd.Categorical(so["matrix"], ["sediment", "water"])
    so["region"] = pd.Categorical(so["region"], ["Nida", "Vidmares"])
    so["season"] = pd.Categorical(so["season"], SEASON_ORDER)
    so["depth_rank"] = so["depth"].map(depth_rank)
    return so.sort_values(["matrix", "depth_rank", "region", "season"])


def phylogenetic_pmoA_tpm(ko: pd.DataFrame, sample_cols: list[str]) -> pd.Series:
    guilds = load_pmoa_guilds(ROOT / "pmoA" / "results" / "lineage_assignment.tsv")
    orf_parts, ko_of = [], {}
    for proj, samples in PROJECT_SAMPLES.items():
        part = parse_orftable_abund(PACK / "marker_orf_abund" / f"{proj}.marker_abund.tsv", samples)
        ko_of.update(part.attrs["ko"])
        orf_parts.append(part)
    orf_ab = pd.concat(orf_parts, axis=1).fillna(0.0).reindex(columns=sample_cols, fill_value=0.0)

    k10944_orfs = [o for o, k in ko_of.items() if "K10944" in k]
    unplaced = sorted(set(k10944_orfs) - set(guilds))

    def sum_orfs(ids: list[str]) -> pd.Series:
        ids = [o for o in ids if o in orf_ab.index]
        if not ids:
            return pd.Series(0.0, index=sample_cols)
        return orf_ab.loc[ids, sample_cols].sum(axis=0)

    raw = {
        "pmoA": sum_orfs([o for o, g in guilds.items() if g == "pmoA"]),
        "amoA": sum_orfs([o for o, g in guilds.items() if g == "amoA"]),
        "alkaneMO": sum_orfs([o for o, g in guilds.items() if g == "alkaneMO"]),
        "unplaced": sum_orfs(unplaced),
    }
    raw_total = sum(raw.values())
    ko_k10944 = ko.loc["pmoA-amoA (K10944)", sample_cols].astype(float)
    scale = np.where(raw_total.to_numpy() > 0, ko_k10944.to_numpy() / raw_total.to_numpy(), 0.0)
    return raw["pmoA"] * pd.Series(scale, index=sample_cols)


def load() -> dict:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    TABDIR.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(PACK / "marker_analysis" / "sample_metadata.csv")
    meta = meta.set_index("sample")
    meta["region_label"] = meta["region"].map(REGION_LABEL)
    meta["panel"] = meta.apply(habitat_panel, axis=1)
    meta["season"] = pd.Categorical(meta["season"], SEASON_ORDER, ordered=True)

    ko = pd.read_csv(PACK / "marker_analysis" / "markers_TPM_by_sample.csv", index_col=0)
    sample_cols = [c for c in ko.columns if c != "process"]
    pmoA_tpm = phylogenetic_pmoA_tpm(ko, sample_cols)

    mcra = ko.loc["mcrA (K00399)", sample_cols].astype(float)
    indicators = pd.DataFrame({"mcrA": mcra, "pmoA": pmoA_tpm})
    indicators.index.name = "sample"
    long = (
        indicators.reset_index()
        .melt(id_vars="sample", var_name="gene", value_name="TPM")
        .merge(meta.reset_index(), on="sample")
    )
    return {
        "ko": ko,
        "meta": meta,
        "sample_cols": sample_cols,
        "pmoA_tpm": pmoA_tpm,
        "long": long,
    }
