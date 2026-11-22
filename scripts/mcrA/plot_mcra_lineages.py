import os, json, heapq, collections
import dendropy

ROOT = r"C:\Users\nojus\OneDrive\Desktop\Bubble project\mcrA"
GRAFT = os.path.join(ROOT, "results", "epa_result.newick")
JPLACE = os.path.join(ROOT, "results", "epa_result.jplace")
NAME_MAP = os.path.join(ROOT, "database", "name_map.tsv")
OUT_FIG = os.path.join(ROOT, "results", "mcrA_lineage_placement.png")
OUT_TSV = os.path.join(ROOT, "results", "lineage_assignment.tsv")
OUT_REF = os.path.join(ROOT, "results", "reference_lineages.tsv")

PROJ2REGION = {"sed_n": "Muddy", "sed_v": "Sandy"}
REGIONS = ["Muddy", "Sandy"]
K_ANCHORS = 5
FAR = 2.0
UNRESOLVED = "Uncultured / unresolved mcrA"

AOM_RULES = [
    ("methanopereden", "ANME-2d (Ca. Methanoperedens)", "AOM"),
    ("anme-1", "ANME-1", "AOM"),
    ("anme_1", "ANME-1", "AOM"),
    ("anme1", "ANME-1", "AOM"),
    ("methanophag", "ANME-1", "AOM"),
    ("anme-2", "ANME-2a-c", "AOM"),
    ("anme_2", "ANME-2a-c", "AOM"),
    ("anme2", "ANME-2a-c", "AOM"),
    ("methanocomeden", "ANME-2a-c", "AOM"),
    ("methanomarin", "ANME-2a-c", "AOM"),
    ("methanogaster", "ANME-2a-c", "AOM"),
    ("anme-3", "ANME-3", "AOM"),
    ("anme3", "ANME-3", "AOM"),
    ("methanovoran", "ANME-3", "AOM"),
    ("anme", "ANME (unspecified)", "AOM"),
]
METH_RULES = [
    (("methanoregula", "methanolinea", "methanospirillum", "methanoculleus",
      "methanogenium", "methanofollis", "methanomicrobium", "methanoplanus",
      "methanosphaerula", "methanocorpusculum", "methanocalculus", "methanolacinia"),
     "Methanomicrobiales", "Methanogen"),
    (("methanosaeta", "methanothrix"), "Methanotrichales (Methanosaeta)", "Methanogen"),
    (("methanosarcin", "methanococcoides", "methanolobus", "methanohalophilus",
      "methanomethylovorans", "methanosalsum", "methanimicrococcus", "methanohalobium"),
     "Methanosarcinales", "Methanogen"),
    (("methanobacterium", "methanobrevibacter", "methanosphaera",
      "methanothermobacter", "methanothermus"), "Methanobacteriales", "Methanogen"),
    (("methanococcus", "methanothermococcus", "methanocaldococcus",
      "methanotorris", "methanofervidicoccus"), "Methanococcales", "Methanogen"),
    (("methanomassiliicoccus", "methanomethylophilus", "methanoplasma",
      "methanogranum"), "Methanomassiliicoccales", "Methanogen"),
    (("methanocella",), "Methanocellales", "Methanogen"),
    (("methanopyrus",), "Methanopyrales", "Methanogen"),
    (("methanonatron",), "Methanonatronarchaeia", "Methanogen"),
    (("methanofastidiosum",), "Methanofastidiosales", "Methanogen"),
]

def lineage_of_name(clean_id, orig):
    s = (str(clean_id) + " " + str(orig)).lower()
    for key, lin, guild in AOM_RULES:
        if key in s:
            return lin, guild
    for keys, lin, guild in METH_RULES:
        if any(k in s for k in keys):
            return lin, guild
    if "methano" in s or "methermicoccus" in s:
        return "Other methanogen", "Methanogen"
    return None, None

def load_name_map(path):
    m = {}
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            next(fh, None)
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2:
                    m[p[0]] = p[1]
    return m

def best_lwr_map(path):
    with open(path) as fh:
        jp = json.load(fh)
    fi = {f: i for i, f in enumerate(jp["fields"])}
    out = {}
    for pl in jp["placements"]:
        name = (pl.get("n") or pl.get("nm"))[0]
        if isinstance(name, list):
            name = name[0]
        best = max(pl["p"], key=lambda x: x[fi["like_weight_ratio"]])
        out[name] = round(float(best[fi["like_weight_ratio"]]), 3)
    return out

def build_adjacency(tree):
    adj = collections.defaultdict(list)
    for node in tree.preorder_node_iter():
        parent = node.parent_node
        if parent is not None:
            w = node.edge.length or 0.0
            w = w if w > 0 else 0.0
            adj[node].append((parent, w))
            adj[parent].append((node, w))
    return adj

def nearest_anchors(start, adj, lineage_of_node, k):
    dist = {start: 0.0}
    pq = [(0.0, 0, start)]
    ctr = 1
    found = []
    while pq and len(found) < k:
        d, _, node = heapq.heappop(pq)
        if d > dist.get(node, float("inf")):
            continue
        info = lineage_of_node.get(node)
        if info is not None:
            found.append((d, info[0], info[1], node.taxon.label if node.taxon else ""))
        for nb, w in adj[node]:
            nd = d + w
            if nd < dist.get(nb, float("inf")):
                dist[nb] = nd
                heapq.heappush(pq, (nd, ctr, nb))
                ctr += 1
    return found

def label_all_by_nearest_anchor(adj, lineage_of_node):

    dist, lab, pq, ctr = {}, {}, [], 0
    for node, info in lineage_of_node.items():
        dist[node] = 0.0
        lab[node] = info
        heapq.heappush(pq, (0.0, ctr, node))
        ctr += 1
    while pq:
        d, _, node = heapq.heappop(pq)
        if d > dist.get(node, float("inf")):
            continue
        base = lab[node]
        for nb, w in adj[node]:
            nd = d + w
            if nd < dist.get(nb, float("inf")):
                dist[nb] = nd
                lab[nb] = base
                heapq.heappush(pq, (nd, ctr, nb))
                ctr += 1
    return lab

def main():
    nmap = load_name_map(NAME_MAP)
    lwr = best_lwr_map(JPLACE)
    tree = dendropy.Tree.get(path=GRAFT, schema="newick", preserve_underscores=True)
    adj = build_adjacency(tree)

    lineage_of_node = {}
    ref_lineage_rows = []
    lineage_ref_count = collections.Counter()
    lineage_guild = {}
    lineage_rep = {}
    query_leaves = []
    for leaf in tree.leaf_node_iter():
        label = leaf.taxon.label if leaf.taxon else None
        if label is None:
            continue
        if label.startswith("mcrA|"):
            query_leaves.append(leaf)
            continue
        lin, guild = lineage_of_name(label, nmap.get(label, ""))
        if lin is not None:
            lineage_of_node[leaf] = (lin, guild)
            lineage_ref_count[lin] += 1
            lineage_guild[lin] = guild
            ref_lineage_rows.append((label, lin, guild))

            if lin not in lineage_rep or ("uncult" not in label.lower()
                                          and "uncult" in lineage_rep[lin].lower()):
                lineage_rep[lin] = label

    print(f"reference anchors (named): {sum(lineage_ref_count.values())} in {len(lineage_ref_count)} lineages")

    node_lineage = label_all_by_nearest_anchor(adj, lineage_of_node)
    lineage_full_count = collections.Counter()
    for leaf in tree.leaf_node_iter():
        lab = leaf.taxon.label if leaf.taxon else None
        if not lab or lab.startswith("mcrA|"):
            continue
        info = node_lineage.get(leaf)
        if info:
            lineage_full_count[info[0]] += 1
    print(f"full reference clades: {sum(lineage_full_count.values())} refs "
          f"across {len(lineage_full_count)} lineages")
    for lin, n in lineage_full_count.most_common():
        print(f"   {lin:38s} {n:5d}  (named anchors: {lineage_ref_count[lin]})  [{lineage_guild[lin]}]")

    q_rows = []
    per_region_conf = collections.defaultdict(collections.Counter)
    per_region_rel = collections.defaultdict(collections.Counter)
    for leaf in query_leaves:
        name = leaf.taxon.label
        anchors = nearest_anchors(leaf, adj, lineage_of_node, K_ANCHORS)
        if not anchors:
            continue
        dist = anchors[0][0]
        lin = anchors[0][1]
        guild = anchors[0][2]
        agree = sum(1 for a in anchors if a[1] == lin)
        conf = round(agree / len(anchors), 2)
        near = nmap.get(anchors[0][3], anchors[0][3]).split("  ")[0]
        if dist <= 0.6:
            rank, related = "genus/order", False
        elif dist <= FAR:
            rank, related = "order", False
        else:
            rank, related = "class-level (related)", True
        proj = name.split("|")[1] if "|" in name else "?"
        region = PROJ2REGION.get(proj, proj)
        (per_region_rel if related else per_region_conf)[region][lin] += 1
        q_rows.append((name, region, lin, guild, "related" if related else "confident",
                       rank, conf, round(dist, 3), lwr.get(name, ""), near))

    print("closest-lineage assignment per region (confident + related):")
    regions = sorted(set(list(per_region_conf) + list(per_region_rel)))
    for reg in regions:
        lins = set(per_region_conf[reg]) | set(per_region_rel[reg])
        parts = [f"{l}:{per_region_conf[reg][l]}+{per_region_rel[reg][l]}r" for l in lins]
        print(f"   {reg}: " + "; ".join(parts))

    q_rows.sort(key=lambda r: (r[1], r[2]))
    with open(OUT_TSV, "w", encoding="utf-8") as of:
        of.write("query\tregion\tclosest_lineage\tguild\tsupport\trank\t"
                 "top5_agreement\tdist_to_anchor\tbest_LWR\tnearest_anchor\n")
        for r in q_rows:
            of.write("\t".join(map(str, r)) + "\n")
    with open(OUT_REF, "w", encoding="utf-8") as of:
        of.write("reference\tlineage\tguild\n")
        for r in sorted(ref_lineage_rows, key=lambda x: (x[1], x[0])):
            of.write("\t".join(r) + "\n")
    print(f"wrote {OUT_TSV} and {OUT_REF}")

    _draw(tree, lineage_rep, lineage_full_count, lineage_guild,
          per_region_conf, per_region_rel)

def _draw(tree, lineage_rep, lineage_ref_count, lineage_guild,
          per_region_conf, per_region_rel):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    from matplotlib.lines import Line2D
    import math

    from matplotlib.colors import to_rgba
    GUILD_COLOR = {"Methanogen": "#1a9850", "AOM": "#d73027", "NA": "#9e9e9e"}
    SMAX = 560.0

    def bubble_size(n, qmax):
        return 45 + SMAX * (n / qmax)

    def tot(reg, lin):
        return per_region_conf[reg][lin] + per_region_rel[reg][lin]

    reps = {lab: lin for lin, lab in lineage_rep.items()}
    keep = set(lineage_rep.values())
    btree = tree.clone(depth=1)
    taxa_keep = [t for t in btree.taxon_namespace if t.label in keep]
    btree.retain_taxa(taxa_keep)

    def height(n):
        if n.is_leaf():
            return 0
        return 1 + max(height(c) for c in n.child_nodes())
    total = height(btree.seed_node)

    leaves = [lf for lf in btree.leaf_node_iter()]
    for i, lf in enumerate(leaves):
        lf._y = i
    ycount = len(leaves)

    def assign(n):
        if n.is_leaf():
            n._x = total
            return n._y
        ys = [assign(c) for c in n.child_nodes()]
        n._y = sum(ys) / len(ys)
        n._x = total - height(n)
        return n._y
    assign(btree.seed_node)

    LABEL_X = total + 0.6
    NCOL_X = LABEL_X + 6.3
    base = LABEL_X + 7.6
    DX = 1.15
    GAP = 1.9
    COL_X = {
        ("Muddy", "conf"): base,
        ("Muddy", "div"): base + DX,
        ("Sandy", "conf"): base + DX + GAP,
        ("Sandy", "div"): base + 2 * DX + GAP,
    }
    region_center = {
        "Muddy": (COL_X[("Muddy", "conf")] + COL_X[("Muddy", "div")]) / 2,
        "Sandy": (COL_X[("Sandy", "conf")] + COL_X[("Sandy", "div")]) / 2,
    }
    right_x = COL_X[("Sandy", "div")]
    fig_h = max(5.0, 0.46 * (ycount + 1.6) + 1.4)
    fig, ax = plt.subplots(figsize=(12.5, fig_h))

    for n in btree.preorder_node_iter():
        if n.parent_node is not None:
            ax.plot([n.parent_node._x, n._x], [n._y, n._y], color="#444", lw=1.0, zorder=1)
        if not n.is_leaf():
            cys = [c._y for c in n.child_nodes()]
            ax.plot([n._x, n._x], [min(cys), max(cys)], color="#444", lw=1.0, zorder=1)

    cmax = max([lineage_ref_count[reps[lf.taxon.label]] for lf in leaves] + [1])

    qmax = max([per_region_conf[r][reps[lf.taxon.label]] for r in REGIONS for lf in leaves]
               + [per_region_rel[r][reps[lf.taxon.label]] for r in REGIONS for lf in leaves]
               + [1])

    for lf in leaves:
        lin = reps[lf.taxon.label]
        guild = lineage_guild[lin]
        col = GUILD_COLOR.get(guild, "#9e9e9e")

        nref = lineage_ref_count[lin]
        hh = 0.16 + 0.34 * (math.log10(nref + 1) / math.log10(cmax + 1))
        ax.add_patch(Polygon([(lf._x, lf._y), (lf._x + 0.55, lf._y - hh),
                              (lf._x + 0.55, lf._y + hh)], closed=True,
                             facecolor=col, edgecolor="#333", lw=0.6, zorder=2))
        ax.text(LABEL_X, lf._y, lin, va="center", ha="left", fontsize=9.5)
        ax.text(NCOL_X, lf._y, f"{nref:,}", va="center", ha="right", fontsize=8, color="#666")

        for r in REGIONS:
            for kind, table in (("conf", per_region_conf), ("div", per_region_rel)):
                v = table[r][lin]
                if v <= 0:
                    continue
                x = COL_X[(r, kind)]
                if kind == "conf":
                    ax.scatter([x], [lf._y], s=bubble_size(v, qmax), facecolor=col,
                               alpha=0.95, edgecolors="#222", linewidths=0.5, zorder=4)
                    tcol = "white"
                else:
                    ax.scatter([x], [lf._y], s=bubble_size(v, qmax),
                               facecolor=to_rgba(col, 0.30), edgecolors=col,
                               linewidths=1.4, linestyle="--", zorder=3)
                    tcol = "#222"
                ax.text(x, lf._y, str(v), va="center", ha="center", fontsize=8,
                        color=tcol, zorder=5, fontweight="bold")

    top = ycount - 0.1
    ax.text(LABEL_X, top + 0.5, "mcrA lineage (collapsed clade)", fontsize=10.5,
            fontweight="bold", va="bottom")
    ax.text(NCOL_X, max(lf._y for lf in leaves) + 0.42, "Reference\nsequences",
            fontsize=7.5, color="#555", ha="center", va="bottom", linespacing=0.9)
    for r in REGIONS:
        ax.text(region_center[r], top + 0.5, r, fontsize=10.5, fontweight="bold",
                ha="center", va="bottom")

    def sw(color, solid, label):
        if solid:
            return Line2D([0], [0], marker='o', color='w', markerfacecolor=color,
                          markeredgecolor="#222", markersize=11, label=label)
        return Line2D([0], [0], marker='o', color='w', markerfacecolor=to_rgba(color, 0.30),
                      markeredgecolor=color, markeredgewidth=1.6, markersize=13, label=label)
    gc = GUILD_COLOR
    combined = [
        sw(gc["Methanogen"], True,  "Methanogen — confident"),
        sw(gc["Methanogen"], False, "Methanogen — closest clade"),
        sw(gc["AOM"], True,  "ANME — confident"),
        sw(gc["AOM"], False, "ANME — closest clade"),
    ]
    ax.set_xlim(-0.5, right_x + 1.4)
    ax.set_ylim(-2.2, ycount + 2.0)
    ax.legend(handles=combined, loc="upper center", frameon=False,
              fontsize=9, ncol=2, columnspacing=1.8, handletextpad=0.6,
              bbox_to_anchor=(0.5, -0.82),
              bbox_transform=matplotlib.transforms.blended_transform_factory(
                  ax.transAxes, ax.transData))
    ax.axis("off")
    plt.subplots_adjust(bottom=0.06, top=0.94)
    saved = None
    for path in (OUT_FIG, OUT_FIG[:-4] + ".pdf"):
        try:
            plt.savefig(path, dpi=220, bbox_inches="tight", pad_inches=0.35)
            saved = saved or path
        except Exception as e:
            print("save failed:", e)
    print(f"wrote {saved}")
    try:
        from shutil import copy2
        supp = r"C:\Users\nojus\OneDrive\Desktop\Bubble project\final_results\Supplementary_images"
        copy2(OUT_FIG, os.path.join(supp, "mcrA_lineage_placement.png"))
        copy2(OUT_FIG[:-4] + ".pdf", os.path.join(supp, "mcrA_lineage_placement.pdf"))
    except Exception as e:
        print("supp copy failed:", e)

if __name__ == "__main__":
    main()
