import os, json, heapq, collections
import dendropy

ROOT = r"C:\Users\nojus\OneDrive\Desktop\Bubble project\pmoA"
GRAFT = os.path.join(ROOT, "results", "epa_result.newick")
JPLACE = os.path.join(ROOT, "results", "epa_result.jplace")
TAX = os.path.join(ROOT, "database", "XmoA_leaf_taxonomy.tsv")
OUT_FIG = os.path.join(ROOT, "results", "pmoA_lineage_placement.png")
OUT_TSV = os.path.join(ROOT, "results", "lineage_assignment.tsv")
OUT_REF = os.path.join(ROOT, "results", "reference_lineages.tsv")

QPREFIX = "pmoA_amoA|"
PROJ2SAMPLE = {
    "sed_n": "Muddy\nsediment",
    "sed_v": "Sandy\nsediment",
    "water_autumn_N_V": "Autumn\nwater",
    "water_spring_vidmares": "Spring\nwater",
    "water_summer": "Summer\nwater",
}
SAMPLES = ["Muddy\nsediment", "Sandy\nsediment", "Autumn\nwater",
           "Spring\nwater", "Summer\nwater"]
K_ANCHORS = 5
GENUS = 0.6
FAR = 2.0

GUILD_COLOR = {"pmoA": "#1a9850", "amoA": "#d73027", "alkaneMO": "#2c7fb8", "NA": "#cccccc"}
GUILD_LONG = {"pmoA": "methanotroph (pmoA)", "amoA": "ammonia oxidiser (amoA)",
              "alkaneMO": "divergent CuMMO (alkane/alkene MO)"}

def classify(lineage):
    L = (lineage or "").lower()
    if "thaumarchaeota" in L:
        return "AOA (Thaumarchaeota) — amoA", "amoA"
    if "nitrosomonadales" in L:
        return "AOB-β (Nitrosomonadales) — amoA", "amoA"
    if "nitrosococcus" in L or "chromatiales" in L:
        return "AOB-γ (Nitrosococcus) — amoA", "amoA"

    if "methylococcales" in L or "methylothermaceae" in L:
        return "Type I MOB (Methylococcales) — pmoA", "pmoA"
    if ("beijerinckiaceae" in L or "methylocystaceae" in L
            or "hyphomicrobiales" in L or "alphaproteobacteria" in L):
        return "Type II MOB (Alphaproteobacteria) — pmoA", "pmoA"
    if "verrucomicrobia" in L or "methylacidiphil" in L:
        return "Verrucomicrobia MOB — pmoA", "pmoA"
    if "nc10" in L or "methylomirabilis" in L:
        return "NC10 (Ca. Methylomirabilis) — pmoA", "pmoA"

    if "actinobacteria" in L:
        return "Actinobacteria CuMMO (alkane MO)", "alkaneMO"
    if "cellvibrionales" in L or "haliea" in L:
        return "Haliea CuMMO (alkene MO)", "alkaneMO"
    return "Unclassified CuMMO", "NA"

def load_tip_lineage(path):
    m = {}
    with open(path, encoding="utf-8") as fh:
        next(fh, None)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if not p or not p[0]:
                continue
            m[f"{p[0]}_XmoA"] = p[3] if len(p) > 3 else ""
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
    tip_lin = load_tip_lineage(TAX)
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
        if label.startswith(QPREFIX):
            query_leaves.append(leaf)
            continue
        lin, guild = classify(tip_lin.get(label, ""))
        lineage_of_node[leaf] = (lin, guild)
        lineage_ref_count[lin] += 1
        lineage_guild[lin] = guild
        ref_lineage_rows.append((label, lin, guild))
        if lin not in lineage_rep:
            lineage_rep[lin] = label

    print(f"reference tips: {sum(lineage_ref_count.values())} in {len(lineage_ref_count)} CuMMO clades")

    node_lineage = label_all_by_nearest_anchor(adj, lineage_of_node)
    lineage_full_count = collections.Counter()
    for leaf in tree.leaf_node_iter():
        lab = leaf.taxon.label if leaf.taxon else None
        if not lab or lab.startswith(QPREFIX):
            continue
        info = node_lineage.get(leaf)
        if info:
            lineage_full_count[info[0]] += 1
    for lin, n in lineage_full_count.most_common():
        print(f"   {lin:42s} {n:4d}  [{lineage_guild[lin]}]")

    q_rows = []
    per_sample_conf = collections.defaultdict(collections.Counter)
    per_sample_rel = collections.defaultdict(collections.Counter)
    dists = []
    guild_tot = collections.Counter()
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
        near = anchors[0][3]
        dists.append(dist)
        guild_tot[guild] += 1
        if dist <= GENUS:
            rank, related = "genus/family", False
        elif dist <= FAR:
            rank, related = "order", False
        else:
            rank, related = "class-level (related)", True
        proj = name.split("|")[1] if "|" in name else "?"
        sample = PROJ2SAMPLE.get(proj, proj)
        (per_sample_rel if related else per_sample_conf)[sample][lin] += 1
        q_rows.append((name, sample.replace("\n", " "), lin, guild, "related" if related else "confident",
                       rank, conf, round(dist, 3), lwr.get(name, ""), near))

    if dists:
        ds = sorted(dists)
        q = lambda p: ds[min(len(ds) - 1, int(p * (len(ds) - 1)))]
        print(f"nearest-anchor distance: min={ds[0]:.3f} 25%={q(.25):.3f} "
              f"median={q(.5):.3f} 75%={q(.75):.3f} max={ds[-1]:.3f}")
    print("guild totals (all fragments):", dict(guild_tot))

    q_rows.sort(key=lambda r: (r[1], r[3], r[2]))
    with open(OUT_TSV, "w", encoding="utf-8") as of:
        of.write("query\tsample\tclosest_lineage\tguild\tsupport\trank\t"
                 "top5_agreement\tdist_to_anchor\tbest_LWR\tnearest_anchor\n")
        for r in q_rows:
            of.write("\t".join(map(str, r)) + "\n")
    with open(OUT_REF, "w", encoding="utf-8") as of:
        of.write("reference\tlineage\tguild\n")
        for r in sorted(ref_lineage_rows, key=lambda x: (x[1], x[0])):
            of.write("\t".join(r) + "\n")
    print(f"wrote {OUT_TSV} and {OUT_REF}")

    _draw(tree, lineage_rep, lineage_full_count, lineage_guild,
          per_sample_conf, per_sample_rel, len(q_rows))

def _draw(tree, lineage_rep, lineage_full_count, lineage_guild,
          per_sample_conf, per_sample_rel, nquery):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    from matplotlib.lines import Line2D
    from matplotlib.colors import to_rgba
    import math

    SMAX = 560.0

    def bubble_size(n, qmax):
        return 45 + SMAX * (n / qmax)

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

    has_div = any(per_sample_rel[s][lin] > 0 for s in SAMPLES for lin in per_sample_rel[s])

    LABEL_X = total + 0.6
    NCOL_X = LABEL_X + 9.5
    base = LABEL_X + 11.2
    DX = 1.15
    COL_X = {}
    sample_center = {}
    x = base
    if has_div:
        GAP = 1.7
        for s in SAMPLES:
            COL_X[(s, "conf")] = x
            COL_X[(s, "div")] = x + DX
            sample_center[s] = x + DX / 2
            x += DX + GAP
        right_x = COL_X[(SAMPLES[-1], "div")]
    else:
        STEP = 2.3
        for s in SAMPLES:
            COL_X[(s, "conf")] = x
            COL_X[(s, "div")] = x
            sample_center[s] = x
            x += STEP
        right_x = COL_X[(SAMPLES[-1], "conf")]
    fig_h = max(5.0, 0.5 * (ycount + 1.6) + 1.6)
    fig, ax = plt.subplots(figsize=(15.5, fig_h))

    for n in btree.preorder_node_iter():
        if n.parent_node is not None:
            ax.plot([n.parent_node._x, n._x], [n._y, n._y], color="#444", lw=1.0, zorder=1)
        if not n.is_leaf():
            cys = [c._y for c in n.child_nodes()]
            ax.plot([n._x, n._x], [min(cys), max(cys)], color="#444", lw=1.0, zorder=1)

    cmax = max([lineage_full_count[reps[lf.taxon.label]] for lf in leaves] + [1])
    qmax = max([per_sample_conf[s][reps[lf.taxon.label]] for s in SAMPLES for lf in leaves]
               + [per_sample_rel[s][reps[lf.taxon.label]] for s in SAMPLES for lf in leaves]
               + [1])

    for lf in leaves:
        lin = reps[lf.taxon.label]
        guild = lineage_guild[lin]
        col = GUILD_COLOR.get(guild, "#9e9e9e")
        nref = lineage_full_count[lin]
        hh = 0.16 + 0.34 * (math.log10(nref + 1) / math.log10(cmax + 1))
        ax.add_patch(Polygon([(lf._x, lf._y), (lf._x + 0.55, lf._y - hh),
                              (lf._x + 0.55, lf._y + hh)], closed=True,
                             facecolor=col, edgecolor="#333", lw=0.6, zorder=2))
        ax.text(LABEL_X, lf._y, lin, va="center", ha="left", fontsize=9.5)
        ax.text(NCOL_X, lf._y, f"{nref}", va="center", ha="right", fontsize=8, color="#666")
        for s in SAMPLES:
            for kind, table in (("conf", per_sample_conf), ("div", per_sample_rel)):
                v = table[s][lin]
                if v <= 0:
                    continue
                xx = COL_X[(s, kind)]
                if kind == "conf":
                    ax.scatter([xx], [lf._y], s=bubble_size(v, qmax), facecolor=col,
                               alpha=0.95, edgecolors="#222", linewidths=0.5, zorder=4)
                    tcol = "white"
                else:
                    ax.scatter([xx], [lf._y], s=bubble_size(v, qmax),
                               facecolor=to_rgba(col, 0.30), edgecolors=col,
                               linewidths=1.4, linestyle="--", zorder=3)
                    tcol = "#222"
                ax.text(xx, lf._y, str(v), va="center", ha="center", fontsize=8,
                        color=tcol, zorder=5, fontweight="bold")

    top = ycount - 0.1
    ax.text(LABEL_X, top + 0.6, "CuMMO lineage (collapsed clade)", fontsize=10.5,
            fontweight="bold", va="bottom")
    ax.text(NCOL_X, max(lf._y for lf in leaves) + 0.58, "Reference\nsequences",
            fontsize=7.5, color="#555", ha="center", va="bottom", linespacing=0.9)
    for s in SAMPLES:
        ax.text(sample_center[s], top + 0.55, s, fontsize=9.5, fontweight="bold",
                ha="center", va="bottom", linespacing=0.95)

    def sw(color, solid, label):
        if solid:
            return Line2D([0], [0], marker='o', color='w', markerfacecolor=color,
                          markeredgecolor="#222", markersize=11, label=label)
        return Line2D([0], [0], marker='o', color='w', markerfacecolor=to_rgba(color, 0.30),
                      markeredgecolor=color, markeredgewidth=1.6, markersize=13, label=label)
    gc = GUILD_COLOR
    GUILD_TITLE = {"pmoA": "Methanotroph (pmoA)", "amoA": "Ammonia oxidiser (amoA)",
                   "alkaneMO": "Divergent CuMMO (alkane/alkene MO)"}

    present = set()
    for s in SAMPLES:
        for lin, c in per_sample_conf[s].items():
            if c > 0:
                present.add((lineage_guild[lin], "conf"))
        for lin, c in per_sample_rel[s].items():
            if c > 0:
                present.add((lineage_guild[lin], "div"))
    has_div = any(k == "div" for _, k in present)
    combined = []
    for g in ("pmoA", "amoA", "alkaneMO"):
        if (g, "conf") in present:
            combined.append(sw(gc[g], True,
                                GUILD_TITLE[g] + (" — confident" if has_div else "")))
        if (g, "div") in present:
            combined.append(sw(gc[g], False, GUILD_TITLE[g] + " — closest clade (divergent)"))
    ax.set_xlim(-0.5, right_x + 1.4)
    ax.set_ylim(-1.8, ycount + 1.8)
    ax.legend(handles=combined, loc="upper center", frameon=False,
              fontsize=9, ncol=min(3, len(combined)),
              columnspacing=1.6, handletextpad=0.6,
              bbox_to_anchor=(0.5, -0.72),
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
        copy2(OUT_FIG, os.path.join(supp, "pmoA_lineage_placement.png"))
        copy2(OUT_FIG[:-4] + ".pdf", os.path.join(supp, "pmoA_lineage_placement.pdf"))
    except Exception as e:
        print("supp copy failed:", e)

if __name__ == "__main__":
    main()
