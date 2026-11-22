import os

PROJ_ROOT = "/mnt/c/Users/nojus/OneDrive/Desktop/Bubble project"
PACK = os.path.join(PROJ_ROOT, "DATA/squeezemeta/FINAL_PCoA_analysis_pack")
FAA_ROOT = os.path.join(PROJ_ROOT, "DATA/squeezemeta/FINAL_PCoA_faa")
OUTDIR = os.path.join(PROJ_ROOT, "pmoA", "queries")

PROJECTS = ["sed_n", "sed_v", "water_autumn_N_V", "water_spring_vidmares", "water_summer"]
KOS = {"K10944": "pmoA_amoA", "K16157": "mmoX"}

def besthit_ids(proj):
    fun3 = os.path.join(PACK, proj, "results", f"08.{proj}.fun3.kegg")
    out = {}
    if not os.path.exists(fun3):
        print(f"[warn] missing {fun3}")
        return out
    with open(fun3, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 2 or not p[0]:
                continue
            for ko, gene in KOS.items():
                if ko in p[1]:
                    out[p[0]] = gene
    return out

def stream_faa(proj, wanted, handles):
    faa = os.path.join(FAA_ROOT, proj, f"03.{proj}.faa")
    mapped = set()
    keep_gene = None
    with open(faa, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                orf = line[1:].split()[0].strip()
                keep_gene = wanted.get(orf)
                if keep_gene is not None:
                    mapped.add(orf)
                    handles[keep_gene].write(f">{keep_gene}|{proj}|{orf}\n")
            elif keep_gene is not None:
                handles[keep_gene].write(line)
    return mapped

def main():
    os.makedirs(OUTDIR, exist_ok=True)
    handles = {g: open(os.path.join(OUTDIR, f"{g}_query_proteins.faa"), "w", encoding="utf-8")
               for g in set(KOS.values())}
    id_rows = []
    totals = {g: {"kegg": 0, "mapped": 0} for g in set(KOS.values())}
    try:
        for proj in PROJECTS:
            wanted = besthit_ids(proj)
            for orf, gene in wanted.items():
                totals[gene]["kegg"] += 1
            mapped = stream_faa(proj, wanted, handles)
            for orf in mapped:
                gene = wanted[orf]
                totals[gene]["mapped"] += 1
                id_rows.append((orf, proj, gene))
            mstr = ", ".join(f"{g}:{sum(1 for o in mapped if wanted[o]==g)}"
                             for g in sorted(set(wanted.values()))) or "none"
            print(f"{proj:22s} kegg={len(wanted):4d}  mapped[{mstr}]")
    finally:
        for h in handles.values():
            h.close()
    with open(os.path.join(OUTDIR, "pmoA_amoA_ORF_ids.tsv"), "w", encoding="utf-8") as f:
        f.write("orf\tproject\tgene\n")
        for r in id_rows:
            f.write("\t".join(r) + "\n")
    print("\n=== totals ===")
    for g, t in totals.items():
        print(f"  {g:10s} KEGG best-hit={t['kegg']:4d}  mapped_to_faa={t['mapped']:4d}  "
              f"blastx_only={t['kegg']-t['mapped']:3d}")
    print(f"\nwrote proteins + id table to {OUTDIR}")

if __name__ == "__main__":
    main()
