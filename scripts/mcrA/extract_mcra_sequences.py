import os

BASE = "/path/to/SqueezeMeta/projects"
PROJECTS = ["sed_n", "sed_v"]

OUT = os.path.join(os.getcwd(), "mcrA_query_nuc.fna")
KO_MCRA = "K00399"

def orf_ids_for(proj):
    kegg = os.path.join(BASE, proj, "results", f"08.{proj}.fun3.kegg")
    wanted = {}
    with open(kegg, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[1] == KO_MCRA:
                wanted[p[0]] = "mcrA"
    return wanted

def stream_extract(proj, wanted, out):
    fna = os.path.join(BASE, proj, "results", f"03.{proj}.fna")
    keep = False
    n = 0
    with open(fna, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                orf = line[1:].split()[0].strip()
                keep = orf in wanted
                if keep:
                    n += 1
                    out.write(f">mcrA|{proj}|{orf}\n")
            elif keep:
                out.write(line)
    return n

if __name__ == "__main__":
    total = 0
    with open(OUT, "w", encoding="utf-8") as out:
        for proj in PROJECTS:
            wanted = orf_ids_for(proj)
            n = stream_extract(proj, wanted, out)
            print(f"{proj}: {len(wanted)} mcrA ORFs in kegg, {n} written")
            total += n
    print(f"TOTAL mcrA nucleotide ORFs written: {total} -> {OUT}")
