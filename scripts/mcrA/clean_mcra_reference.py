import os, re, collections

DB = "/mnt/c/Users/nojus/OneDrive/Desktop/Bubble project/mcrA/database"
SRC = os.path.join(DB, "mcrA4All", "mcrA_April2023.fasta")
OUT_FASTA = os.path.join(DB, "ref_named.fasta")
OUT_MAP = os.path.join(DB, "name_map.tsv")

TRAIL = re.compile(r"\s+\d+\s*bp\s+\w+\s*$", re.IGNORECASE)
NONALNUM = re.compile(r"[^A-Za-z0-9]+")

def clean_name(header):
    name = TRAIL.sub("", header).strip()
    name = NONALNUM.sub("_", name).strip("_")
    return name or "seq"

VALID = set("ACGTRYSWKMBDHVN-")

def norm_seq(s):
    s = (s.upper()
          .replace("U", "T")
          .replace(".", "-")
          .replace("~", "-")
          .replace(" ", ""))
    return "".join(c if c in VALID else "N" for c in s)

def read_fasta(path):
    name, buf = None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:], []
            else:
                buf.append(line)
    if name is not None:
        yield name, "".join(buf)

def main():
    seen = collections.Counter()
    alphabet = collections.Counter()
    widths = set()
    n = 0
    with open(OUT_FASTA, "w", encoding="utf-8") as of, \
         open(OUT_MAP, "w", encoding="utf-8") as mf:
        mf.write("clean_id\toriginal_header\n")
        for header, seq in read_fasta(SRC):
            cid = clean_name(header)
            seen[cid] += 1
            if seen[cid] > 1:
                cid = f"{cid}_{seen[cid]}"
            s = norm_seq(seq)
            widths.add(len(s))
            alphabet.update(s[:200])
            of.write(f">{cid}\n{s}\n")
            mf.write(f"{cid}\t{header}\n")
            n += 1
    print(f"sequences written: {n}")
    print(f"aligned width(s): {sorted(widths)[:5]}{' ...' if len(widths) > 5 else ''}")
    print(f"alphabet sample: {dict(alphabet.most_common())}")
    print(f"-> {OUT_FASTA}")
    print(f"-> {OUT_MAP}")

if __name__ == "__main__":
    main()
