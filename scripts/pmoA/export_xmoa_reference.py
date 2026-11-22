import os

from treesapp.refpkg import ReferencePackage

PKL = os.path.expanduser(
    "~/anaconda3/envs/treesapp/lib/python3.9/site-packages/treesapp/data/XmoA_build.pkl"
)
DB = "/mnt/c/Users/nojus/OneDrive/Desktop/Bubble project/pmoA/database"

def main():
    os.makedirs(DB, exist_ok=True)
    r = ReferencePackage()
    r.f__pkl = PKL
    r.slurp()
    print(f"prefix={r.prefix} molecule={r.molecule} num_seqs={r.num_seqs} model={getattr(r,'sub_model',None)}")

    msa = None
    for attr in ("msa", "f__msa"):
        v = getattr(r, attr, None)
        if isinstance(v, str) and v.strip().startswith(">"):
            msa = v
            break

    if msa is None and hasattr(r, "__dict__"):
        for k, v in r.__dict__.items():
            if isinstance(v, str) and v.count(">") >= 5 and "\n" in v:
                msa = v
                print(f"[msa] recovered from attribute '{k}'")
                break
    if msa:
        with open(os.path.join(DB, "XmoA_ref.faa"), "w") as fh:
            fh.write(msa if msa.endswith("\n") else msa + "\n")
        print(f"wrote XmoA_ref.faa ({msa.count('>')} seqs)")
    else:
        print("[warn] could not locate MSA string; will disband instead")

    tree = None
    for attr in ("tree", "f__tree"):
        v = getattr(r, attr, None)
        if isinstance(v, str) and v.strip().startswith("("):
            tree = v
            break
    if tree:
        with open(os.path.join(DB, "XmoA_ref.nwk"), "w") as fh:
            fh.write(tree if tree.endswith("\n") else tree + "\n")
        print("wrote XmoA_ref.nwk")
    else:
        print("[warn] could not locate tree string; will disband instead")

    with open(os.path.join(DB, "XmoA_model.txt"), "w") as fh:
        fh.write(str(getattr(r, "sub_model", "LG+G4")) + "\n")

    leaves = r.generate_tree_leaf_references_from_refpkg()
    with open(os.path.join(DB, "XmoA_leaf_taxonomy.tsv"), "w") as fh:
        fh.write("leaf_id\taccession\tdescription\tlineage\n")
        for lf in leaves:
            acc = getattr(lf, "description", "")
            desc = getattr(lf, "organism", "") or ""
            lin = getattr(lf, "lineage", "") or ""
            fh.write(f"{lf.number}\t{acc}\t{desc}\t{lin}\n")
    print(f"wrote XmoA_leaf_taxonomy.tsv ({len(leaves)} leaves)")

    if msa is None or tree is None:
        out = os.path.join(DB, "_disband")
        os.makedirs(out, exist_ok=True)
        try:
            r.disband(out)
            print("disband ->", out)
            for root, _, fs in os.walk(out):
                for f in fs:
                    print("   ", os.path.relpath(os.path.join(root, f), out))
        except Exception as e:
            print("disband err:", repr(e))

if __name__ == "__main__":
    main()
