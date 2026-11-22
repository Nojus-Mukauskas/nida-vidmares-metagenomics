import sys
import dendropy

src, dst = sys.argv[1], sys.argv[2]
t = dendropy.Tree.get(path=src, schema="newick", preserve_underscores=True)
t.resolve_polytomies(update_bipartitions=False)
for e in t.preorder_edge_iter():
    if e.length is None:
        e.length = 0.0
t.write(path=dst, schema="newick", unquoted_underscores=True, suppress_rooting=True)
print(f"resolved -> {dst}")
