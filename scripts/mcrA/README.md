# mcrA

SqueezeMeta annotates *mcrA* as K00399. This KO does not distinguish lineages (Methanothrix, Methanobacterium, ANME). Nucleotide sequences of the corresponding ORFs were placed on the mcrA4All reference, and each query was assigned the lineage of the nearest reference tips. Abundance is the K00399 TPM reported by SqueezeMeta. The scripts in this directory perform phylogenetic placement only.

K00399 was not detected in the water projects. Placement was limited to the two sediment projects (`sed_n`, Muddy; `sed_v`, Sandy).

Paths refer to the original local `Bubble project` directory.

## Protocol

1. **Extract query sequences** — `extract_mcra_sequences.py`  
   ORF identifiers annotated as K00399 were taken from `08.*.fun3.kegg`, and the corresponding nucleotide records from `03.*.fna`. mcrA4All is a nucleotide alignment, so the queries are DNA.  
   Output: `mcrA_query_nuc.fna`. Example header: `>mcrA|sed_n|sed_n_461168_292-393`.

2. **Prepare the reference** — `clean_mcra_reference.py`  
   The mcrA4All alignment (April 2023) was reformatted for tree inference. Sequence identifiers were converted to alphanumeric names. Uracil was replaced with thymine, and dots with gap characters.  
   Output: `ref_named.fasta`.

3. **Infer the reference tree** — `build_mcra_reference.sh`  
   Query sequences were excluded from this step. Gap-only columns were removed with trimAl. A starting tree was inferred with FastTree (GTR+GAMMA), and branch lengths were re-estimated with RAxML-NG (GTR+G).  
   EPA-ng requires the alignment (`ref_msa.fasta`) and the tree (`ref_eval.raxml.bestTree`).

4. **Place the query sequences** — `place_mcra_on_tree.sh`  
   Queries were inserted into the reference alignment with MAFFT (`--addfragments --keeplength`, L-INS-i), without adding columns. Each query was placed with EPA-ng and grafted onto the reference tree with gappa.

5. **Assign lineages** — `plot_mcra_lineages.py`  
   Each query was labelled from the nearest mcrA4All reference tips.

```
python scripts/mcrA/extract_mcra_sequences.py
python scripts/mcrA/clean_mcra_reference.py
bash scripts/mcrA/build_mcra_reference.sh
bash scripts/mcrA/place_mcra_on_tree.sh
python scripts/mcrA/plot_mcra_lineages.py
```

Files in the repository: `data/mcrA_sequences.fna`, `data/mcrA_lineages.tsv`.
