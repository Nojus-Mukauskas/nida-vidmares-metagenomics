# pmoA / amoA

SqueezeMeta assigns both *pmoA* and *amoA* to K10944, so the K10944 TPM combines both functions. Query proteins were placed on the TreeSAPP XmoA reference with EPA-ng. Each query was assigned the guild of its closest XmoA reference tip (*pmoA*, *amoA*, or alkane/alkene monooxygenase). Only the *pmoA* guild was retained as methanotrophy.

Proteins (`03.*.faa`) from all five SqueezeMeta projects were analysed. The XmoA reference is a protein alignment and tree.

Paths refer to the original local `Bubble project` directory.

## Protocol

1. **Extract query proteins** — `extract_pmoa_proteins.py`  
   ORFs annotated as K10944 were taken from `08.*.fun3.kegg`, together with K16157 (*mmoX*). K16157 was extracted but excluded from the *pmoA* figures. Sequences were retrieved from `03.*.faa`.  
   Output: `pmoA_amoA_query_proteins.faa`.

2. **Export the XmoA reference** — `export_xmoa_reference.py`  
   The TreeSAPP XmoA alignment, tree and leaf taxonomy were written from the reference package.

3. **Re-estimate the reference tree** — `build_pmoa_reference.sh`  
   Branch lengths were re-estimated with RAxML-NG (`--evaluate`, LG+G4). As for *mcrA*, EPA-ng uses the alignment (`ref_msa.faa`) and the tree (`ref_eval.raxml.bestTree`).

4. **Place the query proteins** — `place_pmoa_on_tree.sh`  
   Placement followed the *mcrA* procedure: MAFFT `--addfragments --keeplength` (L-INS-i; `--anysymbol` for amino-acid sequences), EPA-ng and gappa.

5. **Assign guilds** — `plot_pmoa_lineages.py`  
   Guilds were taken from the taxonomy of the closest XmoA tip. Type I and type II methanotrophs were assigned to *pmoA*; Thaumarchaeota, *Nitrosomonas* and *Nitrosococcus* to *amoA*; the remaining copper monooxygenases to alkane/alkene monooxygenases. The *pmoA* fraction of K10944 was converted to TPM for the figures.

```
python scripts/pmoA/extract_pmoa_proteins.py
python scripts/pmoA/export_xmoa_reference.py
bash scripts/pmoA/build_pmoa_reference.sh
bash scripts/pmoA/place_pmoa_on_tree.sh
python scripts/pmoA/plot_pmoa_lineages.py
```

Files in the repository: `data/pmoA_proteins.faa`, `data/pmoA_lineages.tsv`.
