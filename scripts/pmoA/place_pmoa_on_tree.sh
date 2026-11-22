#!/usr/bin/env bash
set -euo pipefail
export TMPDIR=/tmp

ROOT="/mnt/c/Users/nojus/OneDrive/Desktop/Bubble project/pmoA"
DB="$ROOT/database"
Q="$ROOT/queries/pmoA_amoA_query_proteins.faa"
R="$ROOT/results"
mkdir -p "$R"

[ -f "$Q" ] || { echo "missing $Q"; exit 1; }

mafft --addfragments "$Q" --keeplength --anysymbol --localpair --maxiterate 1000 --thread 4 \
      "$DB/ref_msa.faa" > "$R/combined.aln" 2> "$R/mafft.log"
seqkit grep -nrp '^pmoA_amoA\|' "$R/combined.aln" > "$R/query_aln.faa"
epa-ng --redo \
       --tree "$DB/ref_eval.raxml.bestTree" \
       --ref-msa "$DB/ref_msa.faa" \
       --query "$R/query_aln.faa" \
       --model "$DB/ref_eval.raxml.bestModel" \
       --outdir "$R"
gappa examine graft --jplace-path "$R/epa_result.jplace" --out-dir "$R" --allow-file-overwriting
