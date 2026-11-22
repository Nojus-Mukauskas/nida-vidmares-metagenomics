#!/usr/bin/env bash
set -euo pipefail
export TMPDIR=/tmp

ROOT="/mnt/c/Users/nojus/OneDrive/Desktop/Bubble project/mcrA"
DB="$ROOT/database"
Q="$ROOT/queries/mcrA_query_nuc.fna"
R="$ROOT/results"
mkdir -p "$R"

[ -f "$Q" ] || { echo "missing $Q"; exit 1; }

mafft --addfragments "$Q" --keeplength --localpair --maxiterate 1000 --thread 4 \
      "$DB/ref_msa.fasta" > "$R/combined.aln" 2> "$R/mafft.log"
seqkit grep -nrp '^mcrA\|' "$R/combined.aln" > "$R/query_aln.fasta"
epa-ng --redo \
       --tree "$DB/ref_eval.raxml.bestTree" \
       --ref-msa "$DB/ref_msa.fasta" \
       --query "$R/query_aln.fasta" \
       --model "$DB/ref_eval.raxml.bestModel" \
       --outdir "$R"
gappa examine graft --jplace-path "$R/epa_result.jplace" --out-dir "$R" --allow-file-overwriting
