#!/usr/bin/env bash
set -euo pipefail
export TMPDIR=/tmp

ROOT="/mnt/c/Users/nojus/OneDrive/Desktop/Bubble project/pmoA"
DB="$ROOT/database"
MSA="$DB/ref_msa.faa"
TREE="$DB/ref.nwk"
MODEL="$(cat "$DB/XmoA_model.txt" 2>/dev/null || echo LG+G4)"
MODEL="${MODEL%%[[:space:]]*}"

[ -f "$MSA" ]  || { echo "missing $MSA"; exit 1; }
[ -f "$TREE" ] || { echo "missing $TREE"; exit 1; }

raxml-ng --evaluate \
         --msa "$MSA" \
         --tree "$TREE" \
         --model "$MODEL" \
         --prefix "$DB/ref_eval" \
         --threads 2 \
         --seed 42 \
         --force perf_threads
ls -la "$DB/ref_eval.raxml.bestModel" "$DB/ref_eval.raxml.bestTree"
