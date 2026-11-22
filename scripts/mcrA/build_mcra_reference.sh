#!/usr/bin/env bash
set -euo pipefail
export TMPDIR=/tmp

DB="/mnt/c/Users/nojus/OneDrive/Desktop/Bubble project/mcrA/database"
cd "$DB"

trimal -in ref_named.fasta -out ref_msa.fasta -noallgaps
FastTree -nt -gtr -gamma -log ref_fasttree.log < ref_msa.fasta > ref.tree
python "$(dirname "$0")/make_binary_tree.py" ref.tree ref_bifurcating.tree
raxml-ng --evaluate \
         --msa ref_msa.fasta \
         --tree ref_bifurcating.tree \
         --model GTR+G \
         --prefix ref_eval \
         --threads 4 \
         --force perf_threads \
         --redo
ls -la ref_msa.fasta ref.tree ref_eval.raxml.bestTree ref_eval.raxml.bestModel
