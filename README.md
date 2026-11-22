Scripts, Nextflow configs and SqueezeMeta commands for shotgun metagenomes from two lagoon sites, Nida (Muddy) and Vidmarės (Sandy): sediment and water, three seasons.

Raw reads: NCBI SRA, BioProject [PRJNA1551516](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1551516).

Paths refer to the original analysis directories (local `Bubble project` folder). The SqueezeMeta configuration files also retain the server paths under `/data/squeezemeta/...`.

## Samples

45 Illumina libraries (2 × 150 bp). Repo and pipeline IDs use the original codes; paper figures use season + site (+ depth for water) + replicate:

| Original | Paper | Meaning |
| --- | --- | --- |
| PN1 | SpMu1 | Spring Muddy sediment, replicate 1 |
| PNS2 | SpMuS2 | Spring Muddy surface water, replicate 2 |
| PNB1 | SpMuB1 | Spring Muddy bottom water, replicate 1 |
| PVF1 | SpSaS1 | Spring Sandy surface water (filtered), replicate 1 |
| VV2 | SuSa2 | Summer Sandy sediment, replicate 2 |
| RN3 | AuMu3 | Autumn Muddy sediment, replicate 3 |

- first letter: P spring, V summer, R autumn → Sp / Su / Au  
- second letter: N Nida (Muddy), V Vidmarės (Sandy) → Mu / Sa  
- water: S surface, B bottom, F filtered surface (Vidmarės) → S or B in the paper ID  
- sediment has no depth letter (PN1, VV2, …)

Full mapping for all 45 libraries is in `data/sample_id_map.tsv`. Raw and filtered read pairs and Nonpareil coverage are in Supplementary Table S1 (`Sample` / `Original ID`).

**Excluded from PCoA, genus bars and SqueezeMeta:** SpMuS1 (PNS1) and SpMuS3 (PNS3). SpMuS2 (PNS2) was kept.

| Original | Paper | Raw read pairs | Nonpareil C |
| --- | --- | --- | --- |
| PNS1 | SpMuS1 | 698,578 | 0.35 |
| PNS3 | SpMuS3 | 2,304,027 | 0.39 |
| other 43 libraries | — | ~20.7–43.0 million (median ~27.0 million) | — |

PNS1 contained approximately 40-fold fewer raw read pairs than the median library, and PNS3 approximately 12-fold fewer. Per-sample values are given in Supplementary Table S1.

FASTQ, assemblies, Kraken databases and Nextflow `work/` are not in this repo.

```
config/               detaxizer, taxprofiler, SqueezeMeta
data/                 sample table, ID map, genus list, mcrA / pmoA sequences and placements
scripts/              plots and marker-gene placement
```

## 1. Decontamination

nf-core/detaxizer 1.3.0, two passes, Kraken2 `k2_core_nt`. Kraken2 was run with `--confidence 0.2 --minimum-hit-groups 3`, as set in the `.config` files. These settings differ from the YAML defaults.

**Pass 1 — eukaryotes**, including trimming and host removal. fastp and BBDuk parameters are set in `config/detaxizer.config`:

- fastp: `--cut_right --cut_window_size 4 --correction --detect_adapter_for_pe --trim_poly_g`. Sliding-window mean quality was 20 and the minimum length was 50 (`config/detaxizer_params.yml`)
- BBDuk: `hdist=1`, reference `human_phix_ref.fasta` (T2T-CHM13v2.0 + phiX)

```
nextflow run nf-core/detaxizer -profile mamba \
  -c config/detaxizer.config \
  -params-file config/detaxizer_params.yml
```

**Pass 2 — viruses.** fastp and BBDuk were not repeated:

```
nextflow run nf-core/detaxizer -profile mamba \
  -c config/detaxizer_virus.config \
  -params-file config/detaxizer_virus_params.yml
```

## 2. Coverage

Nonpareil 3.5.5, k-mer mode, k = 24, on virus-filtered **R1** files. Reads shorter than 24 bp were removed first (`seqkit seq -m 24`). Nonpareil terminates when reads shorter than *k* are present.

```
nonpareil -s SAMPLE_R1.fastq -T kmer -k 24 -f fastq -t 12
```

Coverage (C) in Table S1 is the final redundancy value in the `.npo` file, transformed by the coverage factor at an overlap of 50. This corresponds to `Nonpareil.read_data` without the nonlinear least-squares fit.

## 3. Community taxonomy

nf-core/taxprofiler 1.2.6 was run on the same filtered reads. Quality control and host removal were disabled. Kraken2 was run with `--confidence 0.1`. Bracken abundances were estimated at species, genus and family (`-r 150 -t 10`). PCoA and bar plots use genus-level profiles.

```
nextflow run nf-core/taxprofiler -r 1.2.6 -profile mamba \
  -c config/taxprofiler.config \
  -params-file config/taxprofiler_params.yml
```

All 45 libraries were classified. PCoA and bar plots were calculated for the 43 retained libraries, using Bray–Curtis distances on genus-level relative abundance among classified reads. The five SqueezeMeta projects were defined from this ordination by grouping samples with similar communities, rather than assembling each triplicate set separately.

```
Rscript scripts/plot_pcoa.R
Rscript scripts/plot_genus_composition.R
Rscript scripts/plot_methane_genera.R
Rscript scripts/plot_community_methane.R
```

`plot_community_methane.R` produces the community figure used in the paper. Methane-cycling genera are listed in `data/methane_cycling_genera.tsv`.

## 4. Assembly and functional annotation (SqueezeMeta)

SqueezeMeta 1.7.2 was run as coassembly, with read mapping, gene calling and functional annotation, and without binning. Contigs were assembled with MEGAHIT (`meta-sensitive`; minimum length 200 bp). Functional annotation used Diamond blastx in double-pass mode. Input reads were those used for taxprofiler.

Each directory under `config/squeezemeta/` contains one SqueezeMeta project: the command, the stored `SqueezeMeta_conf.pl`, and the sample FASTQ list.

| project | samples |
| --- | --- |
| Muddy sediment (`sed_n`) | PN1–3, VN1–3, RN1–3 |
| Sandy sediment (`sed_v`) | PV1–3, VV1–3, RV1–3 |
| spring Sandy water (`water_spring_vidmares`) | PVF1–3 |
| summer water (`water_summer`) | VVF1–3, VNS1–3, VNB1–3 |
| autumn water and the Muddy spring water libraries not included above (`water_autumn_N_V`) | RNS1–3, RNB1–3, RVF1–3, PNB1–3, PNS2 |

Gene abundances are SqueezeMeta TPM values (`KO.tpm.tsv`). *pmoA* was distinguished from the other K10944 genes by placement on the XmoA tree (section 6).

## 5. mcrA

SqueezeMeta annotates *mcrA* as K00399, which does not distinguish lineages. Nucleotide sequences of those ORFs were placed on mcrA4All. *mcrA* abundance is the K00399 TPM. K00399 was not detected in water.

Protocol: `scripts/mcrA/README.md`.

```
python scripts/mcrA/extract_mcra_sequences.py
python scripts/mcrA/clean_mcra_reference.py
bash scripts/mcrA/build_mcra_reference.sh
bash scripts/mcrA/place_mcra_on_tree.sh
python scripts/mcrA/plot_mcra_lineages.py
```

`data/mcrA_sequences.fna`, `data/mcrA_lineages.tsv`.

## 6. pmoA vs amoA

K10944 comprises *pmoA*, *amoA* and related alkane/alkene monooxygenases. Proteins from all five SqueezeMeta projects were placed on the TreeSAPP XmoA reference with EPA-ng. Each query was assigned the guild of its closest XmoA reference tip (*pmoA*, *amoA*, or alkane monooxygenase). Only the *pmoA* guild was retained as methanotrophy, and that fraction of K10944 was converted to TPM.

Protocol: `scripts/pmoA/README.md`.

```
python scripts/pmoA/extract_pmoa_proteins.py
python scripts/pmoA/export_xmoa_reference.py
bash scripts/pmoA/build_pmoa_reference.sh
bash scripts/pmoA/place_pmoa_on_tree.sh
python scripts/pmoA/plot_pmoa_lineages.py
```

`data/pmoA_proteins.faa`, `data/pmoA_lineages.tsv`.

## 7. mcrA / pmoA abundances

`scripts/get_mcra_and_pmoa_tpm.py` reads SqueezeMeta KO TPM values and calculates phylogenetic *pmoA* TPM. The plotting scripts import it. *mcrA* abundance is K00399 TPM. *pmoA* abundance is the *pmoA*-guild fraction of K10944 after the assignment in section 6.

`plot_mcra_pmoa.py` draws the sediment and water boxplots. Stars are two-sided Welch t-tests on TPM, Holm-adjusted within each comparison: Muddy versus Sandy in each sediment season, and in water Muddy surface versus Sandy and Muddy surface versus Muddy bottom. Spring Muddy surface water is one library and is not tested. The sediment two-way ANOVA (site × season on log10 TPM) was calculated separately.

`plot_methane_heatmap.py` plots methyl-CoM reductase (*mcrA*, *mcrB*, *mcrG*), the methanogen CODH/ACS subunits (*cdhC*, *cdhD*, *cdhE*), hydrogenotrophic markers (*mtd*, *mtrH*), methylotrophic markers (*mtaA*, *mtbB*) and phylogenetic *pmoA*, as log10(TPM + 0.1).

```
python scripts/plot_mcra_pmoa.py
python scripts/plot_methane_heatmap.py
```

## Software

Nonpareil 3.5.5; nf-core/detaxizer 1.3.0; nf-core/taxprofiler 1.2.6; Nextflow 25.10.4; Kraken2 2.1.5; Bracken 3.0.1; SqueezeMeta 1.7.2 (MEGAHIT meta-sensitive); mcrA4All (April 2023); TreeSAPP XmoA; MAFFT, FastTree, RAxML-NG, EPA-ng, gappa; vegan (Bray–Curtis PCoA).
