# Temporal Transcriptomics Identifies Trans-regulation by Multiple lncRNAs in Human iPSCs

Analysis code for the manuscript *"Temporal Transcriptomics Identifies Trans-regulation by Multiple lncRNAs in Human iPSCs"* (Liu et al., Rinn Lab). This repository contains the full pipeline used to go from raw sequencing reads to the gene-level and transcript-level trans-regulatory analyses, the Doxycycline-Responsive Gene Signature (DRGS), the RF-DBN regulatory network, chromatin-accessibility dynamics, and the associated figures.

The study induces 17 individual lncRNAs (via a Tet-On/H11-safe-harbor system) in human iPSCs and profiles high-density time courses of RNA-seq (0–96 h) and ATAC-seq (0–2.5 h) to identify genes and transcripts regulated in *trans*.

## Contents

### 1. Upstream processing (SLURM / Nextflow)

| File | Purpose |
|---|---|
| `run_rnaseq_pipeline.sh` | SLURM/Nextflow driver (`nf-core/rnaseq` v3.14.0) that runs STAR alignment and Salmon pseudo-alignment/quantification on the raw RNA-seq FASTQs against GRCh38.p13/GENCODE v38. |
| `salmon_with_bootstrap.sh` | SLURM array job that re-quantifies every sample directly with `salmon quant` (bias correction + 100 Gibbs inferential-replicate bootstraps) against a GENCODE v43 transcriptome index, producing the per-sample `quant.sf`/inferential-replicate output used by the transcript-level edgeR analysis. |
| `run_atacseq_pipeline.sh` | SLURM/Nextflow driver (`nf-core/atacseq` v2.1.2) that runs the ATAC-seq raw-read processing and consensus peak calling against GRCh38.p13/GENCODE v38. |

### 2. Gene-level RNA-seq analysis

| File | Purpose |
|---|---|
| `rnaseq.Rmd` | Gene-level RNA-seq pipeline: per-line DESeq2 differential expression across the time course, DEG consistency analysis, definition of the 209-gene Doxycycline-Responsive Gene Signature (DRGS), the stemness "ribbon plot," PCA, DEG overlap (UpSet/Venn), the DRGS heatmap, and GO/GSEA enrichment (including meta-analysis against 28 public datasets). |

### 3. Transcript-level (isoform) trans-target analysis 

| File | Purpose |
|---|---|
| `Transcript_level_drgs_pipeline_rmd.md` | End-to-end **transcript-level** pipeline, from raw Salmon quantifications to the published figures: edgeR time-course and inferential-replicate transcript DE, DEG consistency, DRGS heatmap, transgene–transcript correlation filtering (with GFP-control subtraction), curve-deviation target calling, genome-wide target distribution ("Manhattan" plot), the H11 *cis*-window control, and the LINC00847 target heatmap/trajectories. |

### 4. RF-DBN regulatory network

| File | Purpose |
|---|---|
| `RF.Rmd` | Random Forest regression step of the RF-DBN strategy: for each DRGS lncRNA, predicts its own expression from all other genes (plus time) using grouped cross-validation, and extracts `%IncMSE` feature importance to select the top ~2000 candidate genes carried into the DBN step. |
| `dbn_analysis.R` | Main HPC/parallel driver for Dynamic Bayesian Network (DBN) structure learning: partitions each lncRNA's top RF-ranked genes into ~960 overlapping 50-gene feature subspaces, learns a time-respecting DBN (`bnlearn`, hill-climbing with a temporal blacklist) per subspace with 25 bootstrap resamplings, and performs ensemble voting across subspaces/bootstraps to build the robust regulatory network. |
| `dbn_utils.R` | Function library sourced by `dbn_analysis.R`: gene ID conversion, replicate/time-point sample matching, per-time-point differential expression (`limma`), the temporally-constrained DBN builder, lag-correlation analysis, edge scoring/integration, and the bootstrap wrapper. |
| `ensembleDBN.Rmd` | Downstream aggregation and visualization of the DBN output: ensemble-votes bootstrap edges into Gold/Silver/Bronze confidence tiers (Gold = detection frequency of 1.0 across all bootstraps), classifies each lncRNA–gene relationship (downstream target, upstream modulator, synergistic/repressive partner), computes each gene's regulatory "hub degree" (how many lncRNAs target it), and generates the network and hub-gene figures. |

### 5. ATAC-seq / chromatin accessibility

| File | Purpose |
|---|---|
| `atacseq.Rmd` | ATAC-seq analysis: per-condition DESeq2 (LRT) on consensus peak counts across the 0–2.5 h time course, peak filtering and up/down classification, per-condition log2FC heatmaps, and (worked example for HAGLR) linking significant peaks to nearby DEGs by TSS proximity and RNA–ATAC trend concordance. |
| `run_TOBIAS.sh` | SLURM batch script for TF footprinting with TOBIAS on the merged 0 h vs 2.5 h ATAC-seq BAMs: BAM integrity check/merge → `ATACorrect` (Tn5 bias correction) → `ScoreBigwig` (footprint scoring) → `BINDetect` (differential motif footprinting between time points). |
| `TSS_of_DRGs.Rmd` | Builds a BED file of transcription start sites for the 209-gene DRGS set (via `EnsDb.Hsapiens.v86`), used as the region file for the TSS-accessibility profile/heatmap below. |
| `run_deeptools.sh` | Runs `deepTools computeMatrix`/`plotProfile`/`plotHeatmap` on the TOBIAS-corrected 0 h and 2.5 h bigWigs over the DRGS TSS BED file, producing the chromatin-accessibility-at-DRGS-promoters profile and heatmap. |

### 6. Supplementary / reviewer-response analyses

| File | Purpose |
|---|---|
| `baseline_vs_induced_expression.R` | Supplementary figure comparing each transgene's baseline expression (pooled WT+GFP 0 h vs. its own 0 h) to its own maximally-induced timepoint, for all 17 lncRNA lines; also outputs a leakiness ("is 0 h already elevated over baseline") check. |
| `sORFs_miRNAs.Rmd` | Reviewer-response analysis of alternative regulatory mechanisms for the 6 lncRNAs with confirmed trans-targets (DANCR, FENDRR, LINC00667, LINC00847, LNCPRESS1, PNKY): (1) overlap with catalogued Ribo-seq/GENCODE ORFs, (2) miRNA-sponge potential via seed-site scanning and the ENCORI/starBase CLIP database, and (3) sequence homology/complementarity between each lncRNA and its target transcripts. |

## Suggested running order

1. **`run_rnaseq_pipeline.sh` / `salmon_with_bootstrap.sh` / `run_atacseq_pipeline.sh`** — raw-read processing; run these first to generate the STAR/Salmon and ATAC-seq peak outputs consumed by everything below.
2. **`rnaseq.Rmd`** — the gene-level story: DEG calling, the DRGS definition, and its enrichment/consistency across public datasets.
3. **`Transcript_level_drgs_pipeline_rmd.md`** — the transcript/isoform-resolution pipeline; reproduces the trans-target discovery and filtering steps directly from Salmon output.
4. **`RF.Rmd` → `dbn_analysis.R` (+ `dbn_utils.R`) → `ensembleDBN.Rmd`** — the RF-DBN regulatory-network pipeline, in execution order: feature selection, network inference, then ensemble aggregation/visualization.
5. **`atacseq.Rmd` → `run_TOBIAS.sh` → `TSS_of_DRGs.Rmd` → `run_deeptools.sh`** — the chromatin-accessibility side of the analysis: peak-level dynamics, TF footprinting, and the DRGS-promoter accessibility profile.
6. **`baseline_vs_induced_expression.R`** and **`sORFs_miRNAs.Rmd`** — supplementary and reviewer-requested analyses.

