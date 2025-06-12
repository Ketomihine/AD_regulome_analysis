# AD Regulome Analysis

This repository contains scripts and workflows used to process and analyze single cell ATAC-seq and RNA-seq data for studying the Alzheimer's disease (AD) regulome. The code is organised into directories that mirror the different analysis steps. Below is a description of each directory and the files contained within.

## Root files

- `README.org` – Short org‐mode notes on the repository purpose.
- `.gitmodules` – Defines the `epiclust` submodule (empty here).
- `.gitignore` – Standard ignore rules for Python and R projects.

## `integration_linking_modules`

Python, R and shell scripts used for peak–gene linking, module inference and
multi‑omic integration.

<details>
<summary>File listing</summary>

- **aggregate_modules.py** – Aggregate module counts across clusters.
- **auprc.py / auprc_full.py** – Compute precision–recall curves for predicted links.
- **coaccessibility.py / coaccessibility.sh** – Combine gene and peak matrices to
  compute co‑accessibility.
- **concat_bedgraph.sh** – Helper to concatenate bedGraph files.
- **count_peaks_in_fragments.sh** – Count peak overlaps in fragment files.
- **count_peaks_to_h5ad.py** – Convert peak counts into an AnnData object.
- **deg.R** – Differential expression/peak statistics for modules.
- **diff_mod.R** – Differential module analysis.
- **epiclust_generate_interactions.py** – Create candidate peak–gene pairs for
  `epiclust` linking.
- **epiclust_linking.py / epiclust_linking.sh** – Run the `epiclust` linking
  model on pseudobulk data.
- **epiclust_norm_links.py** – Normalise link scores and evaluate performance.
- **epiclust_predict_links.py / epiclust_predict_links.sh** – Predict links from
  training data using logistic regression.
- **filter_tilematrix_to_peaks.py** – Annotate a TileMatrix with peak labels.
- **fit_adgwas.R** – Fit AD GWAS enrichment models to linked peaks.
- **gene_distance.py** – Utility functions to compute distances between peaks and
  genes.
- **gene_estimation.py** – Estimate gene activity from ATAC data with optional
  linking information.
- **generate_ct_bedgraph.sh** – Produce per–cell‑type bigWig/bedGraph tracks.
- **go_tf.R** – GO/TF enrichment analyses on modules.
- **gtf.py** – Load GTF annotations and compute gene length scores.
- **impute.py** – TF‑IDF normalisation and NMF imputation utilities.
- **integrate.py / integrate.sh** – Integrate ATAC and RNA modalities using
  Harmony/BBKNN and transfer labels.
- **linking.py** – Training and prediction framework for distance‑aware linking.
- **load_fragments_from_annot.py / load_fragments_from_annot.sh** – Build
  SnapATAC objects from fragment files using sample annotations.
- **module_pseudobulk.py** – Create module level pseudobulk counts.
- **module_stats.R / module_stats.py** – Summaries and statistics for modules.
- **pseudobulk.py / pseudobulk.sh** – Generate pseudobulk AnnData objects from
  integration results.
- **run_modules.py** – Run `epiclust` module discovery.
- **run_seurat_signac_integration.sh** – Placeholder for Seurat/Signac
  integration job.
- **seurat_signac_integration.R** – Example workflow integrating scRNA and
  scATAC in Seurat/Signac.
- **seurat_signac_multiome.R** – Seurat script for 10x multiome data.
- **seurat_signac_peaks.R** – Peak calling with Signac/Cicero.
- **subtype_enrichment.R** – Differential enrichment of modules across
  sub‑types.
</details>

## `snATAC.processing`

Scripts for preprocessing the single‑nucleus ATAC‑seq dataset and downstream
analyses.

<details>
<summary>File listing</summary>

- **1.scATAC.processing.R** – Main ArchR pipeline for creating the project,
  performing QC, LSI reduction and initial clustering.
- **2.RemoveDoublets.iter1.R** – Remove doublets and rerun clustering (first
  iteration) with extensive QC plots.
- **3.RemoveDoublets.iter2.R** – Second iteration of doublet removal and
  clustering.
- **4.callpeak.and.GREAT_annotation.R** – Call peaks with MACS2, build peak
  matrix and annotate with rGREAT.
- **5.co_accessibility.R** – Compute co‑accessibility loops with ArchR.
- **6.TF.candidate.R** – Identify candidate transcription factor regulators using
  motif enrichment and gene scores.
- **6_2.TF_footprinting.R** – Footprinting analysis on candidate TFs.
- **7.aQTL.calling/** – Scripts for allele‑specific QTL analysis:
  - `1.normalized.R` – Normalise pseudobulk peak counts.
  - `2.hg38Tohg19.sh` – LiftOver peak coordinates to hg19.
  - `3.attach.pos.sh` – Append genomic positions to normalised matrix.
  - `4.bgzip_index.sh` – Sort, bgzip and index matrices.
  - `5.PEER.sh` – Submit PEER factor estimation jobs.
  - `6.generate.Fastqtl.running.script.sh` – Generate FastQTL job scripts.
  - `attach.pos.pl` – Perl helper used by `3.attach.pos.sh`.
  - `get_Fastqtl.script.pl` – Perl script creating FastQTL commands.
  - `run_PEER.R` – PEER factor estimation utility.
- **8.run_differential_02c_NB_wRUV.non_vs_early.R** – Nebula differential peak
  testing (nonAD vs early AD).
- **8_2.run_differential_02c_NB_wRUV.all.R** – Nebula differential testing across
  all pathology groups.
- **9.Cell.composition.R** – Estimate cell‑type composition differences using
  `speckle` and `limma`.
- **get_markerMat.PEC_markers_2.byCluster.re.pl** – Perl script to aggregate
  marker matrices by cluster.
</details>

## `Process.Morabito.etal.data`

Processing steps for the external dataset from Morabito *et al.*

- **1.setup.TSS1.R** – Create ArchR arrows from raw fragments.
- **2.ArchR.clustering.R** – Perform iterative LSI and clustering.
- **3.assign.celltype.R** – Assign cell types to clusters based on marker genes.
- **4_1.epigenome.chromHMM.R** – Annotate peaks with ChromHMM states.
- **4_2.erosion.score.R** – Compute chromatin erosion scores across states.

## `epiclust`

This folder is empty in the repository but is defined as a git submodule
pointing to [kellislab/epiclust](https://github.com/kellislab/epiclust).  It
contains the implementation of the peak‑gene linking and module algorithms.
Clone the submodule to obtain the code.

---

The repository combines R, Python and shell scripts.  Many scripts assume data
paths specific to the authors' environment, so paths may need adaptation before
running.  See the comments within each script for usage details.
