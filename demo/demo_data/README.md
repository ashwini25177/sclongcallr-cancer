# Demo Single-Cell Dataset 🧬

This folder contains a synthetic benchmark dataset modeling a single-cell long-read cancer sequencing experiment (e.g., PacBio MAS-Seq / 10x Chromium).

### Files
1. **`demo_singlecell_tumor.sam`**:
   - Contains 37 full-length spliced reads covering chromosomes `chr17` and `chr1`.
   - Each read has standard 10x single-cell tags:
     - `CB:Z:<cell_barcode>` (Cell Barcode)
     - `UB:Z:<umi>` (Unique Molecular Identifier)
     - `ts:A:+` (Inferred transcript strand)
   - Features simulated biological events:
     - **Linked Heterozygous SNPs**: Loci at `chr17:1000` (C/T) and `chr17:1100` (G/A).
     - **Immune Cells (`Infiltrating_CD8_T`)**: Normal balanced 50:50 diploid expression of Haplotype 1 (`C...G`) and Haplotype 2 (`T...A`).
     - **Tumor Cells (`Malignant_Subclone_A`)**: Severe allelic imbalance / Loss of Heterozygosity (LOH) retaining Haplotype 1 and silencing/deleting Haplotype 2.
     - **Novel Cancer Splice Junction**: Spliced CIGAR `50M200N100M` creating a tumor-specific cryptic exon-exon junction.

2. **`demo_clusters.csv`**:
   - Maps cell barcodes (`CELL_TUMOR_001`, `CELL_IMMUNE_001`, etc.) to phenotypic cell clusters.
