# scLongcallR-Cancer 🧬🔬

**Single-Cell Long-Read RNA-seq Variant Calling, Haplotype Phasing, and Allele-Specific Cancer Genomics**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.0+-ee4c2c.svg)](https://pytorch.org/)
[![Scanpy Compatible](https://img.shields.io/badge/Scanpy-AnnData_1.9+-0052cc.svg)](https://scanpy.readthedocs.io/)
[![Nature Methods](https://img.shields.io/badge/Based_on-Nature_Methods_(2026)-008080.svg)](https://doi.org/10.1038/s41592-026-03045-6)

---

## 🌟 Overview

**`scLongcallR-Cancer`** is a unified computational framework extending the breakthroughs of [longcallR](https://github.com/huangnengCSU/longcallR) (*Huang & Li, Nature Methods 2026*) into the domain of **single-cell long-read cancer transcriptomics** (PacBio MAS-Seq/Iso-Seq and Oxford Nanopore cDNA/dRNA).

While bulk RNA sequencing washes out clonal heterogeneity, single-cell long-read sequencing captures full-length transcripts tagged with Cell Barcodes (`CB`) and Unique Molecular Identifiers (`UB`/`UMI`). `scLongcallR-Cancer` addresses the unique challenges of single-cell cancer sequencing:
1. **Severe single-cell dropout**: Overcomes low read depth via cluster-aware **pseudo-bulk pooling**.
2. **Deep Learning variant calling**: Employs a multi-task **ResNet** trained on 7-channel pileup images to accurately classify genotypes while disentangling true somatic/germline SNPs from **ADAR-mediated A-to-I RNA hyper-editing**.
3. **Graph-based haplotype phasing**: Reconstructs chromosome-scale phase blocks and haplotags individual single-cell reads as $H_1$ or $H_2$.
4. **Cancer Loss of Heterozygosity (LOH)**: Pinpoints monoallelic silencing and two-hit deletions of critical tumor suppressors (*TP53*, *BRCA1/2*, *PTEN*).
5. **Splicing neoantigen discovery**: Identifies unannotated cancer-specific splice junctions on phased mutant transcripts and translates them into in-frame immunogenic neo-epitopes.
6. **Scanpy / AnnData native export**: Directly produces `.h5ad` objects with sparse allelic count layers (`layers['hap1_counts']`, `layers['hap2_counts']`).

---

## 🧠 Is Deep Learning (DL) Used in this Project?

**Yes! Deep Learning is a foundational component of scLongcallR-Cancer.**

### The ResNet-50 Pileup Image Computer Vision Model
Rather than relying on hand-crafted heuristic filters, variant calling is formulated as a multi-channel computer vision problem:
* **Input Tensor**: An aligned pileup image of dimensions **$41\text{ bp} \times \text{depth} \times 7\text{ channels}$** around candidate variant loci:
  1. **Base Channel**: Observed nucleotide ($A, C, G, T$).
  2. **Reference Channel**: Reference genome base duplicated across read depth.
  3. **Insertion Channel**: Flags insertion events at each position.
  4. **Base Quality Channel**: Normalized Phred quality scores.
  5. **Mapping Quality Channel**: Aligner confidence scores ($MapQ$).
  6. **Read Strand Channel**: Orientation of the read ($+1$ or $-1$).
  7. **Transcript Strand Channel**: Spliced transcript direction inferred from splice signals.

* **Multi-Task Classification Heads**:
  - **Zygosity Head**: Predicts $0/0$ (homozygous ref), $0/1$ (heterozygous), $1/1$ (homozygous alt), $1/2$ (biallelic), and **A-to-I RNA Editing**.
  - **Genotype Head**: Predicts diploid nucleotide combinations ($A/A, A/C, \dots, T/T$).
  - **Somatic Subclonal Confidence**: Estimates the likelihood of true somatic variation under copy-number altered tumor purity.

---

## 📐 Architecture & Workflow

```
       Single-Cell Long Reads (BAM with CB:Z & UB:Z tags)
                             │
                             ▼
 ┌─────────────────────────────────────────────────────────┐
 │ 1. Single-Cell Parsing & Quality Filter                 │
 │    • Multi-threaded extraction of CB, UB, & splice CIGAR│
 └───────────────────────────┬─────────────────────────────┘
                             │
                             ▼
 ┌─────────────────────────────────────────────────────────┐
 │ 2. Cluster-Aware Pseudo-Bulk Pooling                    │
 │    • Group single-cell reads by transcriptional cluster │
 │    • Generate high-depth pileup tensors (41 × d × 7)    │
 └───────────────────────────┬─────────────────────────────┘
                             │
                             ▼
 ┌─────────────────────────────────────────────────────────┐
 │ 3. Deep Learning Variant Calling (ResNet)               │
 │    • Predict zygosity, genotype, and filter A-to-I      │
 └───────────────────────────┬─────────────────────────────┘
                             │
                             ▼
 ┌─────────────────────────────────────────────────────────┐
 │ 4. Graph-Based Haplotype Phasing & Read Haplotagging    │
 │    • Seed conflict-free SNP pairs (c1/c2 = 0)           │
 │    • Coordinate ascent to tag reads as H1 or H2         │
 └───────────────────────────┬─────────────────────────────┘
                             │
             ┌───────────────┴───────────────┐
             ▼                               ▼
 ┌───────────────────────────┐   ┌───────────────────────────┐
 │ 5. Cancer LOH Detection   │   │ 6. Splicing Neoantigens   │
 │ • Tumor suppressor LOH    │   │ • Unannotated junctions   │
 │ • Beta-binomial testing   │   │ • In-frame ORF peptide    │
 └───────────┬───────────────┘   └───────────┬───────────────┘
             │                               │
             └───────────────┬───────────────┘
                             ▼
 ┌─────────────────────────────────────────────────────────┐
 │ 7. AnnData (.h5ad) Export for Scanpy / Seurat           │
 │    • layers['hap1_counts'] & layers['hap2_counts']      │
 └─────────────────────────────────────────────────────────┘
```

---

## 🚀 Installation

### 1. From Source (Recommended)
```bash
git clone https://github.com/your-username/sclongcallr-cancer.git
cd sclongcallr-cancer
pip install -e .
```

### 2. Full Dependencies (with PyTorch and Scanpy)
```bash
pip install -e ".[all]"
```

---

## ⚡ Quickstart

### Command-Line Interface (CLI)

Run the end-to-end pipeline on your single-cell tumor BAM file:
```bash
sclongcallr pipeline \
    --bam tumor_singlecell.bam \
    --clusters cell_clusters.csv \
    --outdir results/ \
    --min-depth 6 \
    --min-alt-count 2 \
    --min-loh-depth 10
```

### Python API with Scanpy Integration

```python
import scanpy as sc
from sclongcallr import ScBamParser, PseudoBulkEngine, ScPhasingEngine, detect_loh_events, export_to_anndata

# 1. Parse single-cell long reads
parser = ScBamParser("tumor_mas_seq.bam", min_mapq=10)
reads = list(parser.iter_reads())

# 2. Pseudo-bulk variant calling
engine = PseudoBulkEngine(min_cluster_depth=6)
pools = engine.pool_reads_by_cluster(reads, cell_to_cluster={"CELL_1": "Subclone_A"})
snps = engine.find_candidate_snps(reads)

# 3. Graph-based single-cell phasing
phasing = ScPhasingEngine()
phased_blocks = phasing.phase_variants(snps, reads)

# 4. Detect Cancer Loss of Heterozygosity (LOH)
loh_events = detect_loh_events(phased_blocks, min_reads_per_block=10)

# 5. Export to AnnData
adata = export_to_anndata(phased_blocks, output_path="phased_tumor.h5ad")

# 6. Visualize allelic imbalance on UMAP
sc.pl.umap(adata, color=["cluster", "total_phased_reads"])
```

---

## 🧪 Real-World Cancer Applications

### 1. Loss of Heterozygosity (LOH) in Tumor Suppressors
Identify when malignant subclones delete or epigenetically silence the wild-type copy of *TP53*, *RB1*, or *BRCA1*, leaving only the mutated haplotype expressed.

### 2. Neoantigen Prediction from Novel Splicing
Discover unannotated, cancer-specific splice junctions that occur exclusively on the mutant tumor haplotype. Translate the junction sequences across reading frames to design mRNA cancer vaccines.

### 3. Subclonal Therapy Resistance
Track alternative splicing events that remove targeted epitopes (such as *CD19* exon 2 skipping in leukemia resistant to CAR-T, or *AR-V7* in prostate cancer) to individual tumor subclones before clinical relapse occurs.

---

## 📂 Project Structure

```
sclongcallr-cancer/
├── sclongcallr/                # Core Python package
│   ├── __init__.py
│   ├── cli.py                  # Unified Command-Line Interface
│   ├── bam_parser.py           # Single-cell long read extractor (CB/UB/CIGAR)
│   ├── pseudobulk.py           # Cluster-aware pseudo-bulk pooling engine
│   ├── phasing.py              # Graph phasing & coordinate ascent
│   ├── anndata_export.py       # AnnData (.h5ad) export with hap1/hap2 layers
│   ├── models/                 # Deep Learning Module
│   │   ├── __init__.py
│   │   ├── pileup_image.py     # 7-channel tensor generator (41 x d x 7)
│   │   └── resnet_caller.py    # ResNet-50 variant caller & A-to-I detector
│   └── cancer/                 # Cancer Genomics Module
│       ├── __init__.py
│       ├── aneuploidy.py       # Tumor purity & LOH detection engine
│       └── neoantigen.py       # Splicing neoantigen prediction
├── demo/                       # Demo scripts & test data
│   ├── generate_demo_data.py   # Synthetic single-cell tumor SAM generator
│   └── run_demo.py             # End-to-end verification script
├── tests/                      # Pytest automated test suite
│   ├── test_bam_parser.py
│   ├── test_pileup_and_model.py
│   ├── test_phasing_and_loh.py
│   └── test_neoantigen.py
├── pyproject.toml              # Modern package metadata
├── setup.py                    # Backward compatibility setup
├── requirements.txt            # Dependency list
├── LICENSE                     # MIT License
└── README.md                   # Documentation
```

---

## 🔬 Benchmark & Citation

If you use this software in your research, please cite:

1. **Huang, N., Li, H., et al.** *SNP calling, haplotype phasing and allele-specific analysis with long RNA-seq reads.* **Nature Methods** 23, 1962–1967 (2026). DOI: [10.1038/s41592-026-03045-6](https://doi.org/10.1038/s41592-026-03045-6).
2. **scLongcallR-Cancer**: Single-cell long-read cancer transcriptomics framework (2026).
