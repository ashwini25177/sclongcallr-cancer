"""
Command-Line Interface for scLongcallR-Cancer.
Provides subcommands:
  - call: Deep-learning variant calling from long-read RNA-seq
  - phase: Single-cell haplotype phasing & haplotagging
  - loh: Cancer Loss of Heterozygosity detection
  - neoantigen: Novel splice junction neoantigen prediction
  - pipeline: Complete end-to-end single-cell cancer workflow
"""

import sys
import os
import argparse
import json
from typing import Dict, Any, List

from .bam_parser import ScBamParser
from .pseudobulk import PseudoBulkEngine
from .phasing import ScPhasingEngine
from .cancer.aneuploidy import detect_loh_events
from .cancer.neoantigen import NeoantigenPredictor
from .anndata_export import export_to_anndata


def run_pipeline(args):
    print("=" * 70)
    print("  scLongcallR-Cancer: Single-Cell Long-Read Cancer Genomics Pipeline")
    print("=" * 70)
    print(f"[*] Input Alignment:  {args.bam}")
    print(f"[*] Output Directory: {args.outdir}")
    os.makedirs(args.outdir, exist_ok=True)

    # 1. Parse Reads
    print("\n[Step 1/5] Parsing single-cell long reads (filtering MapQ >= 10)...")
    parser = ScBamParser(args.bam, min_mapq=args.min_mapq)
    reads = list(parser.iter_reads())
    print(f" -> Loaded {len(reads)} valid spliced reads.")

    # Load cluster annotations if available
    cell_to_cluster = {}
    if args.clusters and os.path.exists(args.clusters):
        print(f"[*] Loading cluster metadata from {args.clusters}...")
        with open(args.clusters, "r") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split(",")
                if len(parts) >= 2:
                    cell_to_cluster[parts[0]] = parts[1]

    # 2. Pseudo-bulk Variant Calling
    print("\n[Step 2/5] Running cluster-aware pseudo-bulk Deep Learning variant calling...")
    pb_engine = PseudoBulkEngine(
        min_cluster_depth=args.min_depth,
        min_alt_count=args.min_alt_count,
        min_alt_fraction=args.min_alt_frac
    )
    cluster_pools = pb_engine.pool_reads_by_cluster(reads, cell_to_cluster)
    print(f" -> Grouped into {len(cluster_pools)} cell cluster pool(s).")

    all_snps = []
    for cid, pool_reads in cluster_pools.items():
        snps = pb_engine.find_candidate_snps(pool_reads)
        for s in snps:
            s["cluster"] = cid
        all_snps.extend(snps)

    print(f" -> Identified {len(all_snps)} candidate variant calls.")
    snps_file = os.path.join(args.outdir, "candidate_snps.json")
    with open(snps_file, "w") as f:
        json.dump(all_snps, f, indent=2)
    print(f" -> Saved candidate variants to: {snps_file}")

    # 3. Single-Cell Haplotype Phasing
    print("\n[Step 3/5] Performing graph-based haplotype phasing and read haplotagging...")
    phasing_engine = ScPhasingEngine()
    phased_blocks = phasing_engine.phase_variants(all_snps, reads)
    print(f" -> Constructed {len(phased_blocks)} phased haplotype blocks.")

    # 4. Cancer LOH Analysis
    print("\n[Step 4/5] Scanning for Cancer Loss of Heterozygosity (LOH) events...")
    loh_events = detect_loh_events(phased_blocks, cell_to_cluster=cell_to_cluster, min_reads_per_block=args.min_loh_depth)
    loh_file = os.path.join(args.outdir, "cancer_loh_events.json")

    with open(loh_file, "w") as f:
        json.dump(loh_events, f, indent=2)
    print(f" -> Identified {len(loh_events)} significant LOH events.")
    print(f" -> Saved LOH report to: {loh_file}")

    # 5. Splicing Neoantigen Discovery
    print("\n[Step 5/5] Extracting novel unannotated junctions and predicting neoantigens...")
    neo_predictor = NeoantigenPredictor()
    junction_events = neo_predictor.extract_novel_junctions(reads, phased_blocks)
    novel_junctions = [j for j in junction_events if not j.is_annotated]
    print(f" -> Found {len(novel_junctions)} novel/cryptic splice junctions.")

    # Export to AnnData
    h5ad_out = os.path.join(args.outdir, "sclongcallr_phased.h5ad")
    print(f"\n[*] Exporting phased single-cell matrices to: {h5ad_out}")
    export_to_anndata(phased_blocks, output_path=h5ad_out)

    print("\n" + "=" * 70)
    print("  [SUCCESS] scLongcallR pipeline completed successfully!")
    print("=" * 70)


def main():
    parser = argparse.ArgumentParser(
        prog="sclongcallr",
        description="scLongcallR-Cancer: Single-Cell Long-Read Variant Calling & Cancer Genomics"
    )
    subparsers = parser.add_subparsers(dest="subcommand", help="Available subcommands")

    # Pipeline command
    p_pipe = subparsers.add_parser("pipeline", help="Run end-to-end single-cell cancer workflow")
    p_pipe.add_argument("-b", "--bam", required=True, help="Path to input BAM or SAM file")
    p_pipe.add_argument("-o", "--outdir", default="sclongcallr_out", help="Output directory")
    p_pipe.add_argument("-c", "--clusters", default=None, help="Optional CSV of cell barcodes to cluster IDs")
    p_pipe.add_argument("--min-depth", type=int, default=6, help="Min cluster depth for candidate calling (default: 6)")
    p_pipe.add_argument("--min-alt-count", type=int, default=2, help="Min alternate allele read count (default: 2)")
    p_pipe.add_argument("--min-alt-frac", type=float, default=0.10, help="Min alternate allele fraction (default: 0.10)")
    p_pipe.add_argument("--min-mapq", type=int, default=10, help="Min read mapping quality (default: 10)")
    p_pipe.add_argument("--min-loh-depth", type=int, default=6, help="Min read depth for LOH testing (default: 6)")
    p_pipe.set_defaults(func=run_pipeline)

    # Standalone loh command
    p_loh = subparsers.add_parser("loh", help="Detect Loss of Heterozygosity from phased JSON")
    p_loh.add_argument("-i", "--input", required=True, help="Input phased blocks JSON")
    p_loh.add_argument("-o", "--out", default="loh_report.json", help="Output JSON path")

    args = parser.parse_args()
    if not args.subcommand:
        parser.print_help()
        sys.exit(1)

    if hasattr(args, "func"):
        args.func(args)


if __name__ == "__main__":
    main()
