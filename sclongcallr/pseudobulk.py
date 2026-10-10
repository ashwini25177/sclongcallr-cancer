"""
Pseudo-Bulk Cluster Pooling Engine.
Overcomes single-cell sequencing dropout by pooling single-cell reads 
across transcriptionally-defined cell clusters, calling high-confidence SNPs, 
and projecting genotypes back to single cells.
"""

from collections import defaultdict
from typing import List, Dict, Any, Optional, Tuple
from .bam_parser import ScRead
from .models.resnet_caller import VariantPredictor
from .models.pileup_image import PileupImageGenerator


class PseudoBulkEngine:
    def __init__(
        self,
        min_cluster_depth: int = 6,
        min_alt_count: int = 2,
        min_alt_fraction: float = 0.10,
        window_size: int = 41
    ):
        self.min_cluster_depth = min_cluster_depth
        self.min_alt_count = min_alt_count
        self.min_alt_fraction = min_alt_fraction
        self.window_size = window_size
        self.pileup_gen = PileupImageGenerator(window_size=window_size)
        self.predictor = VariantPredictor()

    def pool_reads_by_cluster(
        self,
        reads: List[ScRead],
        cell_to_cluster: Optional[Dict[str, str]] = None
    ) -> Dict[str, List[ScRead]]:
        """
        Groups reads by cluster ID. If no cluster mapping is provided,
        groups all reads into a single 'all_cells' pseudo-bulk.
        """
        cluster_map = defaultdict(list)
        for r in reads:
            cluster_id = "unassigned"
            if cell_to_cluster and r.cell_barcode:
                cluster_id = cell_to_cluster.get(r.cell_barcode, "unassigned")
            elif not cell_to_cluster:
                cluster_id = "sample_aggregate"
            cluster_map[cluster_id].append(r)
        return cluster_map

    def find_candidate_snps(
        self,
        cluster_reads: List[ScRead],
        ref_genome: Optional[Dict[str, str]] = None
    ) -> List[Dict[str, Any]]:
        """
        Scans pooled reads within a cluster to find loci meeting candidate thresholds.
        """
        # Count base occurrences per (chrom, pos)
        locus_counts = defaultdict(lambda: defaultdict(int))
        reads_by_locus = defaultdict(list)

        for read in cluster_reads:
            # Sample positions covered by read
            for length, op in read.cigar_ops:
                if op in ('M', '=', 'X'):
                    pass  # Read matches reference

        # Build position pileup
        for read in cluster_reads:
            curr_pos = read.pos
            read_idx = 0
            for length, op in read.cigar_ops:
                if op in ('M', '=', 'X'):
                    for i in range(length):
                        pos = curr_pos + i
                        if read_idx + i < len(read.seq):
                            base = read.seq[read_idx + i].upper()
                            if base in ('A', 'C', 'G', 'T'):
                                locus_counts[(read.chrom, pos)][base] += 1
                                reads_by_locus[(read.chrom, pos)].append(read)
                    curr_pos += length
                    read_idx += length
                elif op in ('D', 'N'):
                    curr_pos += length
                elif op in ('I', 'S'):
                    read_idx += length

        candidate_snps = []
        for (chrom, pos), base_counts in locus_counts.items():
            total_depth = sum(base_counts.values())
            if total_depth < self.min_cluster_depth:
                continue

            # Infer reference base if reference sequence provided, else most frequent
            sorted_bases = sorted(base_counts.items(), key=lambda x: x[1], reverse=True)
            ref_base = sorted_bases[0][0]
            if ref_genome and chrom in ref_genome and pos <= len(ref_genome[chrom]):
                ref_base = ref_genome[chrom][pos - 1].upper()

            alt_counts = {b: c for b, c in base_counts.items() if b != ref_base}
            if not alt_counts:
                continue

            top_alt, top_alt_count = max(alt_counts.items(), key=lambda x: x[1])
            alt_fraction = top_alt_count / float(total_depth)

            if top_alt_count >= self.min_alt_count and alt_fraction >= self.min_alt_fraction:
                # Construct flanking reference window
                ref_window = ref_base * self.window_size
                if ref_genome and chrom in ref_genome:
                    half = self.window_size // 2
                    start_w = max(0, pos - 1 - half)
                    end_w = min(len(ref_genome[chrom]), pos + half)
                    extracted = ref_genome[chrom][start_w:end_w].upper()
                    if len(extracted) == self.window_size:
                        ref_window = extracted

                # Convert reads to format for tensor
                formatted_reads = []
                for r in reads_by_locus[(chrom, pos)][:100]:
                    base, q = r.get_base_at_ref_pos(pos)
                    formatted_reads.append({
                        "mapq": r.mapq,
                        "is_reverse": r.is_reverse,
                        "tx_strand": r.tx_strand,
                        "window_bases": (base or ref_base) * self.window_size,
                        "window_quals": [q or 30] * self.window_size,
                        "window_insertions": [False] * self.window_size
                    })

                tensor = self.pileup_gen.create_pileup_tensor(chrom, pos, ref_window, formatted_reads)
                pred = self.predictor.predict(tensor, ref_base, alt_counts, total_depth)

                candidate_snps.append({
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref_base,
                    "alt": top_alt,
                    "depth": total_depth,
                    "alt_count": top_alt_count,
                    "alt_fraction": alt_fraction,
                    "zygosity": pred["zygosity"],
                    "genotype": pred["genotype"],
                    "confidence": pred["confidence"],
                    "is_rna_editing": pred["is_rna_editing"],
                    "dl_used": pred.get("dl_used", False)
                })

        return candidate_snps
