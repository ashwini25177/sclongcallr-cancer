"""
Cancer Aneuploidy and Loss of Heterozygosity (LOH) Engine.
Models subclonal copy number alterations (CNAs), tumor purity, and detects 
Loss of Heterozygosity on tumor suppressor genes directly from phased single-cell RNA reads.
"""

import math
from collections import defaultdict
from typing import List, Dict, Any, Optional, Tuple


try:
    from scipy.stats import binomtest
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

def _exact_binom_test_pvalue(k: int, n: int, p: float = 0.5) -> float:
    """Exact two-sided binomial test using standard library math.comb."""
    if n == 0:
        return 1.0
    p_k = math.comb(n, k)
    total_ways = sum(math.comb(n, x) for x in range(n + 1) if math.comb(n, x) <= p_k)
    return min(1.0, total_ways / (2 ** n))



class CancerAneuploidyModel:
    """
    Likelihood model for allele frequencies under tumor copy number alterations (CNA).
    """
    def __init__(self, tumor_purity: float = 0.80):
        self.tumor_purity = min(max(tumor_purity, 0.01), 1.0)

    def calculate_expected_baf(
        self,
        major_copies: int,
        minor_copies: int,
        mutant_on_major: bool = True
    ) -> float:
        """
        Calculates expected B-Allele Frequency (BAF) in a tumor sample with stromal contamination.
        Normal stromal cells are assumed to be diploid (1 copy A, 1 copy B).
        """
        c_alt = major_copies if mutant_on_major else minor_copies
        c_ref = minor_copies if mutant_on_major else major_copies
        
        tumor_alt = self.tumor_purity * c_alt
        normal_alt = (1.0 - self.tumor_purity) * 1.0
        
        tumor_total = self.tumor_purity * (c_alt + c_ref)
        normal_total = (1.0 - self.tumor_purity) * 2.0
        
        return (tumor_alt + normal_alt) / (tumor_total + normal_total)


def detect_loh_events(
    phased_blocks: List[PhasedBlock],
    cell_to_cluster: Optional[Dict[str, str]] = None,
    min_reads_per_block: int = 6,
    ai_threshold: float = 0.65,
    alpha: float = 0.01
) -> List[Dict[str, Any]]:
    """
    Detects Loss of Heterozygosity (LOH) and significant allelic imbalance events.
    Supports cluster-resolved testing so stromal infiltration does not dilute tumor LOH.
    """
    loh_events = []

    for block in phased_blocks:
        # Group cell counts by cluster
        cluster_counts: Dict[str, Dict[str, int]] = defaultdict(lambda: {"H1": 0, "H2": 0})
        
        for cb, counts in block.cell_haplotype_counts.items():
            cid = cell_to_cluster.get(cb, "sample_aggregate") if cell_to_cluster else "sample_aggregate"
            cluster_counts[cid]["H1"] += counts.get("H1", 0)
            cluster_counts[cid]["H2"] += counts.get("H2", 0)

        for cid, counts in cluster_counts.items():
            h1 = counts["H1"]
            h2 = counts["H2"]
            total_depth = h1 + h2

            if total_depth < min_reads_per_block:
                continue

            ai_ratio = abs(h1 - h2) / float(total_depth)

            if HAS_SCIPY:
                p_val = binomtest(k=h1, n=total_depth, p=0.5, alternative='two-sided').pvalue
            else:
                p_val = _exact_binom_test_pvalue(k=h1, n=total_depth, p=0.5)

            is_significant_loh = (ai_ratio >= ai_threshold) and (p_val < alpha)

            dominant_haplotype = "H1" if h1 >= h2 else "H2"
            repressed_haplotype = "H2" if dominant_haplotype == "H1" else "H1"

            if is_significant_loh:
                loh_events.append({
                    "block_id": block.block_id,
                    "chrom": block.chrom,
                    "start": block.start_pos,
                    "end": block.end_pos,
                    "cluster": cid,
                    "num_snps": len(block.snps),
                    "h1_counts": h1,
                    "h2_counts": h2,
                    "total_depth": total_depth,
                    "ai_ratio": round(ai_ratio, 4),
                    "dominant_haplotype": dominant_haplotype,
                    "repressed_haplotype": repressed_haplotype,
                    "p_value": float(p_val),
                    "status": "LOH_CONFIRMED"
                })

    return loh_events

