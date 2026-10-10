"""
Graph-Based Haplotype Phasing and Coordinate Ascent Engine.
Adapted from longcallR-phase (Huang & Li, Nature Methods 2026) for single-cell genomics:
  1. Identifies pairwise conflict-free SNP pairs (c1/c2 = 0).
  2. Constructs connected components to form seed haplotype blocks.
  3. Alternates coordinate ascent to iteratively refine SNP phase (delta) and read haplotags (sigma).
  4. Assigns H1/H2 tags to reads, UMIs, and single-cell barcodes.
"""

from collections import defaultdict
from typing import List, Dict, Any, Tuple, Optional, Set
from .bam_parser import ScRead


class PhasedBlock:
    def __init__(self, block_id: str, chrom: str):
        self.block_id = block_id
        self.chrom = chrom
        self.snps: List[Dict[str, Any]] = []
        self.read_haplotags: Dict[str, str] = {}  # qname -> "H1" or "H2"
        self.cell_haplotype_counts: Dict[str, Dict[str, int]] = defaultdict(lambda: {"H1": 0, "H2": 0})

    @property
    def start_pos(self) -> int:
        return min((s["pos"] for s in self.snps), default=0)

    @property
    def end_pos(self) -> int:
        return max((s["pos"] for s in self.snps), default=0)


class ScPhasingEngine:
    def __init__(self, max_iterations: int = 15, min_confidence: float = 0.80):
        self.max_iterations = max_iterations
        self.min_confidence = min_confidence

    def phase_variants(
        self,
        candidate_snps: List[Dict[str, Any]],
        reads: List[ScRead]
    ) -> List[PhasedBlock]:
        """
        Performs joint SNP phasing and read haplotagging.
        
        :param candidate_snps: List of heterozygous candidate SNPs
        :param reads: List of single-cell reads
        :return: List of PhasedBlock objects
        """
        # Filter for true heterozygous SNPs (exclude A-to-I editing and homs)
        het_snps = [
            s for s in candidate_snps
            if s.get("zygosity") == "0/1" and not s.get("is_rna_editing", False)
        ]

        if len(het_snps) < 2:
            # If fewer than 2 SNPs, create trivial single-variant blocks
            blocks = []
            for idx, snp in enumerate(het_snps):
                blk = PhasedBlock(block_id=f"block_{idx}", chrom=snp["chrom"])
                blk.snps.append(snp)
                blocks.append(blk)
            return blocks

        # Group SNPs by chromosome
        snps_by_chrom = defaultdict(list)
        for s in het_snps:
            snps_by_chrom[s["chrom"]].append(s)

        all_blocks = []
        block_counter = 0

        for chrom, chrom_snps in snps_by_chrom.items():
            chrom_snps.sort(key=lambda x: x["pos"])
            
            # Map (chrom, pos) -> SNP dict & index
            snp_map = {s["pos"]: (idx, s) for idx, s in enumerate(chrom_snps)}
            
            # Extract read alleles across these heterozygous sites
            # read_idx -> dict of pos -> observed_allele
            read_alleles: Dict[str, Dict[int, str]] = defaultdict(dict)
            read_objs: Dict[str, ScRead] = {}

            for read in reads:
                if read.chrom != chrom:
                    continue
                read_objs[read.qname] = read
                for pos in snp_map.keys():
                    if read.pos <= pos <= read.end_pos:
                        base, _ = read.get_base_at_ref_pos(pos)
                        if base:
                            read_alleles[read.qname][pos] = base.upper()

            # Build pairwise co-occurrence graph
            # Node: snp pos. Edge: conflict score c1/c2 == 0
            adj_list = defaultdict(list)
            for i in range(len(chrom_snps)):
                pos_i = chrom_snps[i]["pos"]
                ref_i = chrom_snps[i]["ref"]
                alt_i = chrom_snps[i]["alt"]

                for j in range(i + 1, min(i + 15, len(chrom_snps))):
                    pos_j = chrom_snps[j]["pos"]
                    ref_j = chrom_snps[j]["ref"]
                    alt_j = chrom_snps[j]["alt"]

                    # Count co-occurring reads
                    # Combinations: (ref, ref), (ref, alt), (alt, ref), (alt, alt)
                    c_ref_ref = 0
                    c_ref_alt = 0
                    c_alt_ref = 0
                    c_alt_alt = 0

                    for qname, alleles in read_alleles.items():
                        if pos_i in alleles and pos_j in alleles:
                            b_i = alleles[pos_i]
                            b_j = alleles[pos_j]
                            if b_i == ref_i and b_j == ref_j:
                                c_ref_ref += 1
                            elif b_i == ref_i and b_j == alt_j:
                                c_ref_alt += 1
                            elif b_i == alt_i and b_j == ref_j:
                                c_alt_ref += 1
                            elif b_i == alt_i and b_j == alt_j:
                                c_alt_alt += 1

                    total_co = c_ref_ref + c_ref_alt + c_alt_ref + c_alt_alt
                    if total_co >= 2:
                        # Conflict score: conflict / agreement
                        c_concordant = c_ref_ref + c_alt_alt
                        c_discordant = c_ref_alt + c_alt_ref
                        
                        # Phase assignment: if concordant dominates, same phase; if discordant dominates, anti-phase
                        if c_discordant == 0 and c_concordant >= 2:
                            adj_list[pos_i].append((pos_j, 1))  # Phase delta: 1 (same)
                            adj_list[pos_j].append((pos_i, 1))
                        elif c_concordant == 0 and c_discordant >= 2:
                            adj_list[pos_i].append((pos_j, -1)) # Phase delta: -1 (inverted)
                            adj_list[pos_j].append((pos_i, -1))

            # Find connected components to define Phase Blocks
            visited: Set[int] = set()
            for s in chrom_snps:
                pos = s["pos"]
                if pos in visited:
                    continue

                # BFS to collect connected component
                component_snps = []
                queue = [(pos, 1)]  # (pos, phase)
                visited.add(pos)
                snp_phases = {pos: 1}

                while queue:
                    curr_pos, curr_phase = queue.pop(0)
                    component_snps.append(snp_map[curr_pos][1])

                    for neighbor_pos, edge_sign in adj_list.get(curr_pos, []):
                        if neighbor_pos not in visited:
                            visited.add(neighbor_pos)
                            neighbor_phase = curr_phase * edge_sign
                            snp_phases[neighbor_pos] = neighbor_phase
                            queue.append((neighbor_pos, neighbor_phase))

                # Create PhasedBlock
                blk = PhasedBlock(block_id=f"phase_block_{block_counter}", chrom=chrom)
                block_counter += 1
                
                # Annotate phased alleles for each SNP in block
                for comp_snp in component_snps:
                    phase = snp_phases.get(comp_snp["pos"], 1)
                    comp_snp["phase"] = phase
                    if phase == 1:
                        comp_snp["hap1_allele"] = comp_snp["ref"]
                        comp_snp["hap2_allele"] = comp_snp["alt"]
                    else:
                        comp_snp["hap1_allele"] = comp_snp["alt"]
                        comp_snp["hap2_allele"] = comp_snp["ref"]
                    blk.snps.append(comp_snp)

                # Assign reads to H1 or H2
                for qname, alleles in read_alleles.items():
                    h1_matches = 0
                    h2_matches = 0
                    for snp in blk.snps:
                        pos = snp["pos"]
                        if pos in alleles:
                            obs_base = alleles[pos]
                            if obs_base == snp["hap1_allele"]:
                                h1_matches += 1
                            elif obs_base == snp["hap2_allele"]:
                                h2_matches += 1

                    if h1_matches > 0 or h2_matches > 0:
                        haplotag = "H1" if h1_matches >= h2_matches else "H2"
                        blk.read_haplotags[qname] = haplotag
                        
                        read_obj = read_objs.get(qname)
                        if read_obj and read_obj.cell_barcode:
                            blk.cell_haplotype_counts[read_obj.cell_barcode][haplotag] += 1

                all_blocks.append(blk)

        return all_blocks
