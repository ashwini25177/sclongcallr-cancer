"""
Cancer Splicing Neoantigen Discovery Engine.
Detects novel unannotated splice junctions linked to phased haplotypes,
translates junction-spanning cDNA sequences into open reading frames (ORFs),
and predicts candidate neo-epitopes for cancer immunotherapies.
"""

from typing import List, Dict, Any, Optional, Tuple, Set
from collections import defaultdict
from ..bam_parser import ScRead
from ..phasing import PhasedBlock

CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}


class SpliceJunctionEvent:
    def __init__(self, chrom: str, donor: int, acceptor: int, is_annotated: bool = False):
        self.chrom = chrom
        self.donor = donor          # 1-based donor position
        self.acceptor = acceptor    # 1-based acceptor position
        self.is_annotated = is_annotated
        self.h1_reads: List[str] = []
        self.h2_reads: List[str] = []
        self.cells: Set[str] = set()

    @property
    def total_reads(self) -> int:
        return len(self.h1_reads) + len(self.h2_reads)

    @property
    def junction_id(self) -> str:
        return f"{self.chrom}:{self.donor}-{self.acceptor}"


class NeoantigenPredictor:
    def __init__(self, epitope_lengths: List[int] = [9, 10]):
        self.epitope_lengths = epitope_lengths

    def translate_dna(self, dna_seq: str) -> str:
        """Translates cDNA sequence into amino acid peptide string."""
        dna_seq = dna_seq.upper().replace("U", "T")
        pep = []
        for i in range(0, len(dna_seq) - 2, 3):
            codon = dna_seq[i:i+3]
            aa = CODON_TABLE.get(codon, "X")
            if aa == '_':  # Stop codon
                break
            pep.append(aa)
        return "".join(pep)

    def extract_novel_junctions(
        self,
        reads: List[ScRead],
        phased_blocks: List[PhasedBlock],
        annotated_junctions: Optional[Set[Tuple[str, int, int]]] = None
    ) -> List[SpliceJunctionEvent]:
        """
        Discovers all splice junctions in reads, tracks their phased haplotypes,
        and marks whether they are novel (unannotated).
        """
        read_to_hap = {}
        for blk in phased_blocks:
            for qname, hap in blk.read_haplotags.items():
                read_to_hap[qname] = hap

        junction_dict: Dict[Tuple[str, int, int], SpliceJunctionEvent] = {}

        for read in reads:
            hap = read_to_hap.get(read.qname, "unphased")
            for donor, acceptor in read.junctions:
                key = (read.chrom, donor, acceptor)
                if key not in junction_dict:
                    is_known = (annotated_junctions is not None) and (key in annotated_junctions)
                    junction_dict[key] = SpliceJunctionEvent(read.chrom, donor, acceptor, is_annotated=is_known)

                event = junction_dict[key]
                if hap == "H1":
                    event.h1_reads.append(read.qname)
                elif hap == "H2":
                    event.h2_reads.append(read.qname)
                if read.cell_barcode:
                    event.cells.add(read.cell_barcode)

        return list(junction_dict.values())

    def predict_neoantigen_peptides(
        self,
        junction_event: SpliceJunctionEvent,
        upstream_flank_dna: str,
        downstream_flank_dna: str
    ) -> List[Dict[str, Any]]:
        """
        Creates junction cDNA sequence: upstream flank + downstream flank,
        translates into 3 reading frames, and finds 9-mer/10-mer peptides spanning the junction.
        """
        junction_cdna = upstream_flank_dna + downstream_flank_dna
        candidates = []

        for frame in range(3):
            sub_seq = junction_cdna[frame:]
            translated_prot = self.translate_dna(sub_seq)

            junction_aa_idx = len(upstream_flank_dna) // 3
            
            for k in self.epitope_lengths:
                # Find all k-mers spanning the junction amino acid index
                start_min = max(0, junction_aa_idx - k + 1)
                start_max = min(junction_aa_idx, len(translated_prot) - k)

                for st in range(start_min, start_max + 1):
                    peptide = translated_prot[st:st+k]
                    if len(peptide) == k and "X" not in peptide:
                        candidates.append({
                            "junction_id": junction_event.junction_id,
                            "chrom": junction_event.chrom,
                            "peptide": peptide,
                            "length": k,
                            "reading_frame": frame,
                            "is_novel_junction": not junction_event.is_annotated,
                            "supporting_h1": len(junction_event.h1_reads),
                            "supporting_h2": len(junction_event.h2_reads),
                            "num_cells": len(junction_event.cells)
                        })

        return candidates
