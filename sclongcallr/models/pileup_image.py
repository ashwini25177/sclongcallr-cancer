"""
Pileup Image Generator for Deep Learning RNA-seq Variant Calling.
Generates a 7-channel tensor for a 41bp genomic window around candidate SNP loci.
Channels:
  1. Base (A, C, G, T encoded as 1..4)
  2. Reference base (encoded 1..4)
  3. Insertion indicator
  4. Base Quality score (normalized Phred)
  5. Mapping Quality score (normalized MapQ)
  6. Read Strand (+1 / -1)
  7. Transcript Strand (+1 / -1)
"""

import numpy as np
from typing import List, Dict, Any, Optional

BASE_MAP = {'A': 1, 'C': 2, 'G': 3, 'T': 4, 'N': 0,
            'a': 1, 'c': 2, 'g': 3, 't': 4, 'n': 0}

class PileupImageGenerator:
    def __init__(self, window_size: int = 41):
        """
        Initialize the pileup generator.
        :param window_size: Total flanking window size (default 41 bp, +/- 20 bp).
        """
        assert window_size % 2 == 1, "Window size must be an odd integer"
        self.window_size = window_size
        self.half_window = window_size // 2

    def create_pileup_tensor(
        self,
        chrom: str,
        pos: int,
        ref_seq: str,
        reads_at_locus: List[Dict[str, Any]],
        max_depth: int = 100
    ) -> np.ndarray:
        """
        Constructs a (window_size, depth, 7) feature tensor.
        
        :param chrom: Chromosome name
        :param pos: 1-based genomic position
        :param ref_seq: Reference sequence string of length >= window_size
        :param reads_at_locus: List of read dicts spanning the window
        :param max_depth: Max depth capping for tensor uniform batching
        :return: np.ndarray of shape (window_size, effective_depth, 7)
        """
        num_reads = min(len(reads_at_locus), max_depth)
        if num_reads == 0:
            # Empty tensor placeholder
            return np.zeros((self.window_size, 1, 7), dtype=np.float32)

        tensor = np.zeros((self.window_size, num_reads, 7), dtype=np.float32)

        # 1. Reference channel (Channel index 1)
        # Convert ref_seq to encoded values
        ref_encoded = np.array([BASE_MAP.get(b, 0) / 4.0 for b in ref_seq[:self.window_size]], dtype=np.float32)
        if len(ref_encoded) < self.window_size:
            padded_ref = np.zeros(self.window_size, dtype=np.float32)
            padded_ref[:len(ref_encoded)] = ref_encoded
            ref_encoded = padded_ref

        for r_idx, read in enumerate(reads_at_locus[:num_reads]):
            # Channel 1: Reference (duplicated across reads)
            tensor[:, r_idx, 1] = ref_encoded

            # Channel 4: Mapping Quality (replicated across 41bp)
            mapq_norm = min(read.get("mapq", 60), 60) / 60.0
            tensor[:, r_idx, 4] = mapq_norm

            # Channel 5: Read Strand (+1 or -1)
            is_reverse = read.get("is_reverse", False)
            tensor[:, r_idx, 5] = -1.0 if is_reverse else 1.0

            # Channel 6: Transcript Strand (+1 or -1)
            tx_strand = read.get("tx_strand", "+")
            tensor[:, r_idx, 6] = 1.0 if tx_strand == "+" else -1.0

            # Fill read bases, insertion, and base quality per window offset
            read_bases = read.get("window_bases", "")  # length should match window_size
            read_quals = read.get("window_quals", [])
            has_ins = read.get("window_insertions", [False] * self.window_size)

            for w_idx in range(min(self.window_size, len(read_bases))):
                base_char = read_bases[w_idx]
                tensor[w_idx, r_idx, 0] = BASE_MAP.get(base_char, 0) / 4.0
                tensor[w_idx, r_idx, 2] = 1.0 if has_ins[w_idx] else 0.0
                if w_idx < len(read_quals):
                    tensor[w_idx, r_idx, 3] = min(read_quals[w_idx], 60) / 60.0

        return tensor
