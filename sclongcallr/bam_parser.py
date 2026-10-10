"""
Single-Cell Long-Read BAM Parser.
Extracts Cell Barcodes (CB:Z), UMIs (UB:Z), alignment CIGAR strings, 
and splice junctions (N CIGAR operations) across individual single cells.
"""

import re
from typing import List, Dict, Tuple, Optional, Any, Iterator

try:
    import pysam
    HAS_PYSAM = True
except ImportError:
    HAS_PYSAM = False


class ScRead:
    """
    Representation of a single-cell long RNA-seq read.
    """
    def __init__(
        self,
        qname: str,
        chrom: str,
        pos: int,
        mapq: int,
        cigar_str: str,
        seq: str,
        qual: List[int],
        is_reverse: bool = False,
        cell_barcode: Optional[str] = None,
        umi: Optional[str] = None,
        tx_strand: str = "+"
    ):
        self.qname = qname
        self.chrom = chrom
        self.pos = pos  # 1-based start position
        self.mapq = mapq
        self.cigar_str = cigar_str
        self.seq = seq
        self.qual = qual
        self.is_reverse = is_reverse
        self.cell_barcode = cell_barcode
        self.umi = umi
        self.tx_strand = tx_strand
        
        # Parsed attributes
        self.cigar_ops = self._parse_cigar(cigar_str)
        self.end_pos = self._calculate_end_pos()
        self.junctions = self._extract_junctions()

    def _parse_cigar(self, cigar_str: str) -> List[Tuple[int, str]]:
        """Parses CIGAR string into list of (length, op) tuples."""
        if not cigar_str or cigar_str == "*":
            return []
        matches = re.findall(r'(\d+)([MIDNSHP=X])', cigar_str)
        return [(int(length), op) for length, op in matches]

    def _calculate_end_pos(self) -> int:
        """Calculates 1-based reference end position of the alignment."""
        curr = self.pos
        for length, op in self.cigar_ops:
            if op in ('M', 'D', 'N', '=', 'X'):
                curr += length
        return curr - 1

    def _extract_junctions(self) -> List[Tuple[int, int]]:
        """
        Extracts splice junctions (intron deletions represented by 'N' in CIGAR).
        Returns list of (donor_pos, acceptor_pos) 1-based coordinates.
        """
        junctions = []
        curr = self.pos
        for length, op in self.cigar_ops:
            if op == 'N':
                donor = curr - 1
                acceptor = curr + length
                junctions.append((donor, acceptor))
                curr += length
            elif op in ('M', 'D', '=', 'X'):
                curr += length
        return junctions

    def get_base_at_ref_pos(self, target_pos: int) -> Tuple[Optional[str], Optional[int]]:
        """
        Retrieves the read base and quality score aligned to a target reference position.
        """
        if target_pos < self.pos or target_pos > self.end_pos:
            return None, None

        curr_ref = self.pos
        curr_read = 0

        for length, op in self.cigar_ops:
            if op in ('M', '=', 'X'):
                if curr_ref <= target_pos < curr_ref + length:
                    offset = target_pos - curr_ref
                    read_idx = curr_read + offset
                    if read_idx < len(self.seq):
                        base = self.seq[read_idx]
                        q = self.qual[read_idx] if read_idx < len(self.qual) else 30
                        return base, q
                    return None, None
                curr_ref += length
                curr_read += length
            elif op == 'I':
                curr_read += length
            elif op in ('D', 'N'):
                if curr_ref <= target_pos < curr_ref + length:
                    # Target position falls inside a deletion or spliced intron
                    return None, None
                curr_ref += length
            elif op == 'S':
                curr_read += length

        return None, None


class ScBamParser:
    """
    Parser for Single-Cell Long-Read BAM/SAM files.
    """
    def __init__(self, file_path: str, min_mapq: int = 10):
        self.file_path = file_path
        self.min_mapq = min_mapq

    def iter_reads(self) -> Iterator[ScRead]:
        """
        Streams ScRead objects from BAM or SAM file.
        """
        if HAS_PYSAM and (self.file_path.endswith(".bam") or self.file_path.endswith(".sam")):
            try:
                with pysam.AlignmentFile(self.file_path, "rb" if self.file_path.endswith(".bam") else "r") as bam:
                    for record in bam:
                        if record.is_unmapped or record.is_secondary or record.is_supplementary:
                            continue
                        if record.mapping_quality < self.min_mapq:
                            continue

                        # Extract tags
                        cb = record.get_tag("CB") if record.has_tag("CB") else None
                        ub = record.get_tag("UB") if record.has_tag("UB") else None
                        ts = record.get_tag("ts") if record.has_tag("ts") else "+"

                        yield ScRead(
                            qname=record.query_name,
                            chrom=record.reference_name,
                            pos=record.reference_start + 1,  # 1-based
                            mapq=record.mapping_quality,
                            cigar_str=record.cigarstring or "*",
                            seq=record.query_sequence or "",
                            qual=list(record.query_qualities) if record.query_qualities is not None else [30] * len(record.query_sequence or ""),
                            is_reverse=record.is_reverse,
                            cell_barcode=cb,
                            umi=ub,
                            tx_strand=str(ts)
                        )
                return
            except Exception:
                pass

        # Native text parser fallback for SAM or text files
        with open(self.file_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("@") or not line.strip():
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 11:
                    continue

                qname = fields[0]
                flag = int(fields[1])
                chrom = fields[2]
                pos = int(fields[3])
                mapq = int(fields[4])
                cigar = fields[5]
                seq = fields[9]
                qual_str = fields[10]

                # Filter unmapped (0x4), secondary (0x100), supplementary (0x800)
                if (flag & 0x4) or (flag & 0x100) or (flag & 0x800):
                    continue
                if mapq < self.min_mapq:
                    continue

                is_reverse = bool(flag & 0x10)
                quals = [ord(c) - 33 for c in qual_str] if qual_str != "*" else [30] * len(seq)

                # Extract optional tags
                cb = None
                ub = None
                ts = "+"
                for tag_field in fields[11:]:
                    if tag_field.startswith("CB:Z:"):
                        cb = tag_field[5:]
                    elif tag_field.startswith("UB:Z:"):
                        ub = tag_field[5:]
                    elif tag_field.startswith("ts:A:"):
                        ts = tag_field[5:]

                yield ScRead(
                    qname=qname,
                    chrom=chrom,
                    pos=pos,
                    mapq=mapq,
                    cigar_str=cigar,
                    seq=seq,
                    qual=quals,
                    is_reverse=is_reverse,
                    cell_barcode=cb,
                    umi=ub,
                    tx_strand=ts
                )
