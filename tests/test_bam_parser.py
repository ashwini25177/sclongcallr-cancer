import pytest
from sclongcallr.bam_parser import ScRead

def test_sc_read_parsing():
    read = ScRead(
        qname="read_test_01",
        chrom="chr17",
        pos=1000,
        mapq=60,
        cigar_str="50M100N50M",
        seq="A" * 100,
        qual=[30] * 100,
        is_reverse=False,
        cell_barcode="AAACCTG-1",
        umi="UMI_01"
    )

    assert read.qname == "read_test_01"
    assert read.cell_barcode == "AAACCTG-1"
    assert read.umi == "UMI_01"
    assert read.pos == 1000
    # End pos: 1000 + 50 + 100 + 50 - 1 = 1199
    assert read.end_pos == 1199

    # Check extracted splice junction: donor = 1000 + 50 - 1 = 1049, acceptor = 1050 + 100 = 1150
    assert len(read.junctions) == 1
    assert read.junctions[0] == (1049, 1150)

    # Check base at ref position
    base, q = read.get_base_at_ref_pos(1010)
    assert base == "A"
    assert q == 30

    # Inside intron (should be None)
    base_intron, _ = read.get_base_at_ref_pos(1100)
    assert base_intron is None
