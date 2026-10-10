import pytest
from sclongcallr.phasing import ScPhasingEngine, PhasedBlock
from sclongcallr.cancer.aneuploidy import detect_loh_events
from sclongcallr.bam_parser import ScRead

def test_phasing_and_loh():
    # Simulate candidate SNPs
    snps = [
        {"chrom": "chr17", "pos": 1000, "ref": "C", "alt": "T", "zygosity": "0/1", "is_rna_editing": False},
        {"chrom": "chr17", "pos": 1100, "ref": "G", "alt": "A", "zygosity": "0/1", "is_rna_editing": False}
    ]

    # Simulate reads:
    # 15 reads H1 (C at 1000, G at 1100)
    # 1 read H2 (T at 1000, A at 1100) -> Strong LOH
    reads = []
    for i in range(15):
        seq = "A" * 99 + "C" + "A" * 99 + "G" + "A" * 100
        reads.append(ScRead(
            qname=f"read_h1_{i}",
            chrom="chr17",
            pos=901,
            mapq=60,
            cigar_str="300M",
            seq=seq,
            qual=[30] * len(seq),
            cell_barcode="TUMOR_CELL_1"
        ))

    # Single H2 read
    seq_h2 = "A" * 99 + "T" + "A" * 99 + "A" + "A" * 100
    reads.append(ScRead(
        qname="read_h2_0",
        chrom="chr17",
        pos=901,
        mapq=60,
        cigar_str="300M",
        seq=seq_h2,
        qual=[30] * len(seq_h2),
        cell_barcode="TUMOR_CELL_1"
    ))

    engine = ScPhasingEngine()
    blocks = engine.phase_variants(snps, reads)

    assert len(blocks) == 1
    blk = blocks[0]
    assert len(blk.snps) == 2

    # Run LOH Detection
    loh_events = detect_loh_events(blocks, min_reads_per_block=5, ai_threshold=0.60)
    assert len(loh_events) == 1
    loh = loh_events[0]
    assert loh["dominant_haplotype"] == "H1"
    assert loh["h1_counts"] == 15
    assert loh["h2_counts"] == 1
    assert loh["status"] == "LOH_CONFIRMED"
