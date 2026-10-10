import pytest
from sclongcallr.cancer.neoantigen import NeoantigenPredictor, SpliceJunctionEvent

def test_neoantigen_translation():
    predictor = NeoantigenPredictor(epitope_lengths=[9])
    
    # Test DNA translation
    dna = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
    # ATG (M) GCC (A) ATT (I) GTA (V) ATG (M) GGC (G) CGC (R) TGA (stop)
    prot = predictor.translate_dna(dna)
    assert prot.startswith("MAIVMGR")

    # Test junction neoantigen prediction
    event = SpliceJunctionEvent(chrom="chr17", donor=5000, acceptor=6000, is_annotated=False)
    event.h1_reads = ["read1", "read2"]
    
    # Upstream exon cDNA + Downstream exon cDNA (needs >= 27 nt for 9-mer peptide)
    upstream = "ATGGCCATTGTAAAA"   # 15 nt -> 5 AA (M A I V K)
    downstream = "ATGGGCCGCTGAAAA" # 15 nt -> 5 AA (M G R _ K)
    
    candidates = predictor.predict_neoantigen_peptides(event, upstream, downstream)
    assert len(candidates) > 0
    assert candidates[0]["is_novel_junction"] is True
    assert candidates[0]["supporting_h1"] == 2

