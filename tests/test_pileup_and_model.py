import pytest
import numpy as np
from sclongcallr.models.pileup_image import PileupImageGenerator
from sclongcallr.models.resnet_caller import VariantPredictor

def test_pileup_tensor_shape():
    generator = PileupImageGenerator(window_size=41)
    ref_seq = "A" * 41
    reads = [
        {
            "mapq": 60,
            "is_reverse": False,
            "tx_strand": "+",
            "window_bases": "A" * 20 + "G" + "A" * 20,
            "window_quals": [30] * 41,
            "window_insertions": [False] * 41
        }
    ]

    tensor = generator.create_pileup_tensor("chr1", 100, ref_seq, reads)
    assert tensor.shape == (41, 1, 7)
    # Channel 1 is reference
    assert tensor[20, 0, 1] > 0
    # Channel 0 is base
    assert tensor[20, 0, 0] > 0

def test_variant_predictor():
    predictor = VariantPredictor()
    tensor = np.zeros((41, 20, 7), dtype=np.float32)

    # 1. Heterozygous variant test
    alt_counts = {"G": 10, "A": 10}
    res = predictor.predict(tensor, ref_base="A", alt_counts=alt_counts, total_depth=20)
    assert res["zygosity"] == "0/1"
    assert res["genotype"] == ("A", "G")
    assert not res["is_rna_editing"]

    # 2. A-to-I editing test
    alt_counts_edit = {"G": 5, "A": 15}
    res_edit = predictor.predict(tensor, ref_base="A", alt_counts=alt_counts_edit, total_depth=20)
    assert res_edit["zygosity"] == "A-to-I"
    assert res_edit["is_rna_editing"]
