"""
Cancer genomics modules for scLongcallR.
Includes:
  - Aneuploidy and Loss of Heterozygosity (LOH) modeling.
  - Splicing neoantigen and novel cancer junction prediction.
"""
from .aneuploidy import CancerAneuploidyModel, detect_loh_events
from .neoantigen import NeoantigenPredictor, SpliceJunctionEvent

__all__ = [
    "CancerAneuploidyModel",
    "detect_loh_events",
    "NeoantigenPredictor",
    "SpliceJunctionEvent"
]
