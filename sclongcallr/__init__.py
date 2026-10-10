"""
scLongcallR-Cancer: Single-Cell Long-Read RNA-seq Variant Calling, 
Haplotype Phasing, and Allele-Specific Cancer Genomics.
"""

__version__ = "0.1.0"
__author__ = "Ashwini"

from .bam_parser import ScBamParser, ScRead
from .pseudobulk import PseudoBulkEngine
from .phasing import ScPhasingEngine, PhasedBlock
from .cancer.aneuploidy import CancerAneuploidyModel, detect_loh_events
from .cancer.neoantigen import NeoantigenPredictor, SpliceJunctionEvent
from .anndata_export import export_to_anndata

__all__ = [
    "ScBamParser",
    "ScRead",
    "PseudoBulkEngine",
    "ScPhasingEngine",
    "PhasedBlock",
    "CancerAneuploidyModel",
    "detect_loh_events",
    "NeoantigenPredictor",
    "SpliceJunctionEvent",
    "export_to_anndata",
]
