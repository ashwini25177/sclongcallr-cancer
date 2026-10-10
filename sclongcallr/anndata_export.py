"""
AnnData Export Engine for Single-Cell Long-Read Transcriptomics.
Converts phased single-cell haplotype counts and alternative splicing events
into standard Scanpy-compatible AnnData (.h5ad) structures.
"""

from typing import List, Dict, Any, Optional
import numpy as np
import pandas as pd
from .phasing import PhasedBlock

try:
    import anndata
    from scipy.sparse import csr_matrix
    HAS_ANNDATA = True
except ImportError:
    HAS_ANNDATA = False


def export_to_anndata(
    phased_blocks: List[PhasedBlock],
    cell_metadata: Optional[Dict[str, Dict[str, Any]]] = None,
    output_path: Optional[str] = None
) -> Any:
    """
    Constructs an AnnData object containing haplotype-resolved counts per cell per block/gene.
    
    Layers included:
      - 'hap1_counts': Counts from Haplotype 1
      - 'hap2_counts': Counts from Haplotype 2
      - 'allelic_ratio': H1 / (H1 + H2)
    """
    # Collect all unique cell barcodes
    all_cells = set()
    for block in phased_blocks:
        all_cells.update(block.cell_haplotype_counts.keys())
    
    if cell_metadata:
        all_cells.update(cell_metadata.keys())

    cell_list = sorted(list(all_cells))
    block_list = [b.block_id for b in phased_blocks]

    n_cells = len(cell_list)
    n_blocks = len(block_list)

    cell_idx_map = {c: i for i, c in enumerate(cell_list)}
    
    # Initialize count matrices
    h1_matrix = np.zeros((n_cells, n_blocks), dtype=np.float32)
    h2_matrix = np.zeros((n_cells, n_blocks), dtype=np.float32)

    for b_idx, block in enumerate(phased_blocks):
        for cb, counts in block.cell_haplotype_counts.items():
            c_idx = cell_idx_map[cb]
            h1_matrix[c_idx, b_idx] = counts.get("H1", 0)
            h2_matrix[c_idx, b_idx] = counts.get("H2", 0)

    total_matrix = h1_matrix + h2_matrix
    with np.errstate(divide='ignore', invalid='ignore'):
        allele_ratio_matrix = np.where(total_matrix > 0, h1_matrix / total_matrix, 0.5)

    # Construct obs DataFrame
    obs_dict = {
        "cell_barcode": cell_list,
        "total_phased_reads": np.sum(total_matrix, axis=1),
        "h1_reads": np.sum(h1_matrix, axis=1),
        "h2_reads": np.sum(h2_matrix, axis=1)
    }
    
    # Add external metadata (e.g. cluster, cell type, subclone)
    if cell_metadata:
        cluster_list = []
        for c in cell_list:
            cluster_list.append(cell_metadata.get(c, {}).get("cluster", "unknown"))
        obs_dict["cluster"] = cluster_list

    obs_df = pd.DataFrame(obs_dict).set_index("cell_barcode")

    # Construct var DataFrame
    var_dict = {
        "block_id": block_list,
        "chrom": [b.chrom for b in phased_blocks],
        "start": [b.start_pos for b in phased_blocks],
        "end": [b.end_pos for b in phased_blocks],
        "num_snps": [len(b.snps) for b in phased_blocks]
    }
    var_df = pd.DataFrame(var_dict).set_index("block_id")

    if HAS_ANNDATA:
        adata = anndata.AnnData(
            X=csr_matrix(total_matrix),
            obs=obs_df,
            var=var_df,
            layers={
                "hap1_counts": csr_matrix(h1_matrix),
                "hap2_counts": csr_matrix(h2_matrix),
                "allele_ratio": allele_ratio_matrix
            }
        )
        if output_path:
            adata.write_h5ad(output_path)
        return adata
    else:
        # Fallback dictionary if anndata is not installed
        data_bundle = {
            "X": total_matrix,
            "obs": obs_df,
            "var": var_df,
            "layers": {
                "hap1_counts": h1_matrix,
                "hap2_counts": h2_matrix,
                "allele_ratio": allele_ratio_matrix
            }
        }
        if output_path:
            # Save as CSVs or numpy zip
            np.savez_compressed(
                output_path.replace(".h5ad", ".npz"),
                total=total_matrix,
                h1=h1_matrix,
                h2=h2_matrix,
                ratio=allele_ratio_matrix
            )
            obs_df.to_csv(output_path.replace(".h5ad", "_obs.csv"))
            var_df.to_csv(output_path.replace(".h5ad", "_var.csv"))
        return data_bundle
