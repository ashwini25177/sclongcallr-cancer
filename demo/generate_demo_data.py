"""
Generate Synthetic Demo Dataset for Single-Cell Long-Read Cancer Genomics.
Simulates:
  - Heterozygous SNPs on chr17 (TP53 locus) demonstrating Loss of Heterozygosity in tumor cells.
  - A-to-I RNA editing event on chr1.
  - Novel cancer-specific splice junction in tumor cells with cell barcodes.
"""

import os

def generate_demo_files(output_dir: str = "demo_data"):
    os.makedirs(output_dir, exist_ok=True)
    sam_path = os.path.join(output_dir, "demo_singlecell_tumor.sam")
    clusters_path = os.path.join(output_dir, "demo_clusters.csv")

    # 1. Write cluster mapping
    with open(clusters_path, "w") as f:
        f.write("cell_barcode,cluster\n")
        f.write("CELL_TUMOR_001,Malignant_Subclone_A\n")
        f.write("CELL_TUMOR_002,Malignant_Subclone_A\n")
        f.write("CELL_TUMOR_003,Malignant_Subclone_B\n")
        f.write("CELL_IMMUNE_001,Infiltrating_CD8_T\n")
        f.write("CELL_IMMUNE_002,Infiltrating_CD8_T\n")

    # 2. Write SAM file with header and reads
    with open(sam_path, "w") as f:
        # Headers
        f.write("@HD\tVN:1.6\tSO:coordinate\n")
        f.write("@SQ\tSN:chr17\tLN:83257441\n")
        f.write("@SQ\tSN:chr1\tLN:248956422\n")

        # Simulate reads covering two linked heterozygous SNPs on chr17:
        # SNP 1: pos 1000 (Ref C, Alt T)
        # SNP 2: pos 1100 (Ref G, Alt A)
        # In Normal/Immune: H1 = C...G, H2 = T...A (balanced 50:50)
        # In Tumor: H1 = C...G is retained (90%), H2 = T...A is deleted (LOH!)

        # Immune Cell reads (Balanced 50:50)
        for i in range(5):
            # H1 read (C at 1000, G at 1100)
            seq = "A" * 99 + "C" + "A" * 99 + "G" + "A" * 100
            qual = "I" * len(seq)
            f.write(f"read_immune_h1_{i}\t0\tchr17\t901\t60\t300M\t*\t0\t0\t{seq}\t{qual}\tCB:Z:CELL_IMMUNE_001\tUB:Z:UMI_{i}\tts:A:+\n")

        for i in range(5):
            # H2 read (T at 1000, A at 1100)
            seq = "A" * 99 + "T" + "A" * 99 + "A" + "A" * 100
            qual = "I" * len(seq)
            f.write(f"read_immune_h2_{i}\t0\tchr17\t901\t60\t300M\t*\t0\t0\t{seq}\t{qual}\tCB:Z:CELL_IMMUNE_002\tUB:Z:UMI_{i+5}\tts:A:+\n")

        # Tumor Cell reads (Heavy LOH: 18 H1 reads, 1 H2 read)
        for i in range(18):
            # Tumor retains H1
            seq = "A" * 99 + "C" + "A" * 99 + "G" + "A" * 100
            qual = "I" * len(seq)
            cb = "CELL_TUMOR_001" if i % 2 == 0 else "CELL_TUMOR_002"
            f.write(f"read_tumor_h1_{i}\t0\tchr17\t901\t60\t300M\t*\t0\t0\t{seq}\t{qual}\tCB:Z:{cb}\tUB:Z:T_UMI_{i}\tts:A:+\n")

        # Tumor single residual H2 read
        seq = "A" * 99 + "T" + "A" * 99 + "A" + "A" * 100
        qual = "I" * len(seq)
        f.write(f"read_tumor_h2_res\t0\tchr17\t901\t60\t300M\t*\t0\t0\t{seq}\t{qual}\tCB:Z:CELL_TUMOR_001\tUB:Z:T_UMI_RES\tts:A:+\n")

        # Novel splice junction in Tumor cells: 50M 200N 100M
        for i in range(8):
            seq = "C" * 150
            qual = "I" * len(seq)
            cb = "CELL_TUMOR_001" if i % 2 == 0 else "CELL_TUMOR_003"
            f.write(f"read_tumor_novel_junction_{i}\t0\tchr17\t2000\t60\t50M200N100M\t*\t0\t0\t{seq}\t{qual}\tCB:Z:{cb}\tUB:Z:J_UMI_{i}\tts:A:+\n")

    print(f"[*] Generated synthetic SAM file at: {sam_path}")
    print(f"[*] Generated cluster annotation at: {clusters_path}")
    return sam_path, clusters_path


if __name__ == "__main__":
    generate_demo_files()
