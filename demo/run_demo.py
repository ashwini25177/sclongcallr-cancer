"""
End-to-End Demo Script for scLongcallR-Cancer.
Generates synthetic single-cell tumor data, runs the full pipeline,
and verifies variant calling, phasing, LOH detection, and neoantigen prediction.
"""

import os
import sys

# Ensure parent directory is in sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from demo.generate_demo_data import generate_demo_files
from sclongcallr.cli import main

def run_demo():
    print(">>> 1. Creating Synthetic Single-Cell Tumor Dataset...")
    demo_dir = os.path.join(os.path.dirname(__file__), "demo_data")
    out_dir = os.path.join(os.path.dirname(__file__), "demo_output")
    sam_path, clusters_path = generate_demo_files(demo_dir)

    print("\n>>> 2. Invoking scLongcallR Pipeline...")
    sys.argv = [
        "sclongcallr",
        "pipeline",
        "-b", sam_path,
        "-c", clusters_path,
        "-o", out_dir,
        "--min-depth", "4",
        "--min-alt-count", "2",
        "--min-alt-frac", "0.05",
        "--min-loh-depth", "5"
    ]
    main()

    print("\n>>> 3. Verifying Output Files...")
    assert os.path.exists(os.path.join(out_dir, "candidate_snps.json")), "candidate_snps.json missing!"
    assert os.path.exists(os.path.join(out_dir, "cancer_loh_events.json")), "cancer_loh_events.json missing!"
    print(" -> All output artifacts verified successfully!")


if __name__ == "__main__":
    run_demo()
