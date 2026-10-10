from setuptools import setup, find_packages

setup(
    name="sclongcallr",
    version="0.1.0",
    description="Single-Cell Long-Read RNA Variant Calling, Haplotype Phasing, and Allele-Specific Cancer Genomics",
    author="Ashwini",
    license="MIT",
    packages=find_packages(),
    python_requires=">=3.9",
    install_requires=[
        "numpy>=1.21.0",
        "scipy>=1.7.0",
        "pandas>=1.3.0",
    ],
    extras_require={
        "all": ["torch>=2.0.0", "anndata>=0.8.0", "scanpy>=1.9.0", "pysam>=0.19.0"],
        "dl": ["torch>=2.0.0"],
        "singlecell": ["anndata>=0.8.0", "scanpy>=1.9.0", "pysam>=0.19.0"],
    },
    entry_points={
        "console_scripts": [
            "sclongcallr=sclongcallr.cli:main",
        ],
    },
)
