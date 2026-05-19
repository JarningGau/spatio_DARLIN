#!/usr/bin/env python3
"""
Generate synthetic 10x Genomics v3 FASTQ files from DARLIN clone barcode output.

This script creates R1 and R2 FASTQ files that mimic 10x v3 sequencing data,
using clone barcodes from DARLIN output and randomly sampled cell barcodes
from the 10x v3 whitelist.

Usage:
    python prepare_matlab_input_BMKS3000.py <sample> <locus> [options]

Arguments:
    sample: Sample name (e.g., 'E14.5_embryo_dox_E10.5')
    locus: Locus code - must be one of ['CA', 'TA', 'RA']
    
Options:
    --called-bc-file: Path to DARLIN output features.tsv file
                      (default: BST_output/{sample}_{locus}/02.Umi2Gene/features.tsv)
    --r1-file: Output path for R1 FASTQ file
               (default: ./slim_fastq/{sample}_{locus}_R1.fastq.gz)
    --r2-file: Output path for R2 FASTQ file
               (default: ./slim_fastq/{sample}_{locus}_R2.fastq.gz)
    --whitelist: Path to 10x v3 cell barcode whitelist
                 (default: ./10xV3_barcodes.txt.gz)
"""

import os
import sys
import gzip
import random
import argparse
import pandas as pd


def get_template(locus: str) -> str:
    """
    Map locus code to template name.
    
    Args:
        locus: Locus code ('CA', 'TA', or 'RA')
        
    Returns:
        Template name string
        
    Raises:
        ValueError: If locus is not one of the valid codes
    """
    if locus == "CA":
        return "cCARLIN"
    elif locus == "TA":
        return "Tigre_2022_v2"
    elif locus == "RA":
        return "Rosa_v2"
    else:
        raise ValueError(f"Invalid locus: {locus}. Locus must be one of ['CA', 'TA', 'RA']")


def gen_umi(length: int) -> str:
    """
    Generate a random UMI sequence.
    
    Args:
        length: Length of UMI sequence
        
    Returns:
        Random UMI sequence string
    """
    bases = ['A', 'T', 'C', 'G']
    umi_sequence = ''.join(random.choice(bases) for _ in range(length))
    return umi_sequence


def gen_base_qual_score(length: int) -> str:
    """
    Generate random quality scores.
    
    Args:
        length: Length of quality score string
        
    Returns:
        Quality score string (Phred+33, range 30-40)
    """
    qual_scores = [chr(random.randint(30, 40) + 33) for _ in range(length)]
    return ''.join(qual_scores)


def gen_10xv3_R1(seq_id: str, cell_barcode: str) -> str:
    """
    Generate a single R1 FASTQ read.
    
    Args:
        seq_id: Sequence identifier
        cell_barcode: 16-base cell barcode
        
    Returns:
        FASTQ formatted read string
        
    Raises:
        ValueError: If cell barcode is not 16 bases long
    """
    if len(cell_barcode) != 16:
        raise ValueError("Cell barcode must be 16 bases long")
    umi_sequence = gen_umi(12)
    sequence = cell_barcode + umi_sequence
    qual_score = gen_base_qual_score(28)
    fastq_read = f"@{seq_id}\n{sequence}\n+\n{qual_score}"
    return fastq_read


def gen_10xv3_R1_fastq(cell_barcodes: list, filename: str) -> None:
    """
    Generate a compressed FASTQ file containing R1 reads for each cell barcode.
    
    Args:
        cell_barcodes: List of cell barcodes, each 16 bases long
        filename: Output filename for the compressed FASTQ file
        
    Raises:
        ValueError: If any cell barcode is not 16 bases long
    """
    prefix = filename.split('/')[-1].replace("_R1.fastq.gz", "")
    with gzip.open(filename, 'wt') as f:
        for i, barcode in enumerate(cell_barcodes):
            if len(barcode) != 16:
                raise ValueError(f"Cell barcode at index {i} must be 16 bases long")
            seq_id = f"{prefix}_{i+1}"
            fastq_read = gen_10xv3_R1(seq_id, barcode)
            f.write(fastq_read + '\n')


def gen_10xv3_R2_fastq(sequences: list, filename: str) -> None:
    """
    Generate a compressed FASTQ file containing R2 reads for each sequence.
    
    Args:
        sequences: List of sequences to be included in the R2 reads
        filename: Output filename for the compressed FASTQ file
    """
    prefix = filename.split('/')[-1].replace("_R2.fastq.gz", "")
    with gzip.open(filename, 'wt') as f:
        for i, sequence in enumerate(sequences):
            seq_id = f"{prefix}_{i+1}"
            qual_score = gen_base_qual_score(len(sequence))
            fastq_read = f"@{seq_id}\n{sequence}\n+\n{qual_score}"
            f.write(fastq_read + '\n')


def gen_10xv3_fastq(locus: str, cell_bc: list, clone_bc: list, R1_file: str, R2_file: str) -> None:
    """
    Generate both R1 and R2 FASTQ files for 10x v3 format.
    
    Args:
        locus: Template name (must start with 'cCARLIN', 'Tigre', or 'Rosa')
        cell_bc: List of cell barcodes
        clone_bc: List of clone barcodes
        R1_file: Output path for R1 FASTQ file
        R2_file: Output path for R2 FASTQ file
        
    Raises:
        ValueError: If locus is not recognized
    """
    if locus.startswith("Tigre"):
        seq_5prime = 'GCTCGGTACCTCGCGAA'
        seq_3prime = 'GTCTTGTCGGTGCCT'
    elif locus.startswith("cCARLIN"):
        seq_5prime = 'GAGCTGTACAAGTAAGCGGC'
        seq_3prime = 'AGAATTCTAACTAGAGCTCGCTGATCAGCCT'
    elif locus.startswith("Rosa"):
        seq_5prime = 'ATGTACAAGTAAAGCGGCC'
        seq_3prime = 'GTCTGCTGTGTGCCT'
    else:
        raise ValueError(f"Invalid locus: {locus}. Locus must start with one of ['cCARLIN', 'Tigre', 'Rosa']")
    
    # Add flanking sequences to clone barcodes
    clone_bc = [seq_5prime + i + seq_3prime for i in clone_bc]
    
    # Create output directory if it doesn't exist
    R2_dir = os.path.dirname(R2_file)
    if R2_dir and not os.path.exists(R2_dir):
        os.makedirs(R2_dir)
    
    # Generate both FASTQ files
    gen_10xv3_R1_fastq(cell_bc, R1_file)
    gen_10xv3_R2_fastq(clone_bc, R2_file)


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description='Generate synthetic 10x v3 FASTQ files from DARLIN clone barcodes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument('sample', type=str, help='Sample name (e.g., E14.5_embryo_dox_E10.5)')
    parser.add_argument('locus', type=str, choices=['CA', 'TA', 'RA'], 
                       help='Locus code: CA, TA, or RA')
    parser.add_argument('--called-bc-file', type=str, default=None,
                       help='Path to DARLIN output features.tsv file')
    parser.add_argument('--r1-file', type=str, default=None,
                       help='Output path for R1 FASTQ file')
    parser.add_argument('--r2-file', type=str, default=None,
                       help='Output path for R2 FASTQ file')
    parser.add_argument('--whitelist', type=str, default='./10xV3_barcodes.txt.gz',
                       help='Path to 10x v3 cell barcode whitelist')
    
    args = parser.parse_args()
    
    # Set default paths if not provided
    if args.called_bc_file is None:
        args.called_bc_file = f"BST_output/{args.sample}_{args.locus}/02.Umi2Gene/features.tsv"
    if args.r1_file is None:
        args.r1_file = f"./slim_fastq/{args.sample}_{args.locus}_R1.fastq.gz"
    if args.r2_file is None:
        args.r2_file = f"./slim_fastq/{args.sample}_{args.locus}_R2.fastq.gz"
    
    # Get template name from locus
    template = get_template(args.locus)
    
    # Load 10x v3 cell barcode whitelist
    print(f"Loading 10x v3 cell barcode whitelist from {args.whitelist}...")
    if not os.path.exists(args.whitelist):
        raise FileNotFoundError(f"Whitelist file not found: {args.whitelist}")
    
    with gzip.open(args.whitelist, 'rt') as f:
        all_cell_barcodes = [line.strip() for line in f]
    
    # Load DARLIN output
    print(f"Loading DARLIN output from {args.called_bc_file}...")
    if not os.path.exists(args.called_bc_file):
        raise FileNotFoundError(f"DARLIN output file not found: {args.called_bc_file}")
    
    data_df = pd.read_csv(args.called_bc_file, sep='\t', header=None, 
                         names=['clone_id', 'clone_bc', 'type'])
    n_lines = len(data_df['clone_id'])
    n_uniq_clone_bcs = len(data_df['clone_id'].unique())
    
    print(f'{n_lines} lines. {n_uniq_clone_bcs} unique clone barcodes')
    
    # Extract clone barcodes
    clone_barcodes = data_df['clone_id'].tolist()
    
    # Randomly sample cell barcodes (one per unique clone)
    if n_uniq_clone_bcs > len(all_cell_barcodes):
        raise ValueError(f"Not enough cell barcodes in whitelist. "
                        f"Need {n_uniq_clone_bcs}, but only {len(all_cell_barcodes)} available.")
    
    sampled_cell_barcodes = random.sample(all_cell_barcodes, n_uniq_clone_bcs)
    
    # Generate 10x v3 FASTQ files
    print(f"Generating R1 FASTQ file: {args.r1_file}")
    print(f"Generating R2 FASTQ file: {args.r2_file}")
    
    gen_10xv3_fastq(
        locus=template,
        cell_bc=sampled_cell_barcodes,
        clone_bc=clone_barcodes,
        R1_file=args.r1_file,
        R2_file=args.r2_file
    )
    
    print("FASTQ file generation completed successfully!")


if __name__ == '__main__':
    main()

