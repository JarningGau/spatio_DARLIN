import os
import sys

sample = sys.argv[1] # 'E14.5_embryo_dox_E10.5'
locus = sys.argv[2] # 'CA'
image_path = sys.argv[3] # '../data/BMKS3000/image'
threads = sys.argv[4]

# This is the config file for BMK pipeline
config = f'''
### Global parameters
## data Data path
FQ1	cutadapt/{sample}_{locus}_R1.trimmed.fastq.gz
FQ2	cutadapt/{sample}_{locus}_R2.trimmed.fastq.gz

## Flu info file Fluorescence decoding file path
FLU	{image_path}/{sample}-HE.txt

## AllheStat.py & CellSplit Tissue identification and cell segmentation
#HE staining image path, whether to identify internal blank in the tissue (1=identify, 0=ignore, default 1)
HE	{image_path}/{sample}-HE.tif
INSIDE	0
GRAY	200
#Whether to do cell segmentation and the fluorescence image path, color channel (select 0, 1, or 2 according to actual channel used)
CellSplit	False
fluorescence	{image_path}/{sample}-FL.tif
fluorescence_channl	   2 
FLGRAY	15
#cells_npy	/path/to/cells/npyfile
#Cell segmentation parameters
#YAML	/path/to/cell_split/parameter/file
#enhance	1

## Ref genome Reference genome version information, STAR index directory, gff/gtf file and features file path
#GenomeVer	xxxx
# INDEX	/mnt/d/resource/star/mm10-2020-A/
# GFF		/mnt/d/resource/star/mm10-2020-A/genes.gtf
FEATURE	./BST_output/{sample}_{locus}/02.Umi2Gene/features.tsv

## output    Output results and output file prefix
OUTDIR	./BST_output/{sample}_{locus}
PREFIX	out

### Local parameters
## fastq2BcUmi    Barcode version type (usually V2 version) and barcode identification threads
BCType	V2
BCThreads	{threads}

## Umi2Gene       SART alignment parameters
#Sjdboverhang	100
#STARThreads	    8

## ENV      Path to python and Rscript; if not provided, use the version in the system environment (to not provide, comment out the following parameters)
##PYTHON	/path/to/python/dir/
##Rscript	/path/to/Rscript/dir/
'''

if not os.path.exists('BST_config'):
    os.makedirs('BST_config')

if not os.path.exists(f'./BST_output/{sample}_{locus}'):
    os.makedirs(f'./BST_output/{sample}_{locus}')

with open(f'BST_config/{sample}_{locus}.config.txt', 'w') as f:
    f.write(config)