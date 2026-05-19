import os
from spatio_darlin.settings import script_dir, QC_dir

##################
## parameters
##################
template=config['template']

if template.startswith('cCARLIN'):
    locus = 'CA'
elif template.startswith('Rosa'):
    locus = 'RA'
elif template.startswith('Tigre'):
    locus = 'TA'
else:
    raise ValueError(f"The template {template} is not supported, please use one of the following: cCARLIN, Rosa, Tigre.")

## IO
SampleList=config['SampleList']
raw_fastq_dir = config['raw_fastq_dir']
image_dir = config['image_dir']
segmentation_dir = config['segmentation_dir']
outdir = config.get('outdir', '.').rstrip('/')
## BSTMatrix parameters
BSTMatrix_threads=config['BSTMatrix']['threads']
## Cutadapt parameters
base_quality_cutoff=config['cutadapt']['base_quality_cutoff']
cutadapt_threads=config['cutadapt']['threads']
## QC parameters
LB_error_rate=config['QC']['LB_error_rate']
major_fraction_threshold_molecule=config['QC']['major_fraction_threshold_molecule']
reads_cutoff=config['QC']['reads_cutoff']
slope_cutoff=config['QC']['slope_cutoff']

##################
## helper functions
##################
def o(*parts):
    return os.path.join(outdir, *parts) if outdir != '.' else os.path.join(*parts)

def get_fastq_file(wildcards, read_end):
    """Get fastq file path, supporting both .fq.gz and .fastq.gz extensions."""
    base_path = f"{raw_fastq_dir}/{wildcards.sample}_{wildcards.locus}_{read_end}"
    # Try .fastq.gz first, then .fq.gz
    fastq_path = base_path + ".fastq.gz"
    fq_path = base_path + ".fq.gz"
    if os.path.exists(fastq_path):
        return fastq_path
    elif os.path.exists(fq_path):
        return fq_path
    else:
        # If neither exists, return .fastq.gz as default (will fail with clear error)
        raise FileNotFoundError(f"Neither {fastq_path} nor {fq_path} exists for sample {wildcards.sample} and locus {wildcards.locus}")

##################
## start the rules
################## 
rule all:
    input:
        expand(o("outs", "{sample}_{locus}", "all.done"), sample=SampleList, locus=locus)

rule extract_DARLIN_barcodes:
    input:
        r1=lambda wildcards: get_fastq_file(wildcards, "R1"),
        r2=lambda wildcards: get_fastq_file(wildcards, "R2")
    output:
        r1=o("cutadapt", "{sample}_{locus}_R1.trimmed.fastq.gz"),
        r2=o("cutadapt", "{sample}_{locus}_R2.trimmed.fastq.gz")
    shell:
        "python {script_dir}/run_cutadapt.py {template} {cutadapt_threads} {raw_fastq_dir} {outdir} {wildcards.sample}_{wildcards.locus} {base_quality_cutoff}"

rule write_BMK_config:
    input:
        r1=o("cutadapt", "{sample}_{locus}_R1.trimmed.fastq.gz"),
        r2=o("cutadapt", "{sample}_{locus}_R2.trimmed.fastq.gz")
    output:
        config=o("BST_config", "{sample}_{locus}.config.txt")
    shell:
        "python {script_dir}/write_BMK_config.py {wildcards.sample} {wildcards.locus} {image_dir} {BSTMatrix_threads} {outdir}"

rule extract_spatial_barcodes:
    input:
        config=o("BST_config", "{sample}_{locus}.config.txt")
    output:
        select_id=o("BST_output", "{sample}_{locus}", "01.fastq2BcUmi", "out.select_id"),
        parsed_barcodes=o("BST_output", "{sample}_{locus}", "01.fastq2BcUmi", "out.bc_umi_read.tsv.id"),
    conda:
        "BST-env"
    shell:
        "BSTMatrix -c {input.config} -s 1"

rule run_DARLIN_QC:
    input:
        select_id=o("BST_output", "{sample}_{locus}", "01.fastq2BcUmi", "out.select_id"),
        parsed_barcodes=o("BST_output", "{sample}_{locus}", "01.fastq2BcUmi", "out.bc_umi_read.tsv.id"),
        r2=o("cutadapt", "{sample}_{locus}_R2.trimmed.fastq.gz")
    output:
        BST_output=directory(o("BST_output", "{sample}_{locus}", "02.Umi2Gene")),
        feat=o("BST_output", "{sample}_{locus}", "02.Umi2Gene", "features.tsv"),
        umi2gene=o("BST_output", "{sample}_{locus}", "02.Umi2Gene", "out.umi_gene.tsv"),
        notebook=o("BST_output", "{sample}_{locus}", "QC_BMKS3000.ipynb")
    shell:
        "papermill {QC_dir}/BMKS3000.ipynb {output.notebook} "
        "-p sample {wildcards.sample} "
        "-p select_id_file {input.select_id} "
        "-p parsed_barcodes_file {input.parsed_barcodes} "
        "-p darlin_reads {input.r2} "
        "-p umi2gene_file {output.umi2gene} "
        "-p features_file {output.feat} "
        "-p LB_error_rate {LB_error_rate} "
        "-p major_fraction_threshold_molecule {major_fraction_threshold_molecule} "
        "-p reads_cutoff {reads_cutoff} "
        "-p slope_cutoff {slope_cutoff} "

rule call_allele:
    input:
        feat=o("BST_output", "{sample}_{locus}", "02.Umi2Gene", "features.tsv")
    output:
        feat_allele=o("BST_output", "{sample}_{locus}", "02.Umi2Gene", "features_allele.tsv"),
        feat_annot=o("BST_output", "{sample}_{locus}", "02.Umi2Gene", "features_annotation.tsv")
    params:
        min_bc_len = 20
    shell:
        "python {script_dir}/annotate_allele.py {locus} {params.min_bc_len} {input.feat} {output.feat_allele} {output.feat_annot}"

rule generate_level_matrix:
    input:
        config=o("BST_config", "{sample}_{locus}.config.txt"),
        feat=o("BST_output", "{sample}_{locus}", "02.Umi2Gene", "features.tsv"),
        umi2gene=o("BST_output", "{sample}_{locus}", "02.Umi2Gene", "out.umi_gene.tsv")
    output:
        level_1_mtx_dir=directory(o("BST_output", "{sample}_{locus}", "05.AllheStat", "level_matrix", "level_1")),
        done=o("BST_output", "{sample}_{locus}", "generate_matrix.done")
    conda:
        "BST-env"
    shell:
        "BSTMatrix -c {input.config} -s 3,4,5 && "
        "touch {output.done}"

rule group_spots_to_cells:
    input:
        level_1_mtx_dir=o("BST_output", "{sample}_{locus}", "05.AllheStat", "level_matrix", "level_1"),
        seg=segmentation_dir+"/{sample}/all_barcode_num.txt",
        pos=segmentation_dir+"/{sample}/barcodes_pos.tsv.gz"
    output:
        mtx=directory(o("BST_output", "{sample}_{locus}", "07.CellSplit")),
        mtx_pos=o("BST_output", "{sample}_{locus}", "07.CellSplit", "barcodes_pos.tsv.gz"),
        mtx_mat=o("BST_output", "{sample}_{locus}", "07.CellSplit", "matrix.mtx.gz"),
        mtx_barcodes=o("BST_output", "{sample}_{locus}", "07.CellSplit", "barcodes.tsv.gz"),
        mtx_features=o("BST_output", "{sample}_{locus}", "07.CellSplit", "features.tsv.gz"),
        done=o("BST_output", "{sample}_{locus}", "group_spots_to_cells.done")
    conda:
        "BST-env"
    shell:
        "mkdir -p {output.mtx} && "
        "python $(dirname $(command -v BSTMatrix))/cell_split/get_mtx.py -i {input.level_1_mtx_dir} -c {input.seg} -o {output.mtx} && "
        "cp {input.pos} {output.mtx} && "
        "touch {output.done}"

rule collect_BST_output:
    input:
        group_spots_to_cells_done=o("BST_output", "{sample}_{locus}", "group_spots_to_cells.done"),
        generate_matrix_done=o("BST_output", "{sample}_{locus}", "generate_matrix.done"),
        feat=o("BST_output", "{sample}_{locus}", "02.Umi2Gene", "features_allele.tsv")
    params:
        level_1_mtx_dir=o("BST_output", "{sample}_{locus}", "05.AllheStat", "level_matrix", "level_1"),
        level_2_mtx_dir=o("BST_output", "{sample}_{locus}", "05.AllheStat", "level_matrix", "level_2"),
        level_3_mtx_dir=o("BST_output", "{sample}_{locus}", "05.AllheStat", "level_matrix", "level_3"),
        level_4_mtx_dir=o("BST_output", "{sample}_{locus}", "05.AllheStat", "level_matrix", "level_4"),
        level_5_mtx_dir=o("BST_output", "{sample}_{locus}", "05.AllheStat", "level_matrix", "level_5"),
        level_6_mtx_dir=o("BST_output", "{sample}_{locus}", "05.AllheStat", "level_matrix", "level_6"),
        level_7_mtx_dir=o("BST_output", "{sample}_{locus}", "05.AllheStat", "level_matrix", "level_7"),
        level_9_mtx_dir=o("BST_output", "{sample}_{locus}", "05.AllheStat", "level_matrix", "level_9"),
        level_18_mtx_dir=o("BST_output", "{sample}_{locus}", "05.AllheStat", "level_matrix", "level_18"),
        cellbin_mtx_dir=o("BST_output", "{sample}_{locus}", "07.CellSplit")
    output:
        outs_dir=directory(o("outs", "{sample}_{locus}")),
        done=o("outs", "{sample}_{locus}", "all.done")
    shell:
        "mkdir -p {output.outs_dir} && "
        "cp -r {params.level_1_mtx_dir} {output.outs_dir}/level_1 && cp {input.feat} {output.outs_dir}/level_1/features_allele.tsv && "
        "cp -r {params.level_2_mtx_dir} {output.outs_dir}/level_2 && cp {input.feat} {output.outs_dir}/level_2/features_allele.tsv && "
        "cp -r {params.level_3_mtx_dir} {output.outs_dir}/level_3 && cp {input.feat} {output.outs_dir}/level_3/features_allele.tsv && "
        "cp -r {params.level_4_mtx_dir} {output.outs_dir}/level_4 && cp {input.feat} {output.outs_dir}/level_4/features_allele.tsv && "
        "cp -r {params.level_5_mtx_dir} {output.outs_dir}/level_5 && cp {input.feat} {output.outs_dir}/level_5/features_allele.tsv && "
        "cp -r {params.level_6_mtx_dir} {output.outs_dir}/level_6 && cp {input.feat} {output.outs_dir}/level_6/features_allele.tsv && "
        "cp -r {params.level_7_mtx_dir} {output.outs_dir}/level_7 && cp {input.feat} {output.outs_dir}/level_7/features_allele.tsv && "
        "cp -r {params.level_9_mtx_dir} {output.outs_dir}/level_9 && cp {input.feat} {output.outs_dir}/level_9/features_allele.tsv && "
        "cp -r {params.level_18_mtx_dir} {output.outs_dir}/level_18 && cp {input.feat} {output.outs_dir}/level_18/features_allele.tsv && "
        "cp -r {params.cellbin_mtx_dir} {output.outs_dir}/cellbin && cp {input.feat} {output.outs_dir}/cellbin/features_allele.tsv && "
        "touch {output.done}"
