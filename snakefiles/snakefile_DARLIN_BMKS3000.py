import os
from darlin import help_functions as hf
from darlin.settings import script_dir, CARLIN_dir, QC_dir
#configfile: "config.yaml"  # command line way to set it: --configfile 'path/to/config'
#workdir: config['data_dir'] # set working directory, a command-line way to set it: --directory 'path/to/your/dir'
# config['data_dir']=str(os.getcwd())

##################
## parameters
##################
cfg_type='sc10xV3'
template=config['template']
base_quality_cutoff=config['cutadapt']['base_quality_cutoff']
cutadapt_cores=config['cutadapt']['cores']

CARLIN_dir=hf.update_CARLIN_dir(CARLIN_dir, config['template'])
print("Updated CARLIN_dir:"+ str(CARLIN_dir))

if template.startswith('cCARLIN'):
    locus = 'CA'
elif template.startswith('Rosa'):
    locus = 'RA'
elif template.startswith('Tigre'):
    locus = 'TA'
else:
    raise ValueError(f"The template {template} is not supported, please use one of the following: cCARLIN, Rosa_v2, Tigre_2022_v2.")

SampleList=config['SampleList']
raw_fastq_dir = config['raw_fastq_dir']
image_dir = config['image_dir']
segmentation_dir = config['segmentation_dir']
kernel=config['python_DARLIN_pipeline']['kernel']
reads_cutoff_denoise=config['python_DARLIN_pipeline']['reads_cutoff_denoise']
distance_relative_threshold=config['python_DARLIN_pipeline']['distance_relative_threshold']
distance_absolute_threshold=config['python_DARLIN_pipeline']['distance_absolute_threshold']
slope_threshold=config['python_DARLIN_pipeline']['slope_threshold']
min_reads_per_allele_group=config['python_DARLIN_pipeline']['min_reads_per_allele_group']
perc_reads_per_allele_group=config['python_DARLIN_pipeline']['perc_reads_per_allele_group']
read_fraction_per_clone_spot_cutoff=config['python_DARLIN_pipeline']['read_fraction_per_clone_spot_cutoff']


##################
## helper functions
##################
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
        expand("outs/{sample}_{locus}/all.done", sample=SampleList, locus=locus)

rule extract_DARLIN_seq:
    input:
        r1=lambda wildcards: get_fastq_file(wildcards, "R1"),
        r2=lambda wildcards: get_fastq_file(wildcards, "R2")
    output:
        r1="cutadapt/{sample}_{locus}_R1.trimmed.fastq.gz",
        r2="cutadapt/{sample}_{locus}_R2.trimmed.fastq.gz"
    conda:
        kernel
    shell:
        "python {script_dir}/run_cutadapt.py {template} {cutadapt_cores} {raw_fastq_dir} ./ {wildcards.sample}_{wildcards.locus} {base_quality_cutoff}"

rule write_BMK_config:
    input:
        r1="cutadapt/{sample}_{locus}_R1.trimmed.fastq.gz",
        r2="cutadapt/{sample}_{locus}_R2.trimmed.fastq.gz"
    output:
        config="BST_config/{sample}_{locus}.config.txt"
    shell:
        "python {script_dir}/write_BMK_config.py {wildcards.sample} {wildcards.locus} {image_dir}"

rule extract_spatial_barcode:
    input:
        config="BST_config/{sample}_{locus}.config.txt"
    output:
        select_id="BST_output/{sample}_{locus}/01.fastq2BcUmi/out.select_id",
        parsed_barcodes="BST_output/{sample}_{locus}/01.fastq2BcUmi/out.bc_umi_read.tsv.id",
    conda:
        "BST-env"
    shell:
        "BSTMatrix -c {input.config} -s 1"

rule run_DARLIN_pipeline:
    input:
        select_id="BST_output/{sample}_{locus}/01.fastq2BcUmi/out.select_id",
        parsed_barcodes="BST_output/{sample}_{locus}/01.fastq2BcUmi/out.bc_umi_read.tsv.id",
        r2="cutadapt/{sample}_{locus}_R2.trimmed.fastq.gz"
    output:
        BST_output=directory("BST_output/{sample}_{locus}/02.Umi2Gene"),
        feat='BST_output/{sample}_{locus}/02.Umi2Gene/features.tsv',
        umi2gene='BST_output/{sample}_{locus}/02.Umi2Gene/out.umi_gene.tsv',
        notebook="BST_output/{sample}_{locus}/QC_BMKS3000.ipynb"
    shell:
        "papermill {QC_dir}/BMKS3000.ipynb -k {kernel} {output.notebook} "
        "-p select_id_file {input.select_id} "
        "-p parsed_barcodes_file {input.parsed_barcodes} "
        "-p features_file {output.feat} "
        "-p umi2gene_file {output.umi2gene} "
        "-p darlin_reads {input.r2} "
        "-p sample {wildcards.sample} "
        "-p output_dir {output.BST_output} "
        "-p reads_cutoff_denoise {reads_cutoff_denoise} "
        "-p distance_relative_threshold {distance_relative_threshold} "
        "-p distance_absolute_threshold {distance_absolute_threshold} "
        "-p slope_threshold {slope_threshold} "
        "-p min_reads_per_allele_group {min_reads_per_allele_group} "
        "-p perc_reads_per_allele_group {perc_reads_per_allele_group} "
        "-p read_fraction_per_clone_spot_cutoff {read_fraction_per_clone_spot_cutoff} "


rule prepare_MATLAB_input:
    input:
        feat='BST_output/{sample}_{locus}/02.Umi2Gene/features.tsv'
    output:
        r1="slim_fastq/{sample}_{locus}_R1.fastq.gz",
        r2="slim_fastq/{sample}_{locus}_R2.fastq.gz"
    shell:
        "python {script_dir}/prepare_matlab_input_BMKS3000.py {wildcards.sample} {wildcards.locus} "
        "--called-bc-file {input.feat} "
        "--r1-file {output.r1} "
        "--r2-file {output.r2} "
        "--whitelist {CARLIN_dir}/cfg/10xV3_barcodes.txt.gz"

rule allele_calling:
    input:
        r1="slim_fastq/{sample}_{locus}_R1.fastq.gz",
        r2="slim_fastq/{sample}_{locus}_R2.fastq.gz"
    output:
        seq="DARLIN/"+config['template']+"_cutoff_override_1/{sample}_{locus}/Actaul_CARLIN_seq.txt",
        anno="DARLIN/"+config['template']+"_cutoff_override_1/{sample}_{locus}/AlleleAnnotations.txt"
    shell:
        """
        read_cutoff_UMI_override=1
        read_cutoff_CB_override=1
        input_dir=$(pwd)/slim_fastq/
        output_dir_path=$(pwd)/DARLIN/{template}_cutoff_override_1
        mkdir -p $output_dir_path/{wildcards.sample}_{wildcards.locus}
        bash {script_dir}/run_CARLIN.sh {CARLIN_dir} $input_dir $output_dir_path {wildcards.sample}_{wildcards.locus} {cfg_type} {template} $read_cutoff_UMI_override $read_cutoff_CB_override
        """

rule update_features:
    input:
        seq="DARLIN/"+config['template']+"_cutoff_override_1/{sample}_{locus}/Actaul_CARLIN_seq.txt",
        anno="DARLIN/"+config['template']+"_cutoff_override_1/{sample}_{locus}/AlleleAnnotations.txt",
        feat='BST_output/{sample}_{locus}/02.Umi2Gene/features.tsv',
    output:
        feat='BST_output/{sample}_{locus}/02.Umi2Gene/features_allele.tsv',
        done='BST_output/{sample}_{locus}/update_features.done'
    shell:
        "python {script_dir}/update_features_BMKS3000.py {input.seq} {input.anno} {input.feat} {output.feat} && "
        "touch {output.done}"

rule generate_level_matrix:
    input:
        config="BST_config/{sample}_{locus}.config.txt",
        feat='BST_output/{sample}_{locus}/02.Umi2Gene/features.tsv',
        umi2gene='BST_output/{sample}_{locus}/02.Umi2Gene/out.umi_gene.tsv'
    output:
        level_1_mtx_dir=directory('BST_output/{sample}_{locus}/05.AllheStat/level_matrix/level_1'),
        done='BST_output/{sample}_{locus}/generate_matrix.done'
    conda:
        "BST-env"
    shell:
        "BSTMatrix -c {input.config} -s 3,4,5 && "
        "touch {output.done}"

rule group_spots_to_cells:
    input:
        level_1_mtx_dir='BST_output/{sample}_{locus}/05.AllheStat/level_matrix/level_1',
        seg=segmentation_dir+"/{sample}/all_barcode_num.txt",
        pos=segmentation_dir+"/{sample}/barcodes_pos.tsv.gz"
    output:
        mtx=directory('BST_output/{sample}_{locus}/07.CellSplit'),
        mtx_pos='BST_output/{sample}_{locus}/07.CellSplit/barcodes_pos.tsv.gz',
        mtx_mat='BST_output/{sample}_{locus}/07.CellSplit/matrix.mtx.gz',
        mtx_barcodes='BST_output/{sample}_{locus}/07.CellSplit/barcodes.tsv.gz',
        mtx_features='BST_output/{sample}_{locus}/07.CellSplit/features.tsv.gz',
        done='BST_output/{sample}_{locus}/group_spots_to_cells.done'
    conda:
        "BST-env"
    shell:
        "mkdir -p {output.mtx} && "
        "python $(dirname $(command -v BSTMatrix))/cell_split/get_mtx.py -i {input.level_1_mtx_dir} -c {input.seg} -o {output.mtx} && "
        "cp {input.pos} {output.mtx} && "
        "touch {output.done}"

rule collect_data:
    input:
        group_spots_to_cells_done = 'BST_output/{sample}_{locus}/group_spots_to_cells.done',
        generate_matrix_done='BST_output/{sample}_{locus}/generate_matrix.done',
        feat = 'BST_output/{sample}_{locus}/02.Umi2Gene/features_allele.tsv'
    params:
        level_1_mtx_dir = "BST_output/{sample}_{locus}/05.AllheStat/level_matrix/level_1",
        level_2_mtx_dir = "BST_output/{sample}_{locus}/05.AllheStat/level_matrix/level_2",
        level_3_mtx_dir = "BST_output/{sample}_{locus}/05.AllheStat/level_matrix/level_3",
        level_4_mtx_dir = "BST_output/{sample}_{locus}/05.AllheStat/level_matrix/level_4",
        level_5_mtx_dir = "BST_output/{sample}_{locus}/05.AllheStat/level_matrix/level_5",
        level_6_mtx_dir = "BST_output/{sample}_{locus}/05.AllheStat/level_matrix/level_6",
        level_7_mtx_dir = "BST_output/{sample}_{locus}/05.AllheStat/level_matrix/level_7",
        level_9_mtx_dir = "BST_output/{sample}_{locus}/05.AllheStat/level_matrix/level_9",
        level_18_mtx_dir = "BST_output/{sample}_{locus}/05.AllheStat/level_matrix/level_18",
        cellbin_mtx_dir = "BST_output/{sample}_{locus}/07.CellSplit"
    output:
        outs_dir = directory('outs/{sample}_{locus}'),
        done='outs/{sample}_{locus}/all.done'
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
