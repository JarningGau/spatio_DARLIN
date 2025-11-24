# BSTMatrix v2.4 Workflow Instructions

## 1. Install BSTMatrix

```bash
cd /path/to/tools
wget http://www.bmkmanu.com/wp-content/uploads/2024/07/BSTMatrix_v2.4.f.1.zip
unzip BSTMatrix_v2.4.f.1.zip
## conda env for BSTMatrix
cd BSTMatrix_v2.4.f.1
conda env create -n BST-env -f environment.yaml

export PATH=/path/to/tools/BSTMatrix_v2.4.f.1:$PATH
```

## 2. Input Data Preparation

1. **Sequencing data**: Paired-end sequencing FASTQ data

2. **Reference genome data**:
   - Genome sequence file
   - GTF file (column 3 must contain exon)
   - GFF file (optional) (column 3 must contain gene, exon)

3. **features.tsv file**: Can be generated from GTF file, reference command:
   ```bash
   perl /path/to/tools/BSTMatrix_v2.4.f.1/tools/features_generate.pl -i xxx.gtf -o features.tsv
   ```

4. **STAR genome index files**: Can be generated from genome sequence file and GTF file, reference command:
   ```bash
   STAR --runThreadN 8 --runMode genomeGenerate --genomeDir star/ --genomeFastaFiles genome.fa --sjdbGTFfile gene.gtf
   ```

5. **Fluorescence decoding file and HE image file**

6. **Fluorescence image file** (at least one of HE and fluorescence images must be provided)

## 3. Configuration File Writing

Configuration file:

```bash
## fq sequencing data file paths, supports .gz format
FQ1     /path/to/read_1.fq.gz
FQ2     /path/to/read_2.fq.gz

## Flu info file - fluorescence decoding file path
FLU     /path/to/flu_info.txt

## AllheStat.py - tissue identification parameters
HE	/path/to/HE.tif               # Bright-field stained image, at least one of HE and fluorescence images must be provided
#INSIDE	1                      # Whether to identify blank areas inside tissue, 0=no, 1=yes
#GRAY   200                    # Tissue image identification grayscale threshold, default is automatic detection

## CellSplit - whether to perform cell segmentation and fluorescence image path, color channel
CellSplit	True                    # Whether to perform cell segmentation analysis, T/True means yes, otherwise no
fluorescence	/path/to/fluorescence.tiff     # Tissue fluorescence image, optional
fluorescence_channl	   0          # Image color channel, default 0
#FLGRAY   15                      # Fluorescence image identification grayscale threshold, default is automatic detection
#cells_npy	/path/to/cells/npyfile     # Existing cell segmentation npy result file, if provided, use this file for analysis

## Cell segmentation parameters
#YAML   /path/to/cell_split/parameter/file  # Cell segmentation parameter file, optional
#enhance        1                     # Choose 1 for dense cells, choose 0 for distinct nuclear particles

## Reference genome STAR index directory and gff/gtf file paths
GenomeVer		xxx                # Genome version information, used in reports
INDEX   /path/to/STAR/index/dir/
GFF     /path/to/ref/gene/gff3/file    # (gtf file can also be used)

## Reference genome features.tsv file path
FEATURE     /path/to/features.tsv

## Output directory and output file prefix
OUTDIR  /path/to/result/dir/
PREFIX  outfile-prefix

### Program parameters

## fastq2BcUmi
BCType			V2       # Barcode version type (generally V2 version)
BCThreads       8         # Number of threads

## Umi2Gene
Sjdboverhang    100        # -sjdboverhang parameter value used in STAR indexing, default 100
STARThreads     8          # STAR alignment thread number

## ENV - paths to python and Rscript, if not provided, use versions in system environment (if not provided, please comment out the following parameters)
PYTHON  /path/to/python/dir/
Rscript   /path/to/Rscript/dir/
```

**Note**: To perform cell segmentation analysis, `CellSplit` in the configuration file must be set to `True`. For plant cell segmentation, fluorescence images are not required. Only provide the bright-field image in the HE parameter and comment out the fluorescence image parameter `fluorescence`.

## 4. Workflow Execution

### 1) Workflow Description

The workflow is divided into 8 steps, as shown below:

- **Step 1**: Run `fastq2BcUmi` to identify barcode and UMI in FASTQ data
- **Step 2**: Run `Umi2Gene` to align reads to the reference genome and obtain gene information for each UMI
- **Step 3**: Run `LinkBcChip` to identify barcode information in fluorescence data and map to chip positions
- **Step 4**: Run `MatrixMake` to obtain gene expression matrix
- **Step 5**: Run `AllheStat` to process HE images
- **Step 6**: Run `cluster.R` for clustering analysis
- **Step 7**: Run `CellSplit` for cell segmentation analysis
- **Step 8**: Run `WebReport` to generate web-based report

### 2) Workflow Parameters

- `-c config.txt`: Data configuration file
- `-s`: Step selection, 0 means run all steps 1-8, or select individual steps to run separately, multiple steps separated by commas

**Note**: When selecting 0, if cell segmentation analysis is to be performed, the `CellSplit` parameter in the configuration file must be set to `True`.

**Note**: For plant cell segmentation, fluorescence images are not required. Only provide the bright-field image in the HE parameter of the configuration file, and set `CellSplit` to `True` to perform analysis.

### 3) Reference Commands

```bash
./BSTMatrix -c config.txt -s 0

./BSTMatrix -c config.txt -s 1,2,3,4,5,6,7,8

./BSTMatrix -c config.txt -s 1,3
```

## 5. Result File Directory Structure

### 1) Directory Structure and Result Description

```
outdir/
├── 01.fastq2BcUmi                          Step 1 - sequence barcode and UMI detection results directory
│   ├── xxx.bc_dist                             Statistics file for different barcode detections
│   ├── xxx.bc_stat                             Statistics file for different barcode detections
│   ├── xxx.bc_umi_read.tsv                     Statistics file for barcode type, corresponding UMI and read counts
│   ├── xxx.bc_umi_read.tsv.id                  Barcode type, corresponding UMI and read ID file
│   ├── xxx.filter                              Reads information file without complete barcode identification
│   ├── xxx.full_stat                           File with read counts and UMI counts corresponding to barcode types
│   ├── xxx.id_map                              ID number correspondence file
│   ├── xxx.qual.stat                           Reads statistics file
│   ├── xxx.select_id                           Reads ID file with complete barcode and UMI identification
│   ├── xxx.stat                                Barcode detection statistics file
│   ├── xxx.umi                                 File with barcode type and UMI corresponding to reads
│   └── xxx.umi_cor.info                        UMI correction information file
│
├── 02.Umi2Gene                            Step 2 - gene expression information results directory
│   ├── xxxAligned.out.bam                     STAR alignment result BAM file
│   ├── xxx.cut0.fq                             Selected read2 sequence file for alignment
│   ├── xxxLog.final.out                        STAR alignment statistics output file
│   ├── xxxLog.out                              STAR alignment log file
│   ├── xxxLog.progress.out                     STAR alignment log file
│   ├── xxx.map2gene                            Reads information file aligned to genes
│   ├── xxxSJ.out.tab                            STAR alignment splice junction information file
│   ├── xxx_STARtmp                             STAR alignment temporary files directory
│   ├── xxx.stat                                Preliminary alignment information statistics file
│   ├── xxx.total.stat                          Alignment statistics file
│   └── xxx.umi_gene.tsv                        File with barcode, corresponding UMI and genes
│
├── 03.LinkBcChip                            Step 3 - barcode spatial localization results directory
│   ├── xxx.barcode_pos.tsv                      File with barcode type corresponding to chip positions
│   ├── xxx.barcode.tsv                          File with chip corresponding to barcode type
│   ├── xxx.barcode_umi.tsv                     File with chip corresponding to barcode position and UMI counts
│   ├── xxx.chip_type                            Chip corresponding to chip type
│   ├── xxx.dup_distinct                         Decoding intermediate file for duplicate chip signals
│   ├── xxx.one_distinct                         Decoding intermediate file for multiple barcodes corresponding to one chip
│   └── xxx.null                                File with unrecognized chip position information
│
├── 04.MatrixMake                           Step 4 - expression matrix results directory
│   ├── xxx.matrix.tsv                           Gene expression matrix file
│   ├── xxx.matrix.tsv.filt                     Filtered file
│   ├── xxx.select.bc_umi_read.tsv              File with barcode corresponding to UMI and read counts
│   ├── xxx.select.umi_gene.tsv                 File with barcode corresponding to UMI and genes
│   ├── xxx.select.umi_gene.tsv.filter          Filtered barcode gene information file
│   ├── xxx.sequencing_saturation.stat           Sequencing saturation statistics file
│   └── xxx.sequencing_saturation.png            Sequencing saturation curve plot
│
├── 05.AllheStat                            Step 5 - tissue expression analysis results directory
│   ├── allhe                                   Tissue region information directory
│   │   ├── he_roi_small.png                       HE tissue region identified PNG image file
│   │   ├── he_roi.tif                             HE tissue region identified TIF image file
│   │   ├── roi_heAuto.json                        Tissue region JSON file
│   │   └── stat.txt                               Tissue region statistics file
│   ├── all_level_stat.txt                       Statistics file for different resolution levels
│   ├── BSTViewer_project                       BSTViewer software input data directory
│   │   ├── cell_split                              Cell segmentation data directory
│   │   ├── cluster                                 Empty
│   │   ├── he_roi_small.png                       HE tissue region identified PNG image file
│   │   ├── he.tif                                 HE stained image file
│   │   ├── imgs                                   Empty
│   │   ├── level_matrix                           Expression matrix directory for different resolution levels
│   │   ├── project_setting.json                    BSTViewer project JSON file
│   │   ├── roi_groups                            Tissue and HE image JSON file directory
│   │   └── subdata                               Expression matrix directory for tissue regions at different resolution levels
│   ├── heAuto_level_matrix                     Expression matrix directory for tissue regions at different resolution levels
│   │   └── subdata                               Expression matrix directory for tissue regions at different resolution levels
│   ├── level_matrix                            Expression matrix directory for chip at different resolution levels
│   │   ├── level_1                                Level 1 expression matrix directory
│   │   ├── level_13                               Level 13 expression matrix directory
│   │   ├── level_2                                Level 2 expression matrix directory
│   │   ├── level_3                                Level 3 expression matrix directory
│   │   ├── level_4                                Level 4 expression matrix directory
│   │   ├── level_5                                Level 5 expression matrix directory
│   │   ├── level_6                                Level 6 expression matrix directory
│   │   └── level_7                                Level 7 expression matrix directory
│   ├── stat.txt                                Tissue region analysis statistics file
│   └── umi_plot                               UMI plotting results directory
│       ├── all_umi_count_small.png                 Chip region UMI-count PNG image
│       ├── all_umi_count.tif                       Chip region UMI-count TIF image
│       ├── roi_umi_count_small.png                 Tissue region UMI-count PNG image
│       ├── roi_umi_count.tif                       Tissue region UMI-count TIF image
│       ├── roi_umi_count_white_small.png            Tissue region white background UMI-count PNG image
│       └── roi_umi_count_white.tif                 Tissue region white background UMI-count TIF image
│
├── 06.Cluster                              Step 6 - clustering analysis results directory
│   ├── L13                                    Level 13 clustering results directory
│   │   ├── cluster.csv                            Clustering results file
│   │   ├── L13_cluster_files                      Clustering HTML appendix files directory
│   │   ├── L13_cluster.html                       Clustering HTML format image file
│   │   ├── L13_cluster.pdf                        Clustering PDF format image file
│   │   ├── L13_cluster.png                        Clustering PNG format image file
│   │   ├── L13_umap_clstr.pdf                     Merged PDF format image file
│   │   ├── L13_umap_clstr.png                     Merged PNG format image file
│   │   ├── L13_umap_files                         UMAP HTML appendix files directory
│   │   ├── L13_umap.html                          UMAP HTML format image file
│   │   ├── L13_umap.pdf                           UMAP PDF format image file
│   │   └── L13_umap.png                           UMAP PNG format image file
│   ├── L3                                     Level 3 clustering results directory
│   │   ├── cluster.csv                            Clustering results file
│   │   ├── L3_cluster_files                       Clustering HTML appendix files directory
│   │   ├── L3_cluster.html                        Clustering HTML format image file
│   │   ├── L3_cluster.pdf                         Clustering PDF format image file
│   │   ├── L3_cluster.png                         Clustering PNG format image file
│   │   ├── L3_umap_clstr.pdf                      Merged PDF format image file
│   │   ├── L3_umap_clstr.png                      Merged PNG format image file
│   │   ├── L3_umap_files                         UMAP HTML appendix files directory
│   │   ├── L3_umap.html                           UMAP HTML format image file
│   │   ├── L3_umap.pdf                            UMAP PDF format image file
│   │   └── L3_umap.png                            UMAP PNG format image file
│   ... ...
│   └── L7                                     Level 7 clustering results directory
│       ├── cluster.csv                            Clustering results file
│       ├── L7_cluster_files                       Clustering HTML appendix files directory
│       ├── L7_cluster.html                        Clustering HTML format image file
│       ├── L7_cluster.pdf                         Clustering PDF format image file
│       ├── L7_cluster.png                         Clustering PNG format image file
│       ├── L7_umap_clstr.pdf                      Merged PDF format image file
│       ├── L7_umap_clstr.png                      Merged PNG format image file
│       ├── L7_umap_files                         UMAP HTML appendix files directory
│       ├── L7_umap.html                           UMAP HTML format image file
│       ├── L7_umap.pdf                            UMAP PDF format image file
│       └── L7_umap.png                            UMAP PNG format image file
│
├── 07.CellSplit                            Step 7 - cell segmentation results directory
│   ├── cell_split_result                       Cell segmentation results directory
│   │   ├── 0_0.npy                                Local cell segmentation result
│   │   ├── 0_0_ori.tif                            Local fluorescence image
│   │   ├── 0_0.tif                                Local fluorescence image after cell identification
│   ... ...
│   │   ├── 9500_9500.npy                          Local cell segmentation result
│   │   ├── 9500_9500_ori.tif                      Local fluorescence image
│   │   ├── 9500_9500.tif                          Local fluorescence image after cell identification
│   │   ├── all_barcode_num.txt                    Cell barcode ID correspondence file [important]
│   │   ├── all_outline.tif                        Fluorescence image with cell nucleus boundaries added
│   │   ├── cell_color.tif                         Identified cell image file
│   │   ├── cellConts.json                         Identified cell JSON file
│   │   ├── cells.npy                              Identified cell NPY file
│   │   ├── colors.npy                             Cell and color correspondence file
│   │   ├── conts.tif                              Cell segmentation tissue boundary information
│   │   ├── fluorescence.tif                       Tissue fluorescence image
│   │   ├── nucleus_color.tif                      Identified cell nucleus image file
│   │   ├── nucleusConts.json                      Identified cell nucleus JSON file
│   │   ├── nucleus.npy                            Identified cell nucleus NPY file
│   │   ├── progress.txt                           Progress percentage file
│   │   └── SegtoBarcode.log                       Log file
│   ├── cluster                                 Clustering results directory
│   │   ├── cell_cluster_color_img.tif             Cell segmentation clustering image without legend TIF image file
│   │   ├── cell_cluster_color_outline_img.tif     Cell segmentation clustering image without legend with white cell boundaries TIF image file
│   │   ├── cell_cluster_with_legend_img.png       Cell clustering image with legend PNG image file
│   │   ├── cell_cluster_with_legend_img_small.png Cell clustering image with legend low-resolution PNG image file
│   │   ├── cell_cluster_with_legend_img.tif       Cell clustering image with legend TIF image file
│   │   ├── cluster.csv                            Clustering results
│   │   ├── cluster_cells_num.csv                  Clustering category cell count statistics file
│   │   ├── clusters_colors.npy                    Clustering category and color correspondence results
│   │   ├── colors.npy                             Cell and color correspondence results
│   │   ├── legend.tif                             Clustering legend
│   │   ├── marker_gene.csv                        Marker gene information file
│   │   ├── object.RDS                             Seurat object results from cell segmentation matrix
│   │   ├── UMAP.pdf                               UMAP clustering results PDF image file
│   │   └── UMAP.png                               UMAP clustering results PNG image file
│   ├── images                                  Cell segmentation related image results directory
│   │   ├── fluorescence_cell_split.png            Fluorescence image cell segmentation results PNG image file
│   │   ├── fluorescence_cell_split_small.png      Fluorescence image cell segmentation results low-resolution PNG image file
│   │   ├── fluorescence_cell_split.tif            Fluorescence image cell segmentation results TIF image file
│   │   ├── fluorescence.png                       Tissue fluorescence PNG image file
│   │   ├── fluorescence_small.png                 Tissue fluorescence low-resolution PNG image file
│   │   ├── fluorescence.tif                       Tissue fluorescence TIF image file
│   │   ├── he_cell_split.png                      Tissue HE stained cell segmentation PNG image file
│   │   ├── he_cell_split_small.png                Tissue HE stained cell segmentation low-resolution PNG image file
│   │   ├── he_cell_split.tif                      Tissue HE stained cell segmentation TIF image file
│   │   └── he_hr.tif                              Tissue HE stained TIF image file
│   └── mtx                                    Cell segmentation matrix results directory
│       ├── barcodes.tsv.gz                        Cell barcode file
│       ├── barcodes_pos.tsv.gz                    Cell barcode position file
│       ├── cells_center.txt                       Cell centroid position file
│       ├── cells_center.tif                       Cell centroid image file
│       ├── features.tsv.gz                        Cell features file
│       ├── matrix.mtx.gz                          Cell matrix file
│       └── stat.xls                               Cell statistics information file
│
├── 08.WebReport                            Step 8 - web-based report results directory
│   ├── src                                     Web-based report src directory
│   ├── xxx.filelist                             Related file information file used for generating web-based report
│   ├── xxx.stat.xls                             Analysis results statistics information file
│   ├── xxx.rs_stat.xls                          Analysis results statistics information file
│   └── xxx.html                                Web-based report file
│
└── xxx                                      Original expression matrix results directory
    ├── barcode_pos.tsv                        Barcode and corresponding chip position file
    ├── barcode.tsv                            Barcode file
    ├── bc_umi_read.tsv.gz                     Barcode corresponding to UMI and read counts file
    ├── features.tsv                           Features file
    ├── matrix.tsv                             Matrix file
    └── umi_gene.tsv.gz                        Barcode corresponding to UMI and gene file
```

**Note:** The following files serve as key inputs to the spatio_DARLIN pipeline for ensuring consistent cell binning between mRNA and lineage barcodes:

- `07.CellSplit/cell_split_result/all_barcode_num.txt`
- `07.CellSplit/mtx/barcodes_pos.tsv.gz`