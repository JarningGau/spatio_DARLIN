cd BMKS3000

## Col1a1 Array (CA)
snakemake -s ../../snakefiles/snakefile_DARLIN_BMKS3000.py \
--configfile config-CA.yaml \
--cores 1 \
--use-conda

## Rosa Array (RA)
snakemake -s ../../snakefiles/snakefile_DARLIN_BMKS3000.py \
--configfile config-RA.yaml \
--cores 1 \
--use-conda

## Tigre Array (TA)
snakemake -s ../../snakefiles/snakefile_DARLIN_BMKS3000.py \
--configfile config-TA.yaml \
--cores 1 \
--use-conda
