## Col1a1 Array (CA)
snakemake -s ../../snakefiles/BMKS3000.smk \
--configfile config-CA.yaml \
--cores 2 \
--use-conda

## Rosa Array (RA)
snakemake -s ../../snakefiles/BMKS3000.smk \
--configfile config-RA.yaml \
--cores 2 \
--use-conda

## Tigre Array (TA)
snakemake -s ../../snakefiles/BMKS3000.smk \
--configfile config-TA.yaml \
--cores 2 \
--use-conda
