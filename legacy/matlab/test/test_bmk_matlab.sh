## Col1a1 Array (CA)
snakemake -s ../../legacy/matlab/snakefiles/BMKS3000_matlab.smk \
--configfile config-CA.yaml \
--cores 2 \
--use-conda

## Rosa Array (RA)
snakemake -s ../../legacy/matlab/snakefiles/BMKS3000_matlab.smk \
--configfile config-RA.yaml \
--cores 2 \
--use-conda

## Tigre Array (TA)
snakemake -s ../../legacy/matlab/snakefiles/BMKS3000_matlab.smk \
--configfile config-TA.yaml \
--cores 2 \
--use-conda
