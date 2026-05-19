# Deprecated MATLAB allele-calling workflow

This directory preserves the pre-1.0.0 Snakemake workflow that calls the MATLAB CARLIN pipeline for allele annotation.

**Requirements (not needed for the default pipeline):**

- MATLAB (CLI `matlab` command)
- [Custom_CARLIN](https://github.com/ShouWenWang-Lab/Custom_CARLIN) cloned as `../CARLIN_pipeline/Custom_CARLIN` relative to the repository root

**Snakemake:**

```bash
cd test/test_BMKS3000
snakemake -s ../../legacy/matlab/snakefiles/BMKS3000_matlab.smk \
  --configfile config-CA.yaml --cores 2 --use-conda
```

Or run `legacy/matlab/test/test_bmk_matlab.sh` from `test/test_BMKS3000` (update paths in that script if needed).

**Default workflow:** use [snakefiles/BMKS3000.smk](../../snakefiles/BMKS3000.smk) with [darlin-core](https://github.com/JarningGau/darlin-core) (see main [README.md](../../README.md)).
