# Changelog

## 1.0.0

- Remove MATLAB from the default workflow; legacy MATLAB pipeline under [legacy/matlab/](../legacy/matlab/)
- Rename Python package `darlin` → `spatio_darlin`
- Replace [darlinpy](https://github.com/JarningGau/darlinpy) with [darlin-core](https://github.com/JarningGau/darlin-core) for allele annotation and cutadapt primers
- Documentation updated for a single Snakemake entry point ([snakefiles/BMKS3000.smk](../snakefiles/BMKS3000.smk))

### Migrating from 0.x

- Reinstall this repository: `pip install -e .` (package name is now `spatio-darlin`)
- Install `darlin-core` from GitHub: `pip install git+https://github.com/JarningGau/darlin-core.git` (not available on PyPI)
- Update imports: `from darlin...` → `from spatio_darlin...`
- Use `snakefiles/BMKS3000.smk` instead of `BMKS3000_matlab.smk` (MATLAB workflow: see [legacy/matlab/README.md](../legacy/matlab/README.md))
