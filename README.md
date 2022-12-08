# Prediction of hosts susceptible to SARS-CoV infection using ACE2 protein sequences

[![DOI](https://zenodo.org/badge/316297169.svg)](https://zenodo.org/badge/latestdoi/316297169)

Code and data used in Mollentze N, Keen D, Munkhbayar U, Biek R, Streicker DG. (2022). *Variation in the ACE2 receptor has limited utility for SARS-CoV-2 host prediction*. eLife 11:e80329. [doi:10.7554/eLife.80329]( https://doi.org/10.7554/eLife.80329).

## Requirements
- Install the [mamba package manager](https://mamba.readthedocs.io/en/latest/installation.html)
- All other requirements can then be installed using:

```
mamba env create -f environment.yml
```

## Repeating published analyses
To repeat all analyses and recreate the figures used in the manuscript, run:
```
conda activate sars_susceptibles
make
```

## Usage notes
1. Most required external data will be downloaded automatically, with two exceptions:
    - **ACE2 sequences**. The `Makefile` contains code to download these, but any sequence updates by NCBI will mean sequences matching the earlier accession number used here cannot be dowloaded (easily). The sequences used here are therefore included in [`data/external/ace2_protein_sequences.fasta`](data/external/ace2_protein_sequences.fasta).
    - **IUCN range data** (required to recreate map figures). Downloads require a log in, see [`data/iucn_range_maps/`](data/iucn_range_maps/) for instructions.

2. This pipeline was designed for workstations, and may need some minor modifications to run on a regular PC:
    - By default, 20 parallel threads will be used. Edit the `--n_threads` argument throughout the [`Makefile`](Makefile) to change this
    - Creating ensemble models may use large amounts of memory. See [`scripts/train_ensemble.R`](scripts/train_ensemble.R) for instructions if this becomes an issue.
