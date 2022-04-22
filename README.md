# Prediction of hosts susceptible to SARS-CoV infection using ACE2 protein sequences

## Usage
```
mamba env create -f environment.yml
conda activate sars_susceptibles
```

## _Development only (MacOS):_
```
export RSTUDIO_WHICH_R=$(which R)
open -a Rstudio
```


## Data definitions

### Shedding
Whether RNA or infectious material was detected in either:
    - nasal swab
    - nasal wash?
    - rectal swab
    - faeces

`shedding_virus`: Infectious virus demonstrated, either by onward transmission or by virus isolation/titration. Marked as negative if these techniques were attempted, but failed (even though e.g. virus isolation may have failed simply because the wrong cell line was used).
`shedding_rna`: RNA detected


## Usage
*Note: This code is fairly resource-intensive, and may require adjustment to run in a non-workstation environment. Models are currently set to use 20 threads, and will require ~25GB RAM. Edit the `n_threads` argument in `Makefile` to reduce this. Training ensemble models requires at least 62GB of RAM - if less is available, see `scripts/train_ensemble.R` for instructions.*