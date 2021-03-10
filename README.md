# Prediction of hosts susceptible to SARS-CoV infection using ACE2 protein sequences

## Usage
```
conda env create -f environment.yml
conda activate sars_susceptibles
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