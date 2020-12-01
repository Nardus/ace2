## Prediction of hosts susceptible to SARS-CoV infection

# Spike protein amino acid alignment
data/calculated/ace2_accessions.txt: data/internal/ace2_metadata.csv
	mkdir -p data/calculated
	cut -d, -f8 $< | awk 'FNR>1{print}' | sort | uniq > $@ 

data/external/ace2_protein_sequences.fasta: data/calculated/ace2_accessions.txt
	mkdir -p data/external
	epost -db protein -input $< -format acc | \
		efetch -format fasta > $@
	awk -f scripts/utils/find_missing_accessions.awk $< $@

data/calculated/ace2_protein_alignment.fasta: data/external/ace2_protein_sequences.fasta
	mafft-einsi --thread 8 $< > $@
	

# Clean metadata
data/calculated/cleaned_infection_data.rds: data/internal/ace2_metadata.csv
	Rscript scripts/prepare_data.R


# Calculate features
data/calculated/features_pairwise_dists.rds: data/calculated/ace2_protein_alignment.fasta \
                                             data/calculated/cleaned_infection_data.rds
	Rscript scripts/calculate_features.R

