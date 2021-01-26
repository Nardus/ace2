## Prediction of hosts susceptible to SARS-CoV infection

.NOTPARALLEL:
.PHONY: all
all: output/plots/performance.png


# ---- Data ----------------------------------------------------------------------------------------
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
	

# ---- Pre-processing ------------------------------------------------------------------------------
# Clean metadata
data/calculated/cleaned_infection_data.rds: data/internal/ace2_metadata.csv
	Rscript scripts/prepare_data.R


# Calculate features
data/calculated/features_pairwise_dists.rds: data/calculated/ace2_protein_alignment.fasta \
                                             data/calculated/cleaned_infection_data.rds
	Rscript scripts/calculate_features.R


# ---- Training ------------------------------------------------------------------------------------
# Output format is "dataset/response_var/feature_set/*"

TRAINING_REQUIREMENTS = data/calculated/cleaned_infection_data.rds \
						data/calculated/features_pairwise_dists.rds

# Feature sets:
# - Categorical AA representation
output/%/aa_categorical/training_results.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* "$(@D)" --aa_categorical --random_seed 21451233

# - Distance AA representation
output/%/aa_distance/training_results.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* "$(@D)" --aa_distance --random_seed 58371313

# - Distance to humans
output/%/distance_to_humans/training_results.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* "$(@D)" --distance_to_humans --random_seed 42943872

# - Binding affinity
output/%/binding_affinity/training_results.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* "$(@D)" --binding_affinity --random_seed 83720341

# - Combined
output/%/combined/training_results.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* "$(@D)" \
		--aa_categorical --aa_distance --distance_to_humans \
		--binding_affinity --binding_affinity \
		--random_seed 73049274


# Enumerate combinations:
FEATURE_SETS = aa_categorical aa_distance distance_to_humans #binding_affinity combined
INFECTION_MODELS = $(patsubst %, output/infection/%/training_results.rds, $(FEATURE_SETS))
SHEDDING_MODELS = $(patsubst %, output/shedding/%/training_results.rds, $(FEATURE_SETS))
TRANSMISSION_MODELS = $(patsubst %, output/transmission/%/training_results.rds, $(FEATURE_SETS))

.PHONY: train
train: $(INFECTION_MODELS) $(SHEDDING_MODELS) $(TRANSMISSION_MODELS)



# ---- Plots ---------------------------------------------------------------------------------------
output/plots/performance.png: output/infection/training_results.rds \
                              output/shedding/training_results.rds \
							  output/transmission/training_results.rds
	Rscript scripts/plotting/plot_performance.R