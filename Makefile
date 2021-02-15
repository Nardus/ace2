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
# TODO: not currently distinguishing by dataset (below should be "all_data")

TRAINING_REQUIREMENTS = data/calculated/cleaned_infection_data.rds \
						data/calculated/features_pairwise_dists.rds

# Try different feature sets:
# - Categorical AA representation
output/all_data/%/aa_categorical/training_results.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) --aa_categorical --random_seed 21451233

# - Distance AA representation
output/all_data/%/aa_distance/training_results.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) --aa_distance --random_seed 58371313

# - Distance to humans
output/all_data/%/distance_to_humans/training_results.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) --distance_to_humans --random_seed 42943872

# - Binding affinity
output/all_data/%/binding_affinity/training_results.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) --binding_affinity --random_seed 83720341

# - Combined
#   TODO: add --binding_affinity
output/all_data/%/combined/training_results.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--aa_categorical --aa_distance --distance_to_humans \
		--random_seed 73049274


# Full model on data from each evidence level:
#  (note that level 1 [natural infection observed] has no negative data, so 
#  can't be included here)
#   TODO: add --binding_affinity to all of these

# - Level 2 (experimental infection)
output/l2_data/%/combined/training_results.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--aa_categorical --aa_distance --distance_to_humans \
		--evidence_min 2 --evidence_max 2 \
		--random_seed 23556284

# - Level 3 (cell culture)
output/l3_data/%/combined/training_results.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--aa_categorical --aa_distance --distance_to_humans \
		--evidence_min 3 --evidence_max 3 \
		--random_seed 43564215

# - Level 1 and 2 (i.e. exclude cell culture, in case it makes things worse)
output/l1+2_data/%/combined/training_results.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--aa_categorical --aa_distance --distance_to_humans \
		--evidence_min 1 --evidence_max 2 \
		--random_seed 23556284


# Enumerate combinations:
# - Feature set models (on all_data only)
RESPONSE_VARS = {"infection/","shedding/","transmission/"}
FEATURE_SETS = {"aa_categorical","aa_distance","distance_to_humans","combined"} # TODO: binding_affinity 

OUT_FOLDERS_1 = $(shell echo $(RESPONSE_VARS)$(FEATURE_SETS))
FEATURE_MODELS = $(patsubst %, output/all_data/%/training_results.rds, $(OUT_FOLDERS_1))

# - Evidence level models (on combined model only)
#   For l3 (cell culture), shedding and transmission do not apply
DATASETS = {"l2_data/","l1+2_data/"}

OUT_FOLDERS_2 = $(shell echo $(DATASETS)$(RESPONSE_VARS))
L1_L2_MODELS = $(patsubst %, output/%combined/training_results.rds, $(OUT_FOLDERS_2))
L3_MODELS = output/l3_data/infection/combined/training_results.rds

.PHONY: train
train: $(FEATURE_MODELS) $(L1_L2_MODELS) $(L3_MODELS)


# ---- Feature selection ---------------------------------------------------------------------------
output/all_data/%/combined/feature_usage.rds: output/all_data/%/combined/training_results.rds
	Rscript scripts/select_features.R $(@D)


# TODO: add binding_affinity here
output/all_data/%/combined+feature_selection_100/training_results.rds: output/all_data/%/combined/feature_usage.rds \
																	   $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--aa_categorical --aa_distance --distance_to_humans \
		--random_seed 38721019 \
		--select_features 100 \
		--feature_importance $<

.PHONY: feature_selection
feature_selection: output/all_data/infection/combined+feature_selection_100/training_results.rds \
				   output/all_data/shedding/combined+feature_selection_100/training_results.rds \
				   output/all_data/transmission/combined+feature_selection_100/training_results.rds


# ---- Plots ---------------------------------------------------------------------------------------
output/plots/performance.png: $(FEATURE_MODELS) $(L2_MODELS) $(L3_MODELS)
	Rscript scripts/plotting/plot_performance.R
	
.PHONY: plots
plots: output/plots/performance.png