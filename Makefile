## Prediction of hosts susceptible to SARS-CoV infection

.NOTPARALLEL:
.PHONY: all
all: output/plots/performance.png


# ---- Data ----------------------------------------------------------------------------------------
# Spike protein amino acid alignment

# - Check that accession csv is up to date
data/internal/ace2_accessions.csv: data/internal/ace2_accessions.xlsx
	$(error ace2_accessions.csv is out of date - re-export the matching excel file)

# - Extract accessions, removing the header row, blank lines, and duplicates:
data/calculated/ace2_accessions.txt: data/internal/ace2_accessions.csv
	mkdir -p data/calculated
	cut -d, -f4 $< | awk 'FNR>1 && NF {print}' | sort | uniq > $@ 

# - Download sequences
data/external/ace2_protein_sequences.fasta: data/calculated/ace2_accessions.txt
	mkdir -p data/external
	epost -db protein -input $< -format acc | \
		efetch -format fasta > $@
	awk -f scripts/utils/find_missing_accessions.awk $< $@

# - Align
data/calculated/ace2_protein_alignment.fasta: data/external/ace2_protein_sequences.fasta
	mafft-einsi --thread 8 $< > $@
	

# ---- Pre-processing ------------------------------------------------------------------------------
# Clean metadata
data/calculated/cleaned_infection_data.rds: data/internal/infection_data.xlsx data/internal/ace2_accessions.xlsx
	Rscript scripts/prepare_data.R


# Calculate features
data/calculated/features_pairwise_dists.rds: data/calculated/ace2_protein_alignment.fasta \
                                             data/calculated/cleaned_infection_data.rds
	Rscript scripts/calculate_features.R


# ---- Training ------------------------------------------------------------------------------------
# Output format is "dataset/response_var/feature_set/*"

TRAINING_REQUIREMENTS = data/calculated/cleaned_infection_data.rds \
						data/calculated/features_pairwise_dists.rds

# - All feature sets
#   TODO: add --binding_affinity
output/all_data/%/combined/training_results.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--aa_categorical --aa_distance --distance_to_humans \
		--random_seed 73049274


# All feature sets on data from each evidence level:
#  (note that level 1 [natural infection observed] has no negative data, so 
#  can't be included separately here)
#   TODO: add --binding_affinity to all of these

# - Level 2 (experimental infection)
output/l2_data/%/combined/training_results.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--aa_categorical --aa_distance --distance_to_humans \
		--evidence_min 2 --evidence_max 2 \
		--random_seed 23556284

# - Level 3-4 (cell culture)
output/l3+4_data/%/combined/training_results.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--aa_categorical --aa_distance --distance_to_humans \
		--evidence_min 3 --evidence_max 4 \
		--random_seed 43564215

# - Level 1 and 2 (i.e. exclude cell culture, in case it makes things worse)
output/l1+2_data/%/combined/training_results.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--aa_categorical --aa_distance --distance_to_humans \
		--evidence_min 1 --evidence_max 2 \
		--random_seed 23556284


# Enumerate combinations:
# - For l3+4 (cell culture), shedding does not apply
DATASETS = {"all_data/","l2_data/","l1+2_data/"}
RESPONSE_VARS = {"infection","shedding"}

OUT_FOLDERS = $(shell echo $(DATASETS)$(RESPONSE_VARS))
L1_L2_MODELS = $(patsubst %, output/%/combined/training_results.rds, $(OUT_FOLDERS))
L3_MODELS = output/l3+4_data/infection/combined/training_results.rds

.PHONY: train_all_features
train_all_features: $(L1_L2_MODELS) $(L3_MODELS)


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

.PHONY: train_feature_selection train
train_feature_selection: output/all_data/infection/combined+feature_selection_100/training_results.rds \
				   		 output/all_data/shedding/combined+feature_selection_100/training_results.rds

train: train_all_features \
       train_feature_selection


# ---- Plots ---------------------------------------------------------------------------------------
output/plots/performance.png: $(FEATURE_MODELS) $(L2_MODELS) $(L3_MODELS)
	Rscript scripts/plotting/plot_performance.R
	
.PHONY: plots
plots: output/plots/performance.png