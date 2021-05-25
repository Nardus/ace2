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
	cut -d, -f4 $< | awk 'FNR>1 && NF {print}' > $@ 

# - Get additional accessions from NCBI Orthologs:
data/calculated/ncbi_accessions.txt: data/internal/NCBI_ACE2_orthologs.csv
	mkdir -p data/calculated
	cut -d, -f7 $< | awk 'FNR>1 && NF {print}' > $@ 

data/calculated/all_accessions.txt: data/calculated/ace2_accessions.txt data/calculated/ncbi_accessions.txt
	cat $^ | sort | uniq > $@ 


# - Download sequences
data/external/ace2_protein_sequences.fasta: data/calculated/all_accessions.txt
	mkdir -p data/external
	epost -db protein -input $< -format acc | \
		efetch -format fasta > $@
	awk -f scripts/utils/find_missing_accessions.awk $< $@

# - Align
data/calculated/ace2_protein_alignment.fasta: data/external/ace2_protein_sequences.fasta
	mkdir -p data/calculated
	mafft-einsi --thread 16 $< > $@


# ---- Binding affinity data -----------------------------------------------------------------------
# From Fischhoff et al:
data/external/binding_affinity.tar.gz: 
	mkdir -p data/external
	curl -l -o $@ "https://zenodo.org/record/4517509/files/ace2-orthologs-dataset.tar.gz"

data/external/haddock_scores/: data/external/binding_affinity.tar.gz
	mkdir -p $@
	tar -xzv -C $@ -f $< "ace2-orthologs-dataset/refined_models/runs/*/ranks.models"

data/calculated/features_haddock_scores.rds: data/external/haddock_scores/ \
											 data/internal/ace2_accessions.xlsx
	mkdir -p data/calculated
	Rscript scripts/extract_haddock_scores.R


# ---- Pre-processing ------------------------------------------------------------------------------
# Clean metadata
data/calculated/cleaned_infection_data.rds: data/internal/infection_data.xlsx data/internal/ace2_accessions.xlsx
	mkdir -p data/calculated
	Rscript scripts/prepare_data.R


# Calculate features
data/calculated/features_pairwise_dists.rds: data/calculated/ace2_protein_alignment.fasta \
                                             data/calculated/cleaned_infection_data.rds
	Rscript scripts/calculate_features.R


# ---- Feature selection ---------------------------------------------------------------------------
# Output format is "dataset/response_var/feature_set/*"

TRAINING_REQUIREMENTS = data/calculated/cleaned_infection_data.rds \
						data/calculated/features_pairwise_dists.rds \
						data/calculated/features_haddock_scores.rds


# Full model (all features)
# - Here using replicates to check whether feature selection is repeatable
output/all_data/%/all_features/training_results.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--aa_categorical --aa_distance --distance_to_humans --binding_affinity \
		--replicates 100 \
		--random_seed 73049274

output/all_data/%/all_features/feature_usage.rds: output/all_data/%/all_features/training_results.rds
	Rscript scripts/select_features.R $(@D)

.PRECIOUS: output/all_data/shedding/all_features/feature_usage.rds \
		   output/all_data/infection/all_features/feature_usage.rds


# Same training as above, but with features reduced
# - last portion of folder name determines number of features kept
output/all_data/infection/feature_selection_%/training_results.rds: output/all_data/infection/all_features/feature_usage.rds \
																	$(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R infection $(@D) \
		--aa_categorical --aa_distance --distance_to_humans --binding_affinity \
		--random_seed 54278762 \
		--select_features $* \
		--feature_importance $<
		
output/all_data/shedding/feature_selection_%/training_results.rds: output/all_data/shedding/all_features/feature_usage.rds \
																   $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R shedding $(@D) \
		--aa_categorical --aa_distance --distance_to_humans --binding_affinity \
		--random_seed 27337413 \
		--select_features $* \
		--feature_importance $<

# Enumerate combinations:
#  - e.g. "output/all_data/infection/feature_selection_10/training_results.rds"
FEATURE_COUNTS = 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
RESPONSE_VARS = infection shedding

FEATURE_MODELS = $(foreach a,$(RESPONSE_VARS), \
					$(foreach b,$(FEATURE_COUNTS), \
						output/all_data/$(a)/feature_selection_$(b)/training_results.rds ))

.PHONY: train_feature_selection
train_feature_selection: $(FEATURE_MODELS)


# ---- Training on data subsets ------------------------------------------------------------------------------------
# As above, output format is "dataset/response_var/feature_set/*"
# - All data models already trained during feature selection

# TODO: reduce number of features below

# All feature sets on data from each evidence level:
#  (note that level 1 [natural infection observed] has no negative data, so 
#  can't be included separately here)

# - Level 2 (experimental infection)
output/l2_data/%/combined/training_results.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--aa_categorical --aa_distance --distance_to_humans --binding_affinity \
		--evidence_min 2 --evidence_max 2 \
		--random_seed 23556284

# - Level 3-4 (cell culture)
output/l3+4_data/%/combined/training_results.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--aa_categorical --aa_distance --distance_to_humans --binding_affinity \
		--evidence_min 3 --evidence_max 4 \
		--random_seed 43564215

# - Level 1 and 2 (i.e. exclude cell culture, in case it makes things worse)
output/l1+2_data/%/combined/training_results.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--aa_categorical --aa_distance --distance_to_humans --binding_affinity \
		--evidence_min 1 --evidence_max 2 \
		--random_seed 23556284


# Enumerate combinations:
# - For l3+4 (cell culture), shedding does not apply
DATASETS = all_data l2_data l1+2_data

OUT_FOLDERS = $(foreach a,$(DATASETS), \
				$(foreach b,$(RESPONSE_VARS), \
					$(a)/$(b) ))

L1_L2_MODELS = $(patsubst %, output/%/combined/training_results.rds, $(OUT_FOLDERS))
L3_MODELS = output/l3+4_data/infection/combined/training_results.rds

.PHONY: train_all_features train
train_all_features: $(L1_L2_MODELS) $(L3_MODELS)

train: train_all_features \
       train_feature_selection

# Fitted models should never be deleted (even when produced as an intermediate file
# for another step):
.PRECIOUS: $(FEATURE_MODELS) $(L1_L2_MODELS) $(L3_MODELS) \
		   output/all_data/infection/all_features/training_results.rds \
		   output/all_data/shedding/all_features/training_results.rds


# ---- Predict other species for which ACE2 sequences are available --------------------------------

output/all_data/infection/feature_selection_2/additional_preds_infection.rds: output/all_data/infection/feature_selection_2/training_results.rds \
																			  output/all_data/shedding/feature_selection_2/training_results.rds \
																			  data/internal/NCBI_ACE2_orthologs.csv
	Rscript scripts/predict_holdout.R


# ---- Plots ---------------------------------------------------------------------------------------
output/plots/feature_selection.png: $(FEATURE_MODELS)
	Rscript scripts/plotting/plot_feature_selection.R

output/plots/performance.png: $(FEATURE_MODELS) $(L2_MODELS) $(L3_MODELS)
	Rscript scripts/plotting/plot_performance.R


# Variable importance
output/plots/varimp_overview_%.png: output/all_data/%/all_features/feature_usage.rds \
									output/all_data/%/all_features/feature_usage_by_iteration.rds
	Rscript scripts/plotting/plot_varimp_overview.R $(word 2,$^) $@


output/plots/varimp_detail_infection.png: output/all_data/infection/all_features/training_results.rds \
										  output/all_data/infection/feature_selection_2/training_results.rds
	Rscript scripts/plotting/plot_varimp_detail.R \
		--full_model $< \
		--best_model echo $(word 2,$^) \
		--output_name $@


output/plots/varimp_infection.png: $(FEATURE_MODELS)
	Rscript scripts/plotting/plot_varimp_infection.R


output/plots/varimp_shedding.png: $(FEATURE_MODELS)
	Rscript scripts/plotting/plot_varimp_shedding.R


.PHONY: plots
plots: output/plots/feature_selection.png \
	   output/plots/performance.png \
	   output/plots/varimp_infection.png \
	   output/plots/varimp_shedding.png
