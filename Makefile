## Prediction of hosts susceptible to SARS-CoV infection

.NOTPARALLEL:
.PHONY: all
all: output/plots/performance.png


# ---- Data ----------------------------------------------------------------------------------------
# Spike protein amino acid sequences
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


# ---- Shapefiles for map figures ------------------------------------------------------------------
# Continent outlines (from https://www.iucnredlist.org/resources/spatialtoolsanddata)
data/external/iucn_base/Land_Masses_and_Ocean_Islands.zip:
	mkdir -p $(@D)
	curl -L -o $@ "https://spatial-data-2020onwards.s3.eu-west-1.amazonaws.com/WebsiteResources/Land_Masses_and_Ocean_Islands.zip"

# IUCN range maps have to downloaded manually (registration needed), so check if data exists and 
# display instructions if needed
data/iucn_range_maps/MAMMALS_TERRESTRIAL_ONLY.zip:
	$(warning $(file < data/iucn_range_maps/README.md))
	$(error Missing file "MAMMALS_TERRESTRIAL_ONLY.zip")

data/iucn_range_maps/MAMMALS_FRESHWATER.zip:
	$(warning $(file < data/iucn_range_maps/README.md))
	$(error Missing file "MAMMALS_FRESHWATER.zip")

data/iucn_range_maps/MAMMALS_MARINE_AND_TERRESTRIAL.zip:
	$(warning $(file < data/iucn_range_maps/README.md))
	$(error Missing file "MAMMALS_MARINE_AND_TERRESTRIAL.zip")

# Extract
data/external/iucn_base/Land_Masses_and_Ocean_Islands.shp: data/external/iucn_base/Land_Masses_and_Ocean_Islands.zip
	unzip -u -d $(@D) $<

data/iucn_range_maps/%.shp: data/iucn_range_maps/%.zip
	unzip -u -d $(@D) $<


# ---- ACE2 alignment -----------------------------------------------------------------------
# Align
data/calculated/ace2_protein_alignment.fasta: data/external/ace2_protein_sequences.fasta
	mkdir -p data/calculated
	mafft-einsi --thread 16 $< > $@

# Phylogeny
data/calculated/gene_tree/ace2_genetree.treefile: data/calculated/ace2_protein_alignment.fasta
	mkdir -p data/calculated/gene_tree
	iqtree -nt 16 -s $< -bb 1000 -pre $(@D)/ace_genetree


# ---- Pre-processing ------------------------------------------------------------------------------
# Retrieve taxonomy
data/calculated/taxonomy.rds: data/internal/NCBI_ACE2_orthologs.csv \
							  data/internal/ace2_accessions.csv
	mkdir -p data/calculated
	Rscript scripts/build_taxonomy_table.R


# Clean metadata
data/calculated/cleaned_infection_data.rds: data/internal/infection_data.xlsx data/internal/ace2_accessions.xlsx
	mkdir -p data/calculated
	Rscript scripts/prepare_data.R


# Calculate features
data/calculated/features_pairwise_dists.rds: data/calculated/ace2_protein_alignment.fasta \
                                             data/calculated/cleaned_infection_data.rds
	Rscript scripts/calculate_features.R


# ---- Training: full model -----------------------------------------------------------------------
# Using all data and all features
# Output format is "dataset/response_var/feature_set/*"

TRAINING_REQUIREMENTS = data/calculated/cleaned_infection_data.rds \
						data/calculated/features_pairwise_dists.rds \
						data/calculated/features_haddock_scores.rds

ALL_FEATURE_SETS =	--aa_categorical \
					--aa_distance \
					--aa_properties \
					--distance_to_humans \
					--distance_to_positive \
					--binding_affinity

output/all_data/%/all_features/predictions.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		$(ALL_FEATURE_SETS) \
		--random_seed 77043274 \
		--n_threads 8


# ---- Training: individual feature sets ----------------------------------------------------------
output/all_data/%/aa_categorical/predictions.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--aa_categorical \
		--random_seed 45323357 \
		--n_threads 10

output/all_data/%/aa_distance/predictions.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--aa_distance \
		--random_seed 58498184 \
		--n_threads 10

output/all_data/%/aa_properties/predictions.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--aa_properties \
		--random_seed 54564253 \
		--n_threads 10

output/all_data/%/distance_to_humans/predictions.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--distance_to_humans \
		--random_seed 75325883 \
		--n_threads 10

output/all_data/%/distance_to_positive/predictions.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--distance_to_positive \
		--random_seed 18681045 \
		--n_threads 10

output/all_data/%/binding_affinity/predictions.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--binding_affinity \
		--random_seed 11386168 \
		--n_threads 10


# ---- Training on data subsets ------------------------------------------------------------------------------------
# As above, output format is "dataset/response_var/feature_set/*"
# - All data models already trained during feature selection

# TODO: reduce number of features below

# All feature sets on data from each evidence level:
#  (note that level 1 [natural infection observed] has no negative data, so 
#  can't be included separately here)

# - Level 2 (experimental infection)
output/l2_data/%/all_features/predictions.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		$(ALL_FEATURE_SETS) \
		--evidence_min 2 --evidence_max 2 \
		--random_seed 23556284 \
		--n_threads 8

# - Level 3-4 (cell culture)
output/l3+4_data/%/all_features/predictions.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		$(ALL_FEATURE_SETS) \
		--evidence_min 3 --evidence_max 4 \
		--random_seed 43564215 \
		--n_threads 8

# - Level 1 and 2 (i.e. exclude cell culture, in case it makes things worse)
output/l1+2_data/%/all_features/predictions.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		$(ALL_FEATURE_SETS) \
		--evidence_min 1 --evidence_max 2 \
		--random_seed 28641685 \
		--n_threads 8


# ---- Enemurate training combinations ------------------------------------------------------------
# Feature subsets
RESPONSE_VARS = infection shedding
FEATURE_SET_NAMES = $(subst --, , $(ALL_FEATURE_SETS))  # Remove leading "--"

FEATURE_FOLDERS =	$(foreach a,$(RESPONSE_VARS), \
						$(foreach b,$(FEATURE_SET_NAMES), \
							$(a)/$(b) ))

FEATURE_MODELS = $(patsubst %, output/all_data/%/predictions.rds, $(FEATURE_FOLDERS))


# Data subsets (includes full data):
# - For l3+4 (cell culture), shedding does not apply
DATASETS = all_data l2_data l1+2_data

DATA_FOLDERS =	$(foreach a,$(DATASETS), \
					$(foreach b,$(RESPONSE_VARS), \
						$(a)/$(b) ))

L1_L2_MODELS = $(patsubst %, output/%/all_features/predictions.rds, $(DATA_FOLDERS))
L3_MODELS = output/l3+4_data/infection/all_features/predictions.rds


# Shortcuts:
.PHONY: train_feature_subsets train_data_subsets train
train_feature_subsets: $(FEATURE_MODELS)
train_data_subsets: $(L1_L2_MODELS) $(L3_MODELS)

train:	train_feature_subsets \
		train_data_subsets

# Fitted models should never be deleted (even when produced as an intermediate file
# for another step):
.PRECIOUS: $(FEATURE_MODELS) $(L1_L2_MODELS) $(L3_MODELS)


# TODO: updates needed below

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


# Maps
output/plots/maps.pdf: data/iucn_range_maps/MAMMALS_TERRESTRIAL_ONLY.shp \
					   data/iucn_range_maps/MAMMALS_FRESHWATER.shp \
					   data/iucn_range_maps/MAMMALS_MARINE_AND_TERRESTRIAL.shp
	$(error Not implemented) # TODO


.PHONY: plots
plots: output/plots/feature_selection.png \
	   output/plots/performance.png \
	   output/plots/varimp_infection.png \
	   output/plots/varimp_shedding.png
