## Prediction of hosts susceptible to SARS-CoV infection

.NOTPARALLEL:
.PHONY: all
all: plots


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

data/calculated/all_haddock_scores.rds: data/external/haddock_scores/
	mkdir -p data/calculated
	Rscript scripts/extract_haddock_scores.R

# Add Huang et al. scores and match taxonomy:
# - Retrieve taxonomy
data/calculated/taxonomy.rds: data/internal/NCBI_ACE2_orthologs.csv \
							  data/internal/ace2_accessions.csv
	mkdir -p data/calculated
	Rscript scripts/build_taxonomy_table.R

# - Get final scores
data/calculated/features_binding_affinity.rds:	data/internal/ace2_accessions.xlsx \
												data/internal/existing_predictions.csv \
												data/calculated/taxonomy.rds \
												data/calculated/all_haddock_scores.rds
	Rscript scripts/process_binding_affinities.R


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
	touch $@

data/iucn_range_maps/%.shp: data/iucn_range_maps/%.zip
	unzip -u -d $(@D) $<
	touch $@


# ---- ACE2 alignment -----------------------------------------------------------------------
# Align
data/calculated/ace2_protein_alignment.fasta: data/external/ace2_protein_sequences.fasta
	mkdir -p data/calculated
	mafft-einsi --thread 16 $< > $@

# Phylogeny
data/calculated/gene_tree/ace_genetree.treefile: data/calculated/ace2_protein_alignment.fasta
	mkdir -p data/calculated/gene_tree
	iqtree -nt 16 -s $< -bb 1000 -pre $(@D)/ace2_genetree


# ---- Pre-processing ------------------------------------------------------------------------------
# Clean metadata
data/calculated/cleaned_infection_data.rds: data/internal/infection_data.xlsx data/internal/ace2_accessions.xlsx
	mkdir -p data/calculated
	Rscript scripts/prepare_data.R


# Calculate features
data/calculated/features_pairwise_dists.rds: data/calculated/ace2_protein_alignment.fasta \
                                             data/calculated/cleaned_infection_data.rds
	Rscript scripts/calculate_features_aa.R

data/calculated/features_phylogeny_eigenvectors.rds: data/internal/timetree_amniota.nwk \
													 data/internal/timetree_mammalia.nwk \
													 data/internal/timetree_aves.nwk \
													 data/calculated/cleaned_infection_data.rds
	Rscript scripts/calculate_features_phylogeny.R

data/calculated/s_binding_alignment_positions.rds: data/calculated/ace2_protein_alignment.fasta \
												   data/internal/ace2_accessions.csv
	Rscript scripts/calculate_features_s_binding.R


# ---- Training: full model -----------------------------------------------------------------------
# Using all data and all ACE2-based features
# Output format is "dataset/response_var/feature_set/*"

TRAINING_REQUIREMENTS = data/calculated/cleaned_infection_data.rds \
						data/calculated/features_pairwise_dists.rds \
						data/calculated/features_binding_affinity.rds

ALL_FEATURE_SETS =	--aa_categorical \
					--aa_distance \
					--aa_properties \
					--distance_to_humans \
					--distance_to_positive \
					--binding_affinity

output/all_data/%/all_features/predictions.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		$(ALL_FEATURE_SETS) \
		--random_seed 77043274


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
		--random_seed 11386168

# Try timetree phylogeny as an alternative to ACE2 sequences:
output/all_data/%/phylogeny/predictions.rds:	$(TRAINING_REQUIREMENTS) \
												data/calculated/features_phylogeny_eigenvectors.rds
	Rscript scripts/train_models.R $* $(@D) \
		--phylogeny \
		--random_seed 34264755 \
		--n_threads 10

# ---- Training: variations on best model ---------------------------------------------------------
# Best ACE2 model restricted to S-binding sites only:
output/all_data/%/aa_distance_s_binding/predictions.rds: $(TRAINING_REQUIREMENTS) \
														 data/calculated/s_binding_alignment_positions.rds
	Rscript scripts/train_models.R $* $(@D) \
		--aa_distance \
		--s_binding_only \
		--random_seed 58498184 \
		--n_threads 20						 

# Phylogeny combined with the best ACE2 model
output/all_data/%/aa_distance_phylogeny/predictions.rds: $(TRAINING_REQUIREMENTS) \
														 data/calculated/features_phylogeny_eigenvectors.rds
	Rscript scripts/train_models.R $* $(@D) \
		--aa_distance \
		--phylogeny \
		--random_seed 34264755 \
		--n_threads 20

output/all_data/%/binding_affinity_phylogeny/predictions.rds:	$(TRAINING_REQUIREMENTS) \
														 		data/calculated/features_phylogeny_eigenvectors.rds
	Rscript scripts/train_models.R $* $(@D) \
		--binding_affinity \
		--phylogeny \
		--random_seed 34264755 \
		--n_threads 20


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
		--random_seed 23556284

# - Level 3-4 (cell culture)
output/l3+4_data/%/all_features/predictions.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		$(ALL_FEATURE_SETS) \
		--evidence_min 3 --evidence_max 4 \
		--random_seed 43564215

# - Level 1 and 2 (i.e. exclude cell culture, in case it makes things worse)
output/l1+2_data/%/all_features/predictions.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		$(ALL_FEATURE_SETS) \
		--evidence_min 1 --evidence_max 2 \
		--random_seed 28641685


# ---- Enemurate training combinations ------------------------------------------------------------
# Feature subsets
RESPONSE_VARS = infection
FEATURE_SET_NAMES = $(subst --, , $(ALL_FEATURE_SETS))  # Remove leading "--"

FEATURE_FOLDERS =	$(foreach a,$(RESPONSE_VARS), \
						$(foreach b,$(FEATURE_SET_NAMES), \
							$(a)/$(b) ))

FEATURE_MODELS = $(patsubst %, output/all_data/%/predictions.rds, $(FEATURE_FOLDERS))
PHYLO_MODELS = $(patsubst %, output/all_data/%/phylogeny/predictions.rds, $(RESPONSE_VARS))


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
train_feature_subsets: $(FEATURE_MODELS) $(PHYLO_MODELS)
train_data_subsets: $(L1_L2_MODELS) $(L3_MODELS)

train:	train_feature_subsets \
		train_data_subsets

# Fitted models should never be deleted (even when produced as an intermediate file
# for another step):
.PRECIOUS: $(FEATURE_MODELS) $(PHYLO_MODELS) $(L1_L2_MODELS) $(L3_MODELS)


# ---- Train an ensemble model ---------------------------------------------------------------------
output/all_data/infection/ensemble/predictions.rds: $(TRAINING_REQUIREMENTS) \
													data/calculated/features_phylogeny_eigenvectors.rds \
													output/all_data/infection/all_features/predictions.rds \
													output/all_data/infection/phylogeny/predictions.rds
	Rscript scripts/train_ensemble.R


# ---- Predict other species for which ACE2 sequences are available --------------------------------
# Phylogeny model predictions will be made simultaneously
output/all_data/infection/all_features/holdout_predictions.rds: output/all_data/infection/all_features/predictions.rds \
																data/calculated/cleaned_infection_data.rds \
																data/internal/NCBI_ACE2_orthologs.csv \
																$(TRAINING_REQUIREMENTS)
	Rscript scripts/predict_holdout.R


# ---- Plots ---------------------------------------------------------------------------------------
# Diagnostic plots
output/plots/performance.png: $(FEATURE_MODELS) $(L1_L2_MODELS) $(L3_MODELS) \
								output/all_data/infection/ensemble/predictions.rds
	Rscript scripts/plotting/plot_performance_diagnostics.R

# Data overview plots
.PHONY: report_data_overview
report_data_overview:	data/internal/infection_data.xlsx \
						data/calculated/cleaned_infection_data.rds \
						data/calculated/cleaned_shedding_data.rds
	Rscript scripts/plotting/report_data_overview_stats.R

output/plots/raw_data_overview.pdf: data/internal/timetree_amniota.nwk \
									data/calculated/cleaned_infection_data.rds \
									data/calculated/cleaned_shedding_data.rds \
									data/calculated/features_pairwise_dists.rds \
									data/internal/NCBI_ACE2_orthologs.csv \
									data/internal/ace2_accessions.csv
	Rscript scripts/plotting/plot_observations_phylogeny.R


output/plots/phylogeny_congruence.pdf:	data/internal/timetree_amniota.nwk \
										data/internal/ace2_accessions.csv \
										data/internal/NCBI_ACE2_orthologs.csv \
										data/calculated/features_pairwise_dists.rds
	Rscript scripts/plotting/plot_phylogeny_congruence.R


# Accuracy
# - Main figure
output/plots/accuracy.pdf:	$(FEATURE_MODELS) \
							output/all_data/infection/ensemble/predictions.rds
	Rscript scripts/plotting/plot_accuracy.R infection $@

# - Supplement (shedding models)
output/plots/accuracy_shedding.pdf:	$(FEATURE_MODELS) \
									output/all_data/shedding/ensemble/predictions.rds
	Rscript scripts/plotting/plot_accuracy.R shedding $@

# - Data quality (evidence level)
output/plots/accuracy_data_subsets.pdf: $(L1_L2_MODELS) $(L3_MODELS)
	Rscript scripts/plotting/plot_accuracy_data_subsets.R


# Variable importance
# - Cluster sites by correlation
output/plots/intermediates/feature_clusters.rds: data/calculated/features_variable_sites.rds
	Rscript scripts/plotting/get_clustered_sites.R

# - Plot
output/plots/varimp_overview.pdf:	output/all_data/infection/all_features/feature_importance.rds \
									output/plots/intermediates/feature_clusters.rds
	Rscript scripts/plotting/plot_varimp_overview.R $< $@


output/plots/site_varimp_supplement.pdf:	output/all_data/infection/aa_distance/feature_importance.rds \
											output/plots/intermediates/feature_clusters.rds
	Rscript scripts/plotting/plot_site_varimp_supplement.R


# Predictions:
# - Comparison to existing predictions
output/plots/intermediates/prediction_dendrogram.rds: 	data/internal/existing_predictions.csv \
														output/all_data/infection/all_features/predictions.rds \
														output/all_data/infection/phylogeny/predictions.rds \
														output/all_data/infection/all_features/holdout_predictions.rds
	Rscript scripts/plotting/get_prediction_dendrogram.R


output/plots/existing_predictions.pdf:	data/internal/timetree_amniota.nwk \
										output/plots/intermediates/prediction_dendrogram.rds \
										output/plots/raw_data_overview.pdf \
										data/internal/existing_predictions.csv \
										output/all_data/infection/all_features/holdout_predictions.rds
	Rscript scripts/plotting/plot_existing_predictions.R


# TODO: add SI plot (plot_prediction_comparison.R)



# - Get taxonomy
# TODO: taxonomy table still needed?
output/plots/intermediates/taxonomy_table.rds: output/all_data/infection/all_features/holdout_predictions.rds
	Rscript scripts/plotting/get_taxonomy.R

# - Prediction overview (phylogeny)
output/plots/holdout_predictions.png: output/all_data/infection/all_features/holdout_predictions.rds \
										output/plots/intermediates/taxonomy_table.rds
	Rscript scripts/plotting/plot_holdout_predictions.R

# - Maps
output/plots/prediction_maps.png:	output/all_data/infection/ensemble/holdout_predictions.rds \
									output/all_data/infection/phylogeny/holdout_predictions.rds \
									data/internal/timetree_mammalia.nwk \
									data/external/iucn_base/Land_Masses_and_Ocean_Islands.shp \
									data/iucn_range_maps/MAMMALS_TERRESTRIAL_ONLY.shp \
									data/iucn_range_maps/MAMMALS_FRESHWATER.shp 
									output/all_data/infection/all_features/holdout_predictions.rds \
									data/calculated/taxonomy.rds
	Rscript scripts/plotting/plot_holdout_maps.R

output/plots/ace2_availability_map_supplement.pdf:	output/all_data/infection/ensemble/holdout_predictions.rds \
													data/external/iucn_base/Land_Masses_and_Ocean_Islands.shp \
													data/iucn_range_maps/MAMMALS_TERRESTRIAL_ONLY.shp \
													data/iucn_range_maps/MAMMALS_FRESHWATER.shp \
													data/calculated/taxonomy.rds
	Rscript scripts/plotting/plot_ace2_availability_map_supplement.R



.PHONY: plots
plots: 	report_data_overview \
		output/plots/raw_data_overview.pdf \
		output/plots/phylogeny_congruence.pdf \
		output/plots/accuracy.pdf \
		output/plots/accuracy_data_subsets.pdf \
		output/plots/varimp_overview.pdf
