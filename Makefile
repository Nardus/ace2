## Prediction of hosts susceptible to SARS-CoV infection

.NOTPARALLEL:
.PHONY: all
all:	train \
		report_data_overview \
		plots


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


# ---- ACE2 alignment ------------------------------------------------------------------------------
# Align
data/calculated/ace2_protein_alignment.fasta: data/external/ace2_protein_sequences.fasta
	mkdir -p data/calculated
	mafft-einsi --thread 16 $< > $@

# Phylogeny
data/calculated/gene_tree/ace2_genetree.treefile: data/calculated/ace2_protein_alignment.fasta
	mkdir -p data/calculated/gene_tree
	iqtree -nt 16 -s $< -bb 1000 -pre $(@D)/ace2_genetree


# ---- Pre-processing ------------------------------------------------------------------------------
# Clean metadata
data/calculated/cleaned_infection_data.rds: data/internal/infection_data.xlsx data/internal/ace2_accessions.xlsx
	mkdir -p data/calculated
	Rscript scripts/prepare_data.R


# Calculate features
data/calculated/features_pairwise_dists.rds: data/calculated/ace2_protein_alignment.fasta \
                                             data/calculated/cleaned_infection_data.rds \
											 data/calculated/taxonomy.rds
	Rscript scripts/calculate_features_aa.R

data/calculated/features_phylogeny_eigenvectors.rds: data/internal/timetree_amniota.nwk \
													 data/internal/timetree_mammalia.nwk \
													 data/internal/timetree_aves.nwk \
													 data/calculated/cleaned_infection_data.rds
	Rscript scripts/calculate_features_phylogeny.R

data/calculated/s_binding_alignment_positions.rds: data/calculated/ace2_protein_alignment.fasta \
												   data/internal/ace2_accessions.csv
	Rscript scripts/calculate_features_s_binding.R


# ---- Training: full model ------------------------------------------------------------------------
# Using all data and all ACE2-based features
# Output format is "dataset/response_var/feature_set/*"

TRAINING_REQUIREMENTS = data/calculated/cleaned_infection_data.rds \
						data/calculated/features_pairwise_dists.rds \
						data/calculated/features_binding_affinity.rds \
						data/calculated/features_phylogeny_eigenvectors.rds

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
		--n_threads 20


# ---- Training: individual feature sets -----------------------------------------------------------
output/all_data/%/aa_categorical/predictions.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--aa_categorical \
		--random_seed 45323357 \
		--n_threads 20

output/all_data/%/aa_distance/predictions.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--aa_distance \
		--random_seed 58498184 \
		--n_threads 20
		
output/all_data/%/_aa_distance2/predictions.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--aa_distance \
		--random_seed 98503201 \
		--n_threads 20

output/all_data/%/aa_properties/predictions.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--aa_properties \
		--random_seed 54564253 \
		--n_threads 20

output/all_data/%/distance_to_humans/predictions.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--distance_to_humans \
		--random_seed 75325883 \
		--n_threads 20

output/all_data/%/distance_to_positive/predictions.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--distance_to_positive \
		--random_seed 18681045 \
		--n_threads 20

output/all_data/%/binding_affinity/predictions.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--binding_affinity \
		--random_seed 11386168 \
		--n_threads 20

# Try timetree phylogeny as an alternative to ACE2 sequences:
output/all_data/%/phylogeny/predictions.rds:	$(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--phylogeny \
		--random_seed 34264755 \
		--n_threads 20

# ---- Training: variations on best models ---------------------------------------------------------
# Best ACE2 model restricted to S-binding sites only:
output/all_data/%/_supplementary_runs/aa_distance_s_binding/predictions.rds:	$(TRAINING_REQUIREMENTS) \
																				data/calculated/s_binding_alignment_positions.rds
	Rscript scripts/train_models.R $* $(@D) \
		--aa_distance \
		--s_binding_only \
		--random_seed 58498184 \
		--n_threads 20						 

# Phylogeny combined with the best ACE2 models
output/all_data/%/aa_distance_phylogeny/predictions.rds: $(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--aa_distance \
		--phylogeny \
		--random_seed 34264755 \
		--n_threads 20

output/all_data/%/binding_affinity_phylogeny/predictions.rds:	$(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--binding_affinity \
		--phylogeny \
		--random_seed 34264755 \
		--n_threads 20

output/all_data/%/all_features_phylogeny/predictions.rds:	$(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		$(ALL_FEATURE_SETS) \
		--phylogeny \
		--random_seed 09524845 \
		--n_threads 20


# ---- Training on data subsets --------------------------------------------------------------------
# As above, output format is "dataset/response_var/feature_set/*"

# All feature sets on data from each evidence level:
#  (note that level 1 [natural infection observed] has no negative data, so 
#  can't be included separately here)

# - Level 2 (experimental infection)
output/l2_data/%/all_features_phylogeny/predictions.rds:	$(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		$(ALL_FEATURE_SETS) \
		--phylogeny \
		--evidence_min 2 --evidence_max 2 \
		--random_seed 23556284 \
		--n_threads 20

# - Level 3-4 (cell culture)
output/l3+4_data/%/all_features_phylogeny/predictions.rds:	$(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		$(ALL_FEATURE_SETS) \
		--phylogeny \
		--evidence_min 3 --evidence_max 4 \
		--random_seed 43564215 \
		--n_threads 20

# - Level 1 and 2 (i.e. exclude cell culture, in case it makes things worse)
output/l1+2_data/%/all_features_phylogeny/predictions.rds:	$(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		$(ALL_FEATURE_SETS) \
		--phylogeny \
		--evidence_min 1 --evidence_max 2 \
		--random_seed 28641685 \
		--n_threads 20


# Without Rhinolophid bats
# - All ACE2 features
output/no_rhinolophids/%/all_features/predictions.rds:	$(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		$(ALL_FEATURE_SETS) \
		--exclude_rhinolophid \
		--random_seed 45124422 \
		--n_threads 20

# - Phylogeny
output/no_rhinolophids/%/phylogeny/predictions.rds:	$(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		--phylogeny \
		--exclude_rhinolophid \
		--random_seed 53241688 \
		--n_threads 20

# - ACE2 + phylogeny
output/no_rhinolophids/%/all_features_phylogeny/predictions.rds:	$(TRAINING_REQUIREMENTS)
	Rscript scripts/train_models.R $* $(@D) \
		$(ALL_FEATURE_SETS) \
		--phylogeny \
		--exclude_rhinolophid \
		--random_seed 97631262 \
		--n_threads 20


# ---- Train ensemble models -----------------------------------------------------------------------
# All ACE2 + phylogeny
output/all_data/%/ensemble_all_features_phylogeny/predictions.rds:	$(TRAINING_REQUIREMENTS) \
																	output/all_data/%/all_features/predictions.rds \
																	output/all_data/%/phylogeny/predictions.rds
	Rscript scripts/train_ensemble.R \
			--m1 output/all_data/$*/all_features \
			--m2 output/all_data/$*/phylogeny \
			--output_path $(@D) \
			--random_seed 10012022

# ACE2 distance and binding affinity
output/all_data/%/ensemble_aa_distance_binding_affinity/predictions.rds:	$(TRAINING_REQUIREMENTS) \
																			output/all_data/%/aa_distance/predictions.rds \
																			output/all_data/%/binding_affinity/predictions.rds
	Rscript scripts/train_ensemble.R \
			--m1 output/all_data/$*/aa_distance \
			--m2 output/all_data/$*/binding_affinity \
			--output_path $(@D) \
			--random_seed 09524845

# ACE2 distance with another version of itself (differing only in random seed)
output/all_data/%/ensemble_aa_distance_self/predictions.rds:	$(TRAINING_REQUIREMENTS) \
																output/all_data/%/aa_distance/predictions.rds \
																output/all_data/%/_aa_distance2/predictions.rds
	Rscript scripts/train_ensemble.R \
			--m1 output/all_data/$*/aa_distance \
			--m2 output/all_data/$*/_aa_distance2 \
			--output_path $(@D) \
			--random_seed 09524845


# ---- Enemurate training combinations -------------------------------------------------------------
# Feature subsets
RESPONSE_VARS = infection shedding
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

L1_L2_MODELS = $(patsubst %, output/%/all_features_phylogeny/predictions.rds, $(DATA_FOLDERS))
L3_MODELS = output/l3+4_data/infection/all_features_phylogeny/predictions.rds

# - Rhinolophids excluded
NO_RHINOLOPHIDS =	output/no_rhinolophids/infection/all_features/predictions.rds \
					output/no_rhinolophids/infection/phylogeny/predictions.rds \
					output/no_rhinolophids/infection/all_features_phylogeny/predictions.rds


# Ensemble models
ENSEMBLE_COMBINATIONS = ensemble_all_features_phylogeny \
						ensemble_aa_distance_binding_affinity \
						ensemble_aa_distance_self

ENSEMBLE_FOLDERS =	$(foreach a,$(RESPONSE_VARS), \
						$(foreach b,$(ENSEMBLE_COMBINATIONS), \
							$(a)/$(b) ))

ENSEMBLE_MODELS = $(patsubst %, output/all_data/%/predictions.rds, $(ENSEMBLE_FOLDERS))


# Other models (used for checks / supplementary plots)
OTHER_MODELS =	output/all_data/infection/aa_distance_phylogeny/predictions.rds \
				output/all_data/infection/_supplementary_runs/aa_distance_s_binding/predictions.rds \
				

# Shortcuts:
.PHONY: train_feature_subsets train_data_subsets train
train_feature_subsets: $(FEATURE_MODELS) $(PHYLO_MODELS)
train_data_subsets: $(L1_L2_MODELS) $(L3_MODELS) $(NO_RHINOLOPHIDS)

train:	train_feature_subsets \
		train_data_subsets \
		$(ENSEMBLE_MODELS) \
		$(OTHER_MODELS)
		

# Fitted models should never be deleted (even when produced as an intermediate file
# for another step):
.PRECIOUS: $(FEATURE_MODELS) $(PHYLO_MODELS) $(L1_L2_MODELS) $(L3_MODELS) $(OTHER_MODELS)


# ---- Predict other species for which ACE2 sequences are available --------------------------------
# Phylogeny and ensemble model predictions will be made simultaneously
output/all_data/infection/aa_distance/holdout_predictions.rds: output/all_data/infection/aa_distance/predictions.rds \
															   output/all_data/infection/ensemble_aa_distance_self/predictions.rds \
															   output/all_data/infection/phylogeny/predictions.rds \
															   data/calculated/cleaned_infection_data.rds \
															   data/internal/NCBI_ACE2_orthologs.csv \
															   data/calculated/taxonomy.rds \
															   $(TRAINING_REQUIREMENTS)
	Rscript scripts/predict_holdout.R


# ---- Plots ---------------------------------------------------------------------------------------
# Diagnostic plots
output/plots/performance.png: $(FEATURE_MODELS) $(L1_L2_MODELS) $(L3_MODELS) $(ENSEMBLE_MODELS)
	Rscript scripts/plotting/plot_performance_diagnostics.R

# Data overview plots
.PHONY: report_data_overview report_distance_metric_correlation
report_data_overview:	data/internal/infection_data.xlsx \
						data/calculated/cleaned_infection_data.rds \
						data/calculated/cleaned_shedding_data.rds
	Rscript scripts/plotting/report_data_overview_stats.R

report_distance_metric_correlation: data/calculated/cleaned_infection_data.rds \
									data/internal/timetree_amniota.nwk\
									data/calculated/features_pairwise_dists.rds \
									data/internal/NCBI_ACE2_orthologs.csv \
									data/internal/ace2_accessions.csv
	Rscript scripts/plotting/report_distance_metric_correlation.R

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
										data/calculated/features_pairwise_dists.rds \
										data/calculated/gene_tree/ace2_genetree.treefile
	Rscript scripts/plotting/plot_phylogeny_congruence.R


# Accuracy
# - Main figure
output/plots/accuracy.pdf:	$(FEATURE_MODELS) \
							$(PHYLO_MODELS) \
							output/all_data/infection/all_features_phylogeny/predictions.rds \
							output/all_data/infection/aa_distance_phylogeny/predictions.rds \
							output/all_data/infection/ensemble_all_features_phylogeny/predictions.rds \
							output/all_data/infection/ensemble_aa_distance_binding_affinity/predictions.rds \
							output/all_data/infection/ensemble_aa_distance_self/predictions.rds 
	Rscript scripts/plotting/plot_accuracy.R infection $@

# - Supplement (S-binding only)
output/plots/accuracy_supplement_sbinding.pdf:	data/calculated/cleaned_infection_data.rds \
												output/all_data/infection/_supplementary_runs/aa_distance_s_binding/predictions.rds
	Rscript scripts/plotting/plot_accuracy_supplement_sbinding.R

# - Supplement (shedding models)
output/plots/accuracy_shedding.pdf:	$(FEATURE_MODELS) \
									output/all_data/shedding/ensemble_all_features_phylogeny/predictions.rds
	Rscript scripts/plotting/plot_accuracy.R shedding $@

# - Supplement (data quality  / evidence level)
output/plots/accuracy_data_subsets.pdf: $(L1_L2_MODELS) $(L3_MODELS)
	Rscript scripts/plotting/plot_accuracy_data_subsets.R

# - Supplement (rhinolophid bats)
output/plots/accuracy_rhinolophid.pdf:	$(NO_RHINOLOPHIDS)
	Rscript scripts/plotting/plot_accuracy_rhinolophid.R $< $@

# - Supplement (non-ACE2 sarbecoviruses)
output/plots/accuracy_non_ace2_all_features.pdf:	output/all_data/infection/all_features/predictions.rds \
													data/calculated/cleaned_infection_data.rds
	Rscript scripts/plotting/plot_accuracy_supplement_non_ace2.R $< $@
	

output/plots/accuracy_non_ace2_phylogeny.pdf:	output/all_data/infection/phylogeny/predictions.rds \
												data/calculated/cleaned_infection_data.rds
	Rscript scripts/plotting/plot_accuracy_supplement_non_ace2.R $< $@


# Variable importance
# - Cluster sites by correlation
output/plots/intermediates/feature_clusters.rds: data/calculated/features_variable_sites.rds
	Rscript scripts/plotting/get_clustered_sites.R

# - Plot
output/plots/varimp_overview.pdf:	output/all_data/infection/all_features/predictions.rds \
									output/all_data/infection/aa_distance/predictions.rds \
									output/plots/intermediates/feature_clusters.rds
	Rscript scripts/plotting/plot_varimp_overview.R


output/plots/varimp_supplement.pdf: output/all_data/infection/all_features/predictions.rds \
									output/plots/intermediates/feature_clusters.rds
	Rscript scripts/plotting/plot_varimp_supplement.R


# Predictions:
# - Comparison to existing predictions
output/plots/intermediates/prediction_dendrogram.rds: 	data/internal/existing_predictions.csv \
														output/all_data/infection/aa_distance/holdout_predictions.rds
	Rscript scripts/plotting/get_prediction_dendrogram.R


output/plots/existing_predictions.pdf:	output/plots/raw_data_overview.pdf \
										data/internal/existing_predictions.csv \
										output/plots/intermediates/prediction_dendrogram.rds \
										output/all_data/infection/aa_distance/predictions.rds \
										output/all_data/infection/phylogeny/predictions.rds \
										output/all_data/infection/ensemble_aa_distance_self/predictions.rds
	Rscript scripts/plotting/plot_existing_predictions.R


# TODO: add SI plot (plot_prediction_comparison.R)
output/plots/existing_predictions_supplement.pdf:	data/internal/timetree_amniota.nwk \
													output/plots/intermediates/prediction_dendrogram.rds \
													output/plots/raw_data_overview.pdf \
													data/calculated/cleaned_infection_data.rds
	Rscript scripts/plotting/plot_existing_predictions_supplement.R


# - Prediction overview (phylogeny)
output/plots/predictions_by_order_supplement.pdf:	output/all_data/infection/aa_distance/holdout_predictions.rds \
													data/calculated/cleaned_infection_data.rds \
													data/calculated/taxonomy.rds \
													data/internal/timetree_amniota.nwk
	Rscript scripts/plotting/plot_predictions_by_order.R

output/plots/predictions_by_family_supplement.pdf:	output/all_data/infection/aa_distance/holdout_predictions.rds \
													data/calculated/cleaned_infection_data.rds \
													data/calculated/taxonomy.rds
	Rscript scripts/plotting/plot_predictions_by_family.R

# - Maps
output/plots/prediction_maps.png:	output/all_data/infection/aa_distance/holdout_predictions.rds \
									data/internal/timetree_mammalia.nwk \
									data/external/iucn_base/Land_Masses_and_Ocean_Islands.shp \
									data/iucn_range_maps/MAMMALS_TERRESTRIAL_ONLY.shp \
									data/iucn_range_maps/MAMMALS_FRESHWATER.shp \
									data/calculated/taxonomy.rds
	Rscript scripts/plotting/plot_prediction_maps.R


output/plots/prediction_maps_by_order.png:	output/all_data/infection/aa_distance/holdout_predictions.rds \
											data/calculated/cleaned_infection_data.rds \
											data/calculated/taxonomy.rds \
											data/iucn_range_maps/MAMMALS_TERRESTRIAL_ONLY.shp
	Rscript scripts/plotting/plot_maps_by_order.R


# - Prioritisation based on scores
output/plots/phylogeny_predictions_supplement.pdf:	output/all_data/infection/aa_distance/holdout_predictions.rds \
													data/internal/timetree_mammalia.nwk \
													data/calculated/cleaned_shedding_data.rds
	Rscript scripts/plotting/plot_phylogeny_predictions_supplement.R


# Supplementary tables
output/si_tables/supplement_training_data.xlsx: data/internal/infection_data.xlsx \
												data/internal/ace2_accessions.csv \
												data/internal/data_citations.bib \
												data/calculated/cleaned_infection_data.rds \
												data/calculated/cleaned_shedding_data.rds \
												data/calculated/taxonomy.rds \
												output/all_data/infection/aa_distance/holdout_predictions.rds \
												data/internal/NCBI_ACE2_orthologs.csv
	Rscript scripts/plotting/make_supplementary_tables.R


# Make all plots
.PHONY: plots
plots: 	report_data_overview \
		report_distance_metric_correlation \
		output/plots/raw_data_overview.pdf \
		output/plots/phylogeny_congruence.pdf \
		output/plots/accuracy.pdf \
		output/plots/accuracy_supplement_sbinding.pdf \
		output/plots/accuracy_shedding.pdf \
		output/plots/accuracy_data_subsets.pdf \
		output/plots/accuracy_rhinolophid.pdf \
		output/plots/accuracy_non_ace2_all_features.pdf \
		output/plots/accuracy_non_ace2_phylogeny.pdf \
		output/plots/varimp_overview.pdf \
		output/plots/varimp_supplement.pdf \
		output/plots/existing_predictions.pdf \
		output/plots/existing_predictions_supplement.pdf \
		output/plots/predictions_by_order_supplement.pdf \
		output/plots/predictions_by_family_supplement.pdf \
		output/plots/prediction_maps.png \
		output/plots/ace2_availability_map_supplement.png \
		output/plots/prediction_maps_by_order.png \
		output/plots/phylogeny_predictions_supplement.pdf \
		output/si_tables/supplement_training_data.xlsx
