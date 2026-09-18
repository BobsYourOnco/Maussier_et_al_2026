#!/bin/bash -l

#Nina's Circadian Clockwork script
NCC_SCRIPT="NCC_tool_v2c_deltaCCD_matchedNull.py"
#1:1 mouse:human orthologue file
ORTHO_FILE="mouse_human_11orthologue_ENSEMBL_mBMAL1toARNTL.tsv"
#Results folder
RESULTS_DIR="results.d"
#Manifest for reference dataset
REFERENCE_MANIFEST="manifest.d/reference_mouse_clock_adrenal_gland.tsv"
#Make output folder
mkdir -p "$RESULTS_DIR"

#Compute delta-CCD for contrasting variables and establish significance without and with matched background significance for 1) MYCN amplification status, 2) harmonised risk, 3) INSS stage, 4) gender, 5) outcome status, and 6) age at diagnosis.
if [[ ! -s "$RESULTS_DIR/reference_mouse_clock_adrenal_gland_12.weighted_rho.csv" ]]; then
  echo "Skipping gene tag 12 for reference_mouse_clock_adrenal_gland: missing reference matrix" >&2
else
  if [[ -s "$RESULTS_DIR/target_clock_all_12.reference_mouse_clock_adrenal_gland.MYCN.delta.tsv" && -s "$RESULTS_DIR/target_clock_all_12.reference_mouse_clock_adrenal_gland.MYCN.delta_draws.tsv" ]]; then
    echo "Skipping delta CCD MYCN: target_clock_all_12.reference_mouse_clock_adrenal_gland; outputs already exist"
  elif [[ -s "manifests.d/target_clock_all_12.tsv" ]]; then
    python "$SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/reference_mouse_clock_adrenal_gland_12.weighted_rho.csv" \
      --l_clock_genes_reference "Bmal1|Clock|Cry1|Cry2|Dbp|Npas2|Nr1d1|Nr1d2|Per1|Per2|Per3|Tef" \
      --reference_species Mouse \
      --input_target_db_file "manifests.d/target_clock_all_12.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species Human \
      --compute_delta_CCD \
      --contrast_variable "mycn_status_unified" \
      --group_a "amplified" \
      --group_b "not_amplified" \
      --delta_n_permutations 5000 \
      --delta_n_bootstraps 2000 \
      --delta_seed 42 \
      --delta_min_samples_per_group 10 \
      --delta_outfile_results "$RESULTS_DIR/target_clock_all_12.reference_mouse_clock_adrenal_gland.MYCN.delta.tsv" \
      --delta_outfile_draws "$RESULTS_DIR/target_clock_all_12.reference_mouse_clock_adrenal_gland.MYCN.delta_draws.tsv"
  fi

  if [[ -s "$RESULTS_DIR/target_clock_all_12.reference_mouse_clock_adrenal_gland.risk.delta.tsv" && -s "$RESULTS_DIR/target_clock_all_12.reference_mouse_clock_adrenal_gland.risk.delta_draws.tsv" ]]; then
    echo "Skipping delta CCD risk: target_clock_all_12.reference_mouse_clock_adrenal_gland; outputs already exist"
  elif [[ -s "manifests.d/target_clock_all_12.tsv" ]]; then
    python "$SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/reference_mouse_clock_adrenal_gland_12.weighted_rho.csv" \
      --l_clock_genes_reference "Bmal1|Clock|Cry1|Cry2|Dbp|Npas2|Nr1d1|Nr1d2|Per1|Per2|Per3|Tef" \
      --reference_species Mouse \
      --input_target_db_file "manifests.d/target_clock_all_12.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species Human \
      --compute_delta_CCD \
      --contrast_variable "risk_group_custom" \
      --group_a "high_risk" \
      --group_b "low_risk" \
      --delta_n_permutations 5000 \
      --delta_n_bootstraps 2000 \
      --delta_seed 42 \
      --delta_min_samples_per_group 10 \
      --delta_outfile_results "$RESULTS_DIR/target_clock_all_12.reference_mouse_clock_adrenal_gland.risk.delta.tsv" \
      --delta_outfile_draws "$RESULTS_DIR/target_clock_all_12.reference_mouse_clock_adrenal_gland.risk.delta_draws.tsv"
  fi

  if [[ -s "$RESULTS_DIR/target_clock_all_12.reference_mouse_clock_adrenal_gland.stage.delta.tsv" && -s "$RESULTS_DIR/target_clock_all_12.reference_mouse_clock_adrenal_gland.stage.delta_draws.tsv" ]]; then
    echo "Skipping delta CCD stage: target_clock_all_12.reference_mouse_clock_adrenal_gland; outputs already exist"
  elif [[ -s "manifests.d/target_clock_all_12.tsv" ]]; then
    python "$SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/reference_mouse_clock_adrenal_gland_12.weighted_rho.csv" \
      --l_clock_genes_reference "Bmal1|Clock|Cry1|Cry2|Dbp|Npas2|Nr1d1|Nr1d2|Per1|Per2|Per3|Tef" \
      --reference_species Mouse \
      --input_target_db_file "manifests.d/target_clock_all_12.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species Human \
      --compute_delta_CCD \
      --contrast_variable "inss_stage_5cat" \
      --group_a "st4" \
      --group_b "st1" \
      --delta_n_permutations 5000 \
      --delta_n_bootstraps 2000 \
      --delta_seed 42 \
      --delta_min_samples_per_group 10 \
      --delta_outfile_results "$RESULTS_DIR/target_clock_all_12.reference_mouse_clock_adrenal_gland.stage.delta.tsv" \
      --delta_outfile_draws "$RESULTS_DIR/target_clock_all_12.reference_mouse_clock_adrenal_gland.stage.delta_draws.tsv"
  fi

  if [[ -s "$RESULTS_DIR/target_clock_all_12.reference_mouse_clock_adrenal_gland.gender.delta.tsv" && -s "$RESULTS_DIR/target_clock_all_12.reference_mouse_clock_adrenal_gland.gender.delta_draws.tsv" ]]; then
    echo "Skipping delta CCD gender: target_clock_all_12.reference_mouse_clock_adrenal_gland; outputs already exist"
  elif [[ -s "manifests.d/target_clock_all_12.tsv" ]]; then
    python "$SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/reference_mouse_clock_adrenal_gland_12.weighted_rho.csv" \
      --l_clock_genes_reference "Bmal1|Clock|Cry1|Cry2|Dbp|Npas2|Nr1d1|Nr1d2|Per1|Per2|Per3|Tef" \
      --reference_species Mouse \
      --input_target_db_file "manifests.d/target_clock_all_12.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species Human \
      --compute_delta_CCD \
      --contrast_variable "sex_unified" \
      --group_a "male" \
      --group_b "female" \
      --delta_n_permutations 5000 \
      --delta_n_bootstraps 2000 \
      --delta_seed 42 \
      --delta_min_samples_per_group 10 \
      --delta_outfile_results "$RESULTS_DIR/target_clock_all_12.reference_mouse_clock_adrenal_gland.gender.delta.tsv" \
      --delta_outfile_draws "$RESULTS_DIR/target_clock_all_12.reference_mouse_clock_adrenal_gland.gender.delta_draws.tsv"
  fi

  if [[ -s "$RESULTS_DIR/target_clock_all_12.reference_mouse_clock_adrenal_gland.alive.delta.tsv" && -s "$RESULTS_DIR/target_clock_all_12.reference_mouse_clock_adrenal_gland.alive.delta_draws.tsv" ]]; then
    echo "Skipping delta CCD alive: target_clock_all_12.reference_mouse_clock_adrenal_gland; outputs already exist"
  elif [[ -s "manifests.d/target_clock_all_12.tsv" ]]; then
    python "$SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/reference_mouse_clock_adrenal_gland_12.weighted_rho.csv" \
      --l_clock_genes_reference "Bmal1|Clock|Cry1|Cry2|Dbp|Npas2|Nr1d1|Nr1d2|Per1|Per2|Per3|Tef" \
      --reference_species Mouse \
      --input_target_db_file "manifests.d/target_clock_all_12.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species Human \
      --compute_delta_CCD \
      --contrast_variable "alive_status" \
      --group_a "dead" \
      --group_b "alive" \
      --delta_n_permutations 5000 \
      --delta_n_bootstraps 2000 \
      --delta_seed 42 \
      --delta_min_samples_per_group 10 \
      --delta_outfile_results "$RESULTS_DIR/target_clock_all_12.reference_mouse_clock_adrenal_gland.alive.delta.tsv" \
      --delta_outfile_draws "$RESULTS_DIR/target_clock_all_12.reference_mouse_clock_adrenal_gland.alive.delta_draws.tsv"
  fi

  if [[ -s "$RESULTS_DIR/target_clock_all_12.reference_mouse_clock_adrenal_gland.age.delta.tsv" && -s "$RESULTS_DIR/target_clock_all_12.reference_mouse_clock_adrenal_gland.age.delta_draws.tsv" ]]; then
    echo "Skipping delta CCD age: target_clock_all_12.reference_mouse_clock_adrenal_gland; outputs already exist"
  elif [[ -s "manifests.d/target_clock_all_12.tsv" ]]; then
    python "$SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/reference_mouse_clock_adrenal_gland_12.weighted_rho.csv" \
      --l_clock_genes_reference "Bmal1|Clock|Cry1|Cry2|Dbp|Npas2|Nr1d1|Nr1d2|Per1|Per2|Per3|Tef" \
      --reference_species Mouse \
      --input_target_db_file "manifests.d/target_clock_all_12.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species Human \
      --compute_delta_CCD \
      --contrast_variable "age_group_18m" \
      --group_a "ge18m" \
      --group_b "lt18m" \
      --delta_n_permutations 5000 \
      --delta_n_bootstraps 2000 \
      --delta_seed 42 \
      --delta_min_samples_per_group 10 \
      --delta_outfile_results "$RESULTS_DIR/target_clock_all_12.reference_mouse_clock_adrenal_gland.age.delta.tsv" \
      --delta_outfile_draws "$RESULTS_DIR/target_clock_all_12.reference_mouse_clock_adrenal_gland.age.delta_draws.tsv"
  fi

  if [[ -s "$RESULTS_DIR/target_clock_age_12.reference_mouse_clock_adrenal_gland.matched.target_full_rho.tsv" && -s "$RESULTS_DIR/target_clock_age_12.reference_mouse_clock_adrenal_gland.matched.background.tsv" && -s "$RESULTS_DIR/target_clock_age_12.reference_mouse_clock_adrenal_gland.matched.CCD.tsv" && -s "$RESULTS_DIR/target_clock_age_12.reference_mouse_clock_adrenal_gland.matched.gene_summary.tsv" && -s "$RESULTS_DIR/target_clock_age_12.reference_mouse_clock_adrenal_gland.matched.gene_sets.tsv" ]]; then
    echo "Skipping matched-null CCD: target_clock_age_12.reference_mouse_clock_adrenal_gland; outputs already exist"
  elif [[ -s "manifests.d/target_clock_age_12.tsv" ]]; then
    python "$SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/reference_mouse_clock_adrenal_gland_12.weighted_rho.csv" \
      --l_clock_genes_reference "Bmal1|Clock|Cry1|Cry2|Dbp|Npas2|Nr1d1|Nr1d2|Per1|Per2|Per3|Tef" \
      --reference_species Mouse \
      --input_target_db_file "manifests.d/target_clock_age_12.tsv" \
      --target_full_rho_info "$RESULTS_DIR/target_clock_age_12.reference_mouse_clock_adrenal_gland.matched.target_full_rho.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species Human \
      --number_randomizations 1000 \
      --background_target_full_rho_info "$RESULTS_DIR/target_clock_age_12.reference_mouse_clock_adrenal_gland.matched.background.tsv" \
      --outfile_results "$RESULTS_DIR/target_clock_age_12.reference_mouse_clock_adrenal_gland.matched.CCD.tsv" \
      --generate_target_correlations \
      --generate_background_target_correlations \
      --compute_CCD \
      --background_matching_mode matched \
      --matching_mean_caliper 0.05 \
      --matching_variance_caliper 0.05 \
      --matching_expand_step 0.05 \
      --matching_max_caliper 1.0 \
      --matching_seed 42 \
      --matched_gene_summary_outfile "$RESULTS_DIR/target_clock_age_12.reference_mouse_clock_adrenal_gland.matched.gene_summary.tsv" \
      --matched_sets_outfile "$RESULTS_DIR/target_clock_age_12.reference_mouse_clock_adrenal_gland.matched.gene_sets.tsv"
  fi

  if [[ -s "$RESULTS_DIR/target_clock_alive_12.reference_mouse_clock_adrenal_gland.matched.target_full_rho.tsv" && -s "$RESULTS_DIR/target_clock_alive_12.reference_mouse_clock_adrenal_gland.matched.background.tsv" && -s "$RESULTS_DIR/target_clock_alive_12.reference_mouse_clock_adrenal_gland.matched.CCD.tsv" && -s "$RESULTS_DIR/target_clock_alive_12.reference_mouse_clock_adrenal_gland.matched.gene_summary.tsv" && -s "$RESULTS_DIR/target_clock_alive_12.reference_mouse_clock_adrenal_gland.matched.gene_sets.tsv" ]]; then
    echo "Skipping matched-null CCD: target_clock_alive_12.reference_mouse_clock_adrenal_gland; outputs already exist"
  elif [[ -s "manifests.d/target_clock_alive_12.tsv" ]]; then
    python "$SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/reference_mouse_clock_adrenal_gland_12.weighted_rho.csv" \
      --l_clock_genes_reference "Bmal1|Clock|Cry1|Cry2|Dbp|Npas2|Nr1d1|Nr1d2|Per1|Per2|Per3|Tef" \
      --reference_species Mouse \
      --input_target_db_file "manifests.d/target_clock_alive_12.tsv" \
      --target_full_rho_info "$RESULTS_DIR/target_clock_alive_12.reference_mouse_clock_adrenal_gland.matched.target_full_rho.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species Human \
      --number_randomizations 1000 \
      --background_target_full_rho_info "$RESULTS_DIR/target_clock_alive_12.reference_mouse_clock_adrenal_gland.matched.background.tsv" \
      --outfile_results "$RESULTS_DIR/target_clock_alive_12.reference_mouse_clock_adrenal_gland.matched.CCD.tsv" \
      --generate_target_correlations \
      --generate_background_target_correlations \
      --compute_CCD \
      --background_matching_mode matched \
      --matching_mean_caliper 0.05 \
      --matching_variance_caliper 0.05 \
      --matching_expand_step 0.05 \
      --matching_max_caliper 1.0 \
      --matching_seed 42 \
      --matched_gene_summary_outfile "$RESULTS_DIR/target_clock_alive_12.reference_mouse_clock_adrenal_gland.matched.gene_summary.tsv" \
      --matched_sets_outfile "$RESULTS_DIR/target_clock_alive_12.reference_mouse_clock_adrenal_gland.matched.gene_sets.tsv"
  fi

  if [[ -s "$RESULTS_DIR/target_clock_gender_12.reference_mouse_clock_adrenal_gland.matched.target_full_rho.tsv" && -s "$RESULTS_DIR/target_clock_gender_12.reference_mouse_clock_adrenal_gland.matched.background.tsv" && -s "$RESULTS_DIR/target_clock_gender_12.reference_mouse_clock_adrenal_gland.matched.CCD.tsv" && -s "$RESULTS_DIR/target_clock_gender_12.reference_mouse_clock_adrenal_gland.matched.gene_summary.tsv" && -s "$RESULTS_DIR/target_clock_gender_12.reference_mouse_clock_adrenal_gland.matched.gene_sets.tsv" ]]; then
    echo "Skipping matched-null CCD: target_clock_gender_12.reference_mouse_clock_adrenal_gland; outputs already exist"
  elif [[ -s "manifests.d/target_clock_gender_12.tsv" ]]; then
    python "$SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/reference_mouse_clock_adrenal_gland_12.weighted_rho.csv" \
      --l_clock_genes_reference "Bmal1|Clock|Cry1|Cry2|Dbp|Npas2|Nr1d1|Nr1d2|Per1|Per2|Per3|Tef" \
      --reference_species Mouse \
      --input_target_db_file "manifests.d/target_clock_gender_12.tsv" \
      --target_full_rho_info "$RESULTS_DIR/target_clock_gender_12.reference_mouse_clock_adrenal_gland.matched.target_full_rho.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species Human \
      --number_randomizations 1000 \
      --background_target_full_rho_info "$RESULTS_DIR/target_clock_gender_12.reference_mouse_clock_adrenal_gland.matched.background.tsv" \
      --outfile_results "$RESULTS_DIR/target_clock_gender_12.reference_mouse_clock_adrenal_gland.matched.CCD.tsv" \
      --generate_target_correlations \
      --generate_background_target_correlations \
      --compute_CCD \
      --background_matching_mode matched \
      --matching_mean_caliper 0.05 \
      --matching_variance_caliper 0.05 \
      --matching_expand_step 0.05 \
      --matching_max_caliper 1.0 \
      --matching_seed 42 \
      --matched_gene_summary_outfile "$RESULTS_DIR/target_clock_gender_12.reference_mouse_clock_adrenal_gland.matched.gene_summary.tsv" \
      --matched_sets_outfile "$RESULTS_DIR/target_clock_gender_12.reference_mouse_clock_adrenal_gland.matched.gene_sets.tsv"
  fi

  if [[ -s "$RESULTS_DIR/target_clock_mycn_12.reference_mouse_clock_adrenal_gland.matched.target_full_rho.tsv" && -s "$RESULTS_DIR/target_clock_mycn_12.reference_mouse_clock_adrenal_gland.matched.background.tsv" && -s "$RESULTS_DIR/target_clock_mycn_12.reference_mouse_clock_adrenal_gland.matched.CCD.tsv" && -s "$RESULTS_DIR/target_clock_mycn_12.reference_mouse_clock_adrenal_gland.matched.gene_summary.tsv" && -s "$RESULTS_DIR/target_clock_mycn_12.reference_mouse_clock_adrenal_gland.matched.gene_sets.tsv" ]]; then
    echo "Skipping matched-null CCD: target_clock_mycn_12.reference_mouse_clock_adrenal_gland; outputs already exist"
  elif [[ -s "manifests.d/target_clock_mycn_12.tsv" ]]; then
    python "$SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/reference_mouse_clock_adrenal_gland_12.weighted_rho.csv" \
      --l_clock_genes_reference "Bmal1|Clock|Cry1|Cry2|Dbp|Npas2|Nr1d1|Nr1d2|Per1|Per2|Per3|Tef" \
      --reference_species Mouse \
      --input_target_db_file "manifests.d/target_clock_mycn_12.tsv" \
      --target_full_rho_info "$RESULTS_DIR/target_clock_mycn_12.reference_mouse_clock_adrenal_gland.matched.target_full_rho.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species Human \
      --number_randomizations 1000 \
      --background_target_full_rho_info "$RESULTS_DIR/target_clock_mycn_12.reference_mouse_clock_adrenal_gland.matched.background.tsv" \
      --outfile_results "$RESULTS_DIR/target_clock_mycn_12.reference_mouse_clock_adrenal_gland.matched.CCD.tsv" \
      --generate_target_correlations \
      --generate_background_target_correlations \
      --compute_CCD \
      --background_matching_mode matched \
      --matching_mean_caliper 0.05 \
      --matching_variance_caliper 0.05 \
      --matching_expand_step 0.05 \
      --matching_max_caliper 1.0 \
      --matching_seed 42 \
      --matched_gene_summary_outfile "$RESULTS_DIR/target_clock_mycn_12.reference_mouse_clock_adrenal_gland.matched.gene_summary.tsv" \
      --matched_sets_outfile "$RESULTS_DIR/target_clock_mycn_12.reference_mouse_clock_adrenal_gland.matched.gene_sets.tsv"
  fi

  if [[ -s "$RESULTS_DIR/target_clock_risk_12.reference_mouse_clock_adrenal_gland.matched.target_full_rho.tsv" && -s "$RESULTS_DIR/target_clock_risk_12.reference_mouse_clock_adrenal_gland.matched.background.tsv" && -s "$RESULTS_DIR/target_clock_risk_12.reference_mouse_clock_adrenal_gland.matched.CCD.tsv" && -s "$RESULTS_DIR/target_clock_risk_12.reference_mouse_clock_adrenal_gland.matched.gene_summary.tsv" && -s "$RESULTS_DIR/target_clock_risk_12.reference_mouse_clock_adrenal_gland.matched.gene_sets.tsv" ]]; then
    echo "Skipping matched-null CCD: target_clock_risk_12.reference_mouse_clock_adrenal_gland; outputs already exist"
  elif [[ -s "manifests.d/target_clock_risk_12.tsv" ]]; then
    python "$SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/reference_mouse_clock_adrenal_gland_12.weighted_rho.csv" \
      --l_clock_genes_reference "Bmal1|Clock|Cry1|Cry2|Dbp|Npas2|Nr1d1|Nr1d2|Per1|Per2|Per3|Tef" \
      --reference_species Mouse \
      --input_target_db_file "manifests.d/target_clock_risk_12.tsv" \
      --target_full_rho_info "$RESULTS_DIR/target_clock_risk_12.reference_mouse_clock_adrenal_gland.matched.target_full_rho.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species Human \
      --number_randomizations 1000 \
      --background_target_full_rho_info "$RESULTS_DIR/target_clock_risk_12.reference_mouse_clock_adrenal_gland.matched.background.tsv" \
      --outfile_results "$RESULTS_DIR/target_clock_risk_12.reference_mouse_clock_adrenal_gland.matched.CCD.tsv" \
      --generate_target_correlations \
      --generate_background_target_correlations \
      --compute_CCD \
      --background_matching_mode matched \
      --matching_mean_caliper 0.05 \
      --matching_variance_caliper 0.05 \
      --matching_expand_step 0.05 \
      --matching_max_caliper 1.0 \
      --matching_seed 42 \
      --matched_gene_summary_outfile "$RESULTS_DIR/target_clock_risk_12.reference_mouse_clock_adrenal_gland.matched.gene_summary.tsv" \
      --matched_sets_outfile "$RESULTS_DIR/target_clock_risk_12.reference_mouse_clock_adrenal_gland.matched.gene_sets.tsv"
  fi

  if [[ -s "$RESULTS_DIR/target_clock_stage_12.reference_mouse_clock_adrenal_gland.matched.target_full_rho.tsv" && -s "$RESULTS_DIR/target_clock_stage_12.reference_mouse_clock_adrenal_gland.matched.background.tsv" && -s "$RESULTS_DIR/target_clock_stage_12.reference_mouse_clock_adrenal_gland.matched.CCD.tsv" && -s "$RESULTS_DIR/target_clock_stage_12.reference_mouse_clock_adrenal_gland.matched.gene_summary.tsv" && -s "$RESULTS_DIR/target_clock_stage_12.reference_mouse_clock_adrenal_gland.matched.gene_sets.tsv" ]]; then
    echo "Skipping matched-null CCD: target_clock_stage_12.reference_mouse_clock_adrenal_gland; outputs already exist"
  elif [[ -s "manifests.d/target_clock_stage_12.tsv" ]]; then
    python "$SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/reference_mouse_clock_adrenal_gland_12.weighted_rho.csv" \
      --l_clock_genes_reference "Bmal1|Clock|Cry1|Cry2|Dbp|Npas2|Nr1d1|Nr1d2|Per1|Per2|Per3|Tef" \
      --reference_species Mouse \
      --input_target_db_file "manifests.d/target_clock_stage_12.tsv" \
      --target_full_rho_info "$RESULTS_DIR/target_clock_stage_12.reference_mouse_clock_adrenal_gland.matched.target_full_rho.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species Human \
      --number_randomizations 1000 \
      --background_target_full_rho_info "$RESULTS_DIR/target_clock_stage_12.reference_mouse_clock_adrenal_gland.matched.background.tsv" \
      --outfile_results "$RESULTS_DIR/target_clock_stage_12.reference_mouse_clock_adrenal_gland.matched.CCD.tsv" \
      --generate_target_correlations \
      --generate_background_target_correlations \
      --compute_CCD \
      --background_matching_mode matched \
      --matching_mean_caliper 0.05 \
      --matching_variance_caliper 0.05 \
      --matching_expand_step 0.05 \
      --matching_max_caliper 1.0 \
      --matching_seed 42 \
      --matched_gene_summary_outfile "$RESULTS_DIR/target_clock_stage_12.reference_mouse_clock_adrenal_gland.matched.gene_summary.tsv" \
      --matched_sets_outfile "$RESULTS_DIR/target_clock_stage_12.reference_mouse_clock_adrenal_gland.matched.gene_sets.tsv"
  fi

fi
