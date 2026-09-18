#!/bin/bash -l

#Leave-one-out script 
LOGO_SCRIPT="leave_one_gene_out_deltaCCD.py"
#Nina's Circadian Clockwork script
NCC_SCRIPT="NCC_tool_v2c_deltaCCD_matchedNull.py"
#1:1 mouse:human orthologue file
ORTHO_FILE="mouse_human_11orthologue_ENSEMBL_mBMAL1toARNTL.tsv"
#Results folder
RESULTS_DIR="results.d"
#Manifest for reference dataset
REFERENCE_MANIFEST="manifest.d/reference_hg10_eAdrnl4v2_KEql10.tsv"
#Make output folder
mkdir -p "$RESULTS_DIR"

#Leave-one-out test for contrasting variables and establish significance for 1) MYCN amplification status, 2) harmonised risk, 3) INSS stage, 4) gender, 5) outcome status, and 6) age at diagnosis.
if [[ ! -s "$RESULTS_DIR/reference_hg10_eAdrnl4v2_KEql10_12.weighted_rho.csv" ]]; then
  echo "Skipping gene tag 12 for reference_hg10_eAdrnl4v2_KEql10: missing reference matrix" >&2
else
  if [[ -d "$RESULTS_DIR/leave_one_out/reference_hg10_eAdrnl4v2_KEql10/12/MYCN" && -n "$(find "$RESULTS_DIR/leave_one_out/reference_hg10_eAdrnl4v2_KEql10/12/MYCN" -type f -size +0c -print -quit 2>/dev/null)" ]]; then
    echo "Skipping leave-one-gene-out MYCN: reference_hg10_eAdrnl4v2_KEql10/12; output directory already has files"
  elif [[ -s "manifests.d/target_clock_all_12.tsv" ]]; then
    mkdir -p "$RESULTS_DIR/leave_one_out/reference_hg10_eAdrnl4v2_KEql10/12/MYCN"
    python "$LOGO_SCRIPT" \
      --ncc_script "$NCC_SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/reference_hg10_eAdrnl4v2_KEql10_12.weighted_rho.csv" \
      --l_clock_genes_reference "ARNTL|CLOCK|CRY1|CRY2|DBP|NPAS2|NR1D1|NR1D2|PER1|PER2|PER3|TEF" \
      --reference_species Human \
      --input_target_db_file "manifests.d/target_clock_all_12.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species Human \
      --contrast_variable "mycn_status_unified" \
      --group_a "amplified" \
      --group_b "not_amplified" \
      --delta_n_permutations 5000 \
      --delta_n_bootstraps 2000 \
      --delta_seed 42 \
      --delta_min_samples_per_group 10 \
      --include_full_panel \
      --save_draws \
      --output_dir "$RESULTS_DIR/leave_one_out/reference_hg10_eAdrnl4v2_KEql10/12/MYCN" \
      --output_prefix "reference_hg10_eAdrnl4v2_KEql10_12_MYCN"
  fi

  if [[ -d "$RESULTS_DIR/leave_one_out/reference_hg10_eAdrnl4v2_KEql10/12/risk" && -n "$(find "$RESULTS_DIR/leave_one_out/reference_hg10_eAdrnl4v2_KEql10/12/risk" -type f -size +0c -print -quit 2>/dev/null)" ]]; then
    echo "Skipping leave-one-gene-out risk: reference_hg10_eAdrnl4v2_KEql10/12; output directory already has files"
  elif [[ -s "manifests.d/target_clock_all_12.tsv" ]]; then
    mkdir -p "$RESULTS_DIR/leave_one_out/reference_hg10_eAdrnl4v2_KEql10/12/risk"
    python "$LOGO_SCRIPT" \
      --ncc_script "$NCC_SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/reference_hg10_eAdrnl4v2_KEql10_12.weighted_rho.csv" \
      --l_clock_genes_reference "ARNTL|CLOCK|CRY1|CRY2|DBP|NPAS2|NR1D1|NR1D2|PER1|PER2|PER3|TEF" \
      --reference_species Human \
      --input_target_db_file "manifests.d/target_clock_all_12.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species Human \
      --contrast_variable "risk_group_custom" \
      --group_a "high_risk" \
      --group_b "low_risk" \
      --delta_n_permutations 5000 \
      --delta_n_bootstraps 2000 \
      --delta_seed 42 \
      --delta_min_samples_per_group 10 \
      --include_full_panel \
      --save_draws \
      --output_dir "$RESULTS_DIR/leave_one_out/reference_hg10_eAdrnl4v2_KEql10/12/risk" \
      --output_prefix "reference_hg10_eAdrnl4v2_KEql10_12_risk"
  fi

  if [[ -d "$RESULTS_DIR/leave_one_out/reference_hg10_eAdrnl4v2_KEql10/12/stage" && -n "$(find "$RESULTS_DIR/leave_one_out/reference_hg10_eAdrnl4v2_KEql10/12/stage" -type f -size +0c -print -quit 2>/dev/null)" ]]; then
    echo "Skipping leave-one-gene-out stage: reference_hg10_eAdrnl4v2_KEql10/12; output directory already has files"
  elif [[ -s "manifests.d/target_clock_all_12.tsv" ]]; then
    mkdir -p "$RESULTS_DIR/leave_one_out/reference_hg10_eAdrnl4v2_KEql10/12/stage"
    python "$LOGO_SCRIPT" \
      --ncc_script "$NCC_SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/reference_hg10_eAdrnl4v2_KEql10_12.weighted_rho.csv" \
      --l_clock_genes_reference "ARNTL|CLOCK|CRY1|CRY2|DBP|NPAS2|NR1D1|NR1D2|PER1|PER2|PER3|TEF" \
      --reference_species Human \
      --input_target_db_file "manifests.d/target_clock_all_12.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species Human \
      --contrast_variable "inss_stage_5cat" \
      --group_a "st4" \
      --group_b "st1" \
      --delta_n_permutations 5000 \
      --delta_n_bootstraps 2000 \
      --delta_seed 42 \
      --delta_min_samples_per_group 10 \
      --include_full_panel \
      --save_draws \
      --output_dir "$RESULTS_DIR/leave_one_out/reference_hg10_eAdrnl4v2_KEql10/12/stage" \
      --output_prefix "reference_hg10_eAdrnl4v2_KEql10_12_stage"
  fi

  if [[ -d "$RESULTS_DIR/leave_one_out/reference_hg10_eAdrnl4v2_KEql10/12/gender" && -n "$(find "$RESULTS_DIR/leave_one_out/reference_hg10_eAdrnl4v2_KEql10/12/gender" -type f -size +0c -print -quit 2>/dev/null)" ]]; then
    echo "Skipping leave-one-gene-out gender: reference_hg10_eAdrnl4v2_KEql10/12; output directory already has files"
  elif [[ -s "manifests.d/target_clock_all_12.tsv" ]]; then
    mkdir -p "$RESULTS_DIR/leave_one_out/reference_hg10_eAdrnl4v2_KEql10/12/gender"
    python "$LOGO_SCRIPT" \
      --ncc_script "$NCC_SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/reference_hg10_eAdrnl4v2_KEql10_12.weighted_rho.csv" \
      --l_clock_genes_reference "ARNTL|CLOCK|CRY1|CRY2|DBP|NPAS2|NR1D1|NR1D2|PER1|PER2|PER3|TEF" \
      --reference_species Human \
      --input_target_db_file "manifests.d/target_clock_all_12.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species Human \
      --contrast_variable "sex_unified" \
      --group_a "male" \
      --group_b "female" \
      --delta_n_permutations 5000 \
      --delta_n_bootstraps 2000 \
      --delta_seed 42 \
      --delta_min_samples_per_group 10 \
      --include_full_panel \
      --save_draws \
      --output_dir "$RESULTS_DIR/leave_one_out/reference_hg10_eAdrnl4v2_KEql10/12/gender" \
      --output_prefix "reference_hg10_eAdrnl4v2_KEql10_12_gender"
  fi

  if [[ -d "$RESULTS_DIR/leave_one_out/reference_hg10_eAdrnl4v2_KEql10/12/alive" && -n "$(find "$RESULTS_DIR/leave_one_out/reference_hg10_eAdrnl4v2_KEql10/12/alive" -type f -size +0c -print -quit 2>/dev/null)" ]]; then
    echo "Skipping leave-one-gene-out alive: reference_hg10_eAdrnl4v2_KEql10/12; output directory already has files"
  elif [[ -s "manifests.d/target_clock_all_12.tsv" ]]; then
    mkdir -p "$RESULTS_DIR/leave_one_out/reference_hg10_eAdrnl4v2_KEql10/12/alive"
    python "$LOGO_SCRIPT" \
      --ncc_script "$NCC_SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/reference_hg10_eAdrnl4v2_KEql10_12.weighted_rho.csv" \
      --l_clock_genes_reference "ARNTL|CLOCK|CRY1|CRY2|DBP|NPAS2|NR1D1|NR1D2|PER1|PER2|PER3|TEF" \
      --reference_species Human \
      --input_target_db_file "manifests.d/target_clock_all_12.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species Human \
      --contrast_variable "alive_status" \
      --group_a "dead" \
      --group_b "alive" \
      --delta_n_permutations 5000 \
      --delta_n_bootstraps 2000 \
      --delta_seed 42 \
      --delta_min_samples_per_group 10 \
      --include_full_panel \
      --save_draws \
      --output_dir "$RESULTS_DIR/leave_one_out/reference_hg10_eAdrnl4v2_KEql10/12/alive" \
      --output_prefix "reference_hg10_eAdrnl4v2_KEql10_12_alive"
  fi

  if [[ -d "$RESULTS_DIR/leave_one_out/reference_hg10_eAdrnl4v2_KEql10/12/age" && -n "$(find "$RESULTS_DIR/leave_one_out/reference_hg10_eAdrnl4v2_KEql10/12/age" -type f -size +0c -print -quit 2>/dev/null)" ]]; then
    echo "Skipping leave-one-gene-out age: reference_hg10_eAdrnl4v2_KEql10/12; output directory already has files"
  elif [[ -s "manifests.d/target_clock_all_12.tsv" ]]; then
    mkdir -p "$RESULTS_DIR/leave_one_out/reference_hg10_eAdrnl4v2_KEql10/12/age"
    python "$LOGO_SCRIPT" \
      --ncc_script "$NCC_SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/reference_hg10_eAdrnl4v2_KEql10_12.weighted_rho.csv" \
      --l_clock_genes_reference "ARNTL|CLOCK|CRY1|CRY2|DBP|NPAS2|NR1D1|NR1D2|PER1|PER2|PER3|TEF" \
      --reference_species Human \
      --input_target_db_file "manifests.d/target_clock_all_12.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species Human \
      --contrast_variable "age_group_18m" \
      --group_a "ge18m" \
      --group_b "lt18m" \
      --delta_n_permutations 5000 \
      --delta_n_bootstraps 2000 \
      --delta_seed 42 \
      --delta_min_samples_per_group 10 \
      --include_full_panel \
      --save_draws \
      --output_dir "$RESULTS_DIR/leave_one_out/reference_hg10_eAdrnl4v2_KEql10/12/age" \
      --output_prefix "reference_hg10_eAdrnl4v2_KEql10_12_age"
  fi

fi

