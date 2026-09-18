#!/bin/bash -l

#Nina's Circadian Clockwork script
NCC_SCRIPT="NCC_tool_v2c_deltaCCD_matchedNull.py"
#1:1 mouse:human orthologue file
ORTHO_FILE="mouse_human_11orthologue_ENSEMBL_mBMAL1toARNTL.tsv"
#Results folder
RESULTS_DIR="results.d"
#Manifest for reference dataset
REFERENCE_MANIFEST="manifest.d/reference_hg10_eAdrnl4v2_KEql10.tsv"
#Number of randomizations
N_RANDOM="1000"
#Make output folder
mkdir -p "$RESULTS_DIR"

#Generate reference correlations
if [[ ! -s "$RESULTS_DIR/reference_hg10_eAdrnl4v2_KEql10_12.weighted_rho.csv" ]]; then
  if [[ -s "$RESULTS_DIR/reference_hg10_eAdrnl4v2_KEql10_12.full_rho.tsv" && -s "$RESULTS_DIR/reference_hg10_eAdrnl4v2_KEql10_12.weighted_rho.csv" ]]; then
    echo "Reference matrix exists for reference_hg10_eAdrnl4v2_KEql10 12"
  else
    python "$NCC_SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_full_rho_info "$RESULTS_DIR/reference_hg10_eAdrnl4v2_KEql10_12.full_rho.tsv" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/reference_hg10_eAdrnl4v2_KEql10_12.weighted_rho.csv" \
      --l_clock_genes_reference "ARNTL|CLOCK|CRY1|CRY2|DBP|NPAS2|NR1D1|NR1D2|PER1|PER2|PER3|TEF" \
      --reference_species Human \
      --input_target_db_file "$REFERENCE_MANIFEST" \
      --target_full_rho_info "$RESULTS_DIR/_tmp.reference_hg10_eAdrnl4v2_KEql10_12.target.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species Human \
      --generate_reference_correlations
  fi
fi

#Compute CCD and significance for 1) age at diagnosis, 2) outcome status, 3) all samples independently of the variable, 4) gender, 5) MYCN amplification status, 6) harmonised risk, and 7) INSS stage.
if [[ ! -s "$RESULTS_DIR/reference_hg10_eAdrnl4v2_KEql10_12.weighted_rho.csv" ]]; then
  echo "Skipping gene tag 12 for reference_hg10_eAdrnl4v2_KEql10: missing reference matrix" >&2
else
  if [[ -s "$RESULTS_DIR/target_clock_age_12.reference_hg10_eAdrnl4v2_KEql10.target_full_rho.tsv" && -s "$RESULTS_DIR/target_clock_age_12.reference_hg10_eAdrnl4v2_KEql10.background.tsv" && -s "$RESULTS_DIR/target_clock_age_12.reference_hg10_eAdrnl4v2_KEql10.CCD.tsv" ]]; then
    echo "Skipping unmatched CCD: target_clock_age_12.reference_hg10_eAdrnl4v2_KEql10; outputs already exist"
  elif [[ -s "manifest.d/target_clock_age_12.tsv" ]]; then
    python "$NCC_SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/reference_hg10_eAdrnl4v2_KEql10_12.weighted_rho.csv" \
      --l_clock_genes_reference "ARNTL|CLOCK|CRY1|CRY2|DBP|NPAS2|NR1D1|NR1D2|PER1|PER2|PER3|TEF" \
      --reference_species Human \
      --input_target_db_file "manifest.d/target_clock_age_12.tsv" \
      --target_full_rho_info "$RESULTS_DIR/target_clock_age_12.reference_hg10_eAdrnl4v2_KEql10.target_full_rho.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species Human \
      --number_randomizations "$N_RANDOM" \
      --background_target_full_rho_info "$RESULTS_DIR/target_clock_age_12.reference_hg10_eAdrnl4v2_KEql10.background.tsv" \
      --outfile_results "$RESULTS_DIR/target_clock_age_12.reference_hg10_eAdrnl4v2_KEql10.CCD.tsv" \
      --generate_target_correlations \
      --generate_background_target_correlations \
      --compute_CCD \
      --background_matching_mode unmatched
  fi

  if [[ -s "$RESULTS_DIR/target_clock_alive_12.reference_hg10_eAdrnl4v2_KEql10.target_full_rho.tsv" && -s "$RESULTS_DIR/target_clock_alive_12.reference_hg10_eAdrnl4v2_KEql10.background.tsv" && -s "$RESULTS_DIR/target_clock_alive_12.reference_hg10_eAdrnl4v2_KEql10.CCD.tsv" ]]; then
    echo "Skipping unmatched CCD: target_clock_alive_12.reference_hg10_eAdrnl4v2_KEql10; outputs already exist"
  elif [[ -s "manifest.d/target_clock_alive_12.tsv" ]]; then
    python "$NCC_SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/reference_hg10_eAdrnl4v2_KEql10_12.weighted_rho.csv" \
      --l_clock_genes_reference "ARNTL|CLOCK|CRY1|CRY2|DBP|NPAS2|NR1D1|NR1D2|PER1|PER2|PER3|TEF" \
      --reference_species Human \
      --input_target_db_file "manifest.d/target_clock_alive_12.tsv" \
      --target_full_rho_info "$RESULTS_DIR/target_clock_alive_12.reference_hg10_eAdrnl4v2_KEql10.target_full_rho.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species Human \
      --number_randomizations "$N_RANDOM" \
      --background_target_full_rho_info "$RESULTS_DIR/target_clock_alive_12.reference_hg10_eAdrnl4v2_KEql10.background.tsv" \
      --outfile_results "$RESULTS_DIR/target_clock_alive_12.reference_hg10_eAdrnl4v2_KEql10.CCD.tsv" \
      --generate_target_correlations \
      --generate_background_target_correlations \
      --compute_CCD \
      --background_matching_mode unmatched
  fi

  if [[ -s "$RESULTS_DIR/target_clock_all_12.reference_hg10_eAdrnl4v2_KEql10.target_full_rho.tsv" && -s "$RESULTS_DIR/target_clock_all_12.reference_hg10_eAdrnl4v2_KEql10.background.tsv" && -s "$RESULTS_DIR/target_clock_all_12.reference_hg10_eAdrnl4v2_KEql10.CCD.tsv" ]]; then
    echo "Skipping unmatched CCD: target_clock_all_12.reference_hg10_eAdrnl4v2_KEql10; outputs already exist"
  elif [[ -s "manifest.d/target_clock_all_12.tsv" ]]; then
    python "$NCC_SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/reference_hg10_eAdrnl4v2_KEql10_12.weighted_rho.csv" \
      --l_clock_genes_reference "ARNTL|CLOCK|CRY1|CRY2|DBP|NPAS2|NR1D1|NR1D2|PER1|PER2|PER3|TEF" \
      --reference_species Human \
      --input_target_db_file "manifest.d/target_clock_all_12.tsv" \
      --target_full_rho_info "$RESULTS_DIR/target_clock_all_12.reference_hg10_eAdrnl4v2_KEql10.target_full_rho.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species Human \
      --number_randomizations "$N_RANDOM" \
      --background_target_full_rho_info "$RESULTS_DIR/target_clock_all_12.reference_hg10_eAdrnl4v2_KEql10.background.tsv" \
      --outfile_results "$RESULTS_DIR/target_clock_all_12.reference_hg10_eAdrnl4v2_KEql10.CCD.tsv" \
      --generate_target_correlations \
      --generate_background_target_correlations \
      --compute_CCD \
      --background_matching_mode unmatched
  fi

  if [[ -s "$RESULTS_DIR/target_clock_gender_12.reference_hg10_eAdrnl4v2_KEql10.target_full_rho.tsv" && -s "$RESULTS_DIR/target_clock_gender_12.reference_hg10_eAdrnl4v2_KEql10.background.tsv" && -s "$RESULTS_DIR/target_clock_gender_12.reference_hg10_eAdrnl4v2_KEql10.CCD.tsv" ]]; then
    echo "Skipping unmatched CCD: target_clock_gender_12.reference_hg10_eAdrnl4v2_KEql10; outputs already exist"
  elif [[ -s "manifest.d/target_clock_gender_12.tsv" ]]; then
    python "$NCC_SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/reference_hg10_eAdrnl4v2_KEql10_12.weighted_rho.csv" \
      --l_clock_genes_reference "ARNTL|CLOCK|CRY1|CRY2|DBP|NPAS2|NR1D1|NR1D2|PER1|PER2|PER3|TEF" \
      --reference_species Human \
      --input_target_db_file "manifest.d/target_clock_gender_12.tsv" \
      --target_full_rho_info "$RESULTS_DIR/target_clock_gender_12.reference_hg10_eAdrnl4v2_KEql10.target_full_rho.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species Human \
      --number_randomizations "$N_RANDOM" \
      --background_target_full_rho_info "$RESULTS_DIR/target_clock_gender_12.reference_hg10_eAdrnl4v2_KEql10.background.tsv" \
      --outfile_results "$RESULTS_DIR/target_clock_gender_12.reference_hg10_eAdrnl4v2_KEql10.CCD.tsv" \
      --generate_target_correlations \
      --generate_background_target_correlations \
      --compute_CCD \
      --background_matching_mode unmatched
  fi

  if [[ -s "$RESULTS_DIR/target_clock_mycn_12.reference_hg10_eAdrnl4v2_KEql10.target_full_rho.tsv" && -s "$RESULTS_DIR/target_clock_mycn_12.reference_hg10_eAdrnl4v2_KEql10.background.tsv" && -s "$RESULTS_DIR/target_clock_mycn_12.reference_hg10_eAdrnl4v2_KEql10.CCD.tsv" ]]; then
    echo "Skipping unmatched CCD: target_clock_mycn_12.reference_hg10_eAdrnl4v2_KEql10; outputs already exist"
  elif [[ -s "manifest.d/target_clock_mycn_12.tsv" ]]; then
    python "$NCC_SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/reference_hg10_eAdrnl4v2_KEql10_12.weighted_rho.csv" \
      --l_clock_genes_reference "ARNTL|CLOCK|CRY1|CRY2|DBP|NPAS2|NR1D1|NR1D2|PER1|PER2|PER3|TEF" \
      --reference_species Human \
      --input_target_db_file "manifest.d/target_clock_mycn_12.tsv" \
      --target_full_rho_info "$RESULTS_DIR/target_clock_mycn_12.reference_hg10_eAdrnl4v2_KEql10.target_full_rho.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species Human \
      --number_randomizations "$N_RANDOM" \
      --background_target_full_rho_info "$RESULTS_DIR/target_clock_mycn_12.reference_hg10_eAdrnl4v2_KEql10.background.tsv" \
      --outfile_results "$RESULTS_DIR/target_clock_mycn_12.reference_hg10_eAdrnl4v2_KEql10.CCD.tsv" \
      --generate_target_correlations \
      --generate_background_target_correlations \
      --compute_CCD \
      --background_matching_mode unmatched
  fi

  if [[ -s "$RESULTS_DIR/target_clock_risk_12.reference_hg10_eAdrnl4v2_KEql10.target_full_rho.tsv" && -s "$RESULTS_DIR/target_clock_risk_12.reference_hg10_eAdrnl4v2_KEql10.background.tsv" && -s "$RESULTS_DIR/target_clock_risk_12.reference_hg10_eAdrnl4v2_KEql10.CCD.tsv" ]]; then
    echo "Skipping unmatched CCD: target_clock_risk_12.reference_hg10_eAdrnl4v2_KEql10; outputs already exist"
  elif [[ -s "manifest.d/target_clock_risk_12.tsv" ]]; then
    python "$NCC_SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/reference_hg10_eAdrnl4v2_KEql10_12.weighted_rho.csv" \
      --l_clock_genes_reference "ARNTL|CLOCK|CRY1|CRY2|DBP|NPAS2|NR1D1|NR1D2|PER1|PER2|PER3|TEF" \
      --reference_species Human \
      --input_target_db_file "manifest.d/target_clock_risk_12.tsv" \
      --target_full_rho_info "$RESULTS_DIR/target_clock_risk_12.reference_hg10_eAdrnl4v2_KEql10.target_full_rho.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species Human \
      --number_randomizations "$N_RANDOM" \
      --background_target_full_rho_info "$RESULTS_DIR/target_clock_risk_12.reference_hg10_eAdrnl4v2_KEql10.background.tsv" \
      --outfile_results "$RESULTS_DIR/target_clock_risk_12.reference_hg10_eAdrnl4v2_KEql10.CCD.tsv" \
      --generate_target_correlations \
      --generate_background_target_correlations \
      --compute_CCD \
      --background_matching_mode unmatched
  fi

  if [[ -s "$RESULTS_DIR/target_clock_stage_12.reference_hg10_eAdrnl4v2_KEql10.target_full_rho.tsv" && -s "$RESULTS_DIR/target_clock_stage_12.reference_hg10_eAdrnl4v2_KEql10.background.tsv" && -s "$RESULTS_DIR/target_clock_stage_12.reference_hg10_eAdrnl4v2_KEql10.CCD.tsv" ]]; then
    echo "Skipping unmatched CCD: target_clock_stage_12.reference_hg10_eAdrnl4v2_KEql10; outputs already exist"
  elif [[ -s "manifest.d/target_clock_stage_12.tsv" ]]; then
    python "$NCC_SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/reference_hg10_eAdrnl4v2_KEql10_12.weighted_rho.csv" \
      --l_clock_genes_reference "ARNTL|CLOCK|CRY1|CRY2|DBP|NPAS2|NR1D1|NR1D2|PER1|PER2|PER3|TEF" \
      --reference_species Human \
      --input_target_db_file "manifest.d/target_clock_stage_12.tsv" \
      --target_full_rho_info "$RESULTS_DIR/target_clock_stage_12.reference_hg10_eAdrnl4v2_KEql10.target_full_rho.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species Human \
      --number_randomizations "$N_RANDOM" \
      --background_target_full_rho_info "$RESULTS_DIR/target_clock_stage_12.reference_hg10_eAdrnl4v2_KEql10.background.tsv" \
      --outfile_results "$RESULTS_DIR/target_clock_stage_12.reference_hg10_eAdrnl4v2_KEql10.CCD.tsv" \
      --generate_target_correlations \
      --generate_background_target_correlations \
      --compute_CCD \
      --background_matching_mode unmatched
  fi

fi
