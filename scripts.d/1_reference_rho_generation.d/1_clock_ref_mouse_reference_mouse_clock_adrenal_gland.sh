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

#Compute reference correlations
if [[ -s "$RESULTS_DIR/reference_mouse_clock_adrenal_gland_12.full_rho.tsv" && -s "$RESULTS_DIR/reference_mouse_clock_adrenal_gland_12.weighted_rho.csv" ]]; then
  echo "Skipping reference correlations: reference_mouse_clock_adrenal_gland 12; outputs already exist"
else
  python "$NCC_SCRIPT" \
    --input_reference_db_file "$REFERENCE_MANIFEST" \
    --reference_full_rho_info "$RESULTS_DIR/reference_mouse_clock_adrenal_gland_12.full_rho.tsv" \
    --reference_weighted_rho_matrix "$RESULTS_DIR/reference_mouse_clock_adrenal_gland_12.weighted_rho.csv" \
    --l_clock_genes_reference "Bmal1|Clock|Cry1|Cry2|Dbp|Npas2|Nr1d1|Nr1d2|Per1|Per2|Per3|Tef" \
    --reference_species Mouse \
    --input_target_db_file "$REFERENCE_MANIFEST" \
    --target_full_rho_info "$RESULTS_DIR/_tmp.reference_mouse_clock_adrenal_gland_12.target.tsv" \
    --reference_target_gene_file "$ORTHO_FILE" \
    --target_species Mouse \
    --generate_reference_correlations
fi
