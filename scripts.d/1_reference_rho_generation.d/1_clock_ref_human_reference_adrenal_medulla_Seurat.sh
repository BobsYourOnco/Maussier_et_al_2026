#!/bin/bash -l

#Nina's Circadian Clockwork script
NCC_SCRIPT="NCC_tool_v2c_deltaCCD_matchedNull.py"
#1:1 mouse:human orthologue file
ORTHO_FILE="mouse_human_11orthologue_ENSEMBL_mBMAL1toARNTL.tsv"
#Results folder
RESULTS_DIR="results.d"
#Manifest for reference dataset
REFERENCE_MANIFEST="manifest.d/reference_adrenal_medulla_Seurat.tsv"
#Make output folder
mkdir -p "$RESULTS_DIR"

#Compute reference correlations
if [[ -s "$RESULTS_DIR/reference_adrenal_medulla_Seurat_12.full_rho.tsv" && -s "$RESULTS_DIR/reference_adrenal_medulla_Seurat_12.weighted_rho.csv" ]]; then
  echo "Skipping reference correlations: reference_adrenal_medulla_Seurat 12; outputs already exist"
else
  python "$NCC_SCRIPT" \
    --input_reference_db_file "$REFERENCE_MANIFEST" \
    --reference_full_rho_info "$RESULTS_DIR/reference_adrenal_medulla_Seurat_12.full_rho.tsv" \
    --reference_weighted_rho_matrix "$RESULTS_DIR/reference_adrenal_medulla_Seurat_12.weighted_rho.csv" \
    --l_clock_genes_reference "ARNTL|CLOCK|CRY1|CRY2|DBP|NPAS2|NR1D1|NR1D2|PER1|PER2|PER3|TEF" \
    --reference_species Human \
    --input_target_db_file "$REFERENCE_MANIFEST" \
    --target_full_rho_info "$RESULTS_DIR/_tmp.reference_adrenal_medulla_Seurat_12.target.tsv" \
    --reference_target_gene_file "$ORTHO_FILE" \
    --target_species Human \
    --generate_reference_correlations
fi

