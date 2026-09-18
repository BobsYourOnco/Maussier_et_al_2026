Purpose
-------
This README describes the Bash workflow used to perform leave-one-gene-out
(LOO) sensitivity analysis for direct Delta-CCD contrasts using a 12-gene
circadian clock panel.

The workflow is reference-agnostic: the reference dataset may be Human or
Mouse, provided that the reference name, species, 12-gene clock list, reference
manifest, and reference-to-target gene mapping are configured consistently.
The neuroblastoma target data are Human.

The purpose of LOO analysis is to test whether a Delta-CCD result depends
disproportionately on any single clock gene. The workflow retains the full
12-gene result and repeats the analysis after omitting each gene in turn.

Analysis overview
-----------------
For each clinical contrast:

    1. Compute the full-panel 12-gene Delta-CCD result.
    2. Omit one clock gene.
    3. Recompute Delta-CCD.
    4. Repeat until each of the 12 genes has been omitted once.
    5. Save the resulting statistics and resampling draws.

This gives 13 analysis states per contrast:

    1 full-panel result
    12 leave-one-gene-out results


Clinical contrasts
------------------
LOO Delta-CCD is performed for six variables:

    MYCN amplification:
        mycn_status_unified
        amplified vs not_amplified

    Harmonized risk:
        risk_group_custom
        high_risk vs low_risk

    INSS stage:
        inss_stage_5cat
        st4 vs st1

    Gender / sex:
        sex_unified
        male vs female

    Outcome status:
        alive_status
        dead vs alive

    Age at diagnosis:
        age_group_18m
        ge18m vs lt18m


Required scripts
----------------
The workflow uses two Python scripts:

    LOGO_SCRIPT
        Leave-one-gene-out wrapper.

        Example:
            leave_one_gene_out_deltaCCD.py

    NCC_SCRIPT
        NCC/CCD script called internally by LOGO_SCRIPT.

        Example:
            NCC_tool_v2c_deltaCCD_matchedNull.py


Required variables
------------------
Define:

    LOGO_SCRIPT
        Path to the leave-one-gene-out wrapper.

    NCC_SCRIPT
        Path to the underlying NCC/CCD analysis script.

    ORTHO_FILE
        Reference-to-target gene mapping file.

    RESULTS_DIR
        Output directory.

    REFERENCE_MANIFEST
        Manifest describing the selected reference dataset.

For a generic wrapper, also define:

    REFERENCE_NAME
    REFERENCE_SPECIES
    REFERENCE_GENES
    TARGET_SPECIES


Reference-specific settings
---------------------------
Example Mouse reference:

    REFERENCE_NAME="reference_mouse_clock"
    REFERENCE_SPECIES="Mouse"
    REFERENCE_GENES="Bmal1|Clock|Cry1|Cry2|Dbp|Npas2|Nr1d1|Nr1d2|Per1|Per2|Per3|Tef"
    TARGET_SPECIES="Human"

Example Human reference:

    REFERENCE_NAME="reference_human_adrenal_gland"
    REFERENCE_SPECIES="Human"
    REFERENCE_GENES="ARNTL|CLOCK|CRY1|CRY2|DBP|NPAS2|NR1D1|NR1D2|PER1|PER2|PER3|TEF"
    TARGET_SPECIES="Human"


Reference matrix requirement
----------------------------
Before any LOO analysis is run, the workflow checks:

    $RESULTS_DIR/${REFERENCE_NAME}_12.weighted_rho.csv

If this file is absent or empty, all LOO analyses are skipped.

Conceptually:

    if weighted reference matrix is missing:
        skip leave-one-gene-out analyses
    else:
        continue


Required target manifest
------------------------
All six contrasts use:

    manifests.d/target_clock_all_12.tsv

This manifest must contain the six clinical variables listed above.


LOO settings
------------
Each analysis uses:

    --delta_n_permutations 5000
    --delta_n_bootstraps 2000
    --delta_seed 42
    --delta_min_samples_per_group 10
    --include_full_panel
    --save_draws

Meaning:

    5000 permutations
        permutation-based significance testing.

    2000 bootstrap resamples
        uncertainty estimation / confidence intervals.

    seed = 42
        reproducibility.

    minimum 10 samples per group
        required for a contrast.

    --include_full_panel
        saves the complete 12-gene result together with LOO results.

    --save_draws
        saves the resampling draws generated during the analysis.


Delta-CCD orientation
---------------------
For each full-panel or LOO state:

    Delta-CCD = CCD(Group A) - CCD(Group B)

Therefore:

    Delta-CCD > 0
        Group A has larger CCD than Group B.

    Delta-CCD < 0
        Group A has smaller CCD than Group B.

Interpretation must always follow the specified Group A / Group B orientation.


Output organization
-------------------
General directory structure:

    $RESULTS_DIR/
        leave_one_out/
            ${REFERENCE_NAME}/
                12/
                    MYCN/
                    risk/
                    stage/
                    gender/
                    alive/
                    age/

Suggested output prefixes:

    ${REFERENCE_NAME}_12_MYCN
    ${REFERENCE_NAME}_12_risk
    ${REFERENCE_NAME}_12_stage
    ${REFERENCE_NAME}_12_gender
    ${REFERENCE_NAME}_12_alive
    ${REFERENCE_NAME}_12_age


Species configuration
---------------------
The reference can be Human or Mouse.

Mouse reference -> Human target:

    --reference_species Mouse
    --target_species Human

Human reference -> Human target:

    --reference_species Human
    --target_species Human

The reference gene list must use symbols appropriate for the reference species.
ORTHO_FILE must support the selected reference-target species pair.


Typical command structure
-------------------------
A generic call follows:

    python "$LOGO_SCRIPT" \
      --ncc_script "$NCC_SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/${REFERENCE_NAME}_12.weighted_rho.csv" \
      --l_clock_genes_reference "$REFERENCE_GENES" \
      --reference_species "$REFERENCE_SPECIES" \
      --input_target_db_file "manifests.d/target_clock_all_12.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species "$TARGET_SPECIES" \
      --contrast_variable "<clinical_variable>" \
      --group_a "<group_a>" \
      --group_b "<group_b>" \
      --delta_n_permutations 5000 \
      --delta_n_bootstraps 2000 \
      --delta_seed 42 \
      --delta_min_samples_per_group 10 \
      --include_full_panel \
      --save_draws \
      --output_dir "$RESULTS_DIR/leave_one_out/${REFERENCE_NAME}/12/<contrast>" \
      --output_prefix "${REFERENCE_NAME}_12_<contrast>"


Skip / resume logic
-------------------
For each contrast, the current Bash script checks whether its output directory
already exists and contains at least one non-empty file.

Conceptually:

    if output directory exists
       AND contains at least one non-empty file:
        skip contrast

    else if target_clock_all_12.tsv exists:
        create output directory
        run leave-one-gene-out analysis

Important: this does not verify that the entire LOO run completed. A partially
completed directory containing even one non-empty file will cause that contrast
to be skipped on the next run.


Recommended safer completion check
----------------------------------
For robust resume behavior, check for one or more specific final summary files
written by LOGO_SCRIPT instead of checking for any file in the directory.

Current logic:

    find <output_dir> -type f -size +0c -print -quit

Safer logic:

    check explicitly for the final LOO summary/result file(s)

This avoids treating an interrupted run as complete.


Workflow summary
----------------
The complete workflow is:

    1. Define the LOO wrapper and NCC/CCD script.
    2. Select a Human or Mouse reference.
    3. Define reference name, species, gene list, manifest, and mapping file.
    4. Confirm that the 12-gene weighted reference matrix exists.
    5. Confirm that manifests.d/target_clock_all_12.tsv exists.
    6. Run LOO Delta-CCD for:
         MYCN
         risk
         stage
         gender
         alive/dead
         age
    7. For each contrast:
         compute the full 12-gene result;
         omit each gene in turn;
         run 5000 permutations;
         run 2000 bootstraps;
         require at least 10 samples per group;
         save draws and result files.


Pseudocode
----------
    REFERENCE_NAME    = selected reference
    REFERENCE_SPECIES = Human or Mouse
    REFERENCE_GENES   = species-appropriate 12-gene clock panel
    TARGET_SPECIES    = Human

    if reference weighted-rho matrix is missing:

        skip all LOO analyses

    else:

        for CONTRAST in MYCN, risk, stage, gender, alive, age:

            OUTPUT_DIR =
                results/leave_one_out/REFERENCE_NAME/12/CONTRAST

            if OUTPUT_DIR exists
               AND contains at least one non-empty file:

                skip CONTRAST

            else if target_clock_all_12 manifest exists:

                create OUTPUT_DIR

                compute full 12-gene Delta-CCD

                for each gene in the 12-gene panel:

                    omit that gene

                    recompute Delta-CCD

                    run 5000 permutations

                    run 2000 bootstraps

                    save results and draws


Interpretation
--------------
The LOO analysis is a robustness analysis, not a separate clinical contrast.

A Delta-CCD result is more robust when:

    the direction remains similar after individual genes are removed;

    the effect magnitude is not driven by one specific gene;

    statistical support is retained across most LOO configurations.

A strong change after removal of one gene suggests that the full-panel result
depends disproportionately on that gene.

The full-panel Delta-CCD result should remain the primary estimate; the LOO
results provide sensitivity evidence.


Troubleshooting
---------------

Reference matrix is missing

    ls -lh "$RESULTS_DIR/${REFERENCE_NAME}_12.weighted_rho.csv"

The Bash -s test requires the file to exist and have size greater than zero.


Target manifest is missing

    ls -lh manifests.d/target_clock_all_12.tsv


LOO analysis is skipped unexpectedly

    find "$RESULTS_DIR/leave_one_out/${REFERENCE_NAME}/12/<contrast>" \
      -type f -size +0c -print

If any non-empty file is present, the current script skips that contrast.


Interrupted run will not resume

Because the current skip logic checks for any non-empty file, a partially
completed directory may be interpreted as complete. Remove/rename the partial
directory or change the Bash logic to check specific final output files.


LOO wrapper cannot be found

    echo "$LOGO_SCRIPT"
    ls -lh "$LOGO_SCRIPT"


NCC/CCD script cannot be found

    echo "$NCC_SCRIPT"
    ls -lh "$NCC_SCRIPT"


Gene symbols are not found

Confirm that REFERENCE_GENES matches REFERENCE_SPECIES and that all 12 genes
are present in the selected reference dataset.


Reference-target mapping fails

    ls -lh "$ORTHO_FILE"

Confirm that the mapping file supports the selected Human-Human or Mouse-Human
configuration.

