Purpose
-------
This Bash workflow generates a 12-gene circadian reference correlation matrix
and then computes clock correlation distance (CCD) and significance for
multiple neuroblastoma clinical groupings using an unmatched random-gene
background.

The workflow is reference-agnostic: the reference dataset may be derived from
either Human or Mouse tissue, provided that the reference name, species,
12-gene clock list, manifest, and orthology mapping are configured consistently.

Typical examples include:

    reference_mouse_clock
    reference_mouse_clock_adrenal_gland
    reference_human_adrenal_gland

The target neuroblastoma data are Human.


Reference-specific settings
---------------------------
For each reference, define or substitute the following values:

    REFERENCE_NAME
        Short reference identifier used in output filenames, for example:

            reference_mouse_clock
            reference_mouse_clock_adrenal_gland
            reference_human_adrenal_gland

    REFERENCE_SPECIES
        Species of the reference expression data:

            Mouse

        or:

            Human

    REFERENCE_GENES
        Pipe-separated 12-gene circadian panel using the gene symbols expected
        for the reference species.

        Example for Mouse:

            Bmal1|Clock|Cry1|Cry2|Dbp|Npas2|Nr1d1|Nr1d2|Per1|Per2|Per3|Tef

        Example for Human:

            ARNTL|CLOCK|CRY1|CRY2|DBP|NPAS2|NR1D1|NR1D2|PER1|PER2|PER3|TEF

The reference species and gene-symbol convention must match the expression
matrix described by the reference manifest.


Required variables
------------------
The Bash environment must define:

    RESULTS_DIR
        Directory where reference and target CCD outputs are written.

    NCC_SCRIPT
        Path to the Python script that generates correlations, randomized
        backgrounds, and CCD statistics.

    REFERENCE_MANIFEST
        Manifest describing the selected reference expression dataset.

    ORTHO_FILE
        Reference-to-target gene mapping file.

        This is required when the reference and target use different species
        or when explicit reference-target gene mapping is needed by NCC_SCRIPT.

    N_RANDOM
        Number of randomizations used for the unmatched background.

In a fully generic wrapper, it is also useful to define:

    REFERENCE_NAME
    REFERENCE_SPECIES
    REFERENCE_GENES


Target species
--------------
For the neuroblastoma analyses described here:

    TARGET_SPECIES=Human

The reference may be either:

    REFERENCE_SPECIES=Mouse

or:

    REFERENCE_SPECIES=Human

If both reference and target are Human, the mapping file should still be
consistent with the expectations of NCC_SCRIPT. If Mouse is used as the
reference, the orthology mapping must correctly connect the mouse clock genes
to their human counterparts.


Required target manifests
-------------------------
The workflow looks for the following manifests under manifest.d/:

    target_clock_age_12.tsv
    target_clock_alive_12.tsv
    target_clock_all_12.tsv
    target_clock_gender_12.tsv
    target_clock_mycn_12.tsv
    target_clock_risk_12.tsv
    target_clock_stage_12.tsv

A target analysis is run only when its corresponding manifest exists and is
non-empty.


Step 1: Generate the reference correlations
-------------------------------------------
The workflow first checks for the weighted reference matrix:

    $RESULTS_DIR/${REFERENCE_NAME}_12.weighted_rho.csv

If the weighted matrix does not exist or is empty, NCC_SCRIPT is run with:

    --generate_reference_correlations

The main reference outputs are:

    ${REFERENCE_NAME}_12.full_rho.tsv
        Full pairwise reference-gene correlation information.

    ${REFERENCE_NAME}_12.weighted_rho.csv
        Weighted reference correlation matrix used in downstream CCD analyses.

A temporary target-format reference file may also be written as:

    _tmp.${REFERENCE_NAME}_12.target.tsv

Reference generation should use:

    --reference_species "$REFERENCE_SPECIES"

and the 12-gene list appropriate for that species:

    --l_clock_genes_reference "$REFERENCE_GENES"

For reference-only correlation generation, the target-side species argument
should correspond to the species represented by the reference manifest, unless
NCC_SCRIPT requires a different configuration.


Conceptually:

    if weighted reference matrix is missing:
        generate reference correlations
    else:
        reuse the existing reference matrix


Step 2: Verify that the reference matrix exists
------------------------------------------------
Before any target CCD analysis is attempted, the workflow checks:

    $RESULTS_DIR/${REFERENCE_NAME}_12.weighted_rho.csv

If this file is absent or empty, all target analyses are skipped.


Step 3: Compute target correlations, unmatched background, and CCD
------------------------------------------------------------------
For each clinical target, the workflow generates three outputs:

    1. target_full_rho.tsv
       Observed target-group correlation information.

    2. background.tsv
       Unmatched random-gene background correlations.

    3. CCD.tsv
       CCD values and associated significance statistics.

Each analysis uses:

    --generate_target_correlations
    --generate_background_target_correlations
    --compute_CCD
    --background_matching_mode unmatched

The number of randomizations is controlled by:

    --number_randomizations "$N_RANDOM"


Generic target output naming
----------------------------
For any selected reference:

    target_clock_<variable>_12.${REFERENCE_NAME}.target_full_rho.tsv
    target_clock_<variable>_12.${REFERENCE_NAME}.background.tsv
    target_clock_<variable>_12.${REFERENCE_NAME}.CCD.tsv

where <variable> is one of:

    age
    alive
    all
    gender
    mycn
    risk
    stage


Clinical analyses
-----------------
The workflow evaluates the following seven target definitions:

1. Age at diagnosis

    Manifest:
        manifest.d/target_clock_age_12.tsv

    Outputs:
        target_clock_age_12.${REFERENCE_NAME}.target_full_rho.tsv
        target_clock_age_12.${REFERENCE_NAME}.background.tsv
        target_clock_age_12.${REFERENCE_NAME}.CCD.tsv


2. Outcome / alive status

    Manifest:
        manifest.d/target_clock_alive_12.tsv

    Outputs:
        target_clock_alive_12.${REFERENCE_NAME}.target_full_rho.tsv
        target_clock_alive_12.${REFERENCE_NAME}.background.tsv
        target_clock_alive_12.${REFERENCE_NAME}.CCD.tsv


3. All samples, independent of clinical grouping

    Manifest:
        manifest.d/target_clock_all_12.tsv

    Outputs:
        target_clock_all_12.${REFERENCE_NAME}.target_full_rho.tsv
        target_clock_all_12.${REFERENCE_NAME}.background.tsv
        target_clock_all_12.${REFERENCE_NAME}.CCD.tsv


4. Gender

    Manifest:
        manifest.d/target_clock_gender_12.tsv

    Outputs:
        target_clock_gender_12.${REFERENCE_NAME}.target_full_rho.tsv
        target_clock_gender_12.${REFERENCE_NAME}.background.tsv
        target_clock_gender_12.${REFERENCE_NAME}.CCD.tsv


5. MYCN amplification status

    Manifest:
        manifest.d/target_clock_mycn_12.tsv

    Outputs:
        target_clock_mycn_12.${REFERENCE_NAME}.target_full_rho.tsv
        target_clock_mycn_12.${REFERENCE_NAME}.background.tsv
        target_clock_mycn_12.${REFERENCE_NAME}.CCD.tsv


6. Harmonized risk

    Manifest:
        manifest.d/target_clock_risk_12.tsv

    Outputs:
        target_clock_risk_12.${REFERENCE_NAME}.target_full_rho.tsv
        target_clock_risk_12.${REFERENCE_NAME}.background.tsv
        target_clock_risk_12.${REFERENCE_NAME}.CCD.tsv


7. INSS stage

    Manifest:
        manifest.d/target_clock_stage_12.tsv

    Outputs:
        target_clock_stage_12.${REFERENCE_NAME}.target_full_rho.tsv
        target_clock_stage_12.${REFERENCE_NAME}.background.tsv
        target_clock_stage_12.${REFERENCE_NAME}.CCD.tsv


Reference-to-target species configuration
-----------------------------------------
The reference-specific analysis should follow this logic:

    reference species = Human or Mouse
    target species    = Human

Examples:

    Human reference -> Human target

        --reference_species Human
        --target_species Human

    Mouse reference -> Human target

        --reference_species Mouse
        --target_species Human

The orthology/reference-target mapping file must be appropriate for the
selected pair.


Gene-symbol conventions
-----------------------
Use gene symbols that match the reference expression dataset.

For Mouse references:

    Bmal1
    Clock
    Cry1
    Cry2
    Dbp
    Npas2
    Nr1d1
    Nr1d2
    Per1
    Per2
    Per3
    Tef

For Human references:

    ARNTL
    CLOCK
    CRY1
    CRY2
    DBP
    NPAS2
    NR1D1
    NR1D2
    PER1
    PER2
    PER3
    TEF

Do not assume that changing only --reference_species is sufficient. The
reference manifest and reference gene list must also correspond to that species.


Skip / resume logic
-------------------
For every target, the workflow first checks whether all three expected outputs
already exist and are non-empty:

    target_full_rho.tsv
    background.tsv
    CCD.tsv

If all three exist, the analysis is skipped.

If one or more outputs are missing, the workflow checks whether the target
manifest exists and is non-empty. If the manifest exists, the complete target
analysis is run again.

Conceptually:

    for each target:
        if target_full_rho exists
           and background exists
           and CCD exists:
               skip analysis

        else if target manifest exists:
               generate target correlations
               generate unmatched random-gene background
               compute CCD

        else:
               do nothing for that target


Background model
----------------
The clinical CCD analyses described here use:

    --background_matching_mode unmatched

Thus, the null/background is generated from unmatched random gene sets rather
than expression-matched or otherwise matched random genes.


Suggested generic Bash variables
--------------------------------
A reference-agnostic wrapper can begin with variables such as:

    REFERENCE_NAME="reference_mouse_clock"
    REFERENCE_SPECIES="Mouse"
    REFERENCE_GENES="Bmal1|Clock|Cry1|Cry2|Dbp|Npas2|Nr1d1|Nr1d2|Per1|Per2|Per3|Tef"
    TARGET_SPECIES="Human"

For a human adrenal reference:

    REFERENCE_NAME="reference_human_adrenal_gland"
    REFERENCE_SPECIES="Human"
    REFERENCE_GENES="ARNTL|CLOCK|CRY1|CRY2|DBP|NPAS2|NR1D1|NR1D2|PER1|PER2|PER3|TEF"
    TARGET_SPECIES="Human"

The remaining workflow can then use these variables instead of hard-coded
reference names, species, and gene symbols.


Workflow summary
----------------
The complete reference-agnostic logic is:

    1. Select a reference dataset.

    2. Set:
         - reference name
         - reference species
         - species-appropriate 12-gene clock list
         - reference manifest
         - reference-target mapping file

    3. Check for the weighted reference correlation matrix.

    4. If it is missing, generate the reference correlations.

    5. Confirm that the weighted reference matrix exists.

    6. For each of the seven target manifests:
         - skip if all expected outputs already exist;
         - otherwise run the analysis if the manifest exists.

    7. For each target analysis:
         - generate observed target correlations;
         - generate an unmatched randomized-gene background;
         - compute CCD and significance.


Pseudocode
----------
    REFERENCE_NAME    = selected reference
    REFERENCE_SPECIES = Human or Mouse
    REFERENCE_GENES   = 12-gene panel matching REFERENCE_SPECIES
    TARGET_SPECIES    = Human
    GENE_SET          = 12

    if ${REFERENCE_NAME}_12.weighted_rho.csv is missing:

        generate reference correlations using:
            REFERENCE_MANIFEST
            REFERENCE_SPECIES
            REFERENCE_GENES

    if reference weighted-rho matrix is still missing:
        stop target CCD processing

    else:
        for TARGET in age, alive, all, gender, mycn, risk, stage:

            if TARGET target_full_rho exists
               AND TARGET background exists
               AND TARGET CCD result exists:

                skip TARGET

            else if TARGET manifest exists:

                generate target correlations using:
                    REFERENCE_SPECIES
                    TARGET_SPECIES=Human
                    ORTHO_FILE

                generate unmatched random-gene background

                compute CCD


Notes
-----
1. The workflow is reference-agnostic only if all reference-specific inputs are
   changed consistently. Reference name, species, manifest, gene symbols, and
   mapping file must describe the same reference dataset.

2. Human and Mouse references can both be compared with Human neuroblastoma
   target data.

3. For a Human reference, use the Human clock-gene symbols expected by the
   reference matrix. For a Mouse reference, use Mouse symbols.

4. A partially completed target analysis is rerun because skipping occurs only
   when all three target output files are present and non-empty.

5. If a target manifest is absent or empty, that target is silently skipped.

6. The "all" target does not represent a clinical contrast; it evaluates the
   complete target sample set independently of the clinical grouping variables.

7. This README describes the unmatched-background workflow. A matched-background
   workflow should be documented separately or explicitly identified by its
   own --background_matching_mode setting.
