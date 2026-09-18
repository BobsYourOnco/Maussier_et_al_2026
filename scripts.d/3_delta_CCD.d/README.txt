Purpose
-------
This README describes the Bash workflow used to compute direct Delta-CCD
contrasts and matched-background CCD significance for a 12-gene circadian
clock panel.

The workflow is reference-agnostic: the reference dataset may be Human or
Mouse, provided that the reference name, reference species, 12-gene clock list,
reference manifest, and reference-to-target gene mapping are configured
consistently.

The target neuroblastoma data are Human.

The workflow has two distinct analysis layers:

    1. Direct Delta-CCD contrasts
       Directly compare CCD between two clinical groups.

    2. Matched-background CCD analyses
       Test each clinical group against matched random-gene backgrounds.

These two layers are complementary and should not be interpreted as the same
statistical test.


Clinical contrasts
------------------
Direct Delta-CCD is computed for six clinical variables:

    1. MYCN amplification status
       amplified vs not_amplified

    2. Harmonized risk
       high_risk vs low_risk

    3. INSS stage
       st4 vs st1

    4. Gender / sex
       male vs female

    5. Outcome status
       dead vs alive

    6. Age at diagnosis
       ge18m vs lt18m


Reference-specific settings
---------------------------
For each reference, define or substitute:

    REFERENCE_NAME
        Short reference identifier used in output filenames.

        Examples:

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
        by the reference dataset.

        Mouse example:

            Bmal1|Clock|Cry1|Cry2|Dbp|Npas2|Nr1d1|Nr1d2|Per1|Per2|Per3|Tef

        Human example:

            ARNTL|CLOCK|CRY1|CRY2|DBP|NPAS2|NR1D1|NR1D2|PER1|PER2|PER3|TEF


Required variables
------------------
The Bash environment must define:

    RESULTS_DIR
        Directory where Delta-CCD and matched-background outputs are written.

    SCRIPT
        Path to the Python NCC/CCD analysis script.

    REFERENCE_MANIFEST
        Manifest describing the selected reference expression dataset.

    ORTHO_FILE
        Reference-to-target gene mapping file.

It is also recommended to define:

    REFERENCE_NAME
    REFERENCE_SPECIES
    REFERENCE_GENES
    TARGET_SPECIES

For these neuroblastoma analyses:

    TARGET_SPECIES=Human


Required manifests
------------------
Direct Delta-CCD uses:

    manifests.d/target_clock_all_12.tsv

This manifest must contain the clinical variables used for the six contrasts.

Matched-background CCD analyses use separate group-specific manifests:

    manifests.d/target_clock_age_12.tsv
    manifests.d/target_clock_alive_12.tsv
    manifests.d/target_clock_gender_12.tsv
    manifests.d/target_clock_mycn_12.tsv
    manifests.d/target_clock_risk_12.tsv
    manifests.d/target_clock_stage_12.tsv


Reference matrix requirement
----------------------------
Before any Delta-CCD or matched-background analysis is run, the workflow checks
for the weighted reference correlation matrix:

    $RESULTS_DIR/${REFERENCE_NAME}_12.weighted_rho.csv

If this file is absent or empty, the analyses are skipped.

Conceptually:

    if weighted reference matrix is missing:
        skip Delta-CCD and matched-background analyses
    else:
        continue


Part 1: Direct Delta-CCD analysis
---------------------------------
Direct Delta-CCD compares the CCD values of two clinical groups.

The general quantity is:

    Delta-CCD = CCD(Group A) - CCD(Group B)

Therefore:

    Delta-CCD > 0
        Group A has a larger CCD than Group B.

    Delta-CCD < 0
        Group A has a smaller CCD than Group B.

The biological interpretation of the sign depends on how CCD is defined in the
analysis pipeline and on the chosen Group A / Group B orientation.


Direct Delta-CCD settings
-------------------------
Each direct contrast uses:

    --compute_delta_CCD
    --delta_n_permutations 5000
    --delta_n_bootstraps 2000
    --delta_seed 42
    --delta_min_samples_per_group 10

Thus, the workflow uses:

    5000 permutations
        for permutation-based significance testing.

    2000 bootstrap resamples
        for uncertainty estimation / confidence intervals.

    seed = 42
        for reproducibility.

    minimum 10 samples per group
        before a contrast is evaluated.


Direct Delta-CCD contrasts
--------------------------

1. MYCN amplification
~~~~~~~~~~~~~~~~~~~~~
Contrast variable:

    mycn_status_unified

Groups:

    Group A = amplified
    Group B = not_amplified

Generic outputs:

    target_clock_all_12.${REFERENCE_NAME}.MYCN.delta.tsv
    target_clock_all_12.${REFERENCE_NAME}.MYCN.delta_draws.tsv


2. Harmonized risk
~~~~~~~~~~~~~~~~~~
Contrast variable:

    risk_group_custom

Groups:

    Group A = high_risk
    Group B = low_risk

Generic outputs:

    target_clock_all_12.${REFERENCE_NAME}.risk.delta.tsv
    target_clock_all_12.${REFERENCE_NAME}.risk.delta_draws.tsv


3. INSS stage
~~~~~~~~~~~~~
Contrast variable:

    inss_stage_5cat

Groups:

    Group A = st4
    Group B = st1

Generic outputs:

    target_clock_all_12.${REFERENCE_NAME}.stage.delta.tsv
    target_clock_all_12.${REFERENCE_NAME}.stage.delta_draws.tsv


4. Gender / sex
~~~~~~~~~~~~~~~
Contrast variable:

    sex_unified

Groups:

    Group A = male
    Group B = female

Generic outputs:

    target_clock_all_12.${REFERENCE_NAME}.gender.delta.tsv
    target_clock_all_12.${REFERENCE_NAME}.gender.delta_draws.tsv


5. Outcome status
~~~~~~~~~~~~~~~~~
Contrast variable:

    alive_status

Groups:

    Group A = dead
    Group B = alive

Generic outputs:

    target_clock_all_12.${REFERENCE_NAME}.alive.delta.tsv
    target_clock_all_12.${REFERENCE_NAME}.alive.delta_draws.tsv


6. Age at diagnosis
~~~~~~~~~~~~~~~~~~~
Contrast variable:

    age_group_18m

Groups:

    Group A = ge18m
    Group B = lt18m

Generic outputs:

    target_clock_all_12.${REFERENCE_NAME}.age.delta.tsv
    target_clock_all_12.${REFERENCE_NAME}.age.delta_draws.tsv


Direct Delta-CCD skip / resume logic
------------------------------------
Each contrast is skipped only when both expected files exist and are non-empty:

    *.delta.tsv
    *.delta_draws.tsv

Conceptually:

    if delta result exists
       AND delta draws exist:

        skip contrast

    else if target_clock_all_12.tsv exists:

        compute Delta-CCD

    else:

        do not run the contrast


Part 2: Matched-background CCD analysis
---------------------------------------
The second analysis layer evaluates CCD significance relative to matched
random-gene backgrounds.

This is not another direct Delta-CCD contrast.

Instead, for each clinical grouping manifest, the workflow:

    1. generates observed target correlations;

    2. generates matched random-gene background correlations;

    3. computes CCD and background-based significance;

    4. records matched-gene diagnostics and the generated matched gene sets.


Matched-background settings
---------------------------
Each matched-background analysis uses:

    --number_randomizations 1000
    --background_matching_mode matched
    --matching_mean_caliper 0.05
    --matching_variance_caliper 0.05
    --matching_expand_step 0.05
    --matching_max_caliper 1.0
    --matching_seed 42

The matched random-gene sets are therefore selected using similarity in gene
expression mean and variance.

Initial matching tolerances:

    mean caliper     = 0.05
    variance caliper = 0.05

If adequate matches are not available, the matching tolerance expands by:

    0.05

up to a maximum caliper of:

    1.0

The matching procedure uses:

    seed = 42

for reproducibility.


Matched-background outputs
--------------------------
For each clinical target, five files are expected:

    1. matched.target_full_rho.tsv
       Observed target correlation information.

    2. matched.background.tsv
       Correlations from matched random-gene sets.

    3. matched.CCD.tsv
       CCD values and matched-background significance.

    4. matched.gene_summary.tsv
       Summary of the genes used in the matching procedure.

    5. matched.gene_sets.tsv
       The actual matched random-gene sets generated during randomization.


Generic output pattern:

    target_clock_<variable>_12.${REFERENCE_NAME}.matched.target_full_rho.tsv
    target_clock_<variable>_12.${REFERENCE_NAME}.matched.background.tsv
    target_clock_<variable>_12.${REFERENCE_NAME}.matched.CCD.tsv
    target_clock_<variable>_12.${REFERENCE_NAME}.matched.gene_summary.tsv
    target_clock_<variable>_12.${REFERENCE_NAME}.matched.gene_sets.tsv


Matched-background clinical targets
-----------------------------------
Matched-background CCD is computed for:

    age
    alive
    gender
    mycn
    risk
    stage

Unlike the direct Delta-CCD analysis, there is no separate "all samples" matched
analysis in this Bash block.


Matched-background skip / resume logic
--------------------------------------
A matched-background analysis is skipped only when all five expected outputs
exist and are non-empty.

Conceptually:

    if target_full_rho exists
       AND matched background exists
       AND CCD result exists
       AND gene summary exists
       AND matched gene sets exist:

        skip matched-background analysis

    else if the target manifest exists:

        generate target correlations
        generate matched random-gene background
        compute CCD
        write matching diagnostics

    else:

        do not run that target


Difference between direct Delta-CCD and matched-background CCD
--------------------------------------------------------------
The two parts of the workflow answer different questions.

Direct Delta-CCD asks:

    Is CCD different between Group A and Group B?

For example:

    Is CCD different between high-risk and low-risk tumors?

Matched-background CCD asks:

    Is the observed CCD for a clinical group unusual relative to matched random
    gene sets?

The matched-background analysis therefore provides supporting evidence about
whether observed clock-gene coordination differs from what would be expected
for genes with similar expression properties.

It does not replace the direct Group A versus Group B Delta-CCD test.


Species configuration
---------------------
The reference can be Human or Mouse.

For a Mouse reference and Human target:

    --reference_species Mouse
    --target_species Human

For a Human reference and Human target:

    --reference_species Human
    --target_species Human

The clock-gene list must use symbols appropriate for the reference species.

The ORTHO_FILE must provide a valid reference-to-target mapping compatible with
the selected species pair.


Suggested generic Bash variables
--------------------------------
Mouse reference example:

    REFERENCE_NAME="reference_mouse_clock"
    REFERENCE_SPECIES="Mouse"
    REFERENCE_GENES="Bmal1|Clock|Cry1|Cry2|Dbp|Npas2|Nr1d1|Nr1d2|Per1|Per2|Per3|Tef"
    TARGET_SPECIES="Human"

Human reference example:

    REFERENCE_NAME="reference_human_adrenal_gland"
    REFERENCE_SPECIES="Human"
    REFERENCE_GENES="ARNTL|CLOCK|CRY1|CRY2|DBP|NPAS2|NR1D1|NR1D2|PER1|PER2|PER3|TEF"
    TARGET_SPECIES="Human"


Recommended naming convention
-----------------------------
Reference matrix:

    ${REFERENCE_NAME}_12.weighted_rho.csv

Direct Delta-CCD:

    target_clock_all_12.${REFERENCE_NAME}.<contrast>.delta.tsv
    target_clock_all_12.${REFERENCE_NAME}.<contrast>.delta_draws.tsv

Matched-background CCD:

    target_clock_<variable>_12.${REFERENCE_NAME}.matched.target_full_rho.tsv
    target_clock_<variable>_12.${REFERENCE_NAME}.matched.background.tsv
    target_clock_<variable>_12.${REFERENCE_NAME}.matched.CCD.tsv
    target_clock_<variable>_12.${REFERENCE_NAME}.matched.gene_summary.tsv
    target_clock_<variable>_12.${REFERENCE_NAME}.matched.gene_sets.tsv


Workflow summary
----------------
The complete workflow is:

    1. Select a Human or Mouse reference.

    2. Define:
         - reference name;
         - reference species;
         - 12-gene clock list;
         - reference manifest;
         - reference-target gene mapping.

    3. Confirm that the weighted reference matrix exists.

    4. Run six direct Delta-CCD contrasts:
         - MYCN;
         - risk;
         - stage;
         - gender;
         - alive/dead;
         - age.

    5. For each direct contrast:
         - compare Group A against Group B;
         - use 5000 permutations;
         - use 2000 bootstraps;
         - require at least 10 samples per group;
         - write Delta-CCD results and resampling draws.

    6. Run matched-background CCD analyses for:
         - age;
         - alive;
         - gender;
         - MYCN;
         - risk;
         - stage.

    7. For each matched analysis:
         - generate target correlations;
         - generate 1000 matched random-gene sets;
         - compute CCD significance;
         - save matched-gene diagnostics.


Pseudocode
----------
    REFERENCE_NAME    = selected reference
    REFERENCE_SPECIES = Human or Mouse
    REFERENCE_GENES   = species-appropriate 12-gene clock panel
    TARGET_SPECIES    = Human

    if reference weighted-rho matrix is missing:

        skip all analyses

    else:

        # Direct Delta-CCD
        for CONTRAST in MYCN, risk, stage, gender, alive, age:

            if delta result exists
               AND delta draws exist:

                skip CONTRAST

            else if target_clock_all_12 manifest exists:

                select contrast variable
                select Group A
                select Group B

                compute Delta-CCD

                run 5000 permutations
                run 2000 bootstraps
                require >= 10 samples/group

                save delta result
                save delta draws


        # Matched-background support
        for TARGET in age, alive, gender, mycn, risk, stage:

            if all five matched output files exist:

                skip TARGET

            else if TARGET manifest exists:

                generate observed target correlations

                generate 1000 matched random-gene backgrounds

                match genes on expression mean and variance

                compute CCD and matched-background significance

                save matched-gene summary
                save matched gene sets


Troubleshooting
---------------

Reference matrix is missing
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Check:

    ls -lh "$RESULTS_DIR/${REFERENCE_NAME}_12.weighted_rho.csv"

The -s Bash test requires the file to exist and have size greater than zero.


Direct Delta-CCD does not run
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Check that:

    manifests.d/target_clock_all_12.tsv

exists and is non-empty.

Also confirm that the required contrast variable and group labels are present in
the manifest.


A matched-background analysis does not run
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Check that the corresponding manifest exists, for example:

    manifests.d/target_clock_risk_12.tsv

or:

    manifests.d/target_clock_stage_12.tsv


Analysis reruns despite previous output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
For direct Delta-CCD, both files must exist:

    delta.tsv
    delta_draws.tsv

For matched-background analysis, all five files must exist:

    matched.target_full_rho.tsv
    matched.background.tsv
    matched.CCD.tsv
    matched.gene_summary.tsv
    matched.gene_sets.tsv

If any required file is absent or empty, that analysis is rerun.


Matched gene sets cannot be generated
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Check whether the expression-matching criteria are too restrictive and inspect:

    matched.gene_summary.tsv
    matched.gene_sets.tsv

The workflow begins with mean and variance calipers of 0.05 and progressively
expands them by 0.05 up to 1.0.


Gene symbols are not found
~~~~~~~~~~~~~~~~~~~~~~~~~~
Confirm that REFERENCE_GENES uses symbols appropriate for REFERENCE_SPECIES and
that all genes are present in the reference dataset.


Reference-target mapping fails
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Check:

    ls -lh "$ORTHO_FILE"

and confirm that the file contains valid mappings for the selected Human-Human
or Mouse-Human reference/target configuration.


Important interpretation note
-----------------------------
Direct Delta-CCD significance and matched-background significance should be
reported separately.

The direct Delta-CCD test provides evidence for a difference between two
clinical groups.

The matched-background analysis evaluates whether each group's observed
clock-gene coordination is unusual relative to matched random genes.

Concordance between these layers can strengthen interpretation, but a
matched-background result should not be described as the p-value for the
direct Delta-CCD contrast.

