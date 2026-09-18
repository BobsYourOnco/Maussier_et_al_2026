Purpose
-------
This README describes the Bash scripts used to generate reference correlation
files for NCC/CCD analysis.

Each Bash block generates correlation outputs for a specific reference dataset
and clock-gene set.

For example, a block configured for:

    reference_human_adrenal_gland

with the 12-gene circadian clock panel generates the reference correlations used
in downstream CCD and Delta-CCD analyses.


Main outputs
------------
For each reference and gene set, two main files are generated:

    <reference_name>_<gene_set>.full_rho.tsv

    <reference_name>_<gene_set>.weighted_rho.csv

For example:

    reference_human_adrenal_gland_12.full_rho.tsv

    reference_human_adrenal_gland_12.weighted_rho.csv

The full_rho file contains the reference correlation information, while the
weighted_rho matrix is used as the reference correlation matrix in downstream
CCD and Delta-CCD calculations.


Required variables
------------------
Before running the Bash script, define the following variables:

    NCC_SCRIPT
        Path to the Python script used to compute reference and target
        correlations.

    REFERENCE_MANIFEST
        Manifest describing the reference expression dataset.

    RESULTS_DIR
        Directory where the generated reference correlation files are written.

    ORTHO_FILE
        Ortholog or reference-target gene-mapping file used to match genes
        between reference and target datasets.


Example
-------
Example variable definitions:

    NCC_SCRIPT="/path/to/ncc_script.py"

    REFERENCE_MANIFEST="/path/to/reference_manifest.tsv"

    RESULTS_DIR="/path/to/results"

    ORTHO_FILE="/path/to/ortholog_file.tsv"


Reference configuration
-----------------------
Each reference-generation block should specify:

    1. Reference name

    2. Gene-set size

    3. Reference gene list

    4. Reference species

    5. Reference manifest

The gene symbols supplied through:

    --l_clock_genes_reference

must match the species and naming convention used in the reference expression
dataset.

For example, a Human reference may use:

    ARNTL|CLOCK|CRY1|CRY2|DBP|NPAS2|NR1D1|NR1D2|PER1|PER2|PER3|TEF

A Mouse reference may use:

    Bmal1|Clock|Cry1|Cry2|Dbp|Npas2|Nr1d1|Nr1d2|Per1|Per2|Per3|Tef


Reference-correlation workflow
------------------------------
The Bash block first checks whether the expected output files already exist and
are non-empty.

Typical logic:

    if full_rho.tsv exists
       AND weighted_rho.csv exists:

        skip reference-correlation generation

    else:

        run NCC_SCRIPT
        generate reference correlations

The Bash test:

    -s

checks that a file exists and has a size greater than zero.


Typical command structure
-------------------------
A reference-generation call generally follows this structure:

    python "$NCC_SCRIPT" \
      --input_reference_db_file "$REFERENCE_MANIFEST" \
      --reference_full_rho_info "$RESULTS_DIR/<reference_name>_<gene_set>.full_rho.tsv" \
      --reference_weighted_rho_matrix "$RESULTS_DIR/<reference_name>_<gene_set>.weighted_rho.csv" \
      --l_clock_genes_reference "<reference_gene_list>" \
      --reference_species "<Human_or_Mouse>" \
      --input_target_db_file "$REFERENCE_MANIFEST" \
      --target_full_rho_info "$RESULTS_DIR/_tmp.<reference_name>_<gene_set>.target.tsv" \
      --reference_target_gene_file "$ORTHO_FILE" \
      --target_species "<Human_or_Mouse>" \
      --generate_reference_correlations


Recommended naming convention
-----------------------------
Use the pattern:

    <reference_name>_<gene_set>.full_rho.tsv

    <reference_name>_<gene_set>.weighted_rho.csv

For example:

    reference_human_adrenal_gland_12.full_rho.tsv

    reference_human_adrenal_gland_12.weighted_rho.csv

This naming scheme makes it easy to identify the reference dataset and gene set
used in downstream CCD and Delta-CCD analyses.


Pseudocode
----------
    REFERENCE = selected reference dataset

    GENE_SET = selected clock-gene set

    if reference full-rho file exists
       AND reference weighted-rho file exists:

        skip generation

    else:

        run NCC_SCRIPT

        read REFERENCE_MANIFEST

        use the species-appropriate reference gene list

        generate full reference correlations

        generate weighted reference correlation matrix


Troubleshooting
---------------

The script reruns even though outputs exist
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Check whether one of the expected files is empty:

    ls -lh \
      "$RESULTS_DIR/reference_human_adrenal_gland_12.full_rho.tsv" \
      "$RESULTS_DIR/reference_human_adrenal_gland_12.weighted_rho.csv"

The Bash -s test requires each file to exist and have a size greater than zero.


Python cannot find the NCC script
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Check that NCC_SCRIPT points to the correct file:

    echo "$NCC_SCRIPT"

    ls -lh "$NCC_SCRIPT"


Gene symbols are not found
~~~~~~~~~~~~~~~~~~~~~~~~~~
Check that the gene symbols supplied through:

    --l_clock_genes_reference

match both the species and the expression manifest.

Human references should use Human gene symbols.

Mouse references should use Mouse gene symbols unless the analysis script
explicitly expects a different naming convention.


Ortholog mapping fails
~~~~~~~~~~~~~~~~~~~~~~
Check that ORTHO_FILE exists:

    ls -lh "$ORTHO_FILE"

Also confirm that it contains mappings for all genes in the reference clock-gene
list.


Reference outputs are missing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Check:

    1. REFERENCE_MANIFEST exists and is non-empty.

    2. RESULTS_DIR exists and is writable.

    3. NCC_SCRIPT completed without an error.

    4. The reference gene symbols are present in the reference expression data.

    5. The reference species is specified correctly.


Workflow summary
----------------
The overall procedure is:

    1. Select a reference dataset and gene set.

    2. Define the reference manifest, output directory, NCC script, and mapping
       file.

    3. Use a species-appropriate clock-gene list.

    4. Check whether the expected reference outputs already exist.

    5. If both outputs are present and non-empty, skip generation.

    6. Otherwise, generate:

         - full reference correlation information

         - weighted reference correlation matrix

    7. Use these files as reference inputs for downstream CCD and Delta-CCD
       analyses.
