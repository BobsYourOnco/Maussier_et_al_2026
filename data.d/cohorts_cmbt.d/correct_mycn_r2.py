#!/usr/bin/env python3

"""
correct_mycn_r2.py

Apply ComBat batch correction to one or more TSV files using
'#mycn_status_unified' as the batch variable, while excluding the
'reporter' column and dropping samples with missing batch labels.

The script preserves:
- all metadata rows starting with '#'
- original row order
- original sample order except for:
    1) removal of 'reporter'
    2) removal of samples with missing batch labels

It also:
- prints metadata variable names
- prints sample counts before and after filtering
- prints batch counts
- fails if any batch has fewer than 2 samples
- fails if ComBat produces NaN or Inf values

Typical use:
    python correct_mycn_r2.py "*R2*.tsv"

Single-file run with explicit output file:
    python correct_mycn_r2.py input.tsv --output corrected.tsv
"""

from pathlib import Path
import argparse
import glob
import warnings

import numpy as np
import pandas as pd
import scanpy as sc


def correct_file(
    file_path: Path,
    batch_row: str = "#mycn_status_unified",
    output_path: str | None = None,
) -> Path:
    """
    Correct a single TSV expression file using ComBat.

    1) What does the method do
    --------------------------
    This method reads one tab-separated gene expression file, separates
    metadata rows (those beginning with '#') from gene expression rows,
    excludes the column named 'reporter', and applies Scanpy/ComBat batch
    correction to the gene expression matrix using the metadata row
    specified by `batch_row` as the batch variable.

    Samples whose batch value is missing after cleaning are excluded
    before correction.

    The method also prints:
    - metadata variable names found in the file
    - number of samples before filtering
    - number of samples excluded due to missing batch labels
    - number of samples retained for correction
    - batch counts per group

    It stops with an error if:
    - all samples are removed
    - fewer than 2 batch groups remain
    - any batch has fewer than 2 samples
    - ComBat produces NaN or Inf values

    The corrected file is reconstructed so that:
    - metadata rows are kept unchanged
    - corrected gene rows are placed below the metadata rows
    - row order is preserved
    - the 'reporter' column is excluded from the output
    - samples with missing batch labels are excluded from the output

    2) Input
    --------
    file_path : Path
        Path to the input TSV file.

    batch_row : str, default="#mycn_status_unified"
        Metadata row to use as the ComBat batch variable.

    output_path : str | None, default=None
        Output filename for the corrected TSV.
        If None, a filename is generated automatically.

    3) Output
    ---------
    Path
        Path to the generated corrected TSV file.
    """
    data = pd.read_csv(file_path, sep="\t", index_col=0)

    # Keep all sample columns except 'reporter'
    sample_cols = [c for c in data.columns if c != "reporter"]

    if len(sample_cols) == 0:
        raise ValueError(
            f"{file_path}: no usable sample columns remain after excluding 'reporter'."
        )

    # Split rows into metadata/covariates and genes
    idx_as_str = data.index.to_series().astype(str)
    covariate_rows = list(data.index[idx_as_str.str.startswith("#")])
    gene_rows = list(data.index[~idx_as_str.str.startswith("#")])

    if len(covariate_rows) == 0:
        raise ValueError(f"{file_path}: no metadata rows starting with '#' were found.")

    if len(gene_rows) == 0:
        raise ValueError(f"{file_path}: no gene rows were found.")

    # Metadata table: rows=samples, columns=metadata fields
    all_covariates = data.loc[covariate_rows, sample_cols].T.copy()

    print(f"\n[INFO] Processing file: {file_path}")
    print(f"[INFO] Metadata variables ({len(all_covariates.columns)}):")
    print("       " + ", ".join(map(str, all_covariates.columns)))
    print(f"[INFO] Samples before filtering: {len(sample_cols)}")

    if batch_row not in all_covariates.columns:
        raise ValueError(f"{file_path}: '{batch_row}' was not found in metadata rows.")

    # Clean batch labels
    batch_variable = (
        all_covariates[batch_row]
        .astype(str)
        .str.strip()
        .replace({"": pd.NA, "NA": pd.NA, "nan": pd.NA, "NaN": pd.NA})
    )

    # Drop samples with missing batch labels
    valid_samples = batch_variable[~batch_variable.isna()].index.tolist()
    dropped_samples = batch_variable[batch_variable.isna()].index.tolist()

    if len(valid_samples) == 0:
        raise ValueError(
            f"{file_path}: all samples were removed because {batch_row} is missing."
        )

    # Restrict everything to valid samples only
    sample_cols = [c for c in sample_cols if c in valid_samples]
    all_covariates = all_covariates.loc[valid_samples].copy()
    batch_variable = batch_variable.loc[valid_samples].copy()

    print(f"[INFO] Batch variable used: {batch_row}")
    print(f"[INFO] Samples excluded for missing batch label: {len(dropped_samples)}")
    if dropped_samples:
        print("       " + ", ".join(dropped_samples))
    print(f"[INFO] Samples retained for correction: {len(sample_cols)}")

    # Require at least 2 groups
    n_groups = batch_variable.nunique()
    if n_groups < 2:
        raise ValueError(
            f"{file_path}: {batch_row} has fewer than 2 groups after filtering; "
            f"ComBat cannot be applied."
        )

    print(f"[INFO] Number of batch groups: {n_groups}")
    print(f"[INFO] Batch levels: {', '.join(map(str, sorted(batch_variable.unique())))}")

    # Print and validate batch counts
    batch_counts = batch_variable.value_counts().sort_index()
    print("[INFO] Batch counts:")
    print(batch_counts.to_string())

    small_batches = batch_counts[batch_counts < 2]
    if not small_batches.empty:
        raise ValueError(
            f"{file_path}: ComBat is unsafe because at least one batch has fewer than "
            f"2 samples: "
            + ", ".join(f"{k}={v}" for k, v in small_batches.items())
        )

    # Gene expression matrix: rows=genes, columns=samples
    gene_data = data.loc[gene_rows, sample_cols].apply(pd.to_numeric, errors="raise")

    # AnnData expects rows=samples, columns=genes
    adata = sc.AnnData(gene_data.T.copy())

    # Align metadata to sample order
    adata.obs = all_covariates.reindex(adata.obs_names).copy()
    adata.obs["batch"] = batch_variable.reindex(adata.obs_names).astype("category")

    # Apply ComBat batch correction
    # Catch runtime warnings so they are visible in output without being silent.
    with warnings.catch_warnings(record=True) as caught_warnings:
        warnings.simplefilter("always")
        sc.pp.combat(adata, key="batch")

    for w in caught_warnings:
        print(f"[WARNING] {type(w.message).__name__}: {w.message}")

    # Explicit post-ComBat sanity check
    nan_count = int(np.isnan(adata.X).sum())
    inf_count = int(np.isinf(adata.X).sum())

    print(f"[INFO] NaN count after ComBat: {nan_count}")
    print(f"[INFO] Inf count after ComBat: {inf_count}")

    if nan_count > 0 or inf_count > 0:
        raise ValueError(
            f"{file_path}: ComBat produced invalid values "
            f"(NaN={nan_count}, Inf={inf_count})."
        )

    # Convert corrected matrix back to rows=genes, columns=samples
    corrected_gene_data = pd.DataFrame(
        adata.X.T,
        index=gene_rows,
        columns=sample_cols,
    )

    # Rebuild final output: unchanged metadata + corrected genes
    corrected_data = pd.concat(
        [
            data.loc[covariate_rows, sample_cols],
            corrected_gene_data.loc[gene_rows, sample_cols],
        ],
        axis=0,
    )

    # Preserve exact row order
    corrected_data = corrected_data.loc[covariate_rows + gene_rows, sample_cols]

    if output_path is None:
        out_path = file_path.with_name(file_path.stem + "_MYCNunified_corrected.tsv")
    else:
        out_path = Path(output_path)

    corrected_data.to_csv(out_path, sep="\t")
    print(f"[INFO] Output written to: {out_path}")

    return out_path


def expand_inputs(inputs: list[str]) -> list[Path]:
    """
    Expand explicit filenames and glob patterns into a unique list of file paths.

    1) What does the method do
    --------------------------
    This method resolves command-line inputs into a de-duplicated list of files.

    2) Input
    --------
    inputs : list[str]
        File paths and/or glob patterns.

    3) Output
    ---------
    list[Path]
        A de-duplicated list of resolved file paths.
    """
    paths = []

    for item in inputs:
        matches = glob.glob(item)
        if matches:
            paths.extend(matches)
        else:
            paths.append(item)

    seen = set()
    unique_paths = []

    for p in paths:
        rp = str(Path(p).resolve())
        if rp not in seen:
            seen.add(rp)
            unique_paths.append(Path(p))

    return unique_paths


def main() -> None:
    """
    Parse command-line arguments and run batch correction on all input files.

    1) What does the method do
    --------------------------
    Reads command-line arguments, expands input files, runs correction,
    and prints success/error messages.

    2) Input
    --------
    Provided through the command line:
    - input files or glob patterns
    - optional --batch-row
    - optional --output for single-file runs

    3) Output
    ---------
    None
        Writes corrected files to disk and prints status to stdout.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Apply ComBat correction to TSV files using "
            "#mycn_status_unified as batch, excluding 'reporter', "
            "dropping samples with missing batch labels, printing batch counts, "
            "and checking for NaN/Inf after ComBat."
        )
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="Input TSV files or glob patterns, e.g. '*R2*.tsv'",
    )
    parser.add_argument(
        "--batch-row",
        default="#mycn_status_unified",
        help=(
            "Metadata row to use as ComBat batch variable "
            "(default: #mycn_status_unified)"
        ),
    )
    parser.add_argument(
        "--output",
        default=None,
        help=(
            "Output filename for the corrected TSV. "
            "Use only when processing a single input file."
        ),
    )
    args = parser.parse_args()

    files = expand_inputs(args.inputs)
    if not files:
        raise SystemExit("No input files matched.")

    if args.output is not None and len(files) != 1:
        raise SystemExit(
            "--output can only be used when exactly one input file is processed."
        )

    for f in files:
        try:
            out = correct_file(
                Path(f),
                batch_row=args.batch_row,
                output_path=args.output,
            )
            print(f"[OK] {f} -> {out}")
        except Exception as e:
            print(f"[ERROR] {f}: {e}")


if __name__ == "__main__":
    main()
