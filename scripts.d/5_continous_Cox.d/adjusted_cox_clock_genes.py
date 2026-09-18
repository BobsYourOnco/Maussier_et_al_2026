#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
adjusted_cox_clock_genes.py

Robust adjusted Cox proportional-hazards workflow for the neuroblastoma
circadian-clock manuscript.

The workflow uses explicit and consistent survival-variable names (`time` and
`event`) throughout model fitting and proportional-hazards diagnostics. It also
writes a survival-variable resolution table for every condition/cohort so that
the source survival variables, event coding, sample counts, and resolution
status are documented transparently.

Workflow overview
-----------------
1. Read log2 and/or MYCN-ComBat-corrected expression data.
2. Resolve cohort-specific overall-survival time and event variables.
3. Harmonize age, MYCN, INSS stage, risk, and sex covariates.
4. Standardize each clock-gene expression value within cohort.
5. Fit univariable and prespecified adjusted Cox models.
6. Test the proportional-hazards assumption for every fitted model.
7. Apply Benjamini-Hochberg correction within prespecified testing families.
8. Combine SEQC and Kocak estimates using fixed-effect inverse-variance meta-analysis.
9. Write detailed, summary, PH-diagnostic, meta-analysis, resolution, warning,
   footnote, and interpretation outputs.

Core outputs
------------
Given --out_prefix adjusted_cox, writes:

  adjusted_cox_clock_genes_detail.tsv/.xlsx
  adjusted_cox_clock_genes_summary.tsv/.xlsx
  adjusted_cox_PH_diagnostics.tsv/.xlsx
  adjusted_cox_meta_analysis.tsv/.xlsx
  adjusted_cox_survival_variable_resolution.tsv
  adjusted_cox_warnings.tsv
  adjusted_cox_footnote.txt
  adjusted_cox_interpretation.md

Input options
-------------
Two input modes are supported.

Mode 1: Direct expression matrices
----------------------------------

python adjusted_cox_clock_genes.py \
  --expression log2:SEQC:/path/SEQC_log2_R2.tsv \
  --expression log2:Kocak:/path/Kocak_log2_R2.tsv \
  --expression mycn_corrected:SEQC:/path/SEQC_cmbt_R2.tsv \
  --expression mycn_corrected:Kocak:/path/Kocak_cmbt_R2.tsv \
  --clinical /path/clinical_metadata.tsv \
  --out_prefix adjusted_cox

Mode 2: Existing NCC/R2 target manifest
---------------------------------------

python adjusted_cox_clock_genes.py \
  --manifest log2:/path/to/target_clock_all_12.tsv \
  --manifest mycn_corrected:/path/to/target_clock_all_12_cmbt.tsv \
  --out_prefix adjusted_cox

Manifest mode expects the manifest used by the NCC workflow, with fields such as
`full_path_input_expression_file`, `separator_of_columns_in_input_expression_file`
and optional row/column exclusion/filter fields. Metadata rows beginning with "#"
inside each R2 file are used as clinical metadata.

Dependencies
------------
Required:
  pandas, numpy, lifelines, openpyxl

Install:
  pip install pandas numpy lifelines openpyxl

Notes
-----
- Primary expression condition is log2.
- mycn_corrected/ComBat is retained as technical sensitivity only.
- asIs/untransformed conditions are ignored by default.
- HRs are per 1 SD increase in expression within cohort.
- Event is always coerced to numeric 0/1 before fitting.
- Survival time must be numeric and positive.
- Adjusted models are not fitted when event counts are insufficient.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
import traceback
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


CLOCK_GENES = [
    "ARNTL", "CLOCK", "CRY1", "CRY2", "DBP", "NPAS2",
    "NR1D1", "NR1D2", "PER1", "PER2", "PER3", "TEF",
]

PRIMARY_COHORTS = {"SEQC", "Kocak"}
ALLOWED_CONDITIONS = {"log2", "mycn_corrected"}

CONDITION_ALIASES = {
    "log2": "log2",
    "log2_normalized": "log2",
    "normalized": "log2",
    "cmbt": "mycn_corrected",
    "combat": "mycn_corrected",
    "mycn_corrected": "mycn_corrected",
    "mycn_combat": "mycn_corrected",
    "mycn-corrected": "mycn_corrected",
    "as_is": "asIs",
    "asis": "asIs",
    "untransformed": "asIs",
    "unnormalized": "asIs",
}

EVENT_TRUE = {"1", "true", "yes", "dead", "deceased", "event", "died", "death", "relapse", "progressed"}
EVENT_FALSE = {"0", "false", "no", "alive", "living", "censored", "no_event", "nonevent", "noevent"}

# Common aliases, including the aliases used by the prior continuous Cox script.
SURVIVAL_ALIASES = {
    "SEQC": {
        "time": ["os_day", "os_days", "os_time", "overall_survival_time", "overall_survival_days"],
        "event": ["os_bin", "os_event", "overall_survival_event", "event", "status"],
    },
    "Kocak": {
        "time": ["overallsurvival", "overall_survival", "os_time", "survival_time"],
        "event": ["osevent", "os_event", "overall_survival_event", "event", "status"],
    },
    "default": {
        "time": [
            "os_time", "overall_survival_time", "overall_survival", "survival_time", "time",
            "os_day", "os_days", "days_to_death_or_last_followup", "follow_up_time",
            "overall_survival_time_days", "OS.time"
        ],
        "event": [
            "os_event", "overall_survival_event", "survival_event", "event", "status",
            "vital_status", "death_event", "dead", "os_bin", "osevent", "OS.event"
        ],
    },
}

CLINICAL_ALIASES = {
    "sample_id": [
        "sample_id", "sample", "sampleid", "sample_ID", "id", "r2_id", "r2_sample",
        "geo_accession", "gsm", "patient_id", "patient", "array"
    ],
    "cohort": ["cohort", "database", "database_name", "dataset", "series"],
    "mycn": [
        "mycn_status_unified", "mycn_status", "MYCN_status", "MYCN", "mycn",
        "mycn_amplification", "mycn_amp", "MYCN_amplification", "mycnamplified"
    ],
    "stage": [
        "inss_stage_5cat", "inss_stage", "INSS", "inss", "stage",
        "tumor_stage", "diagnosis_stage"
    ],
    "risk": [
        "risk_group_custom", "risk_group_harmonised", "risk_group_harmonized",
        "risk_group", "risk", "high_risk", "highrisk", "inrg_risk"
    ],
    "age18": [
        "age_group_18m", "age_ge18m", "age_18m", "age_at_diagnosis_18m",
        "age_over_18_months", "age_binary", "age_group"
    ],
    "age_months": [
        "age_months", "age_at_diagnosis_months", "age_at_diagnosis_in_months",
        "diagnosis_age_months", "age_at_diagnosis", "age"
    ],
    "age_days": ["age_days", "age_at_diagnosis_days", "diagnosis_age_days"],
    "age_years": ["age_years", "age_at_diagnosis_years", "diagnosis_age_years"],
    "sex": ["sex_unified", "sex", "gender", "Sex", "Gender"],
}

GROUP_ALIASES = {
    "mycn_amp": {"amplified", "mycn_amplified", "mycnamp", "mycn_amp", "amp", "a", "1", "yes", "true", "amplification"},
    "mycn_nonamp": {"not_amplified", "non_amplified", "nonamplified", "non-amplified", "notamp", "nonamp", "namp", "na", "0", "no", "false"},
    "high_risk": {"high_risk", "highrisk", "high-risk", "high risk", "hr", "high", "1", "yes", "true"},
    "not_high_risk": {"low_risk", "lowrisk", "low-risk", "low risk", "nonhighrisk", "non_high_risk", "non-high-risk", "nothighrisk", "standardrisk", "standard_risk", "intermediate", "intermediate_risk", "0", "no", "false"},
    "age_ge18": {"ge18m", "ge_18m", "gt18m", "gt_18m", ">=18", ">=18m", "older18m", "1", "yes", "true"},
    "age_lt18": {"lt18m", "lt_18m", "le18m", "le_18m", "<18", "<18m", "under18m", "0", "no", "false"},
    "male": {"male", "m", "1"},
    "female": {"female", "f", "0"},
}


# Parse command-line inputs for expression/manifest modes, clinical-variable overrides, model controls, and outputs.
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the adjusted Cox workflow.
    
    Parameters
    ----------
    None
        Arguments are read from ``sys.argv`` by ``argparse``.
    
    Returns
    -------
    argparse.Namespace
        Parsed input paths, clinical-variable overrides, model settings, and
        output controls used by the rest of the workflow.
    """
    p = argparse.ArgumentParser(description="Robust adjusted Cox models for clock genes.")

    # Direct matrix mode.
    p.add_argument(
        "--expression",
        action="append",
        metavar="CONDITION:COHORT:PATH",
        help="Direct expression matrix. Repeat for each condition/cohort."
    )
    p.add_argument("--clinical", help="Clinical metadata TSV/CSV for direct matrix mode.")

    # Manifest mode.
    p.add_argument(
        "--manifest",
        action="append",
        metavar="CONDITION:PATH",
        help="NCC/R2 target manifest containing expression file paths and embedded metadata rows."
    )
    p.add_argument("--manifest_sep", default="\t")

    # Shared settings.
    p.add_argument("--out_prefix", default="adjusted_cox")
    p.add_argument("--clock_genes", default=",".join(CLOCK_GENES))
    p.add_argument("--include_clock_signature", action="store_true")
    p.add_argument("--ignore_asIs", action="store_true", default=True)
    p.add_argument("--expr_sep", default=None, help="Direct expression separator; auto if omitted.")
    p.add_argument("--clinical_sep", default=None, help="Clinical separator; auto if omitted.")

    # Explicit clinical overrides.
    p.add_argument("--sample_id_col", default=None)
    p.add_argument("--cohort_col", default=None)
    p.add_argument("--time_col", default=None)
    p.add_argument("--event_col", default=None)
    p.add_argument("--mycn_col", default=None)
    p.add_argument("--stage_col", default=None)
    p.add_argument("--risk_col", default=None)
    p.add_argument("--age18_col", default=None)
    p.add_argument("--age_col", default=None)
    p.add_argument("--age_unit", choices=["auto", "days", "months", "years"], default="auto")
    p.add_argument("--sex_col", default=None)

    # Model controls.
    p.add_argument("--include_model_c_sex", action="store_true")
    p.add_argument("--min_samples", type=int, default=30)
    p.add_argument("--min_events", type=int, default=10)
    p.add_argument("--min_events_per_variable", type=float, default=5.0)
    p.add_argument("--collinearity_r", type=float, default=0.95)
    p.add_argument("--max_drop_fraction_warning", type=float, default=0.30)
    p.add_argument("--cox_penalizer", type=float, default=0.0)
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--time_unit_label", default="as_provided", help="Used in the survival-variable resolution table.")

    # Output controls.
    p.add_argument("--write_xlsx", action="store_true", default=True)
    return p.parse_args()


# -----------------------------
# Generic helpers
# -----------------------------

# Normalize arbitrary labels for robust value matching (remove punctuation/case differences).
def norm_token(x: object) -> str:
    """Normalize an arbitrary value for case-insensitive token matching.
    
    Parameters
    ----------
    x : object
        Value to normalize.
    
    Returns
    -------
    str
        Lower-case string containing only alphanumeric characters.
    """
    return re.sub(r"[^a-z0-9]+", "", str(x).strip().lower())


# Normalize column names while retaining underscores between words for alias matching.
def norm_col(x: object) -> str:
    """Normalize a value for robust column-name matching.
    
    Parameters
    ----------
    x : object
        Column name or other label to normalize.
    
    Returns
    -------
    str
        Lower-case label in which non-alphanumeric runs are converted to
        underscores and leading/trailing underscores are removed.
    """
    return re.sub(r"[^a-z0-9]+", "_", str(x).strip().lower()).strip("_")


# Parse a comma-separated gene list and standardize gene symbols to uppercase.
def split_csv(x: str) -> List[str]:
    """Split a comma-separated gene list and standardize entries to upper case.
    
    Parameters
    ----------
    x : str
        Comma-separated string, typically supplied through ``--clock_genes``.
    
    Returns
    -------
    List[str]
        Non-empty, stripped, upper-case entries in their original order.
    """
    return [v.strip().upper() for v in str(x).split(",") if v.strip()]


# Determine whether a text input is tab- or comma-delimited unless explicitly supplied.
def guess_sep(path: Path, explicit: Optional[str]) -> str:
    """Determine the delimiter used by a delimited text file.
    
    Parameters
    ----------
    path : pathlib.Path
        Input file to inspect.
    explicit : Optional[str]
        Explicit delimiter supplied by the user. Escape sequences such as ``\t``
        are decoded before use. If ``None``, the delimiter is inferred.
    
    Returns
    -------
    str
        Delimiter to pass to ``pandas.read_csv``.
    """
    if explicit is not None:
        return explicit.encode().decode("unicode_escape")
    if path.suffix.lower() == ".csv":
        return ","
    try:
        first = path.read_text(errors="ignore").splitlines()[0]
        if first.count(",") > first.count("\t"):
            return ","
    except Exception:
        pass
    return "\t"


# Convert a value to float while treating common missing-value strings as missing.
def safe_float(x: object) -> Optional[float]:
    """Convert a value to float while treating common missing-value strings as missing.
    
    Parameters
    ----------
    x : object
        Value to convert.
    
    Returns
    -------
    Optional[float]
        Numeric value when conversion succeeds; otherwise ``None``.
    """
    if x is None:
        return None
    if isinstance(x, float) and math.isnan(x):
        return None
    s = str(x).strip()
    if s == "" or s.upper() in {"NA", "NAN", "NULL", "NONE"}:
        return None
    try:
        return float(s)
    except Exception:
        return None


# Convert a z statistic to a two-sided standard-normal p value.
def normal_two_sided_p(z: float) -> float:
    """Compute a two-sided standard-normal p value from a z statistic.
    
    Parameters
    ----------
    z : float
        Standard-normal test statistic.
    
    Returns
    -------
    float
        Two-sided p value ``P(|Z| >= |z|)``.
    """
    return float(math.erfc(abs(z) / math.sqrt(2.0)))


# Apply Benjamini-Hochberg FDR adjustment to the valid p values in one testing family.
def bh_adjust(pvalues: Sequence[object]) -> List[Optional[float]]:
    """Apply the Benjamini-Hochberg false-discovery-rate adjustment.
    
    Parameters
    ----------
    pvalues : Sequence[object]
        Sequence of raw p values. Missing, non-numeric, or out-of-range values are
        ignored and retain ``None`` in the returned list.
    
    Returns
    -------
    List[Optional[float]]
        BH-adjusted q values in the same order as the input sequence.
    """
    vals = []
    out: List[Optional[float]] = [None] * len(pvalues)
    for i, p in enumerate(pvalues):
        v = safe_float(p)
        if v is not None and 0 <= v <= 1:
            vals.append((i, v))
    if not vals:
        return out
    vals = sorted(vals, key=lambda z: z[1])
    m = len(vals)
    prev = 1.0
    for rev_rank, (idx, p) in enumerate(reversed(vals), start=1):
        rank = m - rev_rank + 1
        q = min(prev, p * m / rank)
        q = max(0.0, min(1.0, q))
        out[idx] = q
        prev = q
    return out


# Return the first dataframe column matching any candidate alias after normalization.
def first_col(df: pd.DataFrame, candidates: Sequence[str]) -> Optional[str]:
    """Find the first dataframe column matching any candidate alias.
    
    Parameters
    ----------
    df : pandas.DataFrame
        Table whose columns are searched.
    candidates : Sequence[str]
        Candidate column names in priority order.
    
    Returns
    -------
    Optional[str]
        Original dataframe column name for the first match, or ``None`` if no
        candidate is found.
    """
    lookup = {norm_col(c): c for c in df.columns}
    for c in candidates:
        hit = lookup.get(norm_col(c))
        if hit is not None:
            return hit
    return None


# Resolve a requested column, preferring an explicit user override and otherwise using aliases.
def find_col(df: pd.DataFrame, explicit: Optional[str], candidates: Sequence[str]) -> Optional[str]:
    """Resolve a column using an explicit override first, then aliases.
    
    Parameters
    ----------
    df : pandas.DataFrame
        Table whose columns are searched.
    explicit : Optional[str]
        User-supplied column name. Case-insensitive matching is attempted when an
        exact match is not present.
    candidates : Sequence[str]
        Fallback aliases searched in order.
    
    Returns
    -------
    Optional[str]
        Resolved original dataframe column name, or ``None``.
    """
    if explicit:
        if explicit in df.columns:
            return explicit
        lookup = {str(c).strip().lower(): c for c in df.columns}
        hit = lookup.get(str(explicit).strip().lower())
        if hit is not None:
            return hit
    return first_col(df, candidates)


# Map alternative expression-condition labels onto the canonical condition names.
def standardise_condition(x: object) -> str:
    """Map an expression-condition label to the workflow's canonical name.
    
    Parameters
    ----------
    x : object
        Raw condition label such as ``log2``, ``combat`` or ``cmbt``.
    
    Returns
    -------
    str
        Canonical condition label when recognized; otherwise the original label as
        a string.
    """
    key = norm_col(x)
    return CONDITION_ALIASES.get(str(x), CONDITION_ALIASES.get(key, str(x)))


# Harmonize common cohort/database labels to SEQC, Kocak, Versteeg, or NRC.
def standardise_cohort(x: object) -> str:
    """Map dataset/cohort labels to canonical cohort names.
    
    Parameters
    ----------
    x : object
        Raw cohort or dataset identifier.
    
    Returns
    -------
    str
        Canonical cohort name such as ``SEQC`` or ``Kocak`` when recognized;
        otherwise the stripped input string.
    """
    key = norm_token(x)
    if "seqc" in key or "gse62564" in key:
        return "SEQC"
    if "kocak" in key or "gse45547" in key:
        return "Kocak"
    if "versteeg" in key:
        return "Versteeg"
    if "nrc" in key:
        return "NRC"
    return str(x).strip()


# Convert heterogeneous categorical labels into a numeric 1/0 coding.
def canonical_binary(x: object, pos: Iterable[str], neg: Iterable[str]) -> float:
    """Convert heterogeneous binary labels to numeric 1/0 coding.
    
    Parameters
    ----------
    x : object
        Raw value to classify.
    pos : Iterable[str]
        Labels interpreted as the positive class (1).
    neg : Iterable[str]
        Labels interpreted as the negative class (0).
    
    Returns
    -------
    float
        ``1.0`` for a positive label, ``0.0`` for a negative label, or ``NaN``
        when the value cannot be classified.
    """
    key = norm_token(x)
    pos_set = {norm_token(v) for v in pos}
    neg_set = {norm_token(v) for v in neg}
    if key in pos_set:
        return 1.0
    if key in neg_set:
        return 0.0
    return np.nan


# Convert survival-status labels to event=1 and censored/alive=0.
def coerce_event(x: object) -> float:
    """Convert survival-event labels to numeric Cox event coding.
    
    Parameters
    ----------
    x : object
        Raw event/status value.
    
    Returns
    -------
    float
        ``1.0`` for death/event, ``0.0`` for alive/censored, or ``NaN`` when the
        value cannot be resolved.
    """
    key = norm_token(x)
    if key in {norm_token(v) for v in EVENT_TRUE}:
        return 1.0
    if key in {norm_token(v) for v in EVENT_FALSE}:
        return 0.0
    v = safe_float(x)
    if v is None:
        return np.nan
    return 1.0 if v != 0 else 0.0


# Harmonize INSS stage labels before dummy-variable construction.
def canonical_stage_label(x: object) -> str:
    """Standardize heterogeneous INSS stage labels.
    
    Parameters
    ----------
    x : object
        Raw stage value.
    
    Returns
    -------
    str
        Canonical stage label (``stage1``, ``stage2``, ``stage3``, ``stage4``,
        ``stage4s`` or ``NA``). Unrecognized non-missing labels are returned as
        stripped strings.
    """
    key = norm_token(x)
    if key in {"st1", "stage1", "stagei", "inss1", "1", "i"}:
        return "stage1"
    if key in {"st2", "stage2", "stageii", "inss2", "2", "ii", "2a", "2b", "stage2a", "stage2b"}:
        return "stage2"
    if key in {"st3", "stage3", "stageiii", "inss3", "3", "iii"}:
        return "stage3"
    if key in {"st4", "stage4", "stageiv", "inss4", "4", "iv"}:
        return "stage4"
    if key in {"st4s", "stage4s", "4s", "ivs"}:
        return "stage4s"
    if key in {"", "na", "nan", "none", "unknown"}:
        return "NA"
    return str(x).strip()


# Standardize expression within cohort so HRs correspond to a 1-SD expression increase.
def zscore(s: pd.Series) -> pd.Series:
    """Standardize a numeric series to mean 0 and sample standard deviation 1.
    
    Parameters
    ----------
    s : pandas.Series
        Values to standardize.
    
    Returns
    -------
    pandas.Series
        Z-scored values aligned to the input index. If the standard deviation is
        zero or undefined, all returned values are ``NaN``.
    """
    x = pd.to_numeric(s, errors="coerce")
    sd = x.std(ddof=1)
    if pd.isna(sd) or sd == 0:
        return pd.Series(np.nan, index=s.index)
    return (x - x.mean()) / sd


# -----------------------------
# Input readers
# -----------------------------

def read_r2_expression(path: Path, sep: Optional[str], warnings_rows: List[Dict[str, str]]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Read an R2-style expression matrix and embedded metadata.
    
    Parameters
    ----------
    path : pathlib.Path
        Path to the gene-by-sample R2 expression matrix.
    sep : Optional[str]
        Field separator. If omitted, it is inferred from the file.
    warnings_rows : List[Dict[str, str]]
        Mutable warning collector to which data-quality messages are appended.
    
    Returns
    -------
    Tuple[pandas.DataFrame, pandas.DataFrame]
        ``expr``: numeric genes-by-samples expression matrix with upper-case gene
        identifiers; ``obs``: samples-by-metadata table derived from ``#`` rows.
    """
    separator = guess_sep(path, sep)
    df = pd.read_csv(path, sep=separator, header=0, index_col=0, low_memory=False)

    obs = df.loc[df.index.astype(str).str.startswith("#")].copy()
    expr_raw = df.loc[~df.index.astype(str).str.startswith("#")].copy()

    if obs.empty:
        # Simple gene-by-sample matrix without metadata.
        obs = pd.DataFrame(index=df.columns)
    else:
        obs.index = obs.index.astype(str).str.replace("#", "", regex=False)
        obs = obs.T

    expr_num = expr_raw.apply(pd.to_numeric, errors="coerce")
    non_numeric_cols = [c for c in expr_num.columns if expr_num[c].isna().all() and not expr_raw[c].isna().all()]
    if non_numeric_cols:
        warnings_rows.append({
            "severity": "WARNING",
            "topic": "non_numeric_expression_columns_dropped",
            "message": f"{path}: dropped non-numeric columns: {', '.join(map(str, non_numeric_cols))}",
        })
        expr_num = expr_num.drop(columns=non_numeric_cols)
        obs = obs.drop(index=[c for c in non_numeric_cols if c in obs.index], errors="ignore")

    expr_num.index = expr_num.index.astype(str).str.upper()
    expr_num = expr_num.groupby(expr_num.index).mean()

    # Align obs to expression columns.
    obs = obs.reindex(expr_num.columns)
    return expr_num, obs


def read_manifest_inputs(manifest_path: Path, condition: str, args: argparse.Namespace, warnings_rows: List[Dict[str, str]]) -> List[Tuple[str, str, pd.DataFrame, pd.DataFrame, str]]:
    """Load all cohort datasets referenced by one NCC target manifest.
    
    Parameters
    ----------
    manifest_path : pathlib.Path
        Path to the NCC/R2 manifest.
    condition : str
        Canonical expression condition assigned to all rows from this manifest.
    args : argparse.Namespace
        Parsed workflow options, including manifest separator.
    warnings_rows : List[Dict[str, str]]
        Mutable warning collector.
    
    Returns
    -------
    List[Tuple[str, str, pandas.DataFrame, pandas.DataFrame, str]]
        Tuples containing condition, cohort, expression matrix, metadata table,
        and manifest database name for each valid entry.
    """
    cfg = pd.read_csv(manifest_path, sep=args.manifest_sep, header=0, index_col=0, low_memory=False)
    items = []

    for db_name in cfg.index:
        row = cfg.loc[db_name]
        cohort = standardise_cohort(db_name)
        path = Path(str(row["full_path_input_expression_file"])).expanduser()
        if not path.exists():
            warnings_rows.append({
                "severity": "ERROR",
                "topic": "expression_file_missing",
                "message": f"{db_name}: {path} not found",
            })
            continue

        file_sep = str(row.get("separator_of_columns_in_input_expression_file", "\t")).encode().decode("unicode_escape")
        expr, obs = read_r2_expression(path, file_sep, warnings_rows)

        # Respect column exclusions by dropping samples.
        exc = split_optional(row.get("columns_to_exclude_in_input_expression_file"))
        if exc:
            matched = case_insensitive_matches(expr.columns, exc)
            if matched:
                expr = expr.drop(columns=matched)
                obs = obs.drop(index=matched, errors="ignore")

        # Apply base manifest filters.
        expr, obs = apply_manifest_filters(expr, obs, row)

        # Try to add a cohort label to obs if missing.
        obs["_manifest_db_name"] = str(db_name)
        obs["_manifest_cohort"] = cohort
        items.append((condition, cohort, expr, obs, str(db_name)))

    return items


# Parse optional manifest list fields that may be comma- or pipe-separated.
def split_optional(value: object) -> Optional[List[str]]:
    """Parse an optional comma- or pipe-separated manifest field.
    
    Parameters
    ----------
    value : object
        Manifest value that may be missing, scalar, comma-separated, or
        pipe-separated.
    
    Returns
    -------
    Optional[List[str]]
        Stripped non-empty entries, or ``None`` when no value is supplied.
    """
    if value is None or pd.isna(value):
        return None
    s = str(value).strip()
    if not s:
        return None
    return [p.strip() for p in s.replace(",", "|").split("|") if p.strip()]


# Match sample/column labels without case sensitivity while preserving original labels.
def case_insensitive_matches(existing: Sequence[str], targets: Sequence[str]) -> List[str]:
    """Match requested labels to existing labels without case sensitivity.
    
    Parameters
    ----------
    existing : Sequence[str]
        Available labels, typically dataframe column names.
    targets : Sequence[str]
        Labels requested by the manifest or user.
    
    Returns
    -------
    List[str]
        Original labels from ``existing`` that match requested targets.
    """
    lookup = {str(x).strip().lower(): str(x) for x in existing}
    out = []
    for t in targets:
        key = str(t).strip().lower()
        if key in lookup:
            out.append(lookup[key])
    return out


def apply_manifest_filters(expr: pd.DataFrame, obs: pd.DataFrame, row: pd.Series) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Apply sample-selection and exclusion rules from a manifest row.
    
    Parameters
    ----------
    expr : pandas.DataFrame
        Genes-by-samples expression matrix.
    obs : pandas.DataFrame
        Samples-by-metadata table aligned to ``expr`` columns.
    row : pandas.Series
        Manifest row containing optional variable/group/patient filters.
    
    Returns
    -------
    Tuple[pandas.DataFrame, pandas.DataFrame]
        Filtered expression matrix and metadata table containing the same samples.
    """
    variable = row.get("variable_of_interest_in_patients_in_input_expression_file")
    if pd.isna(variable):
        variable = None
    groups_to_select = split_optional(row.get("groups_to_select_in_variable_of_interest_in_patients_in_input_expression_file"))
    groups_to_exclude = split_optional(row.get("groups_to_exclude_in_variable_of_interest_in_patients_in_input_expression_file"))
    patients_to_select = split_optional(row.get("patients_to_select_in_input_expression_file"))

    keep = pd.Series(True, index=obs.index)
    if variable is not None and variable in obs.columns and groups_to_select:
        keep &= obs[variable].isin(groups_to_select)
    if variable is not None and variable in obs.columns and groups_to_exclude:
        keep &= ~obs[variable].isin(groups_to_exclude)
    if patients_to_select:
        keep &= obs.index.isin(patients_to_select)

    obs2 = obs.loc[keep].copy()
    expr2 = expr.loc[:, obs2.index].copy()
    return expr2, obs2


# Parse direct-expression arguments formatted as CONDITION:COHORT:PATH.
def parse_expression_arg(item: str) -> Tuple[str, str, Path]:
    """Parse one direct-expression command-line specification.
    
    Parameters
    ----------
    item : str
        String formatted as ``CONDITION:COHORT:PATH``.
    
    Returns
    -------
    Tuple[str, str, pathlib.Path]
        Canonical condition, canonical cohort, and validated expression-file path.
    
    Raises
    ------
    ValueError
        If the argument does not contain three fields.
    FileNotFoundError
        If the expression file does not exist.
    """
    parts = item.split(":", 2)
    if len(parts) != 3:
        raise ValueError(f"--expression must be CONDITION:COHORT:PATH, got {item}")
    condition = standardise_condition(parts[0])
    cohort = standardise_cohort(parts[1])
    path = Path(parts[2]).expanduser()
    if not path.exists():
        raise FileNotFoundError(path)
    return condition, cohort, path


def read_direct_inputs(args: argparse.Namespace, warnings_rows: List[Dict[str, str]]) -> List[Tuple[str, str, pd.DataFrame, pd.DataFrame, str]]:
    """Load direct expression matrices and join them to clinical metadata.
    
    Parameters
    ----------
    args : argparse.Namespace
        Parsed options containing ``--expression`` entries, clinical metadata path,
        separators, and optional column overrides.
    warnings_rows : List[Dict[str, str]]
        Mutable warning collector.
    
    Returns
    -------
    List[Tuple[str, str, pandas.DataFrame, pandas.DataFrame, str]]
        Condition/cohort/expression/metadata/database tuples for valid direct
        inputs. Returns an empty list when no direct expression inputs were given.
    
    Raises
    ------
    ValueError
        If direct-expression mode is requested without clinical metadata or if a
        sample-ID column cannot be resolved.
    """
    if not args.expression:
        return []
    if not args.clinical:
        raise ValueError("Direct expression mode requires --clinical")

    clin_path = Path(args.clinical).expanduser()
    clin_sep = guess_sep(clin_path, args.clinical_sep)
    clinical = pd.read_csv(clin_path, sep=clin_sep, dtype=str, keep_default_na=False, low_memory=False)

    sample_col = find_col(clinical, args.sample_id_col, CLINICAL_ALIASES["sample_id"])
    cohort_col = find_col(clinical, args.cohort_col, CLINICAL_ALIASES["cohort"])
    if sample_col is None:
        raise ValueError("Could not infer sample ID column in clinical metadata; use --sample_id_col")
    clinical["_sample_id"] = clinical[sample_col].astype(str)
    clinical["_clinical_cohort"] = clinical[cohort_col].map(standardise_cohort) if cohort_col else "NA"
    clinical = clinical.set_index("_sample_id", drop=False)

    items = []
    for expr_arg in args.expression:
        condition, cohort, path = parse_expression_arg(expr_arg)
        expr, obs0 = read_r2_expression(path, args.expr_sep, warnings_rows)
        matched = [s for s in expr.columns if str(s) in clinical.index]
        expr = expr.loc[:, matched].copy()
        obs = clinical.loc[matched].copy()
        if cohort_col:
            obs = obs.loc[obs["_clinical_cohort"].eq(cohort)].copy()
            expr = expr.loc[:, obs.index].copy()
        obs["_manifest_db_name"] = cohort
        obs["_manifest_cohort"] = cohort
        items.append((condition, cohort, expr, obs, cohort))
    return items


def read_all_inputs(args: argparse.Namespace, warnings_rows: List[Dict[str, str]]) -> List[Tuple[str, str, pd.DataFrame, pd.DataFrame, str]]:
    """Collect, standardize, and filter all manifest and direct input datasets.
    
    Parameters
    ----------
    args : argparse.Namespace
        Parsed workflow options.
    warnings_rows : List[Dict[str, str]]
        Mutable warning collector.
    
    Returns
    -------
    List[Tuple[str, str, pandas.DataFrame, pandas.DataFrame, str]]
        Valid primary-cohort datasets restricted to allowed expression conditions.
    """
    items = []
    if args.manifest:
        for item in args.manifest:
            if ":" not in item:
                raise ValueError(f"--manifest must be CONDITION:PATH, got {item}")
            cond_raw, path_str = item.split(":", 1)
            condition = standardise_condition(cond_raw)
            path = Path(path_str).expanduser()
            if not path.exists():
                raise FileNotFoundError(path)
            items.extend(read_manifest_inputs(path, condition, args, warnings_rows))
    items.extend(read_direct_inputs(args, warnings_rows))

    # Ignore asIs/untransformed by default.
    filtered = []
    for condition, cohort, expr, obs, db_name in items:
        if args.ignore_asIs and condition not in ALLOWED_CONDITIONS:
            warnings_rows.append({
                "severity": "INFO",
                "topic": "condition_ignored",
                "message": f"Ignored condition {condition} for {db_name}",
            })
            continue
        if cohort not in PRIMARY_COHORTS:
            warnings_rows.append({
                "severity": "INFO",
                "topic": "cohort_ignored",
                "message": f"Ignored non-primary cohort {cohort} for {db_name}",
            })
            continue
        filtered.append((condition, cohort, expr, obs, db_name))
    return filtered


# -----------------------------
# Variable resolution and covariates
# -----------------------------

def resolve_time_event(obs: pd.DataFrame, cohort: str, db_name: str, args: argparse.Namespace, warnings_rows: List[Dict[str, str]]) -> Tuple[Optional[str], Optional[str], pd.Series, pd.Series, Dict[str, object]]:
    """Resolve and validate survival-time and event variables for one cohort.
    
    Parameters
    ----------
    obs : pandas.DataFrame
        Sample-level clinical metadata.
    cohort : str
        Canonical cohort name used to prioritize cohort-specific aliases.
    db_name : str
        Dataset identifier used in diagnostic messages.
    args : argparse.Namespace
        Parsed options containing optional explicit time/event column overrides and
        the time-unit label.
    warnings_rows : List[Dict[str, str]]
        Mutable warning collector.
    
    Returns
    -------
    Tuple[Optional[str], Optional[str], pandas.Series, pandas.Series, Dict[str, object]]
        Resolved time-column name, event-column name, numeric survival-time series,
        numeric 0/1 event series, and a resolution/audit record.
    """
    time_candidates = list(SURVIVAL_ALIASES.get(cohort, {}).get("time", [])) + SURVIVAL_ALIASES["default"]["time"]
    event_candidates = list(SURVIVAL_ALIASES.get(cohort, {}).get("event", [])) + SURVIVAL_ALIASES["default"]["event"]
    time_col = find_col(obs, args.time_col, time_candidates)
    event_col = find_col(obs, args.event_col, event_candidates)

    if time_col is None or event_col is None:
        warnings_rows.append({
            "severity": "ERROR",
            "topic": "survival_variable_resolution_failed",
            "message": f"{db_name}/{cohort}: time_col={time_col}, event_col={event_col}; available={list(obs.columns)}",
        })
        resolution = {
            "db_name": db_name, "cohort": cohort,
            "requested_time_variable": args.time_col or "auto",
            "resolved_time_variable": time_col or "UNRESOLVED",
            "requested_event_variable": args.event_col or "auto",
            "resolved_event_variable": event_col or "UNRESOLVED",
            "event_coding_rule": "1=death/event,0=censored/alive",
            "time_unit_label": args.time_unit_label,
            "n_start": len(obs), "n_complete": 0, "n_events": 0, "n_censored": 0,
            "missing_fraction": 1.0,
            "status": "failed_resolution",
        }
        return time_col, event_col, pd.Series(np.nan, index=obs.index), pd.Series(np.nan, index=obs.index), resolution

    time = pd.to_numeric(obs[time_col], errors="coerce")
    event = obs[event_col].map(coerce_event)

    valid_time = time.notna() & np.isfinite(time) & (time > 0)
    valid_event = event.notna() & event.isin([0.0, 1.0])
    complete = valid_time & valid_event

    resolution = {
        "db_name": db_name,
        "cohort": cohort,
        "requested_time_variable": args.time_col or "auto",
        "resolved_time_variable": time_col,
        "requested_event_variable": args.event_col or "auto",
        "resolved_event_variable": event_col,
        "event_coding_rule": "1=death/event,0=censored/alive",
        "time_unit_label": args.time_unit_label,
        "n_start": int(len(obs)),
        "n_complete": int(complete.sum()),
        "n_events": int(event.loc[complete].sum()) if complete.any() else 0,
        "n_censored": int((1 - event.loc[complete]).sum()) if complete.any() else 0,
        "missing_fraction": float(1 - complete.mean()) if len(obs) else np.nan,
        "status": "ok" if complete.any() else "no_complete_time_event",
    }

    if not complete.any():
        warnings_rows.append({
            "severity": "ERROR",
            "topic": "no_complete_time_event",
            "message": f"{db_name}/{cohort}: no complete numeric positive time and 0/1 event pairs.",
        })
    return time_col, event_col, time, event, resolution


# Resolve or derive the >=18-month age indicator from binary or continuous age fields.
def infer_age18(obs: pd.DataFrame, args: argparse.Namespace) -> pd.Series:
    """Create a binary age-at-diagnosis >=18-month covariate.
    
    Parameters
    ----------
    obs : pandas.DataFrame
        Sample-level clinical metadata.
    args : argparse.Namespace
        Parsed age-column and age-unit overrides.
    
    Returns
    -------
    pandas.Series
        ``1.0`` for age >=18 months, ``0.0`` for age <18 months, and ``NaN`` when
        age cannot be resolved.
    """
    age18_col = find_col(obs, args.age18_col, CLINICAL_ALIASES["age18"])
    if age18_col:
        return obs[age18_col].map(lambda x: canonical_binary(x, GROUP_ALIASES["age_ge18"], GROUP_ALIASES["age_lt18"]))

    age_col = None
    unit = args.age_unit
    if args.age_col:
        age_col = find_col(obs, args.age_col, [args.age_col])
    elif unit != "auto":
        age_col = find_col(obs, None, CLINICAL_ALIASES[f"age_{unit}"])
    else:
        for u in ["months", "days", "years"]:
            hit = find_col(obs, None, CLINICAL_ALIASES[f"age_{u}"])
            if hit:
                age_col = hit
                unit = u
                break

    if age_col is None:
        return pd.Series(np.nan, index=obs.index)

    age = pd.to_numeric(obs[age_col], errors="coerce")
    if unit == "days":
        months = age / 30.4375
    elif unit == "years":
        months = age * 12.0
    else:
        if unit == "auto":
            med = age.dropna().median()
            if pd.notna(med) and med > 365:
                months = age / 30.4375
            elif pd.notna(med) and med < 25:
                months = age * 12.0
            else:
                months = age
        else:
            months = age
    return (months >= 18.0).astype(float).where(months.notna(), np.nan)


# Harmonize clinical covariates and create INSS-stage dummy variables for Cox adjustment.
def build_covariate_frame(obs: pd.DataFrame, args: argparse.Namespace, warnings_rows: List[Dict[str, str]], db_name: str) -> pd.DataFrame:
    """Construct standardized clinical covariates for Cox adjustment models.
    
    Parameters
    ----------
    obs : pandas.DataFrame
        Sample-level clinical metadata.
    args : argparse.Namespace
        Parsed clinical-column overrides.
    warnings_rows : List[Dict[str, str]]
        Mutable warning collector for unresolved covariates.
    db_name : str
        Dataset identifier included in warnings.
    
    Returns
    -------
    pandas.DataFrame
        Sample-indexed covariate table containing age group, MYCN status, high-risk
        status, sex, and INSS-stage dummy variables when available.
    """
    cov = pd.DataFrame(index=obs.index)

    mycn_col = find_col(obs, args.mycn_col, CLINICAL_ALIASES["mycn"])
    stage_col = find_col(obs, args.stage_col, CLINICAL_ALIASES["stage"])
    risk_col = find_col(obs, args.risk_col, CLINICAL_ALIASES["risk"])
    sex_col = find_col(obs, args.sex_col, CLINICAL_ALIASES["sex"])

    cov["age_group_18m"] = infer_age18(obs, args)

    if mycn_col:
        cov["MYCN_status"] = obs[mycn_col].map(lambda x: canonical_binary(x, GROUP_ALIASES["mycn_amp"], GROUP_ALIASES["mycn_nonamp"]))
    else:
        cov["MYCN_status"] = np.nan
        warnings_rows.append({"severity": "WARNING", "topic": "MYCN_missing", "message": f"{db_name}: MYCN covariate not resolved."})

    if risk_col:
        cov["harmonised_high_risk"] = obs[risk_col].map(lambda x: canonical_binary(x, GROUP_ALIASES["high_risk"], GROUP_ALIASES["not_high_risk"]))
    else:
        cov["harmonised_high_risk"] = np.nan
        warnings_rows.append({"severity": "WARNING", "topic": "risk_missing", "message": f"{db_name}: harmonised risk covariate not resolved."})

    if sex_col:
        cov["sex_male"] = obs[sex_col].map(lambda x: canonical_binary(x, GROUP_ALIASES["male"], GROUP_ALIASES["female"]))
    else:
        cov["sex_male"] = np.nan
        warnings_rows.append({"severity": "INFO", "topic": "sex_missing", "message": f"{db_name}: sex covariate not resolved."})

    if stage_col:
        stage = obs[stage_col].map(canonical_stage_label)
        dummies = pd.get_dummies(stage, prefix="INSS_stage", drop_first=True, dtype=float)
        dummies = dummies.loc[:, [c for c in dummies.columns if not c.endswith("_NA") and "unknown" not in c.lower()]]
        cov = pd.concat([cov, dummies], axis=1)
        if dummies.empty:
            warnings_rows.append({"severity": "WARNING", "topic": "stage_no_dummies", "message": f"{db_name}: stage resolved but no useful dummy columns were generated."})
    else:
        warnings_rows.append({"severity": "WARNING", "topic": "stage_missing", "message": f"{db_name}: INSS stage covariate not resolved."})

    return cov


# Build standardized gene-level predictors and, optionally, a mean clock-expression signature.
def build_expression_features(expr: pd.DataFrame, obs: pd.DataFrame, genes: Sequence[str], include_signature: bool, condition: str, cohort: str, warnings_rows: List[Dict[str, str]]) -> pd.DataFrame:
    """Build standardized clock-gene expression predictors.
    
    Parameters
    ----------
    expr : pandas.DataFrame
        Genes-by-samples expression matrix.
    obs : pandas.DataFrame
        Sample metadata whose index defines the analysis sample order.
    genes : Sequence[str]
        Clock genes to extract.
    include_signature : bool
        If ``True``, also compute the mean standardized clock-expression signature.
    condition : str
        Expression condition used in warning messages.
    cohort : str
        Cohort name used in warning messages.
    warnings_rows : List[Dict[str, str]]
        Mutable warning collector for missing genes.
    
    Returns
    -------
    pandas.DataFrame
        Sample-indexed z-scored expression features, optionally including
        ``CLOCK_SIGNATURE``.
    """
    features = pd.DataFrame(index=obs.index)
    present = []
    for gene in genes:
        g = gene.upper()
        if g not in expr.index:
            warnings_rows.append({"severity": "WARNING", "topic": "gene_missing", "message": f"{condition}/{cohort}: {g} not found."})
            continue
        features[g] = zscore(pd.Series(expr.loc[g, obs.index].astype(float).values, index=obs.index))
        present.append(g)

    if include_signature and present:
        features["CLOCK_SIGNATURE"] = features[present].mean(axis=1, skipna=True)
    return features


# -----------------------------
# Cox fitting
# -----------------------------

def drop_bad_covariates(data: pd.DataFrame, covariates: List[str], event_col: str, args: argparse.Namespace) -> Tuple[List[str], List[str]]:
    """Remove unusable or unstable covariates before Cox fitting.
    
    Parameters
    ----------
    data : pandas.DataFrame
        Complete-case model dataframe containing event status and candidate
        covariates.
    covariates : List[str]
        Candidate adjustment-variable names.
    event_col : str
        Name of the binary event column in ``data``.
    args : argparse.Namespace
        Parsed thresholds, including the collinearity cutoff.
    
    Returns
    -------
    Tuple[List[str], List[str]]
        Covariates retained for fitting and diagnostic notes describing any dropped
        variables.
    """
    notes = []
    kept = []
    for c in covariates:
        if c not in data.columns:
            notes.append(f"missing_covariate:{c}")
            continue
        if data[c].nunique(dropna=True) <= 1:
            notes.append(f"dropped_zero_variance:{c}")
            continue
        kept.append(c)

    # Drop perfectly/highly collinear later covariate in pair.
    if len(kept) > 1:
        corr = data[kept].corr().abs()
        to_drop = set()
        for i, c1 in enumerate(kept):
            for c2 in kept[i+1:]:
                val = corr.loc[c1, c2]
                if pd.notna(val) and val >= args.collinearity_r:
                    to_drop.add(c2)
                    notes.append(f"dropped_collinear:{c2}_with_{c1}_r={val:.3f}")
        kept = [c for c in kept if c not in to_drop]

    # Complete separation-like flag for binary covariates; drop if any event cell is zero.
    final = []
    for c in kept:
        vals = set(pd.Series(data[c]).dropna().unique())
        if vals.issubset({0, 1, 0.0, 1.0}) and len(vals) == 2:
            tab = pd.crosstab(data[c], data[event_col])
            if (tab == 0).any().any():
                notes.append(f"dropped_possible_complete_separation:{c}")
                continue
        final.append(c)
    return final, notes


def fit_cox_model(data_raw: pd.DataFrame, feature: str, covariates: List[str], model_name: str, args: argparse.Namespace) -> Tuple[Dict[str, object], pd.DataFrame]:
    """Fit one Cox proportional-hazards model and run PH diagnostics.
    
    Parameters
    ----------
    data_raw : pandas.DataFrame
        Sample-level dataframe containing ``time``, ``event``, the requested
        expression feature, and candidate clinical covariates.
    feature : str
        Expression feature to model. It is renamed internally to ``expression_z``.
    covariates : List[str]
        Requested adjustment covariates for this model.
    model_name : str
        Human-readable model identifier used for reporting.
    args : argparse.Namespace
        Fitting thresholds, penalizer, and PH-analysis settings.
    
    Returns
    -------
    Tuple[Dict[str, object], pandas.DataFrame]
        Model summary record containing fit status, HR, coefficient, SE, CI, p value,
        sample/event counts and PH summary statistics; plus term-level PH diagnostic
        results when available.
    """
    from lifelines import CoxPHFitter
    from lifelines.statistics import proportional_hazard_test

    n_start = len(data_raw)
    cols = ["time", "event", feature] + covariates
    present_cols = [c for c in cols if c in data_raw.columns]
    data = data_raw[present_cols].replace([np.inf, -np.inf], np.nan).copy()
    data = data.rename(columns={feature: "expression_z"})

    # Initial complete cases for feature + requested covariates.
    data = data.dropna()
    missing_n = n_start - len(data)
    missing_fraction = missing_n / n_start if n_start else np.nan

    # Event/time validation.
    if "time" not in data.columns or "event" not in data.columns:
        return {
            "fit_status": "not_fit_survival_variables_unresolved",
            "model_warnings": "time_or_event_missing_after_resolution",
            "n_start": n_start, "n": len(data), "n_events": 0,
            "n_missing_for_model": missing_n, "missing_fraction_for_model": missing_fraction,
        }, pd.DataFrame()

    data["time"] = pd.to_numeric(data["time"], errors="coerce")
    data["event"] = pd.to_numeric(data["event"], errors="coerce")
    data = data.loc[data["time"].notna() & (data["time"] > 0) & data["event"].isin([0, 1])].copy()

    # Drop bad covariates after complete case filtering.
    requested_covars_present = [c for c in covariates if c in data.columns]
    kept_covars, covar_notes = drop_bad_covariates(data, requested_covars_present, "event", args)
    model_cols = ["time", "event", "expression_z"] + kept_covars
    data = data[model_cols].dropna().copy()

    n = len(data)
    n_events = int(data["event"].sum()) if n else 0
    n_predictors = 1 + len(kept_covars)
    epv = n_events / n_predictors if n_predictors else np.nan

    notes = list(covar_notes)
    fit_status = "not_fit"
    convergence_message = ""
    beta = se = hr = ci_low = ci_high = pvalue = np.nan
    ph_gene_p = ph_global_min_p = np.nan
    ph_df = pd.DataFrame()

    if missing_fraction > args.max_drop_fraction_warning:
        notes.append(f"large_complete_case_loss={missing_fraction:.3f}")

    if n < args.min_samples:
        fit_status = "not_fit_sample_count"
        notes.append(f"n={n}<min_samples={args.min_samples}")
    elif n_events < args.min_events:
        fit_status = "not_fit_event_count"
        notes.append(f"events={n_events}<min_events={args.min_events}")
    elif epv < args.min_events_per_variable:
        fit_status = "not_fit_event_per_variable"
        notes.append(f"EPV={epv:.3f}<min_EPV={args.min_events_per_variable}")
    elif data["expression_z"].nunique(dropna=True) <= 1:
        fit_status = "not_fit_zero_variance_expression"
        notes.append("zero_variance_expression_z")
    else:
        try:
            cph = CoxPHFitter(penalizer=args.cox_penalizer)
            cph.fit(data, duration_col="time", event_col="event", show_progress=False)
            fit_status = "fit_ok"
            row = cph.summary.loc["expression_z"]
            beta = float(row["coef"])
            se = float(row["se(coef)"])
            hr = float(row["exp(coef)"])
            ci_low = float(row["exp(coef) lower 95%"])
            ci_high = float(row["exp(coef) upper 95%"])
            pvalue = float(row["p"])

            try:
                ph = proportional_hazard_test(cph, data, time_transform="rank")
                ph_out = ph.summary.reset_index().rename(columns={"index": "covariate"})
                ph_out["PH_p"] = pd.to_numeric(ph_out["p"], errors="coerce")
                if "expression_z" in set(ph_out["covariate"]):
                    ph_gene_p = safe_float(ph_out.loc[ph_out["covariate"].eq("expression_z"), "PH_p"].iloc[0])
                ph_global_min_p = safe_float(ph_out["PH_p"].min())
                ph_df = ph_out
            except Exception as e:
                notes.append("PH_diagnostics_failed:" + str(e)[:200])
        except Exception as e:
            fit_status = "fit_failed"
            convergence_message = str(e)[:500]
            notes.append("fit_failed")

    rec = {
        "fit_status": fit_status,
        "model_warnings": "|".join(notes) if notes else "none",
        "n_start": n_start,
        "n": n,
        "n_events": n_events,
        "n_censored": int(n - n_events) if n else 0,
        "n_missing_for_model": int(missing_n),
        "missing_fraction_for_model": missing_fraction,
        "n_predictors_final": n_predictors,
        "events_per_variable": epv,
        "covariates_requested": "|".join(covariates) if covariates else "none",
        "covariates_used": "|".join(kept_covars) if kept_covars else "none",
        "hr": hr,
        "beta": beta,
        "se": se,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "pvalue": pvalue,
        "PH_gene_term_p": ph_gene_p,
        "PH_global_min_p": ph_global_min_p,
        "convergence_message": convergence_message,
    }
    return rec, ph_df


# Define the prespecified Cox models: univariable, adjusted A/B, and optional adjusted C.
def model_specs(include_sex: bool) -> Dict[str, List[str]]:
    """Define the prespecified Cox adjustment models.
    
    Parameters
    ----------
    include_sex : bool
        Whether to include model C, which adds sex to model A.
    
    Returns
    -------
    Dict[str, List[str]]
        Mapping from model name to requested covariate tokens. ``STAGE_DUMMIES`` is
        expanded later to the stage dummy columns available in each cohort.
    """
    specs = {
        "univariable_expression": [],
        "adjusted_A_age18_MYCN_stage": ["age_group_18m", "MYCN_status", "STAGE_DUMMIES"],
        "adjusted_B_highrisk_age18": ["harmonised_high_risk", "age_group_18m"],
    }
    if include_sex:
        specs["adjusted_C_age18_MYCN_stage_sex"] = ["age_group_18m", "MYCN_status", "STAGE_DUMMIES", "sex_male"]
    return specs


# Replace the STAGE_DUMMIES placeholder with the actual stage dummy columns available in a cohort.
def expand_stage_token(tokens: List[str], covariate_df: pd.DataFrame) -> List[str]:
    """Replace the stage placeholder with available INSS-stage dummy columns.
    
    Parameters
    ----------
    tokens : List[str]
        Covariate specification containing ordinary column names and optionally
        ``STAGE_DUMMIES``.
    covariate_df : pandas.DataFrame
        Covariate table from which available stage-dummy columns are discovered.
    
    Returns
    -------
    List[str]
        Expanded covariate list ready for model fitting.
    """
    out = []
    for t in tokens:
        if t == "STAGE_DUMMIES":
            out.extend([c for c in covariate_df.columns if c.startswith("INSS_stage_")])
        else:
            out.append(t)
    return out


# Orchestrate input loading, variable resolution, Cox fitting, PH testing, BH correction, and summaries.
def run_all(args: argparse.Namespace) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Execute the complete adjusted Cox workflow across all inputs.
    
    Parameters
    ----------
    args : argparse.Namespace
        Parsed workflow configuration.
    
    Returns
    -------
    Tuple[pandas.DataFrame, pandas.DataFrame, pandas.DataFrame, pandas.DataFrame, pandas.DataFrame, pandas.DataFrame]
        Detail table, summary table, term-level PH diagnostics, fixed-effect
        meta-analysis table, survival-variable resolution table, and warning table.
    """
    warnings_rows: List[Dict[str, str]] = []
    resolution_rows: List[Dict[str, object]] = []
    detail_rows: List[Dict[str, object]] = []
    ph_rows: List[pd.DataFrame] = []

    genes = split_csv(args.clock_genes)
    inputs = read_all_inputs(args, warnings_rows)

    if not inputs:
        warnings_rows.append({"severity": "ERROR", "topic": "no_inputs", "message": "No valid log2/mycn_corrected SEQC/Kocak inputs found."})

    for condition, cohort, expr, obs, db_name in inputs:
        time_col, event_col, time, event, resolution = resolve_time_event(obs, cohort, db_name, args, warnings_rows)
        resolution["condition"] = condition
        resolution_rows.append(resolution)

        if resolution["status"] != "ok":
            continue

        covars = build_covariate_frame(obs, args, warnings_rows, db_name)
        features = build_expression_features(expr, obs, genes, args.include_clock_signature, condition, cohort, warnings_rows)

        data_base = pd.DataFrame(index=obs.index)
        data_base["time"] = time
        data_base["event"] = event
        data_base = pd.concat([data_base, covars, features], axis=1)

        for feature in features.columns:
            for model_name, spec_tokens in model_specs(args.include_model_c_sex).items():
                covar_names = expand_stage_token(spec_tokens, covars)
                rec, ph_df = fit_cox_model(data_base, feature, covar_names, model_name, args)
                rec.update({
                    "condition": condition,
                    "technical_sensitivity_flag": "primary_log2_expression" if condition == "log2" else "MYCN_corrected_expression_technical_sensitivity_only",
                    "cohort": cohort,
                    "db_name": db_name,
                    "feature": feature,
                    "gene_set": "12_clock_genes" if feature != "CLOCK_SIGNATURE" else "clock_signature_secondary",
                    "model": model_name,
                    "resolved_time_variable": time_col,
                    "resolved_event_variable": event_col,
                    "event_coding_rule": "1=death/event,0=censored/alive",
                    "time_unit_label": args.time_unit_label,
                    "HR_interpretation": "per_1_SD_expression_increase",
                })
                detail_rows.append(rec)

                if not ph_df.empty:
                    ph_df = ph_df.copy()
                    ph_df.insert(0, "condition", condition)
                    ph_df.insert(1, "cohort", cohort)
                    ph_df.insert(2, "db_name", db_name)
                    ph_df.insert(3, "feature", feature)
                    ph_df.insert(4, "model", model_name)
                    ph_df["technical_sensitivity_flag"] = rec["technical_sensitivity_flag"]
                    ph_rows.append(ph_df)

    detail = pd.DataFrame(detail_rows)
    ph = pd.concat(ph_rows, ignore_index=True, sort=False) if ph_rows else pd.DataFrame()
    resolution_df = pd.DataFrame(resolution_rows)
    warnings_df = pd.DataFrame(warnings_rows)

    if not detail.empty:
        detail["BH_q_within_condition_cohort_model_gene_set"] = np.nan
        for key, idx in detail.groupby(["condition", "cohort", "model", "gene_set"], dropna=False).groups.items():
            detail.loc[idx, "BH_q_within_condition_cohort_model_gene_set"] = bh_adjust(detail.loc[idx, "pvalue"])
        detail["p_lt_0_05"] = pd.to_numeric(detail["pvalue"], errors="coerce").map(lambda x: "yes" if pd.notna(x) and x < args.alpha else ("no" if pd.notna(x) else "NA"))
        detail["q_lt_0_05"] = pd.to_numeric(detail["BH_q_within_condition_cohort_model_gene_set"], errors="coerce").map(lambda x: "yes" if pd.notna(x) and x < args.alpha else ("no" if pd.notna(x) else "NA"))
        detail["PH_gene_term_p_lt_0_05_possible_violation"] = pd.to_numeric(detail["PH_gene_term_p"], errors="coerce").map(lambda x: "yes" if pd.notna(x) and x < args.alpha else ("no" if pd.notna(x) else "NA"))
        detail["PH_global_min_p_lt_0_05_possible_violation"] = pd.to_numeric(detail["PH_global_min_p"], errors="coerce").map(lambda x: "yes" if pd.notna(x) and x < args.alpha else ("no" if pd.notna(x) else "NA"))

    if not ph.empty and "PH_p" in ph.columns:
        ph["PH_p_lt_0_05_possible_violation"] = pd.to_numeric(ph["PH_p"], errors="coerce").map(lambda x: "yes" if pd.notna(x) and x < args.alpha else ("no" if pd.notna(x) else "NA"))

    summary = summarize_detail(detail, args.alpha)
    meta = fixed_effect_meta(detail, args.alpha)
    return detail, summary, ph, meta, resolution_df, warnings_df


# Collapse detailed gene-model results into cohort/model-level quality and significance summaries.
def summarize_detail(detail: pd.DataFrame, alpha: float) -> pd.DataFrame:
    """Summarize cohort-level Cox fitting and diagnostic results.
    
    Parameters
    ----------
    detail : pandas.DataFrame
        Detailed gene/model-level Cox results.
    alpha : float
        Significance threshold used for p/q and PH flags.
    
    Returns
    -------
    pandas.DataFrame
        Aggregated counts of successful fits, significant associations, PH flags,
        sample/event counts, EPV, and missingness by condition/cohort/model/gene set.
    """
    if detail.empty:
        return pd.DataFrame()
    d = detail.copy()
    d["fit_ok_bool"] = d["fit_status"].eq("fit_ok")
    d["p_sig"] = pd.to_numeric(d["pvalue"], errors="coerce") < alpha
    d["q_sig"] = pd.to_numeric(d["BH_q_within_condition_cohort_model_gene_set"], errors="coerce") < alpha
    d["ph_gene_bad"] = pd.to_numeric(d["PH_gene_term_p"], errors="coerce") < alpha
    d["ph_global_bad"] = pd.to_numeric(d["PH_global_min_p"], errors="coerce") < alpha
    return d.groupby(["condition", "cohort", "model", "gene_set", "technical_sensitivity_flag"], dropna=False).agg(
        n_gene_model_estimates=("feature", "count"),
        n_fit_ok=("fit_ok_bool", "sum"),
        n_p_lt_0_05=("p_sig", "sum"),
        n_q_lt_0_05=("q_sig", "sum"),
        n_gene_PH_possible_violations=("ph_gene_bad", "sum"),
        n_global_PH_possible_violations=("ph_global_bad", "sum"),
        median_n=("n", "median"),
        min_n=("n", "min"),
        median_n_events=("n_events", "median"),
        min_n_events=("n_events", "min"),
        median_events_per_variable=("events_per_variable", "median"),
        median_missing_fraction=("missing_fraction_for_model", "median"),
    ).reset_index()


# Combine fitted SEQC/Kocak gene coefficients by inverse-variance fixed-effect meta-analysis.
def fixed_effect_meta(detail: pd.DataFrame, alpha: float) -> pd.DataFrame:
    """Perform inverse-variance fixed-effect meta-analysis across cohorts.
    
    Parameters
    ----------
    detail : pandas.DataFrame
        Detailed cohort-level Cox results containing fitted beta coefficients and
        standard errors.
    alpha : float
        Significance threshold used for meta-analysis q-value flags.
    
    Returns
    -------
    pandas.DataFrame
        Cross-cohort meta-analysis estimates, confidence intervals, p values, BH q
        values, and cohort-availability indicators for each condition/feature/model.
    """
    if detail.empty:
        return pd.DataFrame()
    d = detail.loc[detail["fit_status"].eq("fit_ok") & detail["cohort"].isin(["SEQC", "Kocak"])].copy()
    d["beta_num"] = pd.to_numeric(d["beta"], errors="coerce")
    d["se_num"] = pd.to_numeric(d["se"], errors="coerce")
    d = d.loc[d["beta_num"].notna() & d["se_num"].notna() & (d["se_num"] > 0)].copy()
    if d.empty:
        return pd.DataFrame()

    records = []
    for key, sub in d.groupby(["condition", "feature", "model", "gene_set"], dropna=False):
        cohorts = sorted(set(sub["cohort"]))
        weights = 1.0 / (sub["se_num"].to_numpy(dtype=float) ** 2)
        betas = sub["beta_num"].to_numpy(dtype=float)
        meta_beta = float(np.sum(weights * betas) / np.sum(weights))
        meta_se = float(math.sqrt(1.0 / np.sum(weights)))
        z = meta_beta / meta_se if meta_se > 0 else np.nan
        p = normal_two_sided_p(z) if np.isfinite(z) else np.nan
        records.append({
            "condition": key[0],
            "technical_sensitivity_flag": "primary_log2_expression" if key[0] == "log2" else "MYCN_corrected_expression_technical_sensitivity_only",
            "feature": key[1],
            "model": key[2],
            "gene_set": key[3],
            "n_cohorts_meta": len(cohorts),
            "cohorts_included": ",".join(cohorts),
            "both_primary_cohorts_available": "yes" if {"SEQC", "Kocak"}.issubset(set(cohorts)) else "no",
            "meta_log_HR": meta_beta,
            "meta_se": meta_se,
            "meta_HR": math.exp(meta_beta),
            "meta_CI95_low": math.exp(meta_beta - 1.96 * meta_se),
            "meta_CI95_high": math.exp(meta_beta + 1.96 * meta_se),
            "meta_z": z,
            "meta_p": p,
        })
    meta = pd.DataFrame(records)
    if meta.empty:
        return meta
    meta["meta_BH_q_within_condition_model_gene_set"] = np.nan
    for key, idx in meta.groupby(["condition", "model", "gene_set"], dropna=False).groups.items():
        meta.loc[idx, "meta_BH_q_within_condition_model_gene_set"] = bh_adjust(meta.loc[idx, "meta_p"])
    meta["meta_q_lt_0_05"] = pd.to_numeric(meta["meta_BH_q_within_condition_model_gene_set"], errors="coerce").map(lambda x: "yes" if pd.notna(x) and x < alpha else ("no" if pd.notna(x) else "NA"))
    return meta


# -----------------------------
# Output
# -----------------------------

# Append an extension without stripping dotted components from the output prefix.
def output_path_with_added_extension(base: Path, extension: str) -> Path:
    """Append an extension without stripping dotted components of an output prefix.
    
    Parameters
    ----------
    base : pathlib.Path
        Base output path, which may itself contain periods.
    extension : str
        Extension to append, including the leading period (for example ``.tsv``).
    
    Returns
    -------
    pathlib.Path
        Output path formed as ``<base.name><extension>`` in the same directory.
    """
    base = Path(base)
    return base.parent / f"{base.name}{extension}"


# Write paired TSV/XLSX outputs while preserving complete dotted prefixes.
def write_table(df: pd.DataFrame, base: Path, write_xlsx: bool = True) -> None:
    """Write a dataframe to TSV and optionally XLSX.
    
    Parameters
    ----------
    df : pandas.DataFrame
        Table to write.
    base : pathlib.Path
        Base output path without the final ``.tsv``/``.xlsx`` extension.
    write_xlsx : bool, default=True
        Whether to also create a formatted Excel workbook.
    
    Returns
    -------
    None
        Files are written to disk as a side effect.
    """
    tsv_path = output_path_with_added_extension(base, ".tsv")
    xlsx_path = output_path_with_added_extension(base, ".xlsx")
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(tsv_path, sep="\t", index=False)
    if write_xlsx:
        try:
            with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
                df.to_excel(writer, index=False, sheet_name="table")
                ws = writer.book["table"]
                ws.freeze_panes = "A2"
                if ws.max_row >= 1 and ws.max_column >= 1:
                    ws.auto_filter.ref = ws.dimensions
                for cells in ws.columns:
                    header = str(cells[0].value or "")
                    ws.column_dimensions[cells[0].column_letter].width = min(max(len(header) + 2, 12), 55)
        except Exception as e:
            sys.stderr.write(f"WARNING: failed to write XLSX {xlsx_path}: {e}\n")


# Write manuscript-ready methodological footnote text describing the Cox analysis.
def write_footnote(path: Path) -> None:
    """Write the manuscript-oriented Cox model footnote.
    
    Parameters
    ----------
    path : pathlib.Path
        Destination text-file path.
    
    Returns
    -------
    None
        The explanatory footnote is written to ``path``.
    """
    text = (
        "Adjusted Cox model footnote. Cox proportional-hazards models were fitted separately for SEQC/GSE62564 "
        "and Kocak/GSE45547 using sample-level survival and clinical annotations. Survival variables were resolved "
        "before fitting and are documented in the survival-variable resolution table. Events were coded as 1 for "
        "death/event and 0 for censored/alive, and survival time was required to be numeric and positive. Clock-gene "
        "expression was standardised within cohort; hazard ratios therefore denote the change in hazard per one "
        "standard deviation increase in expression. Univariable models included expression only. Adjusted model A "
        "included expression, age group (>=18 months), MYCN amplification status and INSS stage covariates. Adjusted "
        "model B included expression, harmonised high-risk status and age group. Adjusted model C additionally included "
        "sex when available and event counts permitted. Covariates with zero variance, possible complete separation or "
        "high collinearity were dropped and reported. Models were not fitted when event counts or event-per-variable "
        "ratios were insufficient. Benjamini-Hochberg q values were computed within each condition, cohort, model and "
        "gene-set family. Proportional-hazards assumptions were assessed using rank-transformed Schoenfeld residual "
        "tests; both the gene-term p value and the minimum p value across model terms are reported. Fixed-effect "
        "meta-analysis used inverse-variance weighting of SEQC and Kocak gene coefficients when both cohort-specific "
        "models were fit. MYCN-corrected/ComBat expression analyses are technical sensitivity analyses and do not "
        "replace the primary log2 expression models."
    )
    path.write_text(text + "\n", encoding="utf-8")


# Write reusable interpretation guidance for significant, attenuated, and caveated Cox results.
def write_interpretation(path: Path) -> None:
    """Write suggested manuscript interpretation and limitation text.
    
    Parameters
    ----------
    path : pathlib.Path
        Destination Markdown-file path.
    
    Returns
    -------
    None
        Interpretation guidance is written to ``path``.
    """
    text = """# Suggested manuscript interpretation

## If adjusted associations remain significant

In adjusted Cox analyses, selected circadian genes retained associations with overall survival after accounting for age group, MYCN status and INSS stage, with hazard ratios reported per one standard deviation increase in expression. These associations should only be described as adjusted support when the corresponding model has `fit_status=fit_ok`, BH q<0.05, and the gene-term proportional-hazards diagnostic does not indicate violation. If the fixed-effect meta-analysis across SEQC and Kocak is also significant and directionally consistent, this may be described as cross-cohort supportive evidence.

## If adjusted associations attenuate

If univariable associations attenuate after adjustment, the survival analysis should be presented as contextual rather than independent prognostic evidence. This would indicate that single-gene clock survival associations are largely captured by established clinical variables such as age, MYCN status, stage and harmonised risk group. The main manuscript inference should remain the Delta-CCD network-level coordination result rather than single-gene prognostic biomarker discovery.

## Limitation wording

Adjusted Cox analyses are limited by cohort-specific annotation differences, complete-case attrition and collinearity among high-risk status, MYCN amplification and stage. Models flagged for low event counts, low event-per-variable ratio, dropped covariates, convergence failure or possible proportional-hazards violation should not be used for claims of independent prognostic value.
"""
    path.write_text(text, encoding="utf-8")


# Main entry point: validate dependencies, run the workflow, and write all outputs.
def main() -> None:
    """Run the command-line adjusted Cox analysis and write all outputs.
    
    Parameters
    ----------
    None
        Configuration is obtained from command-line arguments.
    
    Returns
    -------
    None
        The function writes detail, summary, PH diagnostic, meta-analysis,
        survival-variable resolution, warning, footnote, and interpretation files.
    
    Raises
    ------
    RuntimeError
        If the required ``lifelines`` dependency is unavailable.
    """
    args = parse_args()

    try:
        import lifelines  # noqa
    except Exception as e:
        raise RuntimeError("This script requires lifelines. Install with: pip install lifelines openpyxl") from e

    detail, summary, ph, meta, resolution, warnings_df = run_all(args)

    out = Path(args.out_prefix)
    write_table(detail, out.with_name(out.name + "_clock_genes_detail"), args.write_xlsx)
    write_table(summary, out.with_name(out.name + "_clock_genes_summary"), args.write_xlsx)
    write_table(ph, out.with_name(out.name + "_PH_diagnostics"), args.write_xlsx)
    write_table(meta, out.with_name(out.name + "_meta_analysis"), args.write_xlsx)
    resolution.to_csv(out.with_name(out.name + "_survival_variable_resolution.tsv"), sep="\t", index=False)
    if warnings_df.empty:
        warnings_df = pd.DataFrame([{"severity": "INFO", "topic": "none", "message": "No warnings generated."}])
    warnings_df.to_csv(out.with_name(out.name + "_warnings.tsv"), sep="\t", index=False)
    write_footnote(out.with_name(out.name + "_footnote.txt"))
    write_interpretation(out.with_name(out.name + "_interpretation.md"))

    print(f"Wrote adjusted Cox outputs with prefix: {out}")
    print("Expected tables preserve dotted prefixes, e.g. <prefix>_clock_genes_detail.tsv")


if __name__ == "__main__":
    main()
