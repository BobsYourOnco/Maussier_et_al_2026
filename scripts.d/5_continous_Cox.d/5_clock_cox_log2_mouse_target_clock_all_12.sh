#!/bin/bash -l


FORCE="0"
COX_SCRIPT="adjusted_cox_clock_genes.py"
MANIFEST="manifests.d/target_clock_all_12.tsv"
CONDITION="log2"
OUT_PREFIX="results.d/5_continous_Cox.d/mouse/log2/target_clock_all_12.adjusted_cox"
MAIN_OUTFILE="results.d/5_continous_Cox.d/mouse/log2/target_clock_all_12.adjusted_cox_clock_genes_detail.tsv"
CLOCK_GENES="ARNTL,CLOCK,CRY1,CRY2,DBP,NPAS2,NR1D1,NR1D2,PER1,PER2,PER3,TEF"
MIN_SAMPLES="30"
MIN_EVENTS="10"
MIN_EPV="5.0"
COLLINEARITY_R="0.95"
MAX_DROP_FRAC="0.3"
COX_PENALIZER="0.0"
ALPHA="0.05"
INCLUDE_MODEL_C_SEX="1"
INCLUDE_CLOCK_SIGNATURE="1"
TIME_COL=""
EVENT_COL=""
AGE18_COL=""
AGE_COL=""
AGE_UNIT="auto"
MYCN_COL=""
STAGE_COL=""
RISK_COL=""
SEX_COL=""
TIME_UNIT_LABEL="as_provided"
mkdir -p "$(dirname "$OUT_PREFIX")"
if [[ "$FORCE" != "1" && -s "$MAIN_OUTFILE" ]]; then
  echo "Skipping fixed adjusted Cox; output already exists: $MAIN_OUTFILE"
else
  cmd=(
    python "$COX_SCRIPT"
    --manifest "${CONDITION}:${MANIFEST}"
    --out_prefix "$OUT_PREFIX"
    --clock_genes "$CLOCK_GENES"
    --min_samples "$MIN_SAMPLES"
    --min_events "$MIN_EVENTS"
    --min_events_per_variable "$MIN_EPV"
    --collinearity_r "$COLLINEARITY_R"
    --max_drop_fraction_warning "$MAX_DROP_FRAC"
    --cox_penalizer "$COX_PENALIZER"
    --alpha "$ALPHA"
    --age_unit "$AGE_UNIT"
    --time_unit_label "$TIME_UNIT_LABEL"
  )
  if [[ "$INCLUDE_MODEL_C_SEX" == "1" ]]; then
    cmd+=(--include_model_c_sex)
  fi
  if [[ "$INCLUDE_CLOCK_SIGNATURE" == "1" ]]; then
    cmd+=(--include_clock_signature)
  fi
  if [[ -n "${TIME_COL}" ]]; then
    cmd+=(--time_col "${TIME_COL}")
  fi
  if [[ -n "${EVENT_COL}" ]]; then
    cmd+=(--event_col "${EVENT_COL}")
  fi
  if [[ -n "${AGE18_COL}" ]]; then
    cmd+=(--age18_col "${AGE18_COL}")
  fi
  if [[ -n "${AGE_COL}" ]]; then
    cmd+=(--age_col "${AGE_COL}")
  fi
  if [[ -n "${MYCN_COL}" ]]; then
    cmd+=(--mycn_col "${MYCN_COL}")
  fi
  if [[ -n "${STAGE_COL}" ]]; then
    cmd+=(--stage_col "${STAGE_COL}")
  fi
  if [[ -n "${RISK_COL}" ]]; then
    cmd+=(--risk_col "${RISK_COL}")
  fi
  if [[ -n "${SEX_COL}" ]]; then
    cmd+=(--sex_col "${SEX_COL}")
  fi
  echo "Running: ${cmd[*]}"
  "${cmd[@]}"
fi

