#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Leave-one-gene-out wrapper for direct delta CCD analyses.

This script repeatedly calls NCC_tool_v2c_deltaCCD_matchedNull.py in --compute_delta_CCD mode,
removing one gene at a time from the input gene panel. It is intended for fast robustness panels
such as MYCN-amplified vs non-amplified or high-risk vs low-risk contrasts.
"""

import argparse
import os
import sys
import shlex
import subprocess
from pathlib import Path
from typing import List, Dict

import pandas as pd


def parse_args() -> argparse.Namespace:
	"""
	Parse arguments
	"""
	p = argparse.ArgumentParser(description="Run leave-one-gene-out direct delta CCD robustness analyses")
	p.add_argument("--ncc_script", default="/mnt/data/NCC_tool_v2c_deltaCCD_matchedNull.py",
				   help="Path to the integrated NCC script with direct delta CCD support")
	# Core NCC arguments
	p.add_argument("--input_reference_db_file", required=True,
				   help="Path to the reference database config TSV")
	p.add_argument("--reference_weighted_rho_matrix", required=True,
				   help="Path to the weighted reference rho matrix CSV")
	p.add_argument("--reference_full_rho_info", default="reference_databases_full_rho_info.tsv",
				   help="Optional reference full rho output, only used if --generate_reference_correlations is passed")
	p.add_argument("--l_clock_genes_reference", required=True,
				   help="Gene panel joined by '|' characters")
	p.add_argument("--reference_species", required=True, help="Reference species")
	p.add_argument("--input_target_db_file", required=True,
				   help="Path to the target database config TSV")
	p.add_argument("--reference_target_gene_file", required=True,
				   help="Path to the orthologue mapping table")
	p.add_argument("--target_species", required=True, help="Target species")
	# Direct delta CCD arguments
	p.add_argument("--contrast_variable", required=True,
				   help="Metadata variable used for the direct delta CCD contrast")
	p.add_argument("--group_a", required=True, help="First group label")
	p.add_argument("--group_b", required=True, help="Second group label")
	p.add_argument("--delta_n_permutations", type=int, default=5000,
				   help="Number of label permutations per run")
	p.add_argument("--delta_n_bootstraps", type=int, default=2000,
				   help="Number of bootstrap resamples per run")
	p.add_argument("--delta_seed", type=int, default=42,
				   help="Base seed; each leave-one-gene-out run gets a deterministic offset")
	p.add_argument("--delta_min_samples_per_group", type=int, default=10,
				   help="Minimum samples required in each contrasted group")
	# Wrapper behavior
	p.add_argument("--include_full_panel", action="store_true",
				   help="Also run the full gene panel before leave-one-gene-out runs")
	p.add_argument("--output_dir", required=True,
				   help="Directory for per-run outputs and the aggregated summary")
	p.add_argument("--output_prefix", default="delta_ccd_logo",
				   help="Prefix for generated files")
	p.add_argument("--save_draws", action="store_true",
				   help="Also save permutation/bootstrap draw tables for every run")
	p.add_argument("--python_bin", default=sys.executable,
				   help="Python executable used to call the NCC script")
	# Optional passthroughs to the NCC script
	p.add_argument("--generate_reference_correlations", action="store_true",
				   help="Pass through to NCC script when a fresh reference matrix should be generated")
	p.add_argument("--run_metanalysis_target", action="store_true",
				   help="Pass through to NCC script")

	return p.parse_args()


def sanitize_label(label: str) -> str:
	 """
	Convert an arbitrary label into a filesystem-safe string.
	Characters that are alphanumeric or commonly safe in filenames
	("-", "_", ".") are retained unchanged. All other characters,
	including spaces and punctuation, are replaced with underscores.
	This is useful for constructing output filenames or directory names
	from contrast labels, gene names, reference names, or group labels.
	Example
	-------
	"high risk vs low risk" -> "high_risk_vs_low_risk"
	"MYCN amplified/not amplified" -> "MYCN_amplified_not_amplified"
	"""
	# Store the sanitized characters one at a time.
	safe = []
	for ch in str(label):
		if ch.isalnum() or ch in {"-", "_", "."}:
			safe.append(ch)
		else:
			safe.append("_")
	return "".join(safe)


def build_gene_panels(genes_string: str, include_full_panel: bool) -> List[Dict[str, str]]:
	"""
	Build the full and leave-one-gene-out clock-gene panels.
	The input gene list is expected as a pipe-separated string, for example:
	"Bmal1|Clock|Cry1|Cry2|..."
	If include_full_panel is True, the original complete gene set is included
	first. The function then creates one additional panel for each gene, with
	that gene removed.
	Each returned panel contains:
		panel_name
			Name of the panel, e.g. "full_panel" or "drop_Bmal1".
		dropped_gene
			Gene omitted from the panel. Empty for the full panel.
		genes
			List of genes retained in the panel.
		genes_string
			Pipe-separated version of the retained genes, suitable for passing
			back to the NCC/CCD script.
	Example
	-------
	Input:
		genes_string = "A|B|C"
		include_full_panel = True
	Output panels:
		full_panel -> A|B|C
		drop_A	 -> B|C
		drop_B	 -> A|C
		drop_C	 -> A|B
	"""
	# Split the pipe-separated gene string into an ordered list and ignore empty entries.
	genes = [g for g in genes_string.split("|") if g]
	# Leave-one-gene-out analysis requires a sufficiently large starting panel.
	# With fewer than 3 genes, dropping one gene would leave too little
	# structure for a meaningful correlation-based analysis.
	if len(genes) < 3:
		raise ValueError("Leave-one-gene-out needs at least 3 genes in --l_clock_genes_reference")
	# Store all gene-panel definitions here.
	panels: List[Dict[str, str]] = []
	# Optionally retain the original complete gene panel as a reference result.
	if include_full_panel:
		panels.append({
			"panel_name": "full_panel",
			"dropped_gene": "",
			"genes": genes[:],
			"genes_string": "|".join(genes),
		})
	for gene in genes:
		sub = [g for g in genes if g != gene]
		panels.append({
			"panel_name": f"drop_{gene}",
			"dropped_gene": gene,
			"genes": sub,
			"genes_string": "|".join(sub),
		})
	return panels


def run_one_panel(args: argparse.Namespace, panel: Dict[str, str], run_index: int, output_dir: Path) -> pd.DataFrame:
	"""
	Run the NCC/Delta-CCD analysis for one gene panel.
	A panel can be either:
		- the complete clock-gene panel; or
		- a leave-one-gene-out panel with one gene removed.
	The function:
		1. creates panel-specific output filenames;
		2. builds the command used to call the NCC script;
		3. assigns a reproducible panel-specific random seed;
		4. runs the NCC analysis as a subprocess;
		5. checks that the run completed successfully;
		6. reads the resulting Delta-CCD table;
		7. adds metadata identifying the panel and omitted gene;
		8. returns the annotated result as a pandas DataFrame.
	Parameters
	----------
	args
		Parsed command-line arguments containing NCC script paths,
		reference/target settings, contrast definitions, resampling
		parameters, and output options.
	panel
		Dictionary describing the current gene panel. Expected entries are:
			panel_name
			dropped_gene
			genes
			genes_string
	run_index
		Index of the current panel run. This is added to the base random seed
		so that each panel receives a different but reproducible seed.
	output_dir
		Directory where the panel-specific result and draw files are written.

	Returns
	-------
	pd.DataFrame
		Delta-CCD results for this panel, annotated with the panel name,
		omitted gene, number of retained genes, and retained gene list.
	"""
	# Convert the panel name into a filename-safe label.
	panel_tag = sanitize_label(panel["panel_name"])
	# Define panel-specific output files.
	results_path = output_dir / f"{args.output_prefix}.{panel_tag}.results.tsv"
	draws_path = output_dir / f"{args.output_prefix}.{panel_tag}.draws.tsv"
	# Build the command that will call the underlying NCC/Delta-CCD script.
	cmd = [
		args.python_bin,
		args.ncc_script,
		"--input_reference_db_file", args.input_reference_db_file,
		"--reference_weighted_rho_matrix", args.reference_weighted_rho_matrix,
		"--reference_full_rho_info", args.reference_full_rho_info,
		"--l_clock_genes_reference", panel["genes_string"],
		"--reference_species", args.reference_species,
		"--input_target_db_file", args.input_target_db_file,
		"--reference_target_gene_file", args.reference_target_gene_file,
		"--target_species", args.target_species,
		"--compute_delta_CCD",
		"--contrast_variable", args.contrast_variable,
		"--group_a", args.group_a,
		"--group_b", args.group_b,
		"--delta_n_permutations", str(args.delta_n_permutations),
		"--delta_n_bootstraps", str(args.delta_n_bootstraps),
		"--delta_seed", str(args.delta_seed + run_index),
		"--delta_min_samples_per_group", str(args.delta_min_samples_per_group),
		"--delta_outfile_results", str(results_path),
	]
	# Optionally save permutation/bootstrap draws.
	if args.save_draws:
		cmd.extend(["--delta_outfile_draws", str(draws_path)])
	if args.generate_reference_correlations:
		cmd.append("--generate_reference_correlations")
	if args.run_metanalysis_target:
		cmd.append("--run_metanalysis_target")
	# Print the exact command in a shell-readable form.
	# shlex.quote() safely quotes paths or arguments containing spaces
	# or special shell characters.
	print("Running:", " ".join(shlex.quote(x) for x in cmd))
	# Run the NCC command and capture both standard output and standard error.
	completed = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
	if completed.returncode != 0:
		raise RuntimeError(
			f"NCC run failed for panel {panel['panel_name']}\n"
			f"STDOUT:\n{completed.stdout}\n\nSTDERR:\n{completed.stderr}"
		)
	if not results_path.exists():
		raise FileNotFoundError(f"Expected results file was not created: {results_path}")
	# Add LOO metadata so results from all panels can later be combined
	# into a single summary table.
	df = pd.read_csv(results_path, sep="\t")
	df.insert(0, "panel_name", panel["panel_name"])
	df.insert(1, "dropped_gene", panel["dropped_gene"] if panel["dropped_gene"] else "<none>")
	df.insert(2, "n_genes_in_panel", len(panel["genes"]))
	df.insert(3, "genes_in_panel", "|".join(panel["genes"]))
	return df


def main() -> None:
	"""
	Run
	"""
	args = parse_args()
	output_dir = Path(args.output_dir)
	output_dir.mkdir(parents=True, exist_ok=True)
	#
	if not Path(args.ncc_script).exists():
		raise FileNotFoundError(f"NCC script not found: {args.ncc_script}")
	#
	panels = build_gene_panels(args.l_clock_genes_reference, args.include_full_panel)
	all_results = []
	for idx, panel in enumerate(panels):
		df = run_one_panel(args, panel, idx, output_dir)
		all_results.append(df)
	#
	agg = pd.concat(all_results, axis=0, ignore_index=True)
	#
	# Helpful stability summaries against the full-panel result when present
	if args.include_full_panel and (agg["panel_name"] == "full_panel").any():
		full = agg[agg["panel_name"] == "full_panel"].copy()
		ref_cols = [c for c in [
			"database_name", "delta_ccd", "p_perm_two_sided", "p_perm_one_sided_a_gt_b",
			"delta_ccd_boot_ci_low", "delta_ccd_boot_ci_high"
		] if c in full.columns]
		full = full[["database_name"] + [c for c in ref_cols if c != "database_name"]].drop_duplicates("database_name")
		full = full.rename(columns={
			"delta_ccd": "full_panel_delta_ccd",
			"p_perm_two_sided": "full_panel_p_perm_two_sided",
			"p_perm_one_sided_a_gt_b": "full_panel_p_perm_one_sided_a_gt_b",
			"delta_ccd_boot_ci_low": "full_panel_ci_low",
			"delta_ccd_boot_ci_high": "full_panel_ci_high",
		})
		agg = agg.merge(full, on="database_name", how="left")
		if "delta_ccd" in agg.columns and "full_panel_delta_ccd" in agg.columns:
			agg["delta_ccd_shift_vs_full"] = agg["delta_ccd"] - agg["full_panel_delta_ccd"]
	summary_out = output_dir / f"{args.output_prefix}.summary.tsv"
	agg.to_csv(summary_out, sep="\t", index=False)
	# Concise per-dropped-gene meta summary when meta results exist
	meta_cols = [c for c in [
		"panel_name", "dropped_gene", "meta_delta_ccd", "meta_ci_low", "meta_ci_high", "meta_p_two_sided"
	] if c in agg.columns]
	if set(["panel_name", "dropped_gene"]).issubset(agg.columns) and len(meta_cols) >= 3:
		meta = agg[meta_cols].drop_duplicates().sort_values(["panel_name", "dropped_gene"])
		meta_out = output_dir / f"{args.output_prefix}.meta_summary.tsv"
		meta.to_csv(meta_out, sep="\t", index=False)
	#
	print(f"Wrote aggregated leave-one-gene-out summary: {summary_out}")
	if 'meta_out' in locals():
		print(f"Wrote meta summary: {meta_out}")


if __name__ == "__main__":
	main()
