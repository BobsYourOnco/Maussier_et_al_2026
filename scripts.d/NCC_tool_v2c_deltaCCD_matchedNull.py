#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  NCC Tool.py -Nina's Circadian Clockwork-
#  
#  Copyright 2025 Oscar C: Bedoya Reina & Nina Mausser
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

"""
NCC Tool -Nina's Circadian Clockwork- is a tool to recreate each step in Shilts et al. 2018
paper methods. It also aims at adding new approaches for hypothesis testing.
"""

import pandas as pd
import scanpy as sc
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

import argparse
import math
import os
import warnings
import itertools
import sys

from statsmodels.stats.multitest import multipletests

###############
## Parse arguments ##
###############
def parse_args():
    """
    Method to parse input arguments for the NCC command line tool.
    Input: command line arguments provided by the user when running the script.
    Output: args is an argparse namespace containing the input files, species,
    randomization parameters, switches, and options used throughout the NCC workflow.
    """
    parser = argparse.ArgumentParser(description="Run NCC Tool -Nina's Circadian Clockwork- ")
    # Reference database
    parser.add_argument('--input_reference_db_file', required=True,help='Path to reference table with fields required for the reference')
    parser.add_argument('--reference_full_rho_info', default='reference_databases_full_rho_info.tsv',help='Ouput path to full information .on correlations for databases in "input_reference_db_file"')
    parser.add_argument('--reference_weighted_rho_matrix', default='reference_weighted_rho.csv',help='Ouput path to write weighted correlations for the databases in "input_reference_db_file"')
    parser.add_argument('--l_clock_genes_reference', required=True,default='Bmal1|Clock|Cry1|Cry2|Dbp|Npas2|Nr1d1|Nr1d2|Per1|Per2|Per3|Tef', help='List of genes of interest to compute the correlations, joined by "|"')
    parser.add_argument('--reference_species', required=True,help='Input reference species')
    # Target database
    parser.add_argument('--input_target_db_file', required=True,help="Path to reference table with fields required for the target")
    parser.add_argument('--target_full_rho_info', default='target_databases_full_rho_info.csv',help='Output path to full information on correlations for databases in "input_target_db_file"')
    parser.add_argument('--reference_target_gene_file', required=True,help="File with 1:1 target to reference genes orthologues, one per line -tabular-")
    parser.add_argument('--target_species', required=True,help='Input taget species')
    parser.add_argument('--run_metanalysis_target', action='store_true',help='If pass will run metanalysis in target dataset, and also in the background dataset')
    # Randomization
    parser.add_argument('--number_randomizations', type=int, default=1000,help='Number of randomizations to compute p-values')
    parser.add_argument('--background_target_full_rho_info', default='background_target_full_rho_info.tsv',help='Output file for CCD of randomizations')
    parser.add_argument('--outfile_results', default='CCD_results.tsv',help='Output file for results')
    # Switches
    parser.add_argument('--generate_reference_correlations', action='store_true',help='Pass to generate reference correlations' )
    parser.add_argument('--generate_target_correlations', action='store_true',help='Pass to generate target correlations')
    parser.add_argument('--generate_background_target_correlations', action='store_true',help='Generate background target correlations')
    parser.add_argument('--compute_CCD', action='store_true',help='Pass to compute CCD')
    # Direct between-group delta CCD workflow
    parser.add_argument('--compute_delta_CCD', action='store_true', help='Pass to compute direct between-group delta CCD contrasts with permutation testing and bootstrap confidence intervals')
    parser.add_argument('--contrast_variable', default=None, help='Metadata variable used for direct delta CCD contrasts (for example MYCN_status or risk_group)')
    parser.add_argument('--group_a', default=None, help='First group label for delta CCD contrasts')
    parser.add_argument('--group_b', default=None, help='Second group label for delta CCD contrasts')
    parser.add_argument('--delta_n_permutations', type=int, default=5000, help='Number of label permutations for direct delta CCD contrasts')
    parser.add_argument('--delta_n_bootstraps', type=int, default=2000, help='Number of bootstrap resamples for direct delta CCD confidence intervals')
    parser.add_argument('--delta_seed', type=int, default=42, help='Seed for direct delta CCD permutations and bootstraps')
    parser.add_argument('--delta_min_samples_per_group', type=int, default=10, help='Minimum number of samples required in each contrasted group')
    parser.add_argument('--delta_outfile_results', default='delta_CCD_results.tsv', help='Output TSV for direct delta CCD results')
    parser.add_argument('--delta_outfile_draws', default=None, help='Optional output TSV for direct delta CCD permutation and bootstrap draws')
    # Matched-null robustness workflow for background CCDs
    parser.add_argument('--background_matching_mode', default='unmatched', choices=['unmatched','matched'], help='How to build background random gene sets for CCD significance. "matched" samples gene sets matched to the target gene set on aggregated expression and variance ranks across cohorts')
    parser.add_argument('--matching_mean_caliper', type=float, default=0.05, help='Initial caliper on aggregated mean-expression rank percentiles for matched null gene sampling')
    parser.add_argument('--matching_variance_caliper', type=float, default=0.05, help='Initial caliper on aggregated variance rank percentiles for matched null gene sampling')
    parser.add_argument('--matching_expand_step', type=float, default=0.05, help='Step used to expand expression/variance matching calipers when no matched null candidate is available')
    parser.add_argument('--matching_max_caliper', type=float, default=1.00, help='Maximum expression/variance caliper allowed before falling back to the nearest remaining gene in rank space')
    parser.add_argument('--matching_seed', type=int, default=42, help='Seed for matched null gene sampling')
    parser.add_argument('--matched_gene_summary_outfile', default=None, help='Optional output TSV with per-gene aggregated expression/variance matching summaries used for matched null sampling')
    parser.add_argument('--matched_sets_outfile', default=None, help='Optional output TSV with the sampled matched null gene sets (one column per randomization)')
    #
    return parser.parse_args()
#
args = parse_args()

###################
###################
###  Reference database  ###
###################
###################
#Path to reference table with fields required for the reference
input_reference_db_file = args.input_reference_db_file#'reference.tsv'
#Ouput path to full information on correlations for databases in "input_reference_db_file"
reference_full_rho_info = args.reference_full_rho_info#'reference_databases_full_rho_info.tsv'
#Ouput path to write weighted correlations for the databases in "input_reference_db_file"
reference_weighted_rho_matrix = args.reference_weighted_rho_matrix#'reference_weighted_rho.csv'
#List of genes of interest to compute the correlations
l_clock_genes_reference = args.l_clock_genes_reference#'Bmal1|Clock|Cry1|Cry2|Dbp|Npas2|Nr1d1|Nr1d2|Per1|Per2|Per3|Tef'
#Input reference species
reference_species=args.reference_species#'Mouse'

#################
#################
###  Target database  ###
#################
#################
#Path to reference table with fields required for the target
input_target_db_file = args.input_target_db_file#'target.tsv'
#Output path to full information on correlations for databases in "input_target_db_file"
target_full_rho_info=args.target_full_rho_info#'target_databases_full_rho_info.csv'
#File with 1:1 target to reference genes orthologues, one per line -tabular-
reference_target_gene_file=args.reference_target_gene_file#'mouse_human_11orthologue_ENSEMBL_mBMAL1toARNTL.tsv'
#Input species
target_species=args.target_species#'Human'
#Run meta-analysis in target and bacgrkound databases
run_metanalysis_target=args.run_metanalysis_target#False

#########################
#########################
###  Parameters for randomization  ###
#########################
#########################
#Number of randomizations to compute p-values
number_randomizations = args.number_randomizations#1000
#Output file for CCD of randomizations
background_target_full_rho_info = args.background_target_full_rho_info#'background_target_full_rho_info.tsv'
#Output file for results
outfile_results=args.outfile_results#'CCD_results.tsv'

#############
#############
###  Switches  ###
#############
#############
# Generate reference correlations
generate_reference_correlations = args.generate_reference_correlations#False
# Generate target correlations
generate_target_correlations = args.generate_target_correlations#True
# Generate background target correlations
generate_background_target_correlations= args.generate_background_target_correlations#True
# Compute CCD
compute_CCD = args.compute_CCD#True
# Direct delta CCD contrasts
direct_delta_CCD = args.compute_delta_CCD
contrast_variable = args.contrast_variable
group_a = args.group_a
group_b = args.group_b
delta_n_permutations = args.delta_n_permutations
delta_n_bootstraps = args.delta_n_bootstraps
delta_seed = args.delta_seed
delta_min_samples_per_group = args.delta_min_samples_per_group
delta_outfile_results = args.delta_outfile_results
delta_outfile_draws = args.delta_outfile_draws
# Matched-null robustness settings
background_matching_mode = args.background_matching_mode
matching_mean_caliper = args.matching_mean_caliper
matching_variance_caliper = args.matching_variance_caliper
matching_expand_step = args.matching_expand_step
matching_max_caliper = args.matching_max_caliper
matching_seed = args.matching_seed
matched_gene_summary_outfile = args.matched_gene_summary_outfile
matched_sets_outfile = args.matched_sets_outfile

#####################
#####################
### Preliminar adjustments ###
#####################
#####################
#Make the input genes a list
l_clock_genes_reference = l_clock_genes_reference.split('|')

############
############
### Methods  ###
############
############
##Input file to h5ad
def csvToH5ad(input_exprssn_fl,l_excld_clmns=None,l_excld_rows=None,sep='\t'):
	"""
	Method to generate h5ad from csv files.
	Input: input_exprssn_fl, input expression file, expected to be in the format of R2 (i.e.
	first column==gene name, first raw==sample, and sample metadata in the header
	starting with #). The input file should have sample names in the first row and gene 
	names in the first column .l_excld_clmns is a list of column names to be excluded (if 
	for instance not all columns are samples), l_excld_rows if a list of row names to be 
	excluded (if for instance not all columns are samples). sep is the separator character 
	in the input file.
	Ouput: h5ad file with count in input_exprssn_fl and metadata in obs.
	"""
	######
	#Read file
	######
	pd_input_exprssn_fl = pd.read_csv(input_exprssn_fl,sep=sep,header=0,index_col=0, \
	low_memory=False)#
	if l_excld_clmns is not None:#Adjust for exclude columns
		pd_input_exprssn_fl = pd_input_exprssn_fl.T[~pd_input_exprssn_fl.columns.isin(l_excld_clmns)].T
	if l_excld_rows  is not None:
		pd_input_exprssn_fl = pd_input_exprssn_fl[~pd_input_exprssn_fl.index.isin(l_excld_rows)]
	######
	#Define vars and counts
	######
	pd_obs = pd_input_exprssn_fl[pd_input_exprssn_fl.index.str.startswith('#')]
	pd_counts = pd_input_exprssn_fl[~pd_input_exprssn_fl.index.str.startswith('#')].astype(float)
	#Remove # from observations
	pd_obs.index = pd_obs.index.str.replace('#','')
	######
	#list of samples
	######
	lObs = pd_counts.columns
	#list of variables, i.e. gene names
	lVar = pd_counts.index
	######
	#Create anndata
	######
	adata = sc.AnnData(X=pd_counts.T.to_numpy())
	adata.obs_names=lObs.to_numpy()
	adata.var_names=lVar.to_numpy()
	#Add metada for observations
	adata.obs=pd_obs.T
	#Output
	return(adata)

def cmptPrwsCrrltn(adata,lGnsSbst=None,dGrpSbst=None,dGrpExcld=None,lSmplSbst=None, \
	strict_gene_list=False,target_correlation_matrix=None,rtrn_smpl_counts=False,verbose=True):
	"""
	Method to compute pairwise Spearman's correlations for genes in a adata file.
	Input: adata is a h5ad file with obs_names as samples, obs as group sets, and var_names as
	gene names. lGnsSbst is a subset of genes. dGrpSbst is a  dictionary with the obs variable as key
	and the vaues to select as values. This will be used to subset the groups of interest. dGrpExcld is 
	a  dictionary with the obs variable as key	and the groups to exclude as values. This will be used 
	to exclude groups. lSmplSbst	is a subset of samples.  strict_gene_list is a switch, if True will raise
	an Exception if at least one gene in dGrpSbst is not present in adata. target_correlation_matrix is
	the output file to write the correlation matrix. rtrn_smpl_counts is a switch, if True will return the 
	number of samples analyzed in the final matrix.
	Ouput: df_pairwise_spearmanRho is a dataframe of pairwise Spearman correlations obtain from
	from only the upper corner of the full pairwise correlatio matrix.
	"""
	######
	#Subset the input adata file
	######
	adata_sbst = adata.copy()
	#Subset by genes
	if lGnsSbst is not None:
		sbstVar = adata_sbst.var_names.isin(lGnsSbst)
		absntGns = set(lGnsSbst).difference(set(adata_sbst.var_names))
		if strict_gene_list and absntGns:
			raise Exception('The following genes are absent from adata: %s'%', '.join(sorted(absntGns)))
		else:
			if absntGns:
				warnings.warn('The following genes are absent from adata: %s'%', '.join(sorted(absntGns)))
			adata_sbst = adata_sbst[:,sbstVar]
	#Subset by group - include
	if dGrpSbst is not None:
		assert len(dGrpSbst.keys())==1
		grpKy = list(dGrpSbst.keys())[0]
		grpVl = dGrpSbst[grpKy]
		if type(grpVl) is str:
			grpVl = [grpVl,]
		else:
			assert type(grpVl) is list
		adata_sbst = adata_sbst[adata_sbst.obs[grpKy].isin(grpVl)]
		if set(adata_sbst.obs[grpKy])!=set(grpVl):
			raise Exception('The group(s) %s for key %s is absent from adata: %s'%', '.join(sorted(absntGns)))
	#Subset by group - exclude
	if dGrpExcld is not None:
		assert len(dGrpExcld.keys())==1
		grpKy = list(dGrpExcld.keys())[0]
		grpVl = dGrpExcld[grpKy]
		if type(grpVl) is str:
			grpVl = [grpVl,]
		else:
			assert type(grpVl) is list
		adata_sbst = adata_sbst[~adata_sbst.obs[grpKy].isin(grpVl)]
		print('\t%s samples of the group %s were excluded from adata'%(len(adata_sbst.obs_names), \
		'%s: "%s"'%(grpKy,'","'.join(grpVl))))
	#Subset by sample
	if lSmplSbst is not None:
		absntSmpls = set(lSmplSbst).difference(set(adata_sbst.obs_names))
		if absntSmpls:
			raise Exception('The following samples are absent from adata: %s'%', '.join(sorted(absntSmpls)))
		else:
			adata_sbst = adata_sbst[lSmplSbst]
			print('\t%s samples were selected from adata'%adata_sbst.shape[0])
	######
	#Compute pairwise Spearman correlations
	######
	#Compute Spearman
	df_adata_sbst = adata_sbst.to_df()
	if verbose:
		print('\t%s samples and %s genes were selected in the final computation for correlation matrix'% \
		adata_sbst.shape)
	corr_matrix = df_adata_sbst.corr(method="spearman")
	#output
	if target_correlation_matrix is not None:
		corr_matrix.to_csv(target_correlation_matrix, index=False, header=False)
	#
	if rtrn_smpl_counts:
		return(corr_matrix,adata_sbst.shape[0])
	else:
		return(corr_matrix)

def rtrn2clmn_from_corr_matrix(corr_matrix):
	"""
	Method to return a pandas dataframe non-redundant 2 column matrix from a correlation
	matrix.
	Input: corr_matrix is a square correlation matrix from parwise correlations between genes
	with sorted names in headers and in columns.
	Output: upper_triangle_df is the pandas 2 column matrix with only the upper corner 
	correlation matrix
	"""
	#Test index and headers are the same
	assert (corr_matrix.index==corr_matrix.columns).all()
	corr_matrix.index.name='#Gene_A'
	#Extract the upper triangle and exclude diagonal
	mask = np.triu(np.ones(corr_matrix.shape), k=1).astype(bool)
	upper_triangle_values = corr_matrix.where(mask)
	#Convert the upper triangle to a long-form DataFrame
	upper_triangle_df = upper_triangle_values.stack().reset_index()
	upper_triangle_df.columns = ["#Gene_A", "Gene_B", "Spearman_rho"]
	#output
	return(upper_triangle_df)

def rtrn_corr_matrix_from_2clmn(df_3clmns,value_interest='weighted_rho_(rhoMeta)'):
	"""
	Method to convert a three column pandas dataframe into a correlation matrix.
	Input: df_3clmns is a three column pandas dataframe corresponding to the upper triangle of a 
	square matrix excludng the diagonal. value_interest is the name of the column to add.
	Output: df_3clmns is a square matrix version of the input.
	"""
	#Mirror the upper triangle matrix
	df_3clmns_bottom = df_3clmns.rename(columns={'#Gene_A': 'Gene_B','Gene_B': '#Gene_A'})
	#Combine the original and the mirrored DataFrames
	df_3clmns = pd.concat([df_3clmns,df_3clmns_bottom], ignore_index=True)
	# Pivot to create the square matrix
	df_3clmns = df_3clmns.pivot(index='#Gene_A', columns='Gene_B',values= \
	value_interest)
	# Add diagonal entries with 1
	assert (df_3clmns.columns==df_3clmns.index).all()
	np.fill_diagonal(df_3clmns.values, 1)
	return(df_3clmns)

def wrpr_process_exprssnDB(input_reference_db_file,l_genes_interest,sep='\t',outfile_full_rho_info= \
	None,outfile_weighted_rho_matrix=None,run_metanalysis=True,verbose=True):
	"""
	Wrapper to generate correlation matrices and dataframes
	Input: input_reference_db_file is an input reference file with the columns: "database_name"
	with a database name to work with, "full_path_input_expression_file" with the full path to
	the input database,"separator_of_columns_in_input_expression_file" for the separator of 
	columns in the input file, "columns_to_exclude_in_input_expression_file" for values in the 
	columns that needs to be excluded -multiple values accepted joined by "|"-, 
	"rows_to_exclude_in_input_expression_file" for values that need to be excluded from the 
	rows (using the first column) -multiple values accepted joined by "|"-, 
	"variable_of_interest_in_patients_in_input_expression_file" with the variable in the input file
	that has the groups of interest to include/exclude, 
	"groups_to_select_in_variable_of_interest_in_patients_in_input_expression_file" for groups
	in the variable of interest to be included -multiple values accepted joined by "|"-, 
	"groups_to_exclude_in_variable_of_interest_in_patients_in_input_expression_file" for groups
	on the input expression file to be excluded -multiple values accepted joined by "|"-, 
	"patients_to_select_in_input_expression_file" for	patient names on the rows to be included 
	-multiple values accepted joined by "|"-. 
	l_genes_interest is the list of genes of interest to compute the pairwise correlations.  sep is
	the delimiter for the outfile_full_rho_info in case one is provided. outfile_full_rho_info is an
	output file to write the pairwise correlation in genes in l_genes_interest, in each database in
	input_reference_db_file, and in addition the weighted rho correlations, and the average rho
	correlations. outfile_weighted_rho_matrix is the output file to write the pairwise weighted 
	rho matrix for the genes of interest. run_metanalysis is a switch, if True will run a 
	metanalysis and return df_reference_rho_matrix. Otherwise, it will return df_vctrs_crrltns
	that contains pairwise correlations values for each database independently.
	Output: It returns df_reference_rho_matrix which is the weighted pairwise rho correlations 
	using the metadata approach.	
	"""
	#Read input database
	pd_input_reference_db = pd.read_csv(input_reference_db_file,sep=sep,header=0, \
	index_col=0,low_memory=False)#
	#List all the database names
	lDBs = list(pd_input_reference_db.index)
	#Database vector's storage and db sample size
	df_vctrs_crrltns,df_db_smplSize = None,None
	####
	#Collect information from each database
	####
	for db in lDBs:
		db_info = pd_input_reference_db[pd_input_reference_db.index==db]
		input_exprssn_fl = db_info['full_path_input_expression_file'].iloc[0]
		delimiter = db_info['separator_of_columns_in_input_expression_file'].iloc[0].encode(). \
		decode('unicode_escape')
		#Variable to work with
		vrbl_to_test = db_info["variable_of_interest_in_patients_in_input_expression_file"].iloc[0]
		if pd.isna(vrbl_to_test):
			vrbl_to_test = None
		#Exclude columns
		l_excld_clmns = db_info["columns_to_exclude_in_input_expression_file"].iloc[0]
		if not pd.isna(l_excld_clmns):
			l_excld_clmns = l_excld_clmns.split('|')
		else:
			l_excld_clmns = None
		#Exclude rows
		l_excld_rows = db_info["rows_to_exclude_in_input_expression_file"].iloc[0]
		if not pd.isna(l_excld_rows):
			l_excld_rows = l_excld_rows.split('|')
		else:
			l_excld_rows = None
		#Groups to test
		groups_to_test = db_info \
		["groups_to_select_in_variable_of_interest_in_patients_in_input_expression_file"].iloc[0]
		if not pd.isna(groups_to_test):
			groups_to_test = groups_to_test.split('|')
			try:
				assert vrbl_to_test is not None
			except:
				raise Exception('\tThe input database requires a valid value in "variable_of_interest_in_patients_in_input_expression_file" to continue with the value "%s"'% \
				('groups_to_select_in_variable_of_interest_in_patients_in_input_expression_file:%s'% \
				groups_to_test))
			dGrpSbst={vrbl_to_test:groups_to_test}
		else:
			dGrpSbst = None
		#Groups to exclude
		groups_to_exclude = db_info \
		["groups_to_exclude_in_variable_of_interest_in_patients_in_input_expression_file"]. \
		iloc[0]
		if not pd.isna(groups_to_exclude):
			groups_to_exclude = groups_to_exclude.split('|')
			try:
				assert vrbl_to_test is not None
			except:
				raise Exception('\tThe input database requires a valid value in "variable_of_interest_in_patients_in_input_expression_file" to continue with the value "%s"'% \
				('groups_to_exclude_in_variable_of_interest_in_patients_in_input_expression_file:%s'% \
				groups_to_exclude))
			dGrpExcld = {vrbl_to_test:groups_to_exclude}
		else:
			dGrpExcld = None
		#Patients to include
		ptntsIncld = db_info["patients_to_select_in_input_expression_file"].iloc[0]
		if not pd.isna(ptntsIncld):
			ptntsIncld = ptntsIncld.split('|')
		else:
			ptntsIncld = None
		####
		#Convert to adata file
		adata = csvToH5ad(input_exprssn_fl,l_excld_clmns=l_excld_clmns,l_excld_rows= \
		l_excld_rows,sep=delimiter)
		if verbose:
			print('\t%s data points were obtained from database %s'%(adata.shape[0],db))
		####
		#Compute correlation and create vector of the upper triangle
		corr_matrix,smpl_counts = cmptPrwsCrrltn(adata,lGnsSbst=l_genes_interest,dGrpSbst= \
		dGrpSbst,dGrpExcld=dGrpExcld,lSmplSbst=ptntsIncld,strict_gene_list=True, \
		rtrn_smpl_counts=True,verbose=verbose)
		upper_triangle_df = rtrn2clmn_from_corr_matrix(corr_matrix)
		####
		#Add to merged dataframe
		upper_triangle_df = upper_triangle_df.rename(columns={'Spearman_rho':db})
		pd_n_sample = pd.DataFrame({db:smpl_counts},index=['n_samples'])#number of samples
		if df_vctrs_crrltns is not None:
			df_vctrs_crrltns = pd.merge(df_vctrs_crrltns, upper_triangle_df[['#Gene_A','Gene_B',\
			db]],on=['#Gene_A','Gene_B'], how='left')
			df_db_smplSize = pd.concat([df_db_smplSize,pd_n_sample],axis=1)
		else:
			df_vctrs_crrltns = upper_triangle_df
			df_db_smplSize = pd_n_sample
	####
	#Conduct meta-analysis following Shilts paper
	####
	if run_metanalysis:
		if verbose:
			print('\tConducting meta-analysis of data...')
		#Step 1: Apply Fisher z-transformation
		vctr_crrltns = df_vctrs_crrltns.iloc[:,~(df_vctrs_crrltns.columns.isin(['#Gene_A','Gene_B']))]
		z_values = np.arctanh(vctr_crrltns)
		#Step 2: Compute weights as n_i - 3
		weights = df_db_smplSize-3
		#Step 3: Compute the weighted average of z-values
		weighted_z = np.sum(weights.to_numpy()*z_values,axis=1)/np.sum(weights.to_numpy())
		#Step 4: Apply the inverse Fisher transformation
		weighted_rho = np.tanh(weighted_z)
		#Compute normal average
		average_rho = vctr_crrltns.mean(axis=1)
		#Add all values to the pandas dataframe
		df_vctrs_crrltns['average_rho'] = average_rho
		#Add all values to the pandas dataframe
		df_vctrs_crrltns['weighted_rho_(rhoMeta)'] = weighted_rho
		####
		#Obtain a pairwise matrix from weighted_rho
		####
		if verbose:
			print('\tComputing pairwise reference correlation matrix...')
		df_reference_rho_matrix = df_vctrs_crrltns[['#Gene_A','Gene_B','weighted_rho_(rhoMeta)']]
		df_reference_rho_matrix = rtrn_corr_matrix_from_2clmn(df_reference_rho_matrix)
		df_reference_rho_matrix.index.name=''
		df_reference_rho_matrix.columns.name=''
		###
		#Output
		if outfile_full_rho_info is not None:
			df_vctrs_crrltns.to_csv(outfile_full_rho_info,sep=sep,index=False)
			#Make plot
			outfile_full_rho_info_svg = '%s.svg'%outfile_full_rho_info.replace('.csv','').replace('.tsv','')
			mkHeatMap_fromDF(df_vctrs_crrltns,outfile_full_rho_info_svg)
		if outfile_weighted_rho_matrix is not None:
			df_reference_rho_matrix.to_csv(outfile_weighted_rho_matrix)
			df_reference_rho_matrix_plot = df_reference_rho_matrix.copy()
			np.fill_diagonal(df_reference_rho_matrix_plot.values, np.nan)
			#Make plot			
			plt.figure(figsize=(6, 5))
			sns.heatmap(df_reference_rho_matrix, annot=True, cmap='coolwarm', square=True,  \
			linewidths=0.5, cbar_kws={'label': 'Correlation'})
			plt.title('Pairwise correlations for metanalysis of reference databases')
			plt.tight_layout()
			plt.savefig('%s.svg'%outfile_weighted_rho_matrix.replace('.csv','').replace('.tsv',''))
		if verbose:
			print('Done!...')
		return(df_reference_rho_matrix)
	#
	else:
		###
		#Output
		if outfile_full_rho_info is not None:
			df_vctrs_crrltns.to_csv(outfile_full_rho_info,sep=sep,index=False)
		#
		return(df_vctrs_crrltns)

def mkHeatMap_fromDF(df_vctrs_crrltns,outfile_full_rho_info_svg):
	"""
	Manke an output heatmap plot for a pandas dataframe
	Input: df_vctrs_crrltns the dataframe. outfile_full_rho_info is the output plot file
	"""
	#Identify unique labels
	labels = pd.unique(df_vctrs_crrltns[['#Gene_A','Gene_B']].values.ravel())
	labels.sort()
	#Create Pivot Tables for Each Metric
	matrices = {}# Initialize a dictionary to store matrices
	for metric in df_vctrs_crrltns.columns.difference(['#Gene_A','Gene_B']):
		matrix = pd.DataFrame(np.nan, index=labels, columns=labels)
		for _, row in df_vctrs_crrltns.iterrows():
			a, b = row['#Gene_A'], row['Gene_B']
			value = row[metric]
			matrix.at[a, b] = value
			matrix.at[b, a] = value  # Assuming symmetry
		matrices[metric] = matrix
	#Make plot
	num_metrics = len(matrices)
	fig, axes = plt.subplots(1, num_metrics, figsize=(6 * num_metrics, 5))
	if num_metrics == 1:
		axes = [axes]  # Ensure axes is iterable
	for ax, (metric, matrix) in zip(axes, matrices.items()):
		sns.heatmap(matrix, ax=ax, annot=True, cmap='coolwarm', square=True, linewidths=0.5, cbar_kws={'label': metric})
		ax.set_title(f'Heatmap of {metric}')
		ax.set_xlabel('Clock genes')
		ax.set_ylabel('Clock genes')
	plt.tight_layout()
	plt.savefig(outfile_full_rho_info_svg)

def cmptCCD(sqr_mtrx_db_reference,sqr_mtrx_db_target,scale=False):
	"""
	Method to compute Euclidean distance.
	Input: sqr_mtrx_db_reference is a square matrix for reference. sqr_mtrx_db_target is a square
	matrix for targets. scale is a switch if True will scale the distance following the original 
	implementation of calcCCDSimple (i.e. Spearman correlation by the number of pairs).
	Output: ccd is the euclidean ccd between both excluding nan
	"""
	# Compute the element-wise difference
	diff_matrix = sqr_mtrx_db_reference.values - sqr_mtrx_db_target.values
	# Compute the squared differences, ignoring NaNs
	squared_diff_matrix = np.square(diff_matrix)
	# Sum the squared differences, ignoring NaNs
	sum_squared_diff_matrix = np.nansum(squared_diff_matrix)
	# Compute the Euclidean ccd
	ccd = np.sqrt(sum_squared_diff_matrix)
	#Weight
	if scale:
		nPairs = math.comb(ccd.shape[0],2)
		ccd = ccd/nPairs
	#Output
	return(ccd)

def generate_random_lists_of_genes(full_list_target_genes,number_randomizations,size_of_sample, \
	seed=42):
	"""
	Method to generate random lists of genes for randomization
	Input: full_list_target_genes is a full list of genes of interest to select random list from. 
	number_randomizations is the number of lists of genes to obtain. size_of_sample is the number of
	genes to select. seed is a seed to start the randomization.
	Output: df_rndm_data is a dataframe with lists of random genes.
	"""
	np.random.seed(seed)  # Set the seed for reproducibility
	# Generate a 2D array where each column is a random sample of genes
	# ~ rndm_data = [np.random.choice(full_list_target_genes, size=size_of_sample, replace=True) for _ \
	# ~ in range(number_randomizations)]
	rndm_data = [np.random.permutation(full_list_target_genes)[:size_of_sample] for _ \
	in range(number_randomizations)]
	# Transpose the array to have genes as rows and samples as columns
	df_rndm_data = pd.DataFrame(np.array(rndm_data).T, columns=[f'rand_{i+1}' for i in \
	range(number_randomizations)])
	#
	return(df_rndm_data)


def build_target_gene_matching_summary(input_target_db_file,candidate_genes,sep='\t',min_variance=1e-3,verbose=True):
    """
    Method to generate a per-gene summary for matched-null sampling.
    Input: input_target_db_file is an input target file with the columns required by the NCC
    tool to load each cohort and apply the same cohort-level filters used in the downstream CCD
    analyses. candidate_genes is the full list of target genes to summarize. sep is the delimiter
    in input_target_db_file. min_variance is the minimum within-cohort gene variance required for
    a gene to be retained in the matching summary. verbose is a switch, if True will print the
    number of genes retained per cohort.
    Output: df_summary is a dataframe with one row per gene and aggregated columns summarizing
    mean-expression and variance percentile ranks across cohorts, including mean_rank_pct and
    variance_rank_pct.
    """
    lCohorts = load_target_cohorts_for_delta(input_target_db_file,sep=sep,verbose=verbose)
    lFrames = []
    for db,adata in lCohorts:
        cohort_genes = sorted(set(candidate_genes).intersection(set(adata.var_names)))
        if len(cohort_genes) == 0:
            raise Exception('No candidate genes overlap with cohort %s while building matched-null summaries'%db)
        df_samples_x_genes = adata[:,cohort_genes].to_df().loc[:,cohort_genes]
        gene_mean = df_samples_x_genes.mean(axis=0)
        gene_var = df_samples_x_genes.var(axis=0)
        keep = gene_var > min_variance
        gene_mean = gene_mean.loc[keep]
        gene_var = gene_var.loc[keep]
        df_db = pd.DataFrame({
            f'mean_expr__{db}':gene_mean,
            f'variance__{db}':gene_var,
            f'mean_rank_pct__{db}':gene_mean.rank(method='average',pct=True),
            f'variance_rank_pct__{db}':gene_var.rank(method='average',pct=True)
        })
        if verbose:
            print('\t%s genes retained for matched-null summaries in cohort %s'%(df_db.shape[0],db))
        lFrames.append(df_db)
    df_summary = pd.concat(lFrames,axis=1,join='inner')
    mean_rank_cols = [c for c in df_summary.columns if c.startswith('mean_rank_pct__')]
    variance_rank_cols = [c for c in df_summary.columns if c.startswith('variance_rank_pct__')]
    df_summary['mean_rank_pct'] = df_summary[mean_rank_cols].mean(axis=1)
    df_summary['variance_rank_pct'] = df_summary[variance_rank_cols].mean(axis=1)
    df_summary = df_summary.sort_index()
    return(df_summary)

def _matched_candidate_count(df_matching_summary,target_gene,excluded_genes,mean_caliper=0.05,variance_caliper=0.05):
    """
    Method to count the number of matched candidate genes for one target gene.
    Input: df_matching_summary is a dataframe summarizing the genes available for matched-null
    sampling. target_gene is the gene of interest to match. excluded_genes is a set or list of
    genes that cannot be sampled. mean_caliper and variance_caliper define the maximum absolute
    difference allowed in mean-expression rank and variance rank, respectively.
    Output: n_candidates is the number of genes in df_matching_summary that can be used as
    matched candidates for target_gene under the requested calipers.
    """
    if target_gene not in df_matching_summary.index:
        return(0)
    df_remaining = df_matching_summary.loc[~df_matching_summary.index.isin(excluded_genes)]
    if df_remaining.empty:
        return(0)
    trg = df_matching_summary.loc[target_gene]
    m = (df_remaining['mean_rank_pct'].sub(trg['mean_rank_pct']).abs() <= mean_caliper) & \
    (df_remaining['variance_rank_pct'].sub(trg['variance_rank_pct']).abs() <= variance_caliper)
    return(int(m.sum()))

def _sample_one_matched_gene(df_matching_summary,target_gene,excluded_genes,rng,mean_caliper=0.05,variance_caliper=0.05,expand_step=0.05,max_caliper=1.0):
    """
    Sample a single matched gene for target_gene, excluding the genes in excluded_genes.
    Matching is performed on aggregated mean-expression and variance percentile ranks.
    If no candidate is available inside the current caliper, the caliper is expanded
    iteratively; as a last resort, the nearest remaining gene in 2D rank space is used.
    """
    if target_gene not in df_matching_summary.index:
        raise Exception('The target gene %s is absent from the matched-null summary table'%target_gene)
    df_remaining = df_matching_summary.loc[~df_matching_summary.index.isin(excluded_genes)].copy()
    if df_remaining.empty:
        raise Exception('No remaining genes available for matched-null sampling')
    trg = df_matching_summary.loc[target_gene]
    cur_mean = float(mean_caliper)
    cur_var = float(variance_caliper)
    while True:
        mask = (df_remaining['mean_rank_pct'].sub(trg['mean_rank_pct']).abs() <= cur_mean) & \
        (df_remaining['variance_rank_pct'].sub(trg['variance_rank_pct']).abs() <= cur_var)
        l_candidates = df_remaining.index[mask].tolist()
        if l_candidates:
            return(str(rng.choice(l_candidates)))
        if cur_mean >= max_caliper and cur_var >= max_caliper:
            break
        cur_mean = min(float(max_caliper),cur_mean + float(expand_step))
        cur_var = min(float(max_caliper),cur_var + float(expand_step))
    distances = np.sqrt(np.square(df_remaining['mean_rank_pct'] - trg['mean_rank_pct']) + \
    np.square(df_remaining['variance_rank_pct'] - trg['variance_rank_pct']))
    nearest = distances.min()
    l_candidates = distances.index[distances == nearest].tolist()
    return(str(rng.choice(l_candidates)))

def generate_matched_random_lists_of_genes(target_genes,df_matching_summary,number_randomizations,seed=42,mean_caliper=0.05,variance_caliper=0.05,expand_step=0.05,max_caliper=1.0):
    """
    Generate matched-null random gene sets with the same size as target_genes, matched on
    aggregated expression and variance percentile ranks.
    Output dataframe has one column per randomization and one row per target gene, where the
    index corresponds to the target genes ordered from the most difficult to the easiest to match.
    """
    target_genes = list(target_genes)
    missing = sorted(set(target_genes).difference(df_matching_summary.index))
    if missing:
        raise Exception('The following target genes are absent from the matched-null summary table: %s'%', '.join(missing))
    if len(target_genes) == 0:
        raise Exception('None of the target genes were found in the matched-null summary table')
    s_excluded_base = set(target_genes)
    l_target_order = sorted(target_genes,key=lambda g: (_matched_candidate_count(df_matching_summary,g,s_excluded_base,mean_caliper=mean_caliper,variance_caliper=variance_caliper),g))
    rng = np.random.default_rng(seed)
    d_random_sets = {}
    for rand_idx in range(number_randomizations):
        s_excluded = set(s_excluded_base)
        l_selected = []
        for trg in l_target_order:
            matched_gene = _sample_one_matched_gene(df_matching_summary,trg,s_excluded,rng,mean_caliper=mean_caliper,variance_caliper=variance_caliper,expand_step=expand_step,max_caliper=max_caliper)
            l_selected.append(matched_gene)
            s_excluded.add(matched_gene)
        if len(set(l_selected)) != len(l_selected):
            raise Exception('Matched-null sampling produced duplicate genes within randomization rand_%s'%(rand_idx+1))
        d_random_sets['rand_%s'%(rand_idx+1)] = l_selected
    df_rndm_data = pd.DataFrame(d_random_sets,index=l_target_order)
    return(df_rndm_data)

def rtrn_subset_genes_target(input_target_db_file,full_list_target_genes,sep='\t',min_variance=1e-3):
	"""
	Method to subset a list of genes present in a target database.
	Wrapper to generate correlation matrices and dataframes
	Input: input_target_db_file is an input target file with the columns: "database_name"
	with a database name to work with, "full_path_input_expression_file" with the full path to
	the input database,"separator_of_columns_in_input_expression_file" for the separator of 
	columns in the input file, "columns_to_exclude_in_input_expression_file" for values in the 
	columns that needs to be excluded -multiple values accepted joined by "|"-, 
	"rows_to_exclude_in_input_expression_file" for values that need to be excluded from the 
	rows (using the first column) -multiple values accepted joined by "|"-, 
	"variable_of_interest_in_patients_in_input_expression_file" with the variable in the input file
	that has the groups of interest to include/exclude, 
	"groups_to_select_in_variable_of_interest_in_patients_in_input_expression_file" for groups
	in the variable of interest to be included -multiple values accepted joined by "|"-, 
	"groups_to_exclude_in_variable_of_interest_in_patients_in_input_expression_file" for groups
	on the input expression file to be excluded -multiple values accepted joined by "|"-, 
	"patients_to_select_in_input_expression_file" for	patient names on the rows to be included 
	-multiple values accepted joined by "|"-. full_list_target_genes is the list of genes to subset. 
	sep is 	the delimiter for the outfile_full_rho_info in case one is provided.  min_variance is a 
	value below which genes are going to be filtered.
	Output: full_list_target_genes including only the list of genes present in all databases.
	"""
	pd_input_reference_db = pd.read_csv(input_target_db_file,sep=sep,header=0, \
	index_col=0,low_memory=False)#
	#Work input database
	lDBs = list(pd_input_reference_db.index)
	####
	#Collect information from each database
	####
	for db in lDBs:
		db_info = pd_input_reference_db[pd_input_reference_db.index==db]
		input_exprssn_fl = db_info['full_path_input_expression_file'].iloc[0]
		pd_input_exprssn_fl = pd.read_csv(input_exprssn_fl,sep=sep,header=0,index_col=0, \
		low_memory=False)
		##
		#Exclude columns
		l_excld_clmns = db_info["columns_to_exclude_in_input_expression_file"].iloc[0]
		if not pd.isna(l_excld_clmns):
			l_excld_clmns = l_excld_clmns.split('|')
		else:
			l_excld_clmns = None
		#Exclude rows
		l_excld_rows = db_info["rows_to_exclude_in_input_expression_file"].iloc[0]
		if not pd.isna(l_excld_rows):
			l_excld_rows = l_excld_rows.split('|')
		else:
			l_excld_rows = None
		##
		if l_excld_clmns is not None:#Adjust for exclude columns
			pd_input_exprssn_fl = pd_input_exprssn_fl.T[~pd_input_exprssn_fl.columns.isin(l_excld_clmns)].T
		if l_excld_rows  is not None:
			pd_input_exprssn_fl = pd_input_exprssn_fl[~pd_input_exprssn_fl.index.isin(l_excld_rows)]
		pd_counts = pd_input_exprssn_fl[~pd_input_exprssn_fl.index.str.startswith('#')].astype(float)
		###
		#FIlter for variance > 0
		row_variance = pd_counts.var(axis=1)
		pd_counts = pd_counts[row_variance > min_variance]
		###
		full_list_target_genes = sorted(set(full_list_target_genes).intersection(set(pd_counts.index)))
	return(full_list_target_genes)


def revisit_gnsInterest(input_target_db_file,l_genes_interest,sep='\t',verbose=True,min_variance=1e-3):
	"""
	Wrapper to further filter genes of interest for subsets of genes with valid variance
	Input: input_target_db_file is an input reference file with the columns: "database_name"
	with a database name to work with, "full_path_input_expression_file" with the full path to
	the input database,"separator_of_columns_in_input_expression_file" for the separator of 
	columns in the input file, "columns_to_exclude_in_input_expression_file" for values in the 
	columns that needs to be excluded -multiple values accepted joined by "|"-, 
	"rows_to_exclude_in_input_expression_file" for values that need to be excluded from the 
	rows (using the first column) -multiple values accepted joined by "|"-, 
	"variable_of_interest_in_patients_in_input_expression_file" with the variable in the input file
	that has the groups of interest to include/exclude, 
	"groups_to_select_in_variable_of_interest_in_patients_in_input_expression_file" for groups
	in the variable of interest to be included -multiple values accepted joined by "|"-, 
	"groups_to_exclude_in_variable_of_interest_in_patients_in_input_expression_file" for groups
	on the input expression file to be excluded -multiple values accepted joined by "|"-, 
	"patients_to_select_in_input_expression_file" for	patient names on the rows to be included 
	-multiple values accepted joined by "|"-. 
	l_genes_interest is the list of genes of interest to compute the pairwise correlations.  sep is
	the delimiter for the outfile_full_rho_info in case one is provided. outfile_full_rho_info is an
	output file to write the pairwise correlation in genes in l_genes_interest, in each database in
	input_target_db_file, and in addition the weighted rho correlations, and the average rho
	correlations. outfile_weighted_rho_matrix is the output file to write the pairwise weighted 
	rho matrix for the genes of interest. run_metanalysis is a switch, if True will run a 
	metanalysis and return df_reference_rho_matrix. Otherwise, it will return df_vctrs_crrltns
	that contains pairwise correlations values for each database independently.
	Output: It returns df_reference_rho_matrix which is the weighted pairwise rho correlations 
	using the metadata approach.	
	"""
	#Read input database
	pd_input_reference_db = pd.read_csv(input_target_db_file,sep=sep,header=0, \
	index_col=0,low_memory=False)#
	#List all the database names
	lDBs = list(pd_input_reference_db.index)
	#Dfine set of genes present after filtering
	sGnsIntrstSbst = None
	####
	#Collect information from each database
	####
	for db in lDBs:
		db_info = pd_input_reference_db[pd_input_reference_db.index==db]
		input_exprssn_fl = db_info['full_path_input_expression_file'].iloc[0]
		delimiter = db_info['separator_of_columns_in_input_expression_file'].iloc[0].encode(). \
		decode('unicode_escape')
		#Variable to work with
		vrbl_to_test = db_info["variable_of_interest_in_patients_in_input_expression_file"].iloc[0]
		if pd.isna(vrbl_to_test):
			vrbl_to_test = None
		#Exclude columns
		l_excld_clmns = db_info["columns_to_exclude_in_input_expression_file"].iloc[0]
		if not pd.isna(l_excld_clmns):
			l_excld_clmns = l_excld_clmns.split('|')
		else:
			l_excld_clmns = None
		#Exclude rows
		l_excld_rows = db_info["rows_to_exclude_in_input_expression_file"].iloc[0]
		if not pd.isna(l_excld_rows):
			l_excld_rows = l_excld_rows.split('|')
		else:
			l_excld_rows = None
		#Groups to test
		groups_to_test = db_info \
		["groups_to_select_in_variable_of_interest_in_patients_in_input_expression_file"].iloc[0]
		if not pd.isna(groups_to_test):
			groups_to_test = groups_to_test.split('|')
			try:
				assert vrbl_to_test is not None
			except:
				raise Exception('\tThe input database requires a valid value in "variable_of_interest_in_patients_in_input_expression_file" to continue with the value "%s"'% \
				('groups_to_select_in_variable_of_interest_in_patients_in_input_expression_file:%s'% \
				groups_to_test))
			dGrpSbst={vrbl_to_test:groups_to_test}
		else:
			dGrpSbst = None
		#Groups to exclude
		groups_to_exclude = db_info \
		["groups_to_exclude_in_variable_of_interest_in_patients_in_input_expression_file"]. \
		iloc[0]
		if not pd.isna(groups_to_exclude):
			groups_to_exclude = groups_to_exclude.split('|')
			try:
				assert vrbl_to_test is not None
			except:
				raise Exception('\tThe input database requires a valid value in "variable_of_interest_in_patients_in_input_expression_file" to continue with the value "%s"'% \
				('groups_to_exclude_in_variable_of_interest_in_patients_in_input_expression_file:%s'% \
				groups_to_exclude))
			dGrpExcld = {vrbl_to_test:groups_to_exclude}
		else:
			dGrpExcld = None
		#Patients to include
		ptntsIncld = db_info["patients_to_select_in_input_expression_file"].iloc[0]
		if not pd.isna(ptntsIncld):
			ptntsIncld = ptntsIncld.split('|')
		else:
			ptntsIncld = None
		####
		#Convert to adata file
		adata = csvToH5ad(input_exprssn_fl,l_excld_clmns=l_excld_clmns,l_excld_rows= \
		l_excld_rows,sep=delimiter)
		if verbose:
			print('\t%s data points were obtained from database %s'%(adata.shape[0],db))
		######
		#Subset the input adata file
		######
		adata_sbst = adata.copy()
		#Subset by genes
		if l_genes_interest is not None:
			sbstVar = adata_sbst.var_names.isin(l_genes_interest)
			absntGns = set(l_genes_interest).difference(set(adata_sbst.var_names))
			adata_sbst = adata_sbst[:,sbstVar]
		#Subset by group - include
		if dGrpSbst is not None:
			assert len(dGrpSbst.keys())==1
			grpKy = list(dGrpSbst.keys())[0]
			grpVl = dGrpSbst[grpKy]
			if type(grpVl) is str:
				grpVl = [grpVl,]
			else:
				assert type(grpVl) is list
			adata_sbst = adata_sbst[adata_sbst.obs[grpKy].isin(grpVl)]
			if set(adata_sbst.obs[grpKy])!=set(grpVl):
				raise Exception('The group(s) %s for key %s is absent from adata: %s'%', '.join(sorted(absntGns)))
		#Subset by group - exclude
		if dGrpExcld is not None:
			assert len(dGrpExcld.keys())==1
			grpKy = list(dGrpExcld.keys())[0]
			grpVl = dGrpExcld[grpKy]
			if type(grpVl) is str:
				grpVl = [grpVl,]
			else:
				assert type(grpVl) is list
			adata_sbst = adata_sbst[~adata_sbst.obs[grpKy].isin(grpVl)]
			print('\t%s samples of the group %s were excluded from adata'%(len(adata_sbst.obs_names), \
			'%s: "%s"'%(grpKy,'","'.join(grpVl))))
		#Subset by sample
		if ptntsIncld is not None:
			absntSmpls = set(ptntsIncld).difference(set(adata_sbst.obs_names))
			if absntSmpls:
				raise Exception('The following samples are absent from adata: %s'%', '.join(sorted(absntSmpls)))
			else:
				adata_sbst = adata_sbst[ptntsIncld]
				print('\t%s samples were selected from adata'%adata_sbst.shape[0])
		###
		#Revist genes names
		row_variance = adata_sbst.X.var(axis=0)
		adata_sbst = adata_sbst[:,row_variance > min_variance]
		if sGnsIntrstSbst is None:
			sGnsIntrstSbst=set(adata_sbst.var_names)
		else:
			sGnsIntrstSbst = sGnsIntrstSbst.intersection(set(adata_sbst.var_names))
	#
	return(sorted(sGnsIntrstSbst))


def load_target_cohorts_for_delta(input_target_db_file, sep='\t', verbose=True):
	"""
	Load target cohorts after the same database-level include/exclude filtering used elsewhere
	in the NCC tool, but without collapsing across contrasted groups.
	Output: list of tuples (db_name, adata).
	"""
	pd_input_target_db = pd.read_csv(input_target_db_file,sep=sep,header=0,index_col=0,low_memory=False)
	lDBs = list(pd_input_target_db.index)
	lCohorts = []
	for db in lDBs:
		db_info = pd_input_target_db[pd_input_target_db.index==db]
		input_exprssn_fl = db_info['full_path_input_expression_file'].iloc[0]
		delimiter = db_info['separator_of_columns_in_input_expression_file'].iloc[0].encode().decode('unicode_escape')
		# Variable to work with
		vrbl_to_test = db_info["variable_of_interest_in_patients_in_input_expression_file"].iloc[0]
		if pd.isna(vrbl_to_test):
			vrbl_to_test = None
		# Exclude columns
		l_excld_clmns = db_info["columns_to_exclude_in_input_expression_file"].iloc[0]
		if not pd.isna(l_excld_clmns):
			l_excld_clmns = l_excld_clmns.split('|')
		else:
			l_excld_clmns = None
		# Exclude rows
		l_excld_rows = db_info["rows_to_exclude_in_input_expression_file"].iloc[0]
		if not pd.isna(l_excld_rows):
			l_excld_rows = l_excld_rows.split('|')
		else:
			l_excld_rows = None
		# Groups to include/exclude from the cohort before contrasting
		groups_to_test = db_info["groups_to_select_in_variable_of_interest_in_patients_in_input_expression_file"].iloc[0]
		if not pd.isna(groups_to_test):
			groups_to_test = groups_to_test.split('|')
			if vrbl_to_test is None:
				raise Exception('\tThe input database requires a valid value in "variable_of_interest_in_patients_in_input_expression_file" to continue with "groups_to_select_in_variable_of_interest_in_patients_in_input_expression_file"')
			dGrpSbst={vrbl_to_test:groups_to_test}
		else:
			dGrpSbst = None
		groups_to_exclude = db_info["groups_to_exclude_in_variable_of_interest_in_patients_in_input_expression_file"].iloc[0]
		if not pd.isna(groups_to_exclude):
			groups_to_exclude = groups_to_exclude.split('|')
			if vrbl_to_test is None:
				raise Exception('\tThe input database requires a valid value in "variable_of_interest_in_patients_in_input_expression_file" to continue with "groups_to_exclude_in_variable_of_interest_in_patients_in_input_expression_file"')
			dGrpExcld = {vrbl_to_test:groups_to_exclude}
		else:
			dGrpExcld = None
		ptntsIncld = db_info["patients_to_select_in_input_expression_file"].iloc[0]
		if not pd.isna(ptntsIncld):
			ptntsIncld = ptntsIncld.split('|')
		else:
			ptntsIncld = None
		adata = csvToH5ad(input_exprssn_fl,l_excld_clmns=l_excld_clmns,l_excld_rows=l_excld_rows,sep=delimiter)
		# Apply the same cohort-level filters as wrpr_process_exprssnDB/revisit_gnsInterest
		if dGrpSbst is not None:
			grpKy = list(dGrpSbst.keys())[0]
			grpVl = dGrpSbst[grpKy]
			adata = adata[adata.obs[grpKy].isin(grpVl)]
		if dGrpExcld is not None:
			grpKy = list(dGrpExcld.keys())[0]
			grpVl = dGrpExcld[grpKy]
			adata = adata[~adata.obs[grpKy].isin(grpVl)]
		if ptntsIncld is not None:
			absntSmpls = set(ptntsIncld).difference(set(adata.obs_names))
			if absntSmpls:
				raise Exception('The following samples are absent from adata: %s'%', '.join(sorted(absntSmpls)))
			adata = adata[ptntsIncld]
		if verbose:
			print('\t%s samples loaded for cohort %s after base filters'%(adata.shape[0],db))
		lCohorts.append((db,adata))
	return(lCohorts)

def compute_corr_matrix_from_sample_gene_df(df_samples_x_genes, genes):
	"""
	Compute a pairwise Spearman correlation matrix from a sample x gene dataframe.
	"""
	df_samples_x_genes = df_samples_x_genes.loc[:,genes]
	corr_matrix = df_samples_x_genes.corr(method="spearman")
	return(corr_matrix.loc[genes,genes])

def compute_delta_ccd_for_cohort(adata, genes, contrast_variable, group_a, group_b, sqr_mtrx_db_reference, min_samples_per_group=10, verbose=True):
	"""
	Compute CCD for two contrasted groups and return their direct delta CCD.
	Output dictionary stores the per-group CCD values together with the eligible data needed
	for permutation testing and bootstrapping.
	"""
	if contrast_variable not in adata.obs.columns:
		raise Exception('The metadata variable "%s" is absent from adata'%contrast_variable)
	# Keep only the contrasted groups
	obs = adata.obs.copy()
	obs[contrast_variable] = obs[contrast_variable].astype(str)
	eligible_mask = obs[contrast_variable].isin([str(group_a),str(group_b)])
	adata_eligible = adata[eligible_mask].copy()
	n_group_a = int((adata_eligible.obs[contrast_variable].astype(str)==str(group_a)).sum())
	n_group_b = int((adata_eligible.obs[contrast_variable].astype(str)==str(group_b)).sum())
	if n_group_a < min_samples_per_group or n_group_b < min_samples_per_group:
		raise Exception('Insufficient samples after filtering for %s: %s=%s, %s=%s'%(contrast_variable,group_a,n_group_a,group_b,n_group_b))
	missing_genes = sorted(set(genes).difference(set(adata_eligible.var_names)))
	if missing_genes:
		raise Exception('The following genes are absent from adata: %s'%', '.join(missing_genes))
	df_samples_x_genes = adata_eligible[:,genes].to_df().loc[:,genes]
	group_a_ids = adata_eligible.obs_names[adata_eligible.obs[contrast_variable].astype(str)==str(group_a)].tolist()
	group_b_ids = adata_eligible.obs_names[adata_eligible.obs[contrast_variable].astype(str)==str(group_b)].tolist()
	corr_a = compute_corr_matrix_from_sample_gene_df(df_samples_x_genes.loc[group_a_ids], genes)
	corr_b = compute_corr_matrix_from_sample_gene_df(df_samples_x_genes.loc[group_b_ids], genes)
	mask = np.triu(np.ones(corr_a.shape), k=1).astype(bool)
	upper_a = corr_a.where(mask)
	upper_b = corr_b.where(mask)
	ccd_a = cmptCCD(sqr_mtrx_db_reference,upper_a)
	ccd_b = cmptCCD(sqr_mtrx_db_reference,upper_b)
	delta_ccd = ccd_a - ccd_b
	if verbose:
		print('\t%s=%s samples: %s | %s=%s samples: %s | delta CCD: %s'%(group_a,n_group_a,ccd_a,group_b,n_group_b,ccd_b,delta_ccd))
	return({
		'n_group_a':n_group_a,
		'n_group_b':n_group_b,
		'ccd_group_a':ccd_a,
		'ccd_group_b':ccd_b,
		'delta_ccd':delta_ccd,
		'adata_eligible':adata_eligible,
		'df_samples_x_genes':df_samples_x_genes,
	})

def permutation_test_delta_ccd(df_samples_x_genes, obs, genes, contrast_variable, group_a, group_b, sqr_mtrx_db_reference, observed_delta, n_permutations=5000, seed=42):
	"""
	Permutation test for direct between-group delta CCD.
	"""
	rng = np.random.default_rng(seed)
	labels = obs[contrast_variable].astype(str).to_numpy().copy()
	sample_ids = obs.index.to_numpy()
	l_perm_deltas = []
	for perm_id in range(1,n_permutations+1):
		permuted = labels.copy()
		rng.shuffle(permuted)
		group_a_ids = sample_ids[permuted==str(group_a)]
		group_b_ids = sample_ids[permuted==str(group_b)]
		corr_a = compute_corr_matrix_from_sample_gene_df(df_samples_x_genes.loc[group_a_ids], genes)
		corr_b = compute_corr_matrix_from_sample_gene_df(df_samples_x_genes.loc[group_b_ids], genes)
		mask = np.triu(np.ones(corr_a.shape), k=1).astype(bool)
		delta = cmptCCD(sqr_mtrx_db_reference,corr_a.where(mask)) - cmptCCD(sqr_mtrx_db_reference,corr_b.where(mask))
		l_perm_deltas.append((perm_id,delta))
	df_perm = pd.DataFrame(l_perm_deltas,columns=['perm_id','delta_ccd_perm'])
	r_one_sided = int((df_perm['delta_ccd_perm']>=observed_delta).sum())
	r_two_sided = int((df_perm['delta_ccd_perm'].abs()>=abs(observed_delta)).sum())
	df_perm.attrs['p_one_sided_a_gt_b'] = (r_one_sided + 1) / (len(df_perm) + 1)
	df_perm.attrs['p_two_sided'] = (r_two_sided + 1) / (len(df_perm) + 1)
	return(df_perm)

def bootstrap_delta_ccd(df_samples_x_genes, obs, genes, contrast_variable, group_a, group_b, sqr_mtrx_db_reference, n_bootstraps=2000, seed=42):
	"""
	Bootstrap confidence intervals for direct between-group delta CCD.
	"""
	rng = np.random.default_rng(seed)
	group_a_ids = obs.index[obs[contrast_variable].astype(str)==str(group_a)].to_numpy()
	group_b_ids = obs.index[obs[contrast_variable].astype(str)==str(group_b)].to_numpy()
	l_boot_deltas = []
	for bootstrap_id in range(1,n_bootstraps+1):
		boot_a = rng.choice(group_a_ids,size=len(group_a_ids),replace=True)
		boot_b = rng.choice(group_b_ids,size=len(group_b_ids),replace=True)
		corr_a = compute_corr_matrix_from_sample_gene_df(df_samples_x_genes.loc[boot_a], genes)
		corr_b = compute_corr_matrix_from_sample_gene_df(df_samples_x_genes.loc[boot_b], genes)
		mask = np.triu(np.ones(corr_a.shape), k=1).astype(bool)
		delta = cmptCCD(sqr_mtrx_db_reference,corr_a.where(mask)) - cmptCCD(sqr_mtrx_db_reference,corr_b.where(mask))
		l_boot_deltas.append((bootstrap_id,delta))
	return(pd.DataFrame(l_boot_deltas,columns=['bootstrap_id','delta_ccd_boot']))

def fixed_effect_meta_delta_ccd(l_delta_ccd, l_se):
	"""
	Fixed-effect meta-analysis for per-cohort delta CCD values using bootstrap SEs.
	"""
	v_delta = np.asarray(l_delta_ccd,dtype=float)
	v_se = np.asarray(l_se,dtype=float)
	v_var = np.square(v_se)
	v_w = 1.0 / v_var
	meta_delta = float(np.sum(v_w * v_delta) / np.sum(v_w))
	meta_se = float(np.sqrt(1.0 / np.sum(v_w)))
	meta_z = meta_delta / meta_se if meta_se > 0 else np.nan
	meta_p_two_sided = float(math.erfc(abs(meta_z) / math.sqrt(2.0))) if np.isfinite(meta_z) else np.nan
	return({
		'meta_delta_ccd':meta_delta,
		'meta_se':meta_se,
		'meta_ci_low':meta_delta - 1.96 * meta_se,
		'meta_ci_high':meta_delta + 1.96 * meta_se,
		'meta_z':meta_z,
		'meta_p_two_sided':meta_p_two_sided,
	})

def run_direct_delta_ccd(input_target_db_file, genes, contrast_variable, group_a, group_b, sqr_mtrx_db_reference, sep='\t', n_permutations=5000, n_bootstraps=2000, seed=42, min_samples_per_group=10, outfile_results='delta_CCD_results.tsv', outfile_draws=None, verbose=True):
	"""
	Run direct between-group delta CCD contrasts for all cohorts in input_target_db_file and
	write one cohort-level summary table plus optional permutation/bootstrap draws.
	"""
	if contrast_variable is None or group_a is None or group_b is None:
		raise Exception('Direct delta CCD requires --contrast_variable, --group_a and --group_b')
	lCohorts = load_target_cohorts_for_delta(input_target_db_file,sep=sep,verbose=verbose)
	l_summary = []
	l_draws = []
	l_delta_ccd = []
	l_boot_se = []
	for idx,(db,adata) in enumerate(lCohorts):
		if verbose:
			print('Running direct delta CCD for cohort %s...'%db)
		d_obs = compute_delta_ccd_for_cohort(adata,genes,contrast_variable,group_a,group_b,sqr_mtrx_db_reference,min_samples_per_group=min_samples_per_group,verbose=verbose)
		df_perm = permutation_test_delta_ccd(d_obs['df_samples_x_genes'], d_obs['adata_eligible'].obs.copy(), genes, contrast_variable, group_a, group_b, sqr_mtrx_db_reference, d_obs['delta_ccd'], n_permutations=n_permutations, seed=seed + idx)
		df_boot = bootstrap_delta_ccd(d_obs['df_samples_x_genes'], d_obs['adata_eligible'].obs.copy(), genes, contrast_variable, group_a, group_b, sqr_mtrx_db_reference, n_bootstraps=n_bootstraps, seed=seed + 10000 + idx)
		ci_low,ci_high = np.quantile(df_boot['delta_ccd_boot'],[0.025,0.975])
		boot_se = float(df_boot['delta_ccd_boot'].std(ddof=1))
		if boot_se <= 0 or np.isnan(boot_se):
			boot_se = 1e-8
		l_summary.append({
			'database_name':db,
			'contrast_variable':contrast_variable,
			'group_a':group_a,
			'group_b':group_b,
			'n_group_a':d_obs['n_group_a'],
			'n_group_b':d_obs['n_group_b'],
			'ccd_group_a':d_obs['ccd_group_a'],
			'ccd_group_b':d_obs['ccd_group_b'],
			'delta_ccd':d_obs['delta_ccd'],
			'delta_ccd_boot_ci_low':ci_low,
			'delta_ccd_boot_ci_high':ci_high,
			'delta_ccd_boot_se':boot_se,
			'p_perm_one_sided_a_gt_b':df_perm.attrs['p_one_sided_a_gt_b'],
			'p_perm_two_sided':df_perm.attrs['p_two_sided'],
			'n_permutations':n_permutations,
			'n_bootstraps':n_bootstraps,
			'seed':seed,
		})
		l_delta_ccd.append(float(d_obs['delta_ccd']))
		l_boot_se.append(boot_se)
		if outfile_draws is not None:
			l_draws.append(df_perm.assign(database_name=db,draw_type='permutation').rename(columns={'delta_ccd_perm':'delta_ccd_draw'}))
			l_draws.append(df_boot.assign(database_name=db,draw_type='bootstrap').rename(columns={'delta_ccd_boot':'delta_ccd_draw'}))
	pd_summary = pd.DataFrame(l_summary)
	if len(pd_summary) > 1:
		d_meta = fixed_effect_meta_delta_ccd(l_delta_ccd,l_boot_se)
		for ky,vl in d_meta.items():
			pd_summary.loc[:,ky] = vl
	pd_summary.to_csv(outfile_results,sep='\t',index=False)
	if outfile_draws is not None and len(l_draws):
		pd.concat(l_draws,axis=0,ignore_index=True).to_csv(outfile_draws,sep='\t',index=False)
	return(pd_summary)

#############################
#############################
### Generate  reference target gene pairs ###
#############################
#############################
if reference_species!=target_species:
	print('Generating reference:target gene pairs...')
	try:
		assert os.path.exists(reference_target_gene_file)
	except:
		raise Exception('The input reference:target gene pair file "%s" does not exists!'% \
		reference_target_gene_file)
	#Read correlation matrx from file
	df_target_gene = pd.read_csv(reference_target_gene_file,sep='\t',low_memory=False)
	#Make a dictionary reference to target
	df_target_gene.index=df_target_gene[reference_species]
	l_clock_genes_target = df_target_gene[target_species][l_clock_genes_reference].to_list()
	full_list_target_genes = df_target_gene[target_species].to_list()
	print('Done!...')
else:
	l_clock_genes_target = l_clock_genes_reference
	#Read correlation matrx from file
	df_target_gene = pd.read_csv(reference_target_gene_file,sep='\t',low_memory=False)
	full_list_target_genes = df_target_gene[target_species].to_list()
print('Obtaining list of background genes...')
#0c. Subset list of target genes present in input files
full_list_target_genes = rtrn_subset_genes_target(input_target_db_file, \
full_list_target_genes)
print('\t%s genes will be used to build the background...'%len(full_list_target_genes))
print('Done!...')

##########################
##########################
### Generate  reference correlations ###
##########################
##########################
if generate_reference_correlations:
	print('Generating reference pairwise correlation matrix...')
	#0a. Create reference table
	assert os.path.exists(input_reference_db_file)
	assert os.path.getsize(input_reference_db_file)
	df_reference_rho_matrix = wrpr_process_exprssnDB(input_reference_db_file,l_genes_interest= \
	l_clock_genes_reference,outfile_full_rho_info=reference_full_rho_info,outfile_weighted_rho_matrix = \
	reference_weighted_rho_matrix)
else:
	try:
		assert os.path.exists(reference_weighted_rho_matrix)
	except:
		raise Exception('The input reference rho matrix "%s" does not exists!'% \
		reference_weighted_rho_matrix)
	#Read correlation matrx from file
	print('Reading reference pairwise correlation matrix...')
	df_reference_rho_matrix = pd.read_csv(reference_weighted_rho_matrix,index_col=0, \
	low_memory=False)
	print('Done!...')
#0b. Sort input df_reference_rho_matrix by reference genes
sqr_mtrx_db_reference = df_reference_rho_matrix[l_clock_genes_reference].T[l_clock_genes_reference]
#0c. Extract the upper triangle and exclude diagonal
mask = np.triu(np.ones(sqr_mtrx_db_reference.shape), k=1).astype(bool)
sqr_mtrx_db_reference = sqr_mtrx_db_reference.where(mask)

# Optional direct delta CCD workflow
if direct_delta_CCD:
	print('Running direct delta CCD workflow...')
	run_direct_delta_ccd(input_target_db_file=input_target_db_file, genes=l_clock_genes_target, contrast_variable=contrast_variable, group_a=group_a, group_b=group_b, sqr_mtrx_db_reference=sqr_mtrx_db_reference, sep='\t', n_permutations=delta_n_permutations, n_bootstraps=delta_n_bootstraps, seed=delta_seed, min_samples_per_group=delta_min_samples_per_group, outfile_results=delta_outfile_results, outfile_draws=delta_outfile_draws, verbose=True)
	print('Done!')
	sys.exit(0)

##########################
##########################
###    Generate  target correlations    ###
##########################
##########################
if generate_target_correlations:
	print('Generating target pairwise correlation matrix...')
	#1a. Create target table
	assert os.path.exists(input_target_db_file)
	assert os.path.getsize(input_target_db_file)
	rho_matrix_upper_triangle_target = wrpr_process_exprssnDB(input_target_db_file,l_genes_interest= \
	l_clock_genes_target,outfile_full_rho_info=target_full_rho_info,run_metanalysis= \
	run_metanalysis_target)
	#1b. In case of metanalysis in targets
	if run_metanalysis_target:
		#1ba. Convert reference table in pairwise correlation square matrix
		rho_matrix_upper_triangle_target = rtrn2clmn_from_corr_matrix(rho_matrix_upper_triangle_target)
else:
	try:
		assert os.path.exists(target_full_rho_info)
	except:
		raise Exception('The input target rho matrix "%s" does not exists!'%target_full_rho_info)
	#1c. Read correlation matrx from file
	print('Reading target pairwise correlation matrix...')
	rho_matrix_upper_triangle_target = pd.read_csv(target_full_rho_info,sep='\t',low_memory=False)
	#1d. In case of metanalysis in targets
	if run_metanalysis_target:
		rho_matrix_upper_triangle_target = rho_matrix_upper_triangle_target[['#Gene_A','Gene_B', \
		'weighted_rho_(rhoMeta)']]
		rho_matrix_upper_triangle_target = rho_matrix_upper_triangle_target.rename(columns= \
		{'weighted_rho_(rhoMeta)':'Spearman_rho'})
	#
	print('Done!...')

#################################
#################################
###    Generate  background target correlations    ###
#################################
#################################
if generate_background_target_correlations:
	print('Generating background target pairwise correlation matrix...')
	#1a. Create target table
	assert os.path.exists(input_target_db_file)
	assert os.path.getsize(input_target_db_file)
	#1a2. Revisit list of genes
	print('Revisiting list of background genes...')
	nOrigGns = len(set(full_list_target_genes))
	full_list_target_genes = revisit_gnsInterest(input_target_db_file,full_list_target_genes,verbose=False)
	print('\t%s genes out of %s were selected for the background after further filtering...'%(len(full_list_target_genes),nOrigGns))
	#1b. Create dataframe of random samples
	if background_matching_mode == 'matched':
		print('Building matched-null gene sets using expression and variance ranks...')
		df_matching_summary = build_target_gene_matching_summary(input_target_db_file,full_list_target_genes,sep='\t',verbose=False)
		missing_target_genes = sorted(set(l_clock_genes_target).difference(set(df_matching_summary.index)))
		if missing_target_genes:
			raise Exception('The following target genes are absent from the matched-null summary table: %s'%', '.join(missing_target_genes))
		if matched_gene_summary_outfile is not None:
			df_matching_summary.to_csv(matched_gene_summary_outfile,sep='\t')
		df_rndm_data = generate_matched_random_lists_of_genes(l_clock_genes_target,df_matching_summary,number_randomizations,seed=matching_seed,mean_caliper=matching_mean_caliper,variance_caliper=matching_variance_caliper,expand_step=matching_expand_step,max_caliper=matching_max_caliper)
		if matched_sets_outfile is not None:
			df_rndm_data.to_csv(matched_sets_outfile,sep='\t')
	else:
		df_rndm_data = generate_random_lists_of_genes(full_list_target_genes,number_randomizations, \
		size_of_sample=len(l_clock_genes_target),seed=matching_seed)
	#1c. Set dataframes for random lists of genes
	lRndmztns = df_rndm_data.columns.to_list()
	if not run_metanalysis_target:
		lDBs = [d for d in rho_matrix_upper_triangle_target.columns.to_list() if d not in {'#Gene_A', \
		'Gene_B'}]
		lDBsXlRndmztns = list(itertools.product(lDBs,lRndmztns))
		pd_db_ccd_bckgrnd = pd.DataFrame(columns=['|'.join(p) for p in lDBsXlRndmztns] ,index=['CCD'])
	else:
		pd_db_ccd_bckgrnd = pd.DataFrame(columns=[p for p in lRndmztns] ,index=['CCD'])
	#1d. Iterate over each randomized lists of genes, and initialize lists of 1000 gene interval.
	cnt,fll_randmzd_gns = 0,None
	for rndmztn in lRndmztns:		
		cnt+=1
		if fll_randmzd_gns is None or not cnt%100:
			print('\tComputing CCD for background target randomization %s...'%cnt)	
			#1e. Create a tmp pairwise correlations between al random genes. for the interval cnt, cnt+100
			fll_randmzd_gns = list(set(df_rndm_data[df_rndm_data.columns.to_list()[cnt-1:cnt+99]].to_numpy(). \
			flatten().tolist()))
			rho_matrix_upper_triangle_background_target_100= wrpr_process_exprssnDB(input_target_db_file, \
			l_genes_interest=fll_randmzd_gns,run_metanalysis=run_metanalysis_target,verbose=False)
		#1f. In case of metanalysis in targets
		if run_metanalysis_target:
			#1ga. Return a square matrix from the database of interest
			sqr_mtrx_db_background_target = rho_matrix_upper_triangle_background_target_100. \
			iloc[(rho_matrix_upper_triangle_background_target_100.index.isin(df_rndm_data[rndmztn])), \
			(rho_matrix_upper_triangle_background_target_100.columns.isin(df_rndm_data[rndmztn]))]
			#1gb. Extract the upper triangle and exclude diagonal
			mask = np.triu(np.ones(sqr_mtrx_db_background_target.shape), k=1).astype(bool)
			sqr_mtrx_db_background_target = sqr_mtrx_db_background_target.where(mask)
			#1gc. Compute Euclidean distance - CCD
			eucldn_d = cmptCCD(sqr_mtrx_db_reference,sqr_mtrx_db_background_target)
			#add
			pd_db_ccd_bckgrnd['rand_%s'%cnt] = eucldn_d
		else:
			for db in lDBs:
				rho_matrix_upper_triangle_background_target = rho_matrix_upper_triangle_background_target_100 \
				[(rho_matrix_upper_triangle_background_target_100['#Gene_A'].isin(df_rndm_data[rndmztn].to_list())) & \
				(rho_matrix_upper_triangle_background_target_100['Gene_B'].isin(df_rndm_data[rndmztn].to_list()))]
				#1ga. Return a square matrix from the database of interest
				sqr_mtrx_db_background_target = rtrn_corr_matrix_from_2clmn \
				(rho_matrix_upper_triangle_background_target,value_interest=db)
				#1gb. Extract the upper triangle and exclude diagonal
				mask = np.triu(np.ones(sqr_mtrx_db_background_target.shape), k=1).astype(bool)
				sqr_mtrx_db_background_target = sqr_mtrx_db_background_target.where(mask)
				#1gc. Compute Euclidean distance - CCD
				try:
					eucldn_d = cmptCCD(sqr_mtrx_db_reference,sqr_mtrx_db_background_target)
					pd_db_ccd_bckgrnd['|'.join([db,'rand_%s'%cnt])] = eucldn_d
				except:
					import pdb
					pdb.set_trace()
	#Output
	if background_target_full_rho_info is not None:
		pd_db_ccd_bckgrnd.T.to_csv(background_target_full_rho_info,sep='\t')
	#
else:
	try:
		assert os.path.exists(background_target_full_rho_info)
	except:
		raise Exception('The input background target rho matrix "%s" does not exists!'% \
		background_target_full_rho_info)
	#1c. Read correlation matrx from file
	print('Reading background target pairwise correlation matrix...')
	pd_db_ccd_bckgrnd = pd.read_csv(background_target_full_rho_info,sep='\t',low_memory=False, \
	index_col=0).T
##
print('Done!...')
#
##########################
##########################
###                Compute  CCD               ###
##########################
##########################
if compute_CCD:
	print('Computing CCD...')
	#0. Obtain lists of target databases
	if run_metanalysis_target:
		lDBs = ['Spearman_rho']
	#
	lDBs = sorted(set(rho_matrix_upper_triangle_target.columns).difference({'#Gene_A','Gene_B'}))
	pd_db_ccd = pd.DataFrame(columns=['CCD','CCD_pvalue','CCD_FDR'],index=lDBs)#Start pandas output
	#1. Compute CCD over each database
	for db in lDBs:
		print('\tComputing CCD for target %s...'%db)
		#1a. Return a square matrix from the database of interest
		sqr_mtrx_db_target = rtrn_corr_matrix_from_2clmn(rho_matrix_upper_triangle_target, \
		value_interest=db)
		#1b. Sort by target gene names
		sqr_mtrx_db_target = sqr_mtrx_db_target[l_clock_genes_target].T[l_clock_genes_target]
		#1c. Extract the upper triangle and exclude diagonal
		mask = np.triu(np.ones(sqr_mtrx_db_target.shape), k=1).astype(bool)
		sqr_mtrx_db_target = sqr_mtrx_db_target.where(mask)
		#1d. Compute Euclidean distance - CCD
		eucldn_d = cmptCCD(sqr_mtrx_db_reference,sqr_mtrx_db_target)
		pd_db_ccd.loc[db,'CCD'] = eucldn_d
		###
		#1e. Compute probability
		#Get background
		if run_metanalysis_target:
			assert pd_db_ccd_bckgrnd.shape[1] == number_randomizations
			random_ccds_bckgrnd = pd_db_ccd_bckgrnd.to_numpy().flatten()
		else:
			random_ccds_bckgrnd = pd_db_ccd_bckgrnd.loc[:,pd_db_ccd_bckgrnd.columns.str. \
			startswith('%s|'%db)].to_numpy().flatten()
			try:
				assert len(random_ccds_bckgrnd)
			except:
				raise Exception('The input background target randomization do not contain the database %s!'%db)
		#		
		# one-sided test
		r = np.sum(random_ccds_bckgrnd <= eucldn_d)
		n = len(random_ccds_bckgrnd)
		# Phipson & Smyth correction
		p_value = (r + 1) / (n + 1)
		pd_db_ccd.loc[db,'CCD_pvalue'] = p_value
	###
	#Conduct multiple test correction
	###
	pd_FDR = multipletests(pd_db_ccd['CCD_pvalue'], method='fdr_bh')[1]
	pd_db_ccd.loc[:,'CCD_FDR'] = pd_FDR
	###
	#Write output file
	###
	#Output
	pd_db_ccd.T.to_csv(outfile_results,sep='\t')
	#
	print('Done!')
