#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  preprocess.py
#  
#  Copyright 2025 oscar <oscar@oscar-HP-EliteBook-840-G5>
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

import os
import pandas as pd
import gzip
from io import StringIO

# 0. Get databases
get_databases = False
# 1. Convert codes
convert_codes = True
# 2. Convert tables
convert_tables = True

#####################
#####################
# 0. Get databases
##################### 
#####################
if get_databases:
	#GSE35795	Negoro et al. (2012)	Agilent-014868 4x44K	wt	bladder	12	DD
	os.system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE35nnn/GSE35795/matrix/GSE35795_series_matrix.txt.gz')
	#GSE38625	Geyfman et al. (2012)	Affy Mouse 1.0 ST	wt	skin	13	LD 12:12
	os.system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE38nnn/GSE38625/matrix/GSE38625_series_matrix.txt.gz')
	"""
	Dataset includes samples from anagen and telogen phases, we only used those from anagen.
	For GSE38625_series_matrix.txt.gz:
	GSM946294 	ZT2 whole skin anagen P30
	GSM946295 	ZT6 whole skin anagen P30
	GSM946296 	ZT10 whole skin anagen P30
	GSM946297 	ZT14 whole skin anagen P30
	GSM946298 	ZT18 whole skin anagen P30
	GSM946299 	ZT22 whole skin anagen P30
	GSM946300 	ZT26 whole skin anagen P30
	GSM946301 	ZT30 whole skin anagen P30
	GSM946302 	ZT34 whole skin anagen P30
	GSM946303 	ZT38 whole skin anagen P30
	GSM946304 	ZT42 whole skin anagen P30
	GSM946305 	ZT46 whole skin anagen P30
	GSM946306 	ZT50 whole skin anagen P30
	"""
	#GSE4238	Oster et al. (2006)	Affy Mouse 430 2.0	wt	adrenal gland	24	DD
	os.system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE4nnn/GSE4238/matrix/GSE4238_series_matrix.txt.gz')
	#GSE10644	Hoogerwerf et al. (2008)	Affy Mouse 430 2.0	wt	colon	17	LD 12:12
	os.system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE10nnn/GSE10644/matrix/GSE10644_series_matrix.txt.gz')
	#GSE11923	Hughes et al. (2009)	Affy Mouse 430 2.0	wt	liver	48	DD
	os.system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE11nnn/GSE11923/matrix/GSE11923_series_matrix.txt.gz')
	#GSE54650	Zhang et al. (2014)	Affy Mouse 1.0 ST	wt	brain stem, brown fat	24 each organ	DD
	os.system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE54nnn/GSE54650/matrix/GSE54650_series_matrix.txt.gz')
	"""
	Dataset includes data from 12 organs, we only used data from two.
	For GSE54650_series_matrix.txt.gz:
	-Brown fat-
	GSM1321062 	BFat_CT18
	GSM1321063 	BFat_CT20
	GSM1321064 	BFat_CT22
	GSM1321065 	BFat_CT24
	GSM1321066 	BFat_CT26
	GSM1321067 	BFat_CT28
	GSM1321068 	BFat_CT30
	GSM1321069 	BFat_CT32
	GSM1321070 	BFat_CT34
	GSM1321071 	BFat_CT36
	GSM1321072 	BFat_CT38
	GSM1321073 	BFat_CT40
	GSM1321074 	BFat_CT42
	GSM1321075 	BFat_CT44
	GSM1321076 	BFat_CT46
	GSM1321077 	BFat_CT48
	GSM1321078 	BFat_CT50
	GSM1321079 	BFat_CT52
	GSM1321080 	BFat_CT54
	GSM1321081 	BFat_CT56
	GSM1321082 	BFat_CT58
	GSM1321083 	BFat_CT60
	GSM1321084 	BFat_CT62
	GSM1321085 	BFat_CT64
	-Brain stem-
	GSM1321038 	Bstm_CT18
	GSM1321039 	Bstm_CT20
	GSM1321040 	Bstm_CT22
	GSM1321041 	Bstm_CT24
	GSM1321042 	Bstm_CT26
	GSM1321043 	Bstm_CT28
	GSM1321044 	Bstm_CT30
	GSM1321045 	Bstm_CT32
	GSM1321046 	Bstm_CT34
	GSM1321047 	Bstm_CT36
	GSM1321048 	Bstm_CT38
	GSM1321049 	Bstm_CT40
	GSM1321050 	Bstm_CT42
	GSM1321051 	Bstm_CT44
	GSM1321052 	Bstm_CT46
	GSM1321053 	Bstm_CT48
	GSM1321054 	Bstm_CT50
	GSM1321055 	Bstm_CT52
	GSM1321056 	Bstm_CT54
	GSM1321057 	Bstm_CT56
	GSM1321058 	Bstm_CT58
	GSM1321059 	Bstm_CT60
	GSM1321060 	Bstm_CT62
	GSM1321061 	Bstm_CT64
	"""
	#GSE59396	Haspel et al. (2014)	Illumina MouseRef-8 v2.0	wt	lung	36	LD 12:12
	os.system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59396/matrix/GSE59396_series_matrix.txt.gz')
	"""
	Dataset includes both samples from DD and LD, we only used those from LD.
	For GSE59396_series_matrix.txt.gz:
	GSM1435584	T24L1MR
	GSM1435586	T8L4MR
	GSM1435588	T12L5MR
	GSM1435589	T40L5MR
	GSM1435593	T20L3MR
	GSM1435596	T0L5MR
	GSM1435598	T32L3MR
	GSM1435599	T40L1MR
	GSM1435601	T28L3MR
	GSM1435602	T20L4MR
	GSM1435603	T36L3MR
	GSM1435606	T12L3MR
	GSM1435608	T12L2MR
	GSM1435613	T16L5MR
	GSM1435614	T8L3MR
	GSM1435616	T36L2MR
	GSM1435618	T4L5MR
	GSM1435619	T44L2MR
	GSM1435621	T32L5MR
	GSM1435624	T16L1MR
	GSM1435629	T28L5MR
	GSM1435631	T0L1MR
	GSM1435633	T24L2MR
	GSM1435635	T4L2MR
	GSM1435638	T16L4MR
	GSM1435639	T8L2MR
	GSM1435641	T44L3MR
	GSM1435645	T24L3MR
	GSM1435646	T40L4MR
	GSM1435649	T0L2MR
	GSM1435651	T28L2MR
	GSM1435653	T36L5MR
	GSM1435655	T32L2MR
	GSM1435656	T20L2MR
	GSM1435661	T44L4MR
	GSM1435663	T4L1MR
	"""

#####################
#####################
# 1. Convert codes
##################### 
#####################
if convert_codes:
	oneToOne=False
	oneToMany=True
	#
	##
	#Convert codes
	#Get data from BioMart
	"""
	In https://www.ensembl.org/biomart/martview/de9b7c42eb9e4ea27dfec41848a61870
	Select:
	Dataset
	Mouse genes (GRCm39)
	Filters
	Transcript type: protein_coding 
	Attributes
	Gene name
	AFFY Mouse430 2 probe
	AFFY MoGene 1 0 st v1 probe
	*Saved as "GRCm39_array1.tsv.gz"

	Select:
	Dataset
	Mouse genes (GRCm39)
	Filters
	Transcript type: protein_coding 
	Attributes
	Gene name
	ILLUMINA MouseRef 8 probe
	AGILENT WholeGenome 4x44k v1 probe
	*Saved as "GRCm39_array2.tsv.gz"
	"""
	dArrayDPrbGnm={}#Store folder

	for cmprssDB in ['GRCm39_array1.tsv.gz','GRCm39_array2.tsv.gz']:
		#
		# Replace 'your_file.tsv.gz' with your actual file path
		df = pd.read_csv(cmprssDB, sep='\t', compression='gzip', header=0)
		#Make dictionaries
		for arrayClmn in list(df.columns[1:]):
			if oneToOne:
				# Check for multiple mappings from gene names to arrayClmn
				col1_to_col2 = df.groupby('Gene name')[arrayClmn].nunique()
				col1_consistent = col1_to_col2[col1_to_col2 == 1].index
				# Check for multiple mappings from arrayClmn to gene names
				col2_to_col1 = df.groupby(arrayClmn)['Gene name'].nunique()
				col2_consistent = col2_to_col1[col2_to_col1 == 1].index
				#Filter rows
				filtered_df = df[df['Gene name'].isin(col1_consistent) & df[arrayClmn].isin(col2_consistent)]
				filtered_df = filtered_df[['Gene name',arrayClmn]]
				#Ensuring One-to-One Mappings
				assert (filtered_df.groupby('Gene name')[arrayClmn].nunique()==1).all()
				mapping_dict = dict(zip(filtered_df[arrayClmn], filtered_df['Gene name']))
				####
				#Test genes
				mapping_dict = dict(zip(filtered_df['Gene name'], filtered_df[arrayClmn]))
				[mapping_dict[gn] for gn in ['Bmal1','Clock','Cry1','Cry2','Dbp','Npas2','Nr1d1','Nr1d2','Per1','Per2','Per3','Tef']]
				####
			else:
				assert oneToMany
				# Check for multiple mappings from gene names to arrayClmn
				col2_to_col1_counts = df.groupby(arrayClmn)['Gene name'].nunique()
				# Identify arrayClmn values that map to more than one gene name
				invalid_col2 = col2_to_col1_counts[col2_to_col1_counts > 1].index
				# Filter out rows with arrayClmn values that map to multiple gene name values
				filtered_df = df[~df[arrayClmn].isin(invalid_col2)]
				filtered_df.dropna(subset=['Gene name',arrayClmn], inplace=True)
				#
				mapping_dict = dict()
				for ky,vl in zip(filtered_df[arrayClmn],filtered_df['Gene name']):
					if ky in mapping_dict:
						assert mapping_dict[ky]==vl
					else:
						mapping_dict[ky]=vl
				#Test genes
				# ~ [(prb,mapping_dict[prb]) for prb in mapping_dict if mapping_dict[prb] in ['Bmal1','Clock','Cry1','Cry2','Dbp','Npas2','Nr1d1','Nr1d2','Per1','Per2','Per3','Tef']]
				dArrayDPrbGnm[arrayClmn] = mapping_dict


#####################
#####################
# 2. Convert tables
##################### 
#####################
if convert_tables:
	dPltfrmGEOfl = {'AFFY Mouse430 2 probe':['GSE4238_series_matrix.txt.gz','GSE10644_series_matrix.txt.gz','GSE11923_series_matrix.txt.gz'], \
	'AFFY MoGene 1 0 st v1 probe':['GSE38625_series_matrix.txt.gz','GSE54650_series_matrix.txt.gz'], \
	'ILLUMINA MouseRef 8 probe':['GSE59396_series_matrix.txt.gz'], \
	'AGILENT WholeGenome 4x44k v1 probe':['GSE35795_series_matrix.txt.gz']}
	#
	for pltfm,l_GEOfl in dPltfrmGEOfl.items():
		dPrbGnNm = dArrayDPrbGnm[pltfm]
		#Iterate over each GEOfl
		for GEOfl in l_GEOfl:
			#Open input file			
			valid_lines = []
			with gzip.open(GEOfl, 'rt', encoding='utf-8') as file:
				for line in file:
					num_commas = line.count('\t')
					if line[0]!='!' and num_commas > 2:
						valid_lines.append(line)
			#Process as pandas
			data_io = StringIO('\n'.join(valid_lines))
			df_GEOfl = pd.read_csv(data_io, sep='\t',header=0,index_col=0)
			#Add gene symbols
			df_GEOfl['#gene_symbol'] = df_GEOfl.index.map(dPrbGnNm)
			# Exclude rows where '#gene_symbol' is NaN (i.e., 'probe_id' not in the dictionary)
			df_GEOfl = df_GEOfl.dropna(subset=['#gene_symbol'])
			#Compute expression median
			df_GEOfl['median_expr'] = df_GEOfl.iloc[:,~(df_GEOfl.columns=='#gene_symbol')].median(axis=1)
			# Identify the index of the row with the maximum 'median_expr' for each '#gene_symbol'
			idx = df_GEOfl.groupby('#gene_symbol')['median_expr'].idxmax()
			# Add probe to file
			df_GEOfl['probe_id'] = df_GEOfl.index
			# Select these rows to create the final DataFrame
			df_GEOfl_top_probes = df_GEOfl.loc[idx].reset_index(drop=True)
			df_GEOfl_top_probes.index = df_GEOfl_top_probes['#gene_symbol']
			df_GEOfl_top_probes.drop(columns=['#gene_symbol'], inplace=True)
			#Test
			df_GEOfl_top_probes[df_GEOfl_top_probes.index.isin(['Bmal1','Clock','Cry1','Cry2','Dbp','Npas2','Nr1d1','Nr1d2','Per1','Per2','Per3','Tef'])]
			#Reorganize dataframe and write in csv
			top_probes = df_GEOfl_top_probes.pop('probe_id')
			df_GEOfl_top_probes.insert(0, 'probe_id', top_probes)
			df_GEOfl_top_probes.to_csv(GEOfl.replace('_series_matrix.txt.gz','_filtered_plus_wGeneNames.tsv'),sep='\t')#,compression='gzip')
	#
		
