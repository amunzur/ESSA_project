#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import statsmodels.stats.proportion as smp
import matplotlib.pyplot as plt
from scipy.stats import binom
import scipy.stats as stats

# =============================================================================
# Define functions
# =============================================================================
def ctDNA_fraction_het_SNPs(cnv_path, THRESHOLD_copy_status, THRESHOLD_het_SNPs_count, THRESHOLD_het_loss_count, THRESHOLD_chrom_count): 
	'''
	Given a df of copy numbers, calculate the ctDNA% based on the heterozygous SNPs and generate a df with copy number status. 
	cnvs: copy number df
	THRESHOLD_copy_status = -1: subset to samples with this CN change, usually a LOH (== -1)
	THRESHOLD_het_SNPs_count = 4: the minimum number of SNPs a sample should have in order to do ctDNA% calculation based on het SNPs. (> 4)
	THRESHOLD_het_loss_count = 1: number of genes where a sample has a het loss (> 1)
	THRESHOLD_chrom_count = 1: the minimum number of distinct chroms where the sample should have a het loss (> 1)
	'''
	
	# Assign a copy status to the samples based on the log ratio
	cnvs_main = pd.read_csv(cnv_path, index_col=None, sep='\t', dtype='str')
	cnvs = cnvs_main.apply(pd.to_numeric, errors='ignore')
	cnvs = cnvs.rename(columns={'Gene': 'GENE', 'Chrom':'CHROMOSOME', 'Start':'START', 'End':'END', 'Coverage logratio':"Log_ratio", 'Sample':"Sample_ID"})

	# Add a new copy_status column based on the log ratio values and SNP allele fractions
	conditions = [
		(cnvs["Log_ratio"] < -0.7), 
		(cnvs['Log_ratio'] <= -0.20) & (cnvs['Median heterozygous SNP allele fraction'] >= 0.60), 
		(cnvs['Log_ratio'] > -0.3) & (cnvs['Log_ratio'] < 0.3), 
		(cnvs['Log_ratio'] >= 0.20) & (cnvs['Median heterozygous SNP allele fraction'] >= 0.60), 
		(cnvs['Log_ratio'] > 0.7)]

	values = [-2, -1, 0, 1, 2] # corresponding copy numbers to the conditions above
	cnvs["Copy_status"] = np.select(conditions, values)
	cnvs['Copy_status'] = cnvs['Copy_status'].fillna(0) # anything else that doesn't fit the criteria above gets a 0.

	# sex chromosomes get special treatment for CN calculation. We don't use SNP allele fraction for them. 
	conditions = [
		(cnvs.loc[cnvs.CHROMOSOME.isin(["chrX", "chrY"]), 'Log_ratio'] < -0.7), 
		(cnvs.loc[cnvs.CHROMOSOME.isin(["chrX", "chrY"]), 'Log_ratio'] <= -0.20), 
		(cnvs.loc[cnvs.CHROMOSOME.isin(["chrX", "chrY"]), 'Log_ratio'] > -0.3) & (cnvs.loc[cnvs.CHROMOSOME.isin(["chrX", "chrY"]), 'Log_ratio'] < 0.3), 
		(cnvs.loc[cnvs.CHROMOSOME.isin(["chrX", "chrY"]), 'Log_ratio'] >= 0.20), 
		(cnvs.loc[cnvs.CHROMOSOME.isin(["chrX", "chrY"]), 'Log_ratio'] > 0.7)]

	values = [-2, -1, 0, 1, 2] # corresponding copy numbers to the conditions above
	cnvs.loc[cnvs.CHROMOSOME.isin(["chrX", "chrY"]), "Copy_status"] = np.select(conditions, values)
	cnvs.loc[cnvs.CHROMOSOME.isin(["chrX", "chrY"]),'Copy_status'] = cnvs['Copy_status'].fillna(0) # anything else that doesn't fit the criteria above gets a 0.

	cnvs['Patient_ID'] = cnvs["Sample_ID"].str.split("-cfDNA|-WBC", expand = True)[0] # patient id from sample id

	# Make a copy of the df to compute ctDNA% on applicable samples
	cnvs_ctdna = cnvs.copy()
	cnvs_ctdna = cnvs_ctdna[cnvs_ctdna['Copy_status'] == THRESHOLD_copy_status]
	cnvs_ctdna['Hetz_loss_count'] = cnvs_ctdna['GENE'].groupby(cnvs_ctdna['Sample_ID']).transform('count') # count the number of genes that passed the criteria above for each patient
	cnvs_ctdna_copy = cnvs_ctdna.copy() # make a copy to return and do a merge later on
	cnvs_ctdna['Chrom_count'] = cnvs_ctdna['CHROMOSOME'].groupby(cnvs_ctdna['Sample_ID']).transform('count')
	cnvs_ctdna['Sample_median_gene_SNP_VAF_all_genes'] = cnvs_ctdna['Median heterozygous SNP allele fraction'].groupby(cnvs_ctdna['Sample_ID']).transform('median') # Get "median of medians" per sample_ID

	# Filter and calculate ctDNA% based on SNPs 
	foo=(cnvs_ctdna["Number of heterozygous SNPs"] > THRESHOLD_het_SNPs_count) & (cnvs_ctdna["Hetz_loss_count"] > THRESHOLD_het_loss_count) & (cnvs_ctdna['Chrom_count'] > THRESHOLD_chrom_count)
	cnvs_ctdna = cnvs_ctdna[foo] # filtering
	cnvs_ctdna['ctDNA_fraction_SNP_estimate'] = 2 - 1/cnvs_ctdna['Sample_median_gene_SNP_VAF_all_genes'] #Calculate tumour content from HSAFs in sLOH regions

	# merge the chosen columns back to the original df
	cnvs_ctdna = cnvs_ctdna[['ctDNA_fraction_SNP_estimate', 'Hetz_loss_count', 'Sample_ID']].drop_duplicates()
	cnvs = pd.merge(cnvs, cnvs_ctdna, how='left', on=['Sample_ID'])

	# move the Patient_ID col right after the Sample_ID col, pop it off and put it back to the correct position
	x = cnvs.pop("Patient_ID")
	cnvs.insert(1, 'Patient_ID', x)

	x = cnvs_ctdna_copy.pop("Patient_ID")
	cnvs_ctdna_copy.insert(1, 'Patient_ID', x)

	return([cnvs, cnvs_ctdna_copy])

def generate_SNV_table(snv_path, cnvs, sample_labels, keep_only_called_variants, keep_only_coding_variants): 
	'''
	# Clean up and modify the SNVs table.
	sample_labels: Unique labels taken from the cnvs file. This ensures that we have the same list of samples in cnv and snv tables. 
	'''
	snvs_main = pd.read_csv(snv_path, index_col = None, sep = '\t', dtype='str')
	snvs = snvs_main.apply(pd.to_numeric, errors = 'ignore')
	snvs = snvs.drop(labels = "CHROM.1", axis = 1)

	snvs = pd.melt(snvs, id_vars = ['CHROM', 'POSITION', 'REF', 'ALT', 'GENE', 'EFFECT', 'NOTES']) # wide to long format
	snvs = snvs.rename(columns = {'value': 'Allele_frequency','variable': 'Sample_ID'}, inplace = False)
	snvs = pd.merge(snvs, sample_labels, how = 'left', on = 'Sample_ID') # so far no filtering, all muts included

	# Subset to called variants only or keep everything
	if keep_only_called_variants == True: 
		snvs = snvs[snvs['Allele_frequency'].str.contains("\*")]
	else: 
		snvs.loc[snvs['Allele_frequency'].str.contains('\*'), 'Independently_detected'] = True
		snvs['Independently_detected'] = snvs['Independently_detected'].fillna(False)

	snvs[['Mutation_type', 'Protein_annotation']] = snvs['EFFECT'].str.split(" ", expand = True).iloc[:, 0:2]
	snvs['Protein_annotation'] = snvs['Protein_annotation'].str.replace('site', 'Splice site', regex = False)
	snvs['Protein_annotation'] = snvs['Protein_annotation'].str.replace(',', '').replace('.', '')

	# snvs['Mutation_type'] = snvs['EFFECT'].str.split(" ", expand=True)[0]
	snvs['Mutation_type'] = snvs['Mutation_type'].str.replace('Splice','Splice site', regex=False)

	# Boolean column for whether mutation is coding or noncoding
	snvs['Coding'] = snvs['EFFECT'].str.contains('frameshift|missense|synonymous|splice|stopgain|nonsense', case = False)

	# Sometimes a gene name is repeated twice with a semicolon in between which causes issues with merging with the cnvs table.
	snvs['GENE'] = snvs['GENE'].str.replace("(;.*)", "", regex = True) # Replace it with the gene name
	
	if keep_only_coding_variants == True: 
		snvs = snvs[snvs['Coding'] == True] # subset to coding variants
	else: 
		pass

	snvs['Putatively_deleterious'] = snvs['EFFECT'].str.contains('frameshift|splice|stopgain|nonsense', case = False)
	snvs['Read_depth'] = snvs['Allele_frequency'].str.split(" ", expand = True)[1].str.extract('.*\((.*)\).*').apply(pd.to_numeric, errors='ignore') # remove paranthesis and convert the col to numeric using apply
	snvs['Allele_frequency'] = pd.to_numeric(snvs['Allele_frequency'].str.split("%", expand = True)[0])*0.01 # percentage to numeric
	snvs['Mutant_reads'] = snvs['Allele_frequency'] * snvs['Read_depth'].round()
	snvs['Patient_ID'] = snvs["Sample_ID"].str.split("-cfDNA|-WBC", expand = True)[0] # patient id from sample id

	# Add CNV information to SNVs - nice to have some data redundancy so you don't have to merge files constantly
	snvs = pd.merge(snvs, cnvs[['Sample_ID', 'GENE', 'Copy_status', 'Log_ratio']], how = 'left', on = ['Sample_ID', 'GENE'])
	snvs['Lower_95_CI_VAF'], snvs['Upper_95_CI_VAF'] = smp.proportion_confint(snvs['Mutant_reads'], snvs['Read_depth'], alpha = 0.05, method = 'normal')
	snvs['Mutation_VAF_STD'] = binom.std(snvs['Read_depth'], snvs['Allele_frequency'])
	snvs['Mutation_VAF_STD'] = snvs['Mutation_VAF_STD'] / snvs['Read_depth']

	snvs['Adj_allele_frequency'] = snvs.apply(lambda row: get_adj_maf(row['Allele_frequency'], row['Read_depth']), axis = 1) # also add a new col for the adjusted allele frequency

	return(snvs)

def ctDNA_fraction_mutations(snvs, THRESHOLD_log_ratio, THRESHOLD_read_depth, keep_only_coding_variants): 
	'''
	# Filter the df with all snvs to retain only the applicable ones to compute ctDNA fraction
	THRESHOLD_log_ratio: only keep the .... with log ratio less than (0.2)
	THRESHOLD_read_depth: only keep the variants with depth more than (40)
	'''
	snv_subset = snvs.copy()
	snv_subset2 = snvs.copy()
	
	# We have the option to compute ctDNA% based on coding variants only. This may or may not make any difference based on what the muts with the highest vaf are.
	if keep_only_coding_variants == True: 
		snv_subset = snv_subset[snv_subset['Coding'] == True] # subset to coding variants
		snv_subset2 = snv_subset2[snv_subset2['Coding'] == True] # subset to coding variants
	else: 
		pass
	
	# This df contains ctDNA fraction estimates using mutations on the X chromosome, considering only genes with a copy status of -1, 1 or 0.
	snv_subset2 = snv_subset2[(snv_subset2['CHROM'] == "chrX") & 
							(snv_subset2['Read_depth'] > THRESHOLD_read_depth) & 
							(snv_subset2.Copy_status.isin([-1, 0])) & 
							(snv_subset2['Independently_detected'] == True)].reset_index(drop = True)

	snv_subset2["Max_Sample_ID_adj_VAF"] = snv_subset2.groupby('Sample_ID')["Adj_allele_frequency"].transform('max') # for each sample the vaf of the mutation with the highest vaf
	snv_subset2 = snv_subset2[snv_subset2["Adj_allele_frequency"] == snv_subset2["Max_Sample_ID_adj_VAF"]] # only keep the muts with the highest vaf
	snv_subset2['row_max'] = snv_subset2[['Read_depth']].max(axis = 1) # helps break tie and subset if there are more than one muts with the same adjusted vaf
	snv_subset2 = snv_subset2.loc[snv_subset2.groupby('Sample_ID')['row_max'].idxmax()] # only keep the muts with the highest vaf

	snv_subset2['Mutation_ctDNA_fraction'] = snv_subset2['Max_Sample_ID_adj_VAF'] # ctDNA% estimate using mutation calls in the sex chromosomes

	############################################################################################################
	
	snv_subset["Max_Sample_ID_adj_VAF"] = snv_subset.groupby('Sample_ID')["Adj_allele_frequency"].transform('max') # for each sample the vaf of the mutation with the highest vaf
	snv_subset = snv_subset[snv_subset["Adj_allele_frequency"] == snv_subset["Max_Sample_ID_adj_VAF"]] # only keep the muts with the highest vaf
	snv_subset['row_max'] = snv_subset[['Read_depth']].max(axis = 1) # helps break tie and subset if there are more than one muts with the same adjusted vaf
	snv_subset = snv_subset.loc[snv_subset.groupby('Sample_ID')['row_max'].idxmax()] # only keep the muts with the highest vaf

	# Continue filtering snv_subset to choose the right mutations for TF calculation.	
	snv_subset = snv_subset[(snv_subset['Log_ratio'] < THRESHOLD_log_ratio) & (snv_subset['Read_depth'] > THRESHOLD_read_depth) & (snv_subset['Independently_detected'] == True) & (~snv_subset['CHROM'].isin(['chrX', 'chrY']))].reset_index(drop = True)

	snv_subset['Mutation_ctDNA_fraction'] = 2/(1+1/(snv_subset['Max_Sample_ID_adj_VAF'])) # ctDNA% estimate using mutation calls
	snv_subset = snv_subset[snv_subset['Log_ratio'].notna()]

	############################################################################################################

	# At this point we have two dfs with TF calculated based on two sets of mutations. We consider the sex chromosome calculations only if the somatic mutation calculations aren't available. 
	# Both dfs are returned to be processed later on. 
	# snv_subset: somatic mutations 
	# snv_subset2: X chromosome mutations

	return([snv_subset, snv_subset2])

def add_TF_method(snv_subset): 
	'''
	Add a new col to the snv_subset dataframe to indicate how the TF estimate was computed. 
	Assuming some merge was done on the snv_subset so that it has both mutation based and SNP based TF estimates. 
	'''

	col1         = "Mutation_ctDNA_fraction"
	col2         = "ctDNA_fraction_SNP_estimate"

	conditions   = [pd.isna(snv_subset[col1]) & pd.notna(snv_subset[col2]), # Mutation estimate unavailable: choose SNP estimate
					pd.isna(snv_subset[col2]) & pd.notna(snv_subset[col1]), # SNP estimate unavailable: choose mutation estimate
					pd.isna(snv_subset[col1]) & pd.isna(snv_subset[col2]), # both NA, leave NA
					pd.notna(snv_subset[col1]) & pd.notna(snv_subset[col2])] # neither is NA: choose mutation

	choices      = [ "SNP", "Mutation", "NA", "Mutation"]
	snv_subset["TF_method"] = np.select(conditions, choices)

	return(snv_subset)

def add_final_TF(snv_subset): 
	'''
	Add a new col to the snv_subset dataframe to indicate the final TF estimate. 
	Same criteria as add_TF_method(snv_subset) is used.
	'''

	col1 = "Mutation_ctDNA_fraction"
	col2 = "ctDNA_fraction_SNP_estimate"
	
	snv_subset["Final_TF"] = (snv_subset[col1].fillna(0) + snv_subset[col2].fillna(0)).replace(0, np.nan)

	return(snv_subset)

#Adjust for 95% quantile outlier (when using max somatic VAF to calculate ctDNA%)
def get_adj_maf(maf, depth):
    alt = maf*depth
    for p in np.arange(0,1,0.005):
        dist = stats.binom(depth, p)
        if dist.cdf(alt) < 0.95:
            maf = p; break;
    return maf;


# =============================================================================
# Load raw data
# =============================================================================

base_path = '/groups/wyattgrp/users/amunzur/essa/data'
output_path = '/groups/wyattgrp/users/amunzur/essa/output_data'

# Calculate ctDNA fraction from heterozygous SNPs and generate a processed copy number changes and log ratio table
cnv_path = os.path.join(base_path, "gene_cna.tsv")
cnvs = ctDNA_fraction_het_SNPs(cnv_path = cnv_path, THRESHOLD_copy_status = -1, THRESHOLD_het_SNPs_count = 4, THRESHOLD_het_loss_count = 1, THRESHOLD_chrom_count = 1)[0]
cnvs_ctdna_copy = ctDNA_fraction_het_SNPs(cnv_path = cnv_path, THRESHOLD_copy_status = -1, THRESHOLD_het_SNPs_count = 4, THRESHOLD_het_loss_count = 1, THRESHOLD_chrom_count = 1)[1]

# Generate SNVs table 
sample_labels = cnvs[['Sample_ID']].drop_duplicates()
snv_path = os.path.join(base_path, "essa_mutations.tsv")
snvs = generate_SNV_table(snv_path = snv_path, cnvs = cnvs, sample_labels = sample_labels, keep_only_called_variants = False, keep_only_coding_variants = False)

# Select the appropriate snvs to calculate tumor content based on mutations, we just computed the SNP based on ctDNA% estimates based on the cnvs file
snv_subset = ctDNA_fraction_mutations(snvs = snvs, THRESHOLD_log_ratio = 0.2, THRESHOLD_read_depth = 40, keep_only_coding_variants = True)[0] # TF based on somatic muts
snv_subset2 = ctDNA_fraction_mutations(snvs = snvs, THRESHOLD_log_ratio = 0.2, THRESHOLD_read_depth = 40, keep_only_coding_variants = True)[1] # TF based on sex chromosomes

snv_subset = pd.merge(snv_subset, cnvs[["Sample_ID", "Patient_ID", "ctDNA_fraction_SNP_estimate"]].drop_duplicates(), how = "outer", on = ["Sample_ID", 'Patient_ID']) # add SNP estimate, this also helps understand which samples have an NA for a TF estimate
snv_subset = pd.merge(snv_subset, cnvs_ctdna_copy[['Hetz_loss_count', 'Sample_ID', 'Patient_ID']].drop_duplicates(), how='outer', on = ['Sample_ID', "Patient_ID"]) # add number of genes with het loss
# At this point now we have information about which samples have a TF estimate based on somatic mutations. We will make use of snv_subset2 to compute a TF for samples without a TF estimate based on somatic mutations. 

idx = snv_subset2["Sample_ID"].isin(snv_subset[snv_subset['Mutation_ctDNA_fraction'].isnull()].Sample_ID) # index of sample IDs in snv_subset2 that doesn't have a somatic mutation based TF estimate 
snv_subset2 = snv_subset2[idx] # subset to samples whose TF will be calculated based on sex chr calculations

# Final df that contains the TF estimates calculated by somatic and sex chr estimates. 
snv_subset = pd.concat([snv_subset, snv_subset2]).drop_duplicates(subset = "Sample_ID", keep='last').reset_index(drop = True) # slap the two dfs together, snv_subset contains TF based on somatic mutations, snv_subset2 contains TF based on sex chr mutations 

snvs = pd.merge(snvs, snv_subset[['Max_Sample_ID_adj_VAF', 'Mutation_ctDNA_fraction', 'Sample_ID']].drop_duplicates(), how='left', on='Sample_ID') # Merge snv_subset and snv to add ctDNA% estimates based on mutation counts to the snvs table

# subset to chosen columns in this order
reordered_columns = ['Sample_ID','Patient_ID','CHROM', 'POSITION', 'REF', 'ALT',
       'GENE', 'Protein_annotation', 'Mutation_type', 'NOTES', 'Allele_frequency', 'Independently_detected', 'Coding', 
       'Putatively_deleterious', 'Read_depth', 'Mutant_reads', 'Copy_status',
       'Log_ratio', 'Lower_95_CI_VAF', 'Upper_95_CI_VAF', 'Max_Sample_ID_adj_VAF', 'Mutation_ctDNA_fraction']

snvs = snvs[reordered_columns]
snv_subset = snv_subset[reordered_columns + ['ctDNA_fraction_SNP_estimate']]

snv_subset = add_TF_method(snv_subset) # add a col to indicate 
snv_subset = add_final_TF(snv_subset)

# Clean up the TF estimates, add a new col to indicate how the TF was computed: mutation, SNP or manual. 

# Save both tsv and csv files for ease of access
snvs.to_csv(os.path.join(output_path, "all_snvs_essa.tsv"), sep='\t', index=False, na_rep='NA')
snv_subset.to_csv(os.path.join(output_path, "ctdna_frac_variants.tsv"), sep='\t', index=False, na_rep='NA')
cnvs.to_csv(os.path.join(output_path, "all_cnvs_essa.tsv"), sep='\t', index=False, na_rep='NA')

snvs.to_csv(os.path.join(output_path, "all_snvs_essa.csv"), index=False, na_rep='NA')
snv_subset.to_csv(os.path.join(output_path, "ctdna_frac_variants.csv"), index=False, na_rep='NA')
cnvs.to_csv(os.path.join(output_path, "all_cnvs_essa.csv"), index=False, na_rep='NA')