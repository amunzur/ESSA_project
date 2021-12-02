#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import statsmodels.stats.proportion as smp
import matplotlib.pyplot as plt
from scipy.stats import binom
import scipy.stats as stats

base_path = '/groups/wyattgrp/users/amunzur/essa/data'
output_path = '/groups/wyattgrp/users/amunzur/essa/output_data'
utilities_path = '/groups/wyattgrp/users/amunzur/essa/scripts/UTILITIES_data_consolidation.py'
exec(open(utilities_path).read()) # source functions 

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

snv_subset = add_TF_method(snv_subset) # add a col to indicate how the final TF estimate was computed
snv_subset = add_final_TF(snv_subset) # decide on the final TF value

cnvs = add_AR_abs_CN(cnvs, snv_subset) # once we have th TF_final, we can compute the abs AR cN value for AR gain and amplification cases

# Clean up the TF estimates, add a new col to indicate how the TF was computed: mutation, SNP or manual. 

# Save both tsv and csv files for ease of access
snvs.to_csv(os.path.join(output_path, "all_snvs_essa.tsv"), sep='\t', index=False, na_rep='NA')
snv_subset.to_csv(os.path.join(output_path, "ctdna_frac_variants.tsv"), sep='\t', index=False, na_rep='NA')
cnvs.to_csv(os.path.join(output_path, "all_cnvs_essa.tsv"), sep='\t', index=False, na_rep='NA')

snvs.to_csv(os.path.join(output_path, "all_snvs_essa.csv"), index=False, na_rep='NA')
snv_subset.to_csv(os.path.join(output_path, "ctdna_frac_variants.csv"), index=False, na_rep='NA')
cnvs.to_csv(os.path.join(output_path, "all_cnvs_essa.csv"), index=False, na_rep='NA')