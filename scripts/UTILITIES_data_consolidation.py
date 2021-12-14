def reindex_column(df, colname, new_position):
	'''
	Move a col to a new place
	'''
	x = df.pop(colname)
	df.insert(new_position, colname, x)

	return(df)


def ctDNA_fraction_het_SNPs(cnv_path, THRESHOLD_copy_status, THRESHOLD_het_SNPs_count, THRESHOLD_het_loss_count, THRESHOLD_chrom_count): 
	'''
	Given a of copy numbers, calculate the ctDNA% based on the heterozygous SNPs and generate a df with copy number status. 
	cnvs: copy number df
	THRESHOLD_copy_status = -1: subset to samples with this CN change, usually a LOH (== -1)
	THRESHOLD_het_SNPs_count = 4: the minimum number of SNPs a sample should have in order to do ctDNA% calculation based on het SNPs. (> 4)
	THRESHOLD_het_loss_count = 1: number of genes where a sample has a het loss (> 1)
	THRESHOLD_chrom_count = 1: the minimum number of distinct chroms where the sample should have a het loss (> 1)
	'''
	
	# Assign a copy status to the samples based on the log ratio
	cnvs_main = pd.read_csv(cnv_path, index_col=None, sep='\t', dtype='str')
	cnvs = cnvs_main.apply(pd.to_numeric, errors='ignore')
	cnvs = cnvs.rename(columns={'Gene': 'GENE', 'Chrom':'CHROMOSOME', 'Start':'START', 'End':'END', 'Coverage logratio':"Log_ratio", 'Sample':"Sample_ID", "Median heterozygous SNP allele fraction": "msf"})

	# Fix some naming issues, this sample should have been C1D1
	cnvs.loc[cnvs["Sample_ID"] == "ESSA-04-002-cfDNA-1A-C1D28-2020Dec08", "Sample_ID"] = "ESSA-04-002-cfDNA-1A-C1D1-2020Dec08"

	# Add a new copy_status column based on the log ratio values and SNP allele fractions
	conditions = [
		(cnvs['Log_ratio'] <= -1.0),
		(cnvs['Log_ratio'] <= -0.15) & (cnvs['msf'] >= 0.6),
		(cnvs['Log_ratio'] <= -0.3),
		(cnvs['Log_ratio'] >= 0.15) & (cnvs['msf'] >= 0.6),
		(cnvs['Log_ratio'] >= 0.3),
		(cnvs['Log_ratio'] >= 0.7)]

	values = [-2, -1, -1, 1, 1, 2] # corresponding copy numbers to the conditions above
	cnvs["Copy_status"] = np.select(conditions, values, default = 0) # if no conditions met, gets a value of 0. 

	# sex chromosomes get special treatment for CN calculation. We don't use SNP allele fraction for them. 
	conditions = [
		(cnvs.loc[cnvs.CHROMOSOME.isin(["chrX", "chrY"]), 'Log_ratio'] < -0.7), 
		(-0.7 <= cnvs.loc[cnvs.CHROMOSOME.isin(["chrX", "chrY"]), 'Log_ratio']) & (cnvs.loc[cnvs.CHROMOSOME.isin(["chrX", "chrY"]), 'Log_ratio'] < -0.3), 
		(-0.3 <= cnvs.loc[cnvs.CHROMOSOME.isin(["chrX", "chrY"]), 'Log_ratio']) & (cnvs.loc[cnvs.CHROMOSOME.isin(["chrX", "chrY"]), 'Log_ratio'] < 0.3), 
		(0.3 <= cnvs.loc[cnvs.CHROMOSOME.isin(["chrX", "chrY"]), 'Log_ratio']) & (cnvs.loc[cnvs.CHROMOSOME.isin(["chrX", "chrY"]), 'Log_ratio'] < 0.7), 
		(0.7 < cnvs.loc[cnvs.CHROMOSOME.isin(["chrX", "chrY"]), 'Log_ratio'])]

	values = [-2, -1, 0, 1, 2] # corresponding copy numbers to the conditions above
	cnvs.loc[cnvs.CHROMOSOME.isin(["chrX", "chrY"]), "Copy_status"] = np.select(conditions, values)
	# cnvs.loc[cnvs.CHROMOSOME.isin(["chrX", "chrY"]),'Copy_status'] = cnvs['Copy_status'] # anything felse that doesn't fit the criteria above gets a 0.

	cnvs['Patient_ID'] = cnvs["Sample_ID"].str.split("-cfDNA|-WBC", expand = True)[0] # patient id from sample id
	cnvs['Cycle'] = cnvs["Sample_ID"].str.replace('-1A', '', regex = True).str.split("-", expand = True)[4] # cycle 

	# Make a copy of the df to compute ctDNA% on applicable samples
	cnvs_ctdna = cnvs.copy()
	cnvs_ctdna = cnvs_ctdna[cnvs_ctdna['Copy_status'] == THRESHOLD_copy_status]
	cnvs_ctdna['Hetz_loss_count'] = cnvs_ctdna['GENE'].groupby(cnvs_ctdna['Sample_ID']).transform('count') # count the number of genes that passed the criteria above for each patient
	cnvs_ctdna_copy = cnvs_ctdna.copy() # make a copy to return and do a merge later on
	cnvs_ctdna['Chrom_count'] = cnvs_ctdna['CHROMOSOME'].groupby(cnvs_ctdna['Sample_ID']).transform('count')
	cnvs_ctdna['Sample_median_gene_SNP_VAF_all_genes'] = cnvs_ctdna['msf'].groupby(cnvs_ctdna['Sample_ID']).transform('median') # Get "median of medians" per sample_ID

	# Filter and calculate ctDNA% based on SNPs 
	foo=(cnvs_ctdna["Number of heterozygous SNPs"] > THRESHOLD_het_SNPs_count) & (cnvs_ctdna["Hetz_loss_count"] > THRESHOLD_het_loss_count) & (cnvs_ctdna['Chrom_count'] > THRESHOLD_chrom_count)
	cnvs_ctdna = cnvs_ctdna[foo] # filtering
	cnvs_ctdna['ctDNA_fraction_SNP_estimate'] = 2 - 1/cnvs_ctdna['Sample_median_gene_SNP_VAF_all_genes'] #Calculate tumour content from HSAFs in sLOH regions

	# merge the chosen columns back to the original df
	cnvs_ctdna = cnvs_ctdna[['ctDNA_fraction_SNP_estimate', 'Hetz_loss_count', 'Sample_ID']].drop_duplicates()
	cnvs = pd.merge(cnvs, cnvs_ctdna, how='left', on=['Sample_ID'])

	# move the Patient_ID col right after the Sample_ID col, pop it off and put it back to the correct position
	cnvs = reindex_column(cnvs, "Patient_ID", 1)
	cnvs = reindex_column(cnvs, "Cycle", 2)
	cnvs_ctdna_copy = reindex_column(cnvs_ctdna_copy, "Patient_ID", 1)
	cnvs_ctdna_copy = reindex_column(cnvs_ctdna_copy, "Cycle", 2)

	# Fix some naming issues again, these don't get the right cycle because their naming isnt as expected
	cnvs.loc[cnvs["Sample_ID"] == "ESSA-01-006-cfDNA-2021Jan26", "Cycle"] = "C1D1"
	cnvs.loc[cnvs["Sample_ID"] == "ESSA-09-005-cfDNA-2-C1D1-2021Jan26", "Cycle"] = "C1D1"

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

	# This sample is mislabelled, it is actually C1D1.
	snvs.loc[snvs["Sample_ID"] == "ESSA-04-002-cfDNA-1A-C1D28-2020Dec08", "Sample_ID"] = "ESSA-04-002-cfDNA-1A-C1D1-2020Dec08"

	# Subset to called variants only or keep everything
	if keep_only_called_variants == True: 
		snvs = snvs[snvs['Allele_frequency'].str.contains("\*")]
		snvs['Independently_detected'] = True
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
	# snvs['Cycle'] = snvs["Sample_ID"].str.replace('-1A', '', regex=False).str.split("-", expand = True)[4] # cycle from sample id 

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

	# Fill the na TF estimates in the mutation column with values from the SNP column
	snv_subset["Final_TF"] = snv_subset[col1].fillna(snv_subset[col2])
	
	return(snv_subset)

#Adjust for 95% quantile outlier (when using max somatic VAF to calculate ctDNA%)
def get_adj_maf(maf, depth):
    alt = maf*depth
    for p in np.arange(0,1,0.005):
        dist = stats.binom(depth, p)
        if dist.cdf(alt) < 0.95:
            maf = p; break;
    return maf;

def compute_AR_abs_CN(log_ratio_value, TF_value): 
	'''
	Simple formula to compute absolute AR copy number that can be applied to a df rowwise with an apply function.
	'''
	x = (2**log_ratio_value - 1)/TF_value + 1

	return(x)


def add_AR_abs_CN(cnvs, snv_subset): 
	'''
	For chosen cases, compute the absolute CN for AR and add a a new column. Only applied to cases with AR gain (CN = 1) or amplification (CN = 2)
	snv_subset should be processed such that it has the final TF estimates. We do a left merge to carry this info over to the cnvs df and use it in the AR CN calculation.
	'''
	
	cnvs = pd.merge(cnvs, snv_subset[["Sample_ID", "Patient_ID", "Final_TF"]], how = "left", on = ["Sample_ID", "Patient_ID"])
	cnvs["AR_abs_CN"] = cnvs.apply(lambda x: compute_AR_abs_CN(x['Log_ratio'], x['Final_TF']), axis = 1) # add a new col with the AR CN
	
	# Add NA to cases in the AR_abs_CN column if the AR CN isn't 1 or 2
	cnvs["AR_abs_CN"] = np.where(
		((cnvs["GENE"] == "AR") & (~cnvs['Copy_status'].isin([1, 2])) | (cnvs["GENE"] != "AR")),
		np.nan, 
		cnvs["AR_abs_CN"])

	return(cnvs)