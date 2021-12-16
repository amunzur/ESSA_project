import pandas as pd
import string

PATH_snp = "/groups/wyattgrp/users/amunzur/essa/data/hetz_snps.vcf" # path to snp calls, raw and unmelted
PATH_cna = "/groups/wyattgrp/users/amunzur/essa/data/gene_cna.csv" # path to the cna calls, raw and unmelted

snp = pd.read_csv(PATH_snp, sep = '\t')
snp = pd.melt(snp, id_vars = ['CHROM','POSITION','REF','ALT','NOTES'])

snp = snp.rename(columns = {'variable':'Sample_ID'})
snp[["Read_depth", "Mutant_reads"]] = snp["value"].str.split(pat=":", expand=True)[[0, 1]]
snp['VAF'] = snp['Mutant reads'] / snp['Read depth']
snp['Patient_ID'] = snp['Sample_ID'].str.split('-', expand = True)[[0, 1, 2]].str.join("-")
snp['Sample_type'] = snp['Sample_ID'].str.split('_').str[2].str.strip(string.digits)
snp.loc[snp['Sample_ID'].str.contains('WBC'), 'Sample type'] = 'gDNA'

gl = snp.copy()
gl = gl[gl['Sample type'].isin(['NL','gDNA','WBC'])]
gl = gl[gl['value'].str.contains('\*')]
gl['Called'] = True

snp = snp[~snp['Sample type'].isin(['gDNA','WBC'])]

gl = gl[['CHROM','POSITION','REF','ALT','Patient ID','Called']]
snp = snp.merge(gl, on = ['CHROM','POSITION','REF','ALT','Patient ID'])

def label_genes(snp, gene):
    snp['GENE'] = None
    for index, row in gene.iterrows():
        start = row['START']-1000
        end = row['END']+1000
        chrom = row['CHROMOSOME']
        snp.loc[(snp['CHROM']==chrom)&(snp['POSITION']<=end)&(snp['POSITION']>=start), 'GENE'] = row['GENE']
    return snp;

cn = pd.read_csv(PATH_cna)

cn = cn[['GENE','CHROMOSOME','START','END']].drop_duplicates()

snp = label_genes(snp,cn)

snp.to_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/hetz_snps/M1RP_snps_melted.tsv', sep = '\t')
