def upload_data(PATH_cn, PATH_muts, PATH_TF, PATH_CRPC_2020): 
    '''
    Upload the CN, mutations and the TF files to pandas data frames considering the file types.
    '''

    if PATH_cn.split(".")[1] == "tsv": 
        df_cn = pd.read_csv(PATH_cn, sep = "\t")
        df_muts = pd.read_csv(PATH_muts, sep = "\t")
        df_TF = pd.read_csv(PATH_TF, sep = "\t")
        crpc2020 = pd.read_csv(PATH_CRPC_2020, sep = "\t")

    elif PATH_cn.split(".")[1] == "csv":
        df_cn = pd.read_csv(PATH_cn)
        df_muts = pd.read_csv(PATH_muts)
        df_TF = pd.read_csv(PATH_TF)
        crpc2020 = pd.read_csv(PATH_CRPC_2020)

    return([df_cn, df_muts, df_TF, crpc2020])

def process_df_TF(df_TF, df_cn):
    '''
    Modify the tumor fraction df to choose the needed columns and convert TF to percentage, and sort by sample ID 
    '''
    df_TF = pd.read_csv(PATH_TF)
    
    df_TF = df_TF[["Sample_ID", "Patient_ID", "Final_TF"]]
    df_TF = pd.merge(df_TF, df_cn[["Sample_ID", "Cycle"]].drop_duplicates(), how = "left", on = "Sample_ID") # Add the cycle info we can divide the df into two
    df_TF = df_TF.sort_values(by = "Sample_ID") # sort by sample ID
    df_TF["Final_TF"] = (df_TF["Final_TF"].fillna(0)*100).round(0).astype(int) # convert to percentage

    df_TF_1 = df_TF[df_TF.Cycle == "C1D1"].reset_index(drop = True)
    df_TF_28 = df_TF[df_TF.Cycle == "C1D28"].reset_index(drop = True)        

    return([df_TF_1, df_TF_28])

def process_df_cn(df_cn, df_TF, subset_to_crpc2020, crpc2020): 
    
    df_cn = df_cn[["Sample_ID", "Patient_ID", "Cycle", "GENE", "Copy_status"]] # Choose the cols we need 
    if subset_to_crpc2020 == True: df_cn = df_cn[df_cn["GENE"].isin(crpc2020["Gene name"])] # only keep the genes in the CRPC2020

    # Map the copy number status to a dictionary to assign colors
    cn_dict = {-2:'#3f60ac', -1:'#9cc5e9', 0:'#e6e7e8', 1:'#f59496', 2:'#ee2d24'}
    df_cn['Color'] = df_cn['Copy_status'].map(cn_dict)

    # Map the cycle colors to a dictionary to assign colors 
    cn_dict = {"C1D1":'Red', "C1D28":'Black'}
    df_cn['Color_cycle'] = df_cn['Cycle'].map(cn_dict)

    df_cn_1 = df_cn[df_cn.Cycle == "C1D1"].reset_index(drop = True)
    df_cn_28 = df_cn[df_cn.Cycle == "C1D28"].reset_index(drop = True)        

    return([df_cn_1, df_cn_28])

def process_df_muts(df_muts, df_cn, mutations_to_filter_out, filter_out_WBC, subset_to_crpc2020, subset_to_df_cn_genes, crpc2020):
    '''
    # Modify the df_muts by choosing columns and subsetting the list of mutations 
    df_muts = df with the mutations
    df_cn = df with the copy number changes
    mutations_to_filter_out = pandas series, with the mutations types to filter out
    filter_out_WBC = boolean, whether to keep or filter out the WBC mutations
    subset_to_crpc2020 = boolean, whether to subset to the genes in the panel
    subset_to_df_cn_genes = boolean, whether to subset to genes only in the df_cn
    crpc2020 = the df that contains the genes in the CRPC2020 panel
    '''
    
    if subset_to_crpc2020 == True: df_muts = df_muts[df_muts["GENE"].isin(crpc2020["Gene name"])] # only keep the genes in the CRPC2020
    if subset_to_df_cn_genes == True: df_muts = df_muts[df_muts["GENE"].isin(df_cn["GENE"])] # only keep the genes where we have a copy number estimate for
    if filter_out_WBC == True: df_muts = df_muts[~df_muts.Sample_ID.str.contains("WBC", case = False)].reset_index(drop = True) # filter out WBC samples
    df_muts = df_muts[df_muts.Mutation_type.ne(mutations_to_filter_out)].reset_index(drop = True) # filter out some mutation types 

    df_muts = df_muts[df_muts.Independently_detected.eq(True)].reset_index(drop = True) # only keep the independently detected mutations
    df_muts = df_muts[["Sample_ID", "Patient_ID", "CHROM", "POSITION", "GENE", "Mutation_type"]]
    df_muts['Mutation_count'] = df_muts["Sample_ID"].groupby(df_muts['Sample_ID']).transform('count') # number of mutations each patient has

    # Make sure both dfs have the complete list of patients, if one patient has no mutations add 0 to mutation counts  
    df_muts = pd.merge(df_muts, df_cn["Sample_ID"].drop_duplicates(), how = "right", on = "Sample_ID").fillna(0)      

    # Map the mutation types to colors
    mut_dict = {'Missense':'#79B443', 
                'Frameshift':'#BD4398', 
                'Non-frameshift': '#66FFFF',
                'Splice site':'#FFC907',
                'Stopgain':'#FFC907'}

    # '3-UTR':'#8c69ff', 
    #             '5-UTR':'#FF5733'
    
    df_muts['Color'] = df_muts['Mutation_type'].map(mut_dict)

    df_muts = pd.merge(df_muts, df_cn[["Sample_ID", "Cycle"]].drop_duplicates(), how = "left", on = "Sample_ID") # Add the cycle info we can divide the df into two
    df_muts['Patient_ID'] = df_muts["Sample_ID"].str.split("-cfDNA|-WBC", expand = True)[0] # patient id from sample id

    df_muts_1 = df_muts[df_muts.Cycle == "C1D1"].reset_index(drop = True)
    df_muts_28 = df_muts[df_muts.Cycle == "C1D28"].reset_index(drop = True)

    return([df_muts_1, df_muts_28])

def plot_CN(df_cn, subplot, bar_height, bar_width, offset):
   
    patients = df_cn["Patient_ID"].unique().tolist()
    genes = sorted(df_cn['GENE'].unique().tolist(), reverse = True) # sort alphabetically

    for patient in patients:
        bottom = offset 
        for gene in genes:
            row = df_cn.loc[(df_cn['GENE'] == gene) & (df_cn['Patient_ID'] == patient)]
            color = row['Color'].values[0]

            subplot.bar(patient, bar_height, bottom = bottom, color = color, zorder = 10, width = bar_width * 1.2)

            bottom += 1

def plot_mutations(df_muts, subplot, mut_size):

    df_muts = df_muts.loc[df_muts['GENE'] != 0]

    patients = df_muts["Patient_ID"].unique().tolist()
    genes = sorted(df_muts['GENE'].unique().tolist(), reverse = True) # sort alphabetically

    gene_pos = {genes[i]: list(range(0,len(genes)))[i] for i in range(len(genes))}
    patient_pos = {patients[i]: list(range(0,len(patients)))[i] for i in range(len(patients))}

    for i, row in df_muts.iterrows():
        patient = row['Patient_ID']
        mut_type = row['Mutation_type']
        gene = row['GENE']
        color = row['Color']
        marker_type = "s"
            
        # check if there is another mutation in the same sample/gene combination
        x = [patient, gene] == df_muts[["Patient_ID", "GENE"]]
        if sum(x.all(axis = 1)) == 1: # ONE mutation only in the same sample and gene combination        
            subplot.scatter(x = patient_pos[patient], y = gene_pos[gene], c = color, s = mut_size, marker = marker_type, zorder = 100)

        elif sum(x.all(axis = 1)) == 2: # TWO mutations in the same sample and gene combination
            subplot.scatter(x = patient_pos[patient] + 0.08, y = gene_pos[gene] + 0.08, c = color, s = mut_size, marker = marker_type, zorder = 100)
            subplot.scatter(x = patient_pos[patient] - 0.04, y = gene_pos[gene] - 0.08, c = color, s = mut_size, marker = marker_type, zorder = 100)
        
        elif sum(x.all(axis = 1)) == 3: # THREE mutations 
            subplot.scatter(x = patient_pos[patient], y = gene_pos[gene], c = "black", s = mut_size, marker = "^", zorder = 100)

def AES_heatmap(df_cn, name_subplot, subplot):

    patients = df_cn["Patient_ID"].unique().tolist()
    genes = sorted(df_cn['GENE'].unique().tolist(), reverse = True) # sort alphabetically

    if name_subplot.endswith("1"): 
        subplot.tick_params(labelrotation = 90, direction = "out", length = 0)
        subplot.set_yticks(list(range(0, len(genes))))
        subplot.set_yticklabels(genes, rotation = 0, ha = "right")
        subplot.set_xticks(list(range(0, len(patients))))
        subplot.set_xticklabels(patients, ha = "center")
    else:
        subplot.tick_params(labelrotation = 90, direction = "out")
        subplot.set_yticks([])
        subplot.set_yticklabels([])
        subplot.set_xticks(list(range(0, len(patients))))
        subplot.set_xticklabels(patients, ha = "center")
        subplot.spines['top'].set_visible(False)
        subplot.spines['right'].set_visible(False)
        subplot.spines['left'].set_visible(False)

    subplot.spines['top'].set_visible(False)
    subplot.spines['right'].set_visible(False)
    subplot.spines['bottom'].set_visible(False)
    subplot.spines['left'].set_visible(False)


def AES_mutcounts(df_muts, name_subplot, subplot):

    subplot.set_xticks([])
    subplot.tick_params(which = 'both', length = 0)

    subplot.spines['top'].set_visible(False)
    subplot.spines['right'].set_visible(False)
    subplot.spines['bottom'].set_visible(False)
    subplot.spines['left'].set_visible(False)

    subplot.grid(zorder = 0, linestyle='--', linewidth = 0.5)

    if name_subplot.endswith("1"): 
        subplot.set_xticks([])
        subplot.set_yticks([0, 10, 20])
        subplot.set_yticklabels(["0", "10", "20"])
        subplot.set_ylabel("Mutation \n count", labelpad = 20, rotation = 0, va = 'center')
        subplot.grid(zorder = 0, linestyle='--', linewidth = 0.5)
        
    else:
        subplot.set_yticks([])
        subplot.set_yticklabels([])

def AES_TF(df_TF, name_subplot, subplot):

    subplot.set_xticks([])
    subplot.tick_params(which = 'both', length = 0)

    subplot.spines['top'].set_visible(False)
    subplot.spines['bottom'].set_visible(False)
    subplot.spines['right'].set_visible(False)
    subplot.spines['left'].set_visible(False)

    subplot.grid(zorder = 0, linestyle='--', linewidth = 0.5)

    if name_subplot.endswith("1"): 
        subplot.set_yticks([0, 50, 100])
        subplot.set_yticklabels(["0", "50", "100"])
        subplot.set_ylabel("ctDNA %", labelpad=17, rotation = 0, va = 'center')
    else:
        subplot.set_yticks([0, 50, 100])
        subplot.set_yticklabels([])
        subplot.tick_params(which = 'both', length = 0)

def AES_legend(subplot):

    subplot.spines['top'].set_visible(False)
    subplot.spines['bottom'].set_visible(False)
    subplot.spines['right'].set_visible(False)
    subplot.spines['left'].set_visible(False)
