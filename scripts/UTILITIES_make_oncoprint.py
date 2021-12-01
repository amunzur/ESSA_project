def upload_data(PATH_cn, PATH_muts, PATH_TF, PATH_CRPC_2020, PATH_gene_list): 
    '''
    Upload the CN, mutations and the TF files to pandas data frames considering the file types.
    '''
    if PATH_cn.split(".")[1] == "tsv": 
        df_cn = pd.read_csv(PATH_cn, sep = "\t")
        df_muts = pd.read_csv(PATH_muts, sep = "\t")
        df_TF = pd.read_csv(PATH_TF, sep = "\t")
        crpc2020 = pd.read_csv(PATH_CRPC_2020, sep = "\t")
        gene_list = pd.read_csv(PATH_gene_list, sep = "\t")

    elif PATH_cn.split(".")[1] == "csv":
        df_cn = pd.read_csv(PATH_cn)
        df_muts = pd.read_csv(PATH_muts)
        df_TF = pd.read_csv(PATH_TF)
        crpc2020 = pd.read_csv(PATH_CRPC_2020)
        gene_list = pd.read_csv(PATH_gene_list)

    return([df_cn, df_muts, df_TF, crpc2020, gene_list])

def reorder_dfs(df_cn_main, df_muts):
    '''
    Change the order of genes in the df_cn and df_muts to put the genes with highest number of changes to the top. 
    '''
    # df_cn_main = df_cn.copy
    df_cn = df_cn_main.copy()
    df_cn = df_cn[["GENE", "Copy_status"]]
    df_cn = df_cn[df_cn["Copy_status"] != 0].reset_index(drop = True) # filter out the genes with CN 0
    del df_cn["Copy_status"]
    df_cn["Gene_counts_cn"] = df_cn["GENE"].groupby(df_cn["GENE"]).transform("count")
    df_cn = df_cn.sort_values(by = "Gene_counts_cn", ascending = False).reset_index(drop = True) 

    df_muts = df_muts[["GENE"]]
    df_muts = df_muts[df_muts["GENE"] != 0].reset_index(drop = True)
    df_muts["Gene_counts_muts"] = df_muts["GENE"].groupby(df_muts["GENE"]).transform("count")
    df_muts = df_muts.sort_values(by = "Gene_counts_muts", ascending = False).reset_index(drop = True) 

    combined = pd.merge(df_cn, df_muts, how = "outer", on = ["GENE"])
    combined = combined.fillna(0).drop_duplicates().reset_index(drop = True) 
    combined["Total_alterations"] = combined["Gene_counts_cn"] + combined["Gene_counts_muts"]
    combined = combined.sort_values(by = "Total_alterations", ascending = False).reset_index(drop = True) 

    # combined df has the right order of genes we are interested in, now we sort the other dfs in the same order
    # df_cn = pd.merge(combined.GENE, df_cn, how = "left", on = "GENE") # now genes in correct order
    df_cn_main = pd.merge(combined, df_cn_main, how = "left", on = "GENE").reset_index(drop = True) 

    return(df_cn_main)

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

def process_df_cn(df_cn, df_TF, subset_to_crpc2020, crpc2020, subset_to_gene_list, gene_list): 
    
    df_cn = df_cn[["Sample_ID", "Patient_ID", "Cycle", "GENE", "Copy_status"]] # Choose the cols we need 
    if subset_to_crpc2020 == True: df_cn = df_cn[df_cn["GENE"].isin(crpc2020["Gene name"])] # only keep the genes in the CRPC2020
    if subset_to_gene_list == True: 
        df_cn = pd.merge(gene_list, df_cn, how = "left", on = "GENE")
        df_cn = df_cn.dropna(subset=["Sample_ID"]) # the merge above introduces artifacts if the gene isnt in df_muts, drop those. 

    # Map the copy number status to a dictionary to assign colors
    cn_dict = {-2:'#3f60ac', -1:'#9cc5e9', 0:'#e6e7e8', 1:'#f59496', 2:'#ee2d24'}
    df_cn['Color'] = df_cn['Copy_status'].map(cn_dict)

    # Map the cycle colors to a dictionary to assign colors 
    cn_dict = {"C1D1":'Red', "C1D28":'Black'}
    df_cn['Color_cycle'] = df_cn['Cycle'].map(cn_dict)

    df_cn = df_cn.iloc[::-1] # Reverse the order of genes so that they appear correctly on the heatmap

    return(df_cn)

def process_df_muts(df_muts, df_cn, mutations_to_filter_out, filter_out_WBC, subset_to_crpc2020, subset_to_df_cn_genes, crpc2020, subset_to_gene_list, gene_list):
    '''
    # Modify the df_muts by choosing columns and subsetting the list of mutations 
    df_muts = df with the mutations
    df_cn = df with the copy number changes
    mutations_to_filter_out = pandas series, with the mutations types to filter out
    filter_out_WBC = boolean, whether to keep or filter out the WBC mutations
    subset_to_crpc2020 = boolean, whether to subset to the genes in the panel
    subset_to_df_cn_genes = boolean, whether to subset to genes only in the df_cn
    crpc2020 = the df that contains the genes in the CRPC2020 panel\
    gene_list = gene list to subset and annotate the genes with the pathways they belong to
    '''
    
    if subset_to_crpc2020 == True: df_muts = df_muts[df_muts["GENE"].isin(crpc2020["Gene name"])] # only keep the genes in the CRPC2020
    if subset_to_df_cn_genes == True: df_muts = df_muts[df_muts["GENE"].isin(df_cn["GENE"])] # only keep the genes where we have a copy number estimate for
    if subset_to_gene_list == True: 
        df_muts = pd.merge(gene_list, df_muts, how = "left", on = "GENE").reset_index(drop = True) 
        df_muts = df_muts.dropna(subset=["Sample_ID"]) # the merge above introduces artifacts if the gene isnt in df_muts, drop those. 

    if filter_out_WBC == True: df_muts = df_muts[~df_muts.Sample_ID.str.contains("WBC", case = False)].reset_index(drop = True) # filter out WBC samples
    df_muts = df_muts[~df_muts.Mutation_type.isin(mutations_to_filter_out)].reset_index(drop = True) # filter out some mutation types 

    df_muts = df_muts[df_muts.Independently_detected.eq(True)].reset_index(drop = True) # only keep the independently detected mutations
    df_muts = df_muts[["Sample_ID", "Patient_ID", "CHROM", "POSITION", "GENE", "PATHWAY", "Mutation_type"]]
    
    df_muts['Mutation_count'] = df_muts["Sample_ID"].groupby(df_muts['Sample_ID']).transform('count') # number of mutations each patient has
    df_muts['Color'] = df_muts['Mutation_type'].map(mut_dict)
    df_muts['Patient_ID'] = df_muts["Sample_ID"].str.split("-cfDNA|-WBC", expand = True)[0] # patient id from sample id

    df_muts = pd.merge(df_muts, df_cn[["Sample_ID", "Cycle"]].drop_duplicates(), how = "left", on = "Sample_ID") # Add the cycle info we can divide the df into two

    return(df_muts)

def plot_CN(df_cn, subplot, bar_height, bar_width, offset):
   
    patients = df_cn["Patient_ID"].unique().tolist()
    genes = df_cn['GENE'].unique().tolist() # sort alphabetically

    for patient in patients:
        bottom = offset 
        for gene in genes:
            row = df_cn.loc[(df_cn['GENE'] == gene) & (df_cn['Patient_ID'] == patient)]
            color = row['Color'].values[0]

            subplot.bar(patient, bar_height, bottom = bottom, color = color, zorder = 10, width = bar_width * 1.2)

            bottom += 1

def plot_mutations(df_muts, df_cn, subplot, mut_size, mut_dict):
    '''
    Plot mutations based on the order of patients and genes from the df_cn data frame
    '''

    df_muts = df_muts.loc[df_muts['GENE'] != 0]

    patients = df_cn["Patient_ID"].unique().tolist()
    genes = df_cn['GENE'].unique().tolist()

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
    genes = df_cn['GENE'].unique().tolist()
    subplot.tick_params(which = 'both', length = 0)

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

    # subplot.grid(zorder = 0, linestyle='--', linewidth = 0.5)

    if name_subplot.endswith("1"): 
        subplot.set_xticks([])
        subplot.set_yticks([0, 10, 20])
        subplot.set_yticklabels(["0", "10", "20"])
        subplot.set_ylabel("Mutation \n count", labelpad = 20, rotation = 0, va = 'center')
        
    else:
        subplot.set_yticks([])
        subplot.set_yticklabels([])
        subplot.grid(zorder = 0, linestyle='--', linewidth = 0.5)

def AES_TF(df_TF, name_subplot, subplot):

    subplot.set_xticks([])
    subplot.tick_params(which = 'both', length = 0)

    subplot.spines['top'].set_visible(False)
    subplot.spines['bottom'].set_visible(False)
    subplot.spines['right'].set_visible(False)
    subplot.spines['left'].set_visible(False)

    # subplot.grid(zorder = 0, linestyle='--', linewidth = 0.5)

    if name_subplot.endswith("1"): 
        subplot.set_yticks([0, 50, 100])
        subplot.set_yticklabels(["0", "50", "100"])
        subplot.set_ylabel("ctDNA %", labelpad=17, rotation = 0, va = 'center')
    else:
        subplot.set_yticks([0, 50, 100])
        subplot.set_yticklabels([])
        subplot.tick_params(which = 'both', length = 0)

def AES_legend(subplot, subplot_heatmap, y_position, mut_dict):

    subplot.spines['top'].set_visible(False)
    subplot.spines['bottom'].set_visible(False)
    subplot.spines['right'].set_visible(False)
    subplot.spines['left'].set_visible(False)
    subplot.tick_params(which = 'both', length = 0)
    subplot.set_xticklabels([])
    subplot.set_yticklabels([])

    # move the legend axes down a bit based on the location of the heatmap
    pos1 = subplot_heatmap.get_position()
    pos2 = subplot.get_position()
    
    points4 = pos1.get_points()
    points5 = pos2.get_points()
    
    points5[1][1] = points4[0][1] - 0.08
    pos2.set_points(points5)
    subplot.set_position(pos2)

    legend_cn_dict = {"Deep deletion":'#3f60ac', "Deletion":'#9cc5e9', "Neutral":'#e6e7e8', "Gain":'#f59496', "Amplification":'#ee2d24'}

    mut_dict_shape = {'Somatic':'s', '>2 mutations': '^'}
    mut_dict_shape_color = {'Somatic':'#B0B0B0', '>2 mutations': 'black'}

    # legend1
    handles_cn = []
    for key in legend_cn_dict:
        handle = mpatches.Patch(color = legend_cn_dict.get(key), label = key)
        handles_cn.append(handle)

    # legend2
    handles_muts = []
    for key in mut_dict:   
        handle = mpatches.Patch(color = mut_dict.get(key), label = key)
        handles_muts.append(handle)

    # legend3
    handles_mut_shapes = []
    for key in mut_dict_shape:    
        handle = Line2D([0], [0], linestyle = "none", marker = mut_dict_shape.get(key), label = key, markerfacecolor = mut_dict.get(key), color = mut_dict_shape_color.get(key), markersize=5)
        handles_mut_shapes.append(handle)

    legend1 = fig.legend(handles = handles_cn, bbox_to_anchor=(0.33, y_position - 0.45), frameon=False, title = "Copy number variants", title_fontsize = 12)
    legend2 = fig.legend(handles = handles_muts, bbox_to_anchor=(0.56, y_position - 0.45), frameon=False, title = "Mutations", title_fontsize = 12)
    legend3 = fig.legend(handles = handles_mut_shapes, bbox_to_anchor=(0.46, y_position - 0.53), frameon=False)

    # align the legend titlesshapes_dict = {}
    legend1._legend_box.align = "left"
    legend2.get_title().set_position((-40, 0))

def AES_add_annots(subplot):

    pos_array = subplot.get_position().get_points()
    xpos = np.subtract(pos_array[1, 1], pos_array[1, 0])/2
    ypos = 0.70

    subplot.annotate('More text in axes ax2', 
             (xpos, ypos), # these are the coordinates to position the label
             color = 'black',
             size = 10)







