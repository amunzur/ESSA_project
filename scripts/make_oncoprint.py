#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  7 16:43:00 2021

@author: amunzur
"""

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D


########################
# DEFINE VARIABLES
########################

# PATH_sample_ids = "/groups/wyattgrp/users/amunzur/onboarding/data/m1rp_patient_sample_IDs.tsv"
PATH_cn = "/groups/wyattgrp/users/amunzur/essa/output_data/all_cnvs_essa.csv"
PATH_muts = "/groups/wyattgrp/users/amunzur/essa/output_data/all_snvs_essa.csv"
PATH_TF = "/groups/wyattgrp/users/amunzur/essa/output_data/ctdna_frac_variants.csv"
PATH_CRPC_2020 = "/groups/wyattgrp/users/amunzur/useful_data/Panel_CRPC2020.csv"
PATH_figure = "/groups/wyattgrp/users/amunzur/essa/figures/oncoprint.pdf"
PATH_utilities_file = "/groups/wyattgrp/users/amunzur/essa/scripts/UTILITIES_make_oncoprint.py"

exec(open(PATH_utilities_file).read()) # source functions 

mutations_to_filter_out = pd.Series(["Intronic", "Intergenic", "Synonymous", "Upstream", "3'-UTR", "5'-UTR"])

fig_width = 8
fig_height = 16

########################
# MODIFY FILES
########################

[df_cn, df_muts, df_TF, crpc2020] = upload_data(PATH_cn, PATH_muts, PATH_TF, PATH_CRPC_2020)
[df_TF_1, df_TF_28] = process_df_TF(df_TF = df_TF, df_cn = df_cn)
[df_cn_1, df_cn_28] = process_df_cn(df_cn = df_cn, df_TF = df_TF, subset_to_crpc2020 = True, crpc2020 = crpc2020)
[df_muts_1, df_muts_28] = process_df_muts(df_muts = df_muts, df_cn = df_cn, mutations_to_filter_out = mutations_to_filter_out, filter_out_WBC = True, subset_to_crpc2020 = True, subset_to_df_cn_genes = True, crpc2020 = crpc2020)

########################
# PREPARE TO PLOT
########################
bar_height = 0.9
bar_width = 0.7

# offset = -(bar_height/3)
offset = -0.4

fig = plt.figure(figsize=(fig_width, fig_height))
gs  = gridspec.GridSpec(nrows = 4, ncols = 2, height_ratios = [0.7, 0.7, 20, 2], width_ratios = [1, 1], wspace = 0.08, hspace = 0.02)
# gs.update(wspace=0.015, hspace=0.05)# set the spacing between axes. 

sub_top_1 = fig.add_subplot(gs[0,0]) # ctDNA
sub_top_28 = fig.add_subplot(gs[0,1]) # ctDNA

sub_mutcount_1 = fig.add_subplot(gs[1,0]) # mut count
sub_mutcount_28 = fig.add_subplot(gs[1,1]) # mut count

sub_bottom_1 = fig.add_subplot(gs[2,0]) # heatmap
sub_bottom_28 = fig.add_subplot(gs[2,1]) # heatmap

sub_legend = fig.add_subplot(gs[3,:])

########################
# PLOTTING
########################
# Plot ctDNA fraction in the top subplot
sub_top_1.bar(df_TF_1["Sample_ID"], df_TF_1["Final_TF"], color = "#202020", zorder = 3)
sub_top_28.bar(df_TF_28["Sample_ID"], df_TF_28["Final_TF"], color = "#202020", zorder = 3)

# plot the mutation count 
sub_mutcount_1.bar(df_muts_1["Sample_ID"], df_muts_1["Mutation_count"], color = "#202020", zorder = 3)
sub_mutcount_28.bar(df_muts_28["Sample_ID"], df_muts_28["Mutation_count"], color = "#202020", zorder = 3)

plot_CN(df_cn_1, sub_bottom_1, bar_height, bar_width, offset)
plot_CN(df_cn_28, sub_bottom_28, bar_height, bar_width, offset)

plot_mutations(df_muts_1, sub_bottom_1, 5)
plot_mutations(df_muts_28, sub_bottom_28, 5)
        
########################
# STYLING
########################
AES_TF(df_TF_1, "sub_top_1", sub_top_1)
AES_TF(df_TF_28, "sub_top_28", sub_top_28)

AES_mutcounts(df_muts_1, "sub_mutcount_1", sub_mutcount_1)
AES_mutcounts(df_muts_28, "sub_mutcount_28", sub_mutcount_28)

AES_heatmap(df_cn_1, "sub_bottom_1", sub_bottom_1)
AES_heatmap(df_cn_28, "sub_bottom_28", sub_bottom_28)

AES_legend(sub_legend)



y_position = 0.72

    
########################
# LEGEND
########################
    
legend_cn_dict = {"Deep deletion":'#3f60ac', "Deletion":'#9cc5e9', "Neutral":'#e6e7e8', "Gain":'#f59496', "Amplification":'#ee2d24'}

mut_dict = {'Missense':'#79B443', 
            'Frameshift':'#BD4398', 
            'Non-frameshift': '#66FFFF',
            'Splice site':'#FFC907',
            'Stopgain':'#FFC907'}

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

fig.tight_layout()


fig.savefig(PATH_figure, bbox_extra_artists=(), bbox_inches='tight')














# set up font size 
font = {'family' : 'normal', 'weight' : 'normal', 'size'   : 4}
plt.rc('font', **font)
# plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.1)
plt.subplots_adjust(wspace=None, hspace=None)
fig.savefig(PATH_figure, bbox_extra_artists=(), bbox_inches='tight', pad_inches = 0)








