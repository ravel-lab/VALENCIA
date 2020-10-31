#!/usr/bin/env python3

#importing packages to be used
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import numpy as np
import sys
import matplotlib
import seaborn as sns
import random
import matplotlib.patches as mpatches
from skbio.diversity.alpha import shannon


taxa_color_scheme = {'Lactobacillus_crispatus':'#ff0000','Gardnerella_vaginalis':'#20b2aa','g_Lactobacillus':'#eef06c',
                     'Lactobacillus_iners':'#ff8c00','Lactobacillus_gasseri':'#7fff00','Lactobacillus_jensenii':'#333333',
                     'Enterococcus_faecalis':'#bbffff','Raoultella_planticola':'#ffc2ff','g_Peptoniphilus':'#CCCC00',
                     'Sneathia_sanguinegens':'#c3cdce','Atopobium_vaginae':'#0000cd','g_Atopobium':'#0000cd',
                     'Lacotbacillus_helveticus':'#00ccff','BVAB3':'#3cb371','g_Anaerococcus':'#87cefa','g_Gardnerella':'#20b2aa',
                     'Megasphaera_sp_type_1':'#0000ff','Streptococcus_agalactiae':'#ffc0cb','g_Megasphaera':'#008B45',
                     'Streptococcus_oralis':'#ff66ff','Prevotella_bivia':'#bfbfbf','Aerococcus_christensenii':'#bebebe',
                     'Anaerococcus_tetradius':'#87cefa','Gemella':'#daa520','Prevotella_genogroup_1':'#b0b0b0',
                     'Lactobacillus_vaginalis':'#ffffff','other':'#808080','Bifidobacterium_longum':'#c1ffc1',
                     'Bifidobacterium_breve':'#c1ffc1','Eggerthella':'#DE7710','Mycoplasma_hominis':'#10DE4E',
                     'Porphyromonas_bennonis':'#DE4310','Eubacterium_saphenum':'#8F10DE','Fusobacterium_nucleatum':'#CD853F',
                     'Fusobacterium_gonidiaformans':'#CD853F','Streptococcus_equinus':'#ffc0cb','Peptostreptococcus_anaerobius':'#DEDB10',
                     'Arcanobacterium_phocae':'#8c10de','Bacteroides_uniformis':'#de1058','Ureaplasma_parvum':'#9999ff',
                     'Peptoniphilus_harei':'#CCCC00','Mobiluncus_mulieris':'#f08080','Megasphaera_sp._type_2':'#008B45',
                     'BVAB1':'#b31900','BVAB2':'#3bb16f','g_Escherichia.Shigella':'#12456b','Peptoniphilus_lacrimalis':'#d2b48c',
                     'Veillonella_montpellierensis':'#ff8c69','Prevotella_genogroup_3':'#b36200','Parvimonas_micra':'#cdcd00',
                     'Corynebacterium_accolens':'#ffff00','Finegoldia_magna':'#800080','Prevotella_genogroup_2':'#b0b0b0',
                     'Staphylococcus_epidermidis':'#ffffff','Prevotella_timonensis':'#bfbfbf','g_Streptococcus':'#ffc0cb','g_Bifidobacterium':'#c1ffc1','g_Enterococcus':'#bbffff'
                     ,'g_Staphylococcus':'#800080','g_Finegoldia':'#800080','g_Prevotella':'#bfbfbf','g_Sneathia':'#d1e8eb','g_Aerococcus':'#e8e1ba'
                     ,'g_Leptotrichia':'#7ca386','g_Veillonella':'#499996','g_Dialister':'#997dbd','g_Corynebacterium_1':'#ffff00'
                     ,'g_Varibaculum':'#094717','g_Delftia':'#090c47','g_Corynebacterium_1':'#ffff00'}

def get_color(taxa):

    chars = '0123456789ABCDEF'

    if taxa in taxa_color_scheme:
        
        taxa_color = taxa_color_scheme[taxa]

    else:
        taxa_color_scheme[taxa] = '#'+''.join(random.sample(chars,6))
        taxa_color = taxa_color_scheme[taxa]

    return taxa_color


matplotlib.rc('font', serif='Helvetica Neue')
matplotlib.rc('axes',linewidth=2)

#reading the data in 
data = pd.read_csv(sys.argv[1],sep=",")

CSTs = ['I','II','III','V','IV-A','IV-B','IV-C']

bar_width = 1

plot_count=0
data_all_rel = data[data.columns[6:]].div(data['read_count'],axis=0)
data_all = pd.concat([data[data.columns[0:6]],data_all_rel],axis=1)


legend_entries = {}

stacked_fig, (stacked_axs,similarity_axs) = plt.subplots(1,2, figsize=(24,10), facecolor='w', edgecolor='k',gridspec_kw={'width_ratios': [6, 1]})
stacked_fig.subplots_adjust(left=0.18,right=0.98,bottom=0.1,top=0.9, wspace=0.2,hspace=0.05)

data_sub = data_all

data_sub_micro = data_sub[data_sub.columns[6:]]

top_columns = data_sub_micro.sum(axis=0).sort_values(ascending=False).head(15).index

data_sub_micro = data_sub_micro[top_columns]

data_sub_micro['g_Enterococcus'] = data_sub['g_Enterococcus']
data_sub_micro['g_Staphylococcus'] = data_sub['g_Staphylococcus']
data_sub_micro['BVAB1'] = data_sub['BVAB1']

data_sub_micro['other'] = data_sub_micro.apply(lambda row: 1.0 - row.sum(), axis=1)

data_sub = pd.concat([data_sub[data_sub.columns[0:6]],data_sub_micro],axis=1,sort=False)

data_sub['bottom_count'] = pd.Series([0.0 for x in range(len(data_sub.index))], index=data_sub.index)

data_sub['x_values'] = data_sub.index
for taxa in range(6,len(top_columns)+9):

    taxa = data_sub.columns[taxa]

    print(taxa)

    taxa_color = get_color(taxa)

    if taxa not in legend_entries:

        legend_entries[taxa] = taxa_color_scheme[taxa]
    
    stacked_axs.bar(data_sub['x_values'],data_sub[taxa],width=bar_width,bottom=data_sub['bottom_count'],color=taxa_color,clip_on=False)
    
    data_sub['bottom_count'] = data_sub['bottom_count'] + data_sub[taxa]

stacked_axs.set_ylim([0.0,1.0])
stacked_axs.tick_params(width=2)
stacked_axs.set_xticks([1250,3100,3900,5350,7250,8350,9050,11000,13000], minor=False)
stacked_axs.set_yticks([0,0.2,0.4,0.6,0.8,1.0],minor=False)
stacked_axs.set_yticklabels(['0','0.2','0.4','0.6','0.8','1.0'],minor=False,fontsize=20)
stacked_axs.set_xticklabels(['I-A','I-B','II','III-A','III-B','V','IV-A','IV-B','IV-C'],fontsize=20)
stacked_axs.set_xlabel('CSTs assigned by VALENCIA',fontsize=24)
stacked_axs.set_ylabel('Relative abundance',fontsize=24)

stacked_axs.yaxis.grid(color='gray')
stacked_axs.set_axisbelow(True)

taxa_text_legend = {'Lactobacillus_iners':'L. iners','Lactobacillus_crispatus':'L. crispatus','Gardnerella_vaginalis':'G. vaginalis','Lactobacillus_jensenii':'L. jensenii',
                    'Atopobium_vaginae':"A. vaginae","Lactobacillus_gasseri":'L. gasseri','BVAB1':'BVAB1','g_Streptococcus':'Streptococcus','Sneathia_sanguinegens':'S. sanguinegens',
                    'g_Megasphaera':'Megasphaera','g_Fastidiosipila':'Fastidiosipila','g_Bifidobacterium':'Bifidobacterium','g_Anaerococcus':'Anaerococcus','g_Peptoniphilus':'Peptoniphilus',
                    'Prevotella_timonensis':'P. timonensis','g_Enterococcus':'Enterococcus','g_Staphylococcus':'Staphylococcus','other':'other'}

patch_list = []
for taxa in legend_entries:
    taxa_text = taxa_text_legend[taxa]

    if taxa_text == "BVAB1" or taxa_text == "other":
        data_key = mpatches.Patch(color=legend_entries[taxa],label="%s" %(taxa_text))
    
    else:    
        data_key = mpatches.Patch(color=legend_entries[taxa],label="$\it{%s}$" %(taxa_text))
    
    patch_list.append(data_key)

stacked_fig.legend(handles=patch_list,loc=6 ,ncol=1,fontsize=16)



CST_color_scheme = {'I-A':'#ff6868','I-B':'#ffd4da','II':'#b4ff68','III-A':'#ffbc6b','III-B':'#e4a67b','IV-A':'#c1adec','IV-B':'#91a8ed',
                    'IV-C0':'#989898','IV-C1':'#ffc0cb','IV-C2':'#a8e5e5','IV-C3':'#9acc9a','IV-C4':'#800080','V':'#ffff71'}

CSTs = ['I-A','I-B','II','III-A','III-B','V','IV-A','IV-B','IV-C0','IV-C1','IV-C2','IV-C3','IV-C4']

#calculating shannon diversity
data['shannon'] = data.apply(lambda y: shannon(list(y)[6:205]),axis=1)

#building the plot
    
#creating x axis location variables
loc=12

for CST in CSTs:
    boxprops = dict(linewidth=1, color="k")
    medianprops = dict(linewidth=1,color="k")    

    box = similarity_axs.boxplot(x=data[data['subCST'] == CST].shannon,positions=[loc],notch=True,widths=[0.5],patch_artist=True,boxprops=boxprops,medianprops=medianprops,vert=False)
    
    patch = box['boxes']
    for patch in box['boxes']:

        patch.set_facecolor(CST_color_scheme[CST])
    loc = loc - 1

similarity_axs.set_yticks([0,1,2,3,4,5,6,7,8,9,10,11,12])
CSTs.reverse()
similarity_axs.set_yticklabels(CSTs,fontsize=20,verticalalignment='center')
#similarity_axs.xaxis.tick_right()
similarity_axs.set_xticks([0,1,2,3,4,5])
similarity_axs.set_xticklabels(['0','1','2','3','4','5'],fontsize=16)
similarity_axs.set_xlabel("Shannon diversity",fontsize=24)
similarity_axs.xaxis.grid(color='#d4d4d4')
similarity_axs.set_axisbelow(True)

stacked_axs.text(x=-200,y=1.02,s="a",fontsize=44,clip_on=False)
similarity_axs.text(x=-3,y=12.5,s="b",fontsize=44,clip_on=False)

stacked_fig.savefig('all_barplots.pdf')









