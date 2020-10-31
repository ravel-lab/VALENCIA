#!/usr/bin/env python3

#importing packages to be used
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import sys
import matplotlib
import matplotlib.patches as mpatches



matplotlib.rc('text', usetex = True)
matplotlib.rc('font', serif='Helvetica Neue')
matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica Neue']})
matplotlib.rc('lines', linewidth=1.5, color='black')
matplotlib.rcParams['axes.linewidth'] = 1.5


#reading in the data and manipulating
data_all=pd.read_csv(sys.argv[1],sep=",")

subject_data = pd.pivot_table(data_all,index=['PID','Race','age'],columns='subCST',values='sampleID',aggfunc=len,margins=True)

subject_data = subject_data[subject_data.columns[0:13]].div(subject_data['All'],axis=0)

#reseting index and dropping last row which contains summary data
subject_data = subject_data.reset_index().iloc[:-1]
subject_data = subject_data.drop_duplicates(subset='PID',keep='first')

####summaryizing by race

race_data = subject_data.groupby(['Race']).sum()
race_data = race_data.div(race_data.sum(axis=1),axis=0).reset_index()

print(race_data)

#color scheme
CST_color_scheme = {'I-A':'#ff6868','I-B':'#ffd4da','II':'#b4ff68','III-A':'#ffbc6b','III-B':'#e4a67b','IV-A':'#c1adec','IV-B':'#91a8ed',
                    'IV-C0':'#989898','IV-C1':'#ffc0cb','IV-C2':'#a8e5e5','IV-C3':'#9acc9a','IV-C4':'#800080','V':'#ffff71'}


stacked_fig, stacked_axs = plt.subplots(1,1, figsize=(10,6), facecolor='w', edgecolor='k')
stacked_fig.subplots_adjust(left=0.3,right=0.98,bottom=0.2,top=0.9,hspace = 0.2, wspace=0.2)

bar_width = 0.15


race_data['bottom_count'] = pd.Series([0.0 for x in range(len(race_data.index))], index=race_data.index)

cst_list = ['I-A','I-B','II','III-A','III-B','V','IV-A','IV-B','IV-C0','IV-C1','IV-C2','IV-C3','IV-C4']

race_data.insert(1,"loc",[0.2,0.4,0.6,0.8,1.0])


yloc = 0.95

for cst in cst_list:
    print(cst)

    stacked_axs.bar(race_data['loc'],race_data[cst],width=bar_width,bottom=race_data['bottom_count'],color=CST_color_scheme[cst],clip_on=False,zorder=2)
    race_data['bottom_count'] = race_data['bottom_count'] + race_data[cst]

    #making the legend
    rect = mpatches.Rectangle((-0.09,yloc),width=0.06,height=0.05,linewidth=1,edgecolor='none',facecolor=CST_color_scheme[cst],clip_on=False)
    stacked_axs.add_patch(rect)
    stacked_axs.text(-.025,yloc+0.025,s=r'%s' %(cst),fontsize=10,ha='left',va='center')
    yloc=yloc-0.05


stacked_axs.set_xlim(0.1,1.1)
stacked_axs.set_ylim(0,1)

stacked_axs.set_xticks(list(race_data['loc']), minor=False,)

stacked_axs.set_xticklabels(['Asian','Black','Hispanic','Other','White'],fontsize=14)
stacked_axs.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0],minor=False)
stacked_axs.set_yticklabels([0.0,0.2,0.4,0.6,0.8,1.0],fontsize=14)

stacked_axs.xaxis.set_tick_params(width=1.5)
stacked_axs.yaxis.set_tick_params(width=1.5)

stacked_axs.axvline(x=-0.1,ymin=0.725,ymax=0.975,c='black',clip_on=False)
stacked_axs.axvline(x=-0.1,ymin=0.375,ymax=0.675,c='black',clip_on=False)
stacked_axs.text(-.125,0.85,s=r'\textit{Lactobacillus}-dom.',fontsize=10,ha='center',va='center',rotation=90)
stacked_axs.text(-.125,0.525,s=r'Non-\textit{Lactobacillus}-dom.',fontsize=10,ha='center',va='center',rotation=90)

stacked_fig.savefig('race_cst.pdf')
