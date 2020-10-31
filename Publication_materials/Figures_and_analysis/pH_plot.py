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

#Setting CST color scheme
CST_color_scheme = {'I-A':'#ff6868','I-B':'#ffd4da','II':'#b4ff68','III-A':'#ffbc6b','III-B':'#e4a67b','IV-A':'#c1adec','IV-B':'#91a8ed',
                    'IV-C0':'#989898','IV-C1':'#ffc0cb','IV-C2':'#a8e5e5','IV-C3':'#9acc9a','IV-C4':'#800080','V':'#ffff71'}

pH_color_scheme = {'lowest':'#f59b42','low':'#f7ed6c','high':'#d7f59a','highest':'#2b9e38'}


data_all = pd.read_csv(sys.argv[1])

stacked_fig, stacked_axs = plt.subplots(1,1, figsize=(12,6), facecolor='w', edgecolor='k')
stacked_fig.subplots_adjust(left=0.2,right=0.98,bottom=0.2,top=0.9,hspace = 0.2, wspace=0.2)

bar_width = 0.15

pH_list = ['lowest','low','high','highest']
data_all['bottom_count'] = pd.Series([0.0 for x in range(len(data_all.index))], index=data_all.index)

for pH_cat in pH_list:
    print(pH_cat)

    stacked_axs.bar(data_all['loc'],data_all[pH_cat],width=bar_width,bottom=data_all['bottom_count'],color=pH_color_scheme[pH_cat],clip_on=False,zorder=2)
    data_all['bottom_count'] = data_all['bottom_count'] + data_all[pH_cat]

stacked_axs.set_xlim(0,4)
stacked_axs.set_ylim(0,1)

stacked_axs.set_xticks(list(data_all['loc']), minor=False,)

stacked_axs.set_xticklabels(['I-A','I-B','','III-A','III-B','','','','IV-C0','IV-C1','IV-C2','IV-C3','IV-C4'],fontsize=12)
stacked_axs.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0],minor=False)
stacked_axs.set_yticklabels([0.0,0.2,0.4,0.6,0.8,1.0],fontsize=12)
stacked_axs.set_ylabel('Proportion of samples in pH category',fontsize=20)

stacked_axs.yaxis.grid(color='black',zorder=3,linewidth=1.5)

stacked_axs.text(0.3,-0.125,s=r'I',fontsize=20,ha='center')
stacked_axs.text(0.8,-0.125,s=r'II',fontsize=20,ha='center')
stacked_axs.text(1.3,-0.125,s=r'III',fontsize=20,ha='center')
stacked_axs.text(1.8,-0.125,s=r'V',fontsize=20,ha='center')
stacked_axs.text(2.2,-0.125,s=r'IV-A',fontsize=20,ha='center')
stacked_axs.text(2.6,-0.125,s=r'IV-B',fontsize=20,ha='center')
stacked_axs.text(3.4,-0.125,s=r'IV-C',fontsize=20,ha='center')

stacked_axs.text(1,1.02,s=r'\textit{Lactobacillus}-dominated',fontsize=20,ha='center')
stacked_axs.text(3,1.02,s=r'Non-\textit{Lactobacillus}-dominated',fontsize=20,ha='center')

#stacked_axs.axvspan(0, 0.6, facecolor='#989898',alpha=0.1,zorder=4)
#stacked_axs.axvspan(0.6, 1, facecolor='#989898',alpha=0.1,zorder=4)
#stacked_axs.axvspan(1, 1.6, facecolor='#989898',alpha=0.1,zorder=4)
#stacked_axs.axvspan(1.6, 2, facecolor='#989898',alpha=0.1,zorder=4)
#stacked_axs.axvspan(2, 2.4, facecolor='#989898',alpha=0.1,zorder=4)
#stacked_axs.axvspan(2.4, 2.8, facecolor='#989898',alpha=0.1,zorder=4)
#stacked_axs.axvspan(2.8, 4, facecolor='#989898',alpha=0.1,zorder=4)

stacked_axs.xaxis.set_tick_params(width=1.5)
stacked_axs.yaxis.set_tick_params(width=1.5)

len_tick = -0.05
stacked_axs.axvline(x=0.6,ymin=len_tick,ymax=1,c='black',clip_on=False)
stacked_axs.axvline(x=1,ymin=len_tick,ymax=1,c='black',clip_on=False)
stacked_axs.axvline(x=1.6,ymin=len_tick,ymax=1,c='black',clip_on=False)
stacked_axs.axvline(x=2,ymin=len_tick,ymax=1,c='black',clip_on=False)
stacked_axs.axvline(x=2.4,ymin=len_tick,ymax=1,c='black',clip_on=False)
stacked_axs.axvline(x=2.8,ymin=len_tick,ymax=1,c='black',clip_on=False)


#making the legend
highest_rect = mpatches.Rectangle((-0.5,.9),width=.15,height=.1,linewidth=1,edgecolor='none',facecolor='#2b9e38',clip_on=False)
high_rect = mpatches.Rectangle((-0.5,.8),width=.15,height=.1,linewidth=1,edgecolor='none',facecolor='#d7f59a',clip_on=False)
low_rect = mpatches.Rectangle((-0.5,.7),width=.15,height=.1,linewidth=1,edgecolor='none',facecolor='#f7ed6c',clip_on=False)
lowest_rect = mpatches.Rectangle((-0.5,.6),width=.15,height=.1,linewidth=1,edgecolor='none',facecolor='#f59b42',clip_on=False)

stacked_axs.add_patch(highest_rect)
stacked_axs.add_patch(high_rect)
stacked_axs.add_patch(low_rect)
stacked_axs.add_patch(lowest_rect)

stacked_axs.text(-.55,0.95,s=r'$\geq$5.5',fontsize=14,ha='right',va='center')
stacked_axs.text(-.55,0.85,s=r'\underline{5.0}-\underline{5.5}',fontsize=14,ha='right',va='center')
stacked_axs.text(-.55,0.75,s=r'4.6-5.0',fontsize=14,ha='right',va='center')
stacked_axs.text(-.55,0.65,s=r'$\leq$4.5',fontsize=14,ha='right',va='center')


stacked_fig.savefig('pH_barplot.pdf')

