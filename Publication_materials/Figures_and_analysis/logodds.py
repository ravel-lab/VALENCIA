#!/usr/bin/python3

#loading in required python modules
from Bio import SeqIO
import sys
import matplotlib.pyplot as plt
from statistics import median
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib
import random

data = pd.read_csv(sys.argv[1],sep=",")

cst_colors = ['#ff6868','#ffd4da','#b4ff68','#ffbc6b','#e4a67b','#c1adec','#91a8ed','#989898','#ffc0cb','#a8e5e5','#9acc9a','#800080','#ffff71']

stacked_fig, stacked_axs = plt.subplots(1,1, figsize=(12,8), facecolor='w', edgecolor='k')
stacked_fig.subplots_adjust(left=0.1,right=0.9,bottom=0.1,top=0.9,hspace = 0.2, wspace=0.2)

stacked_axs.scatter(x=data['x_loc'],y=data['Estimate'],c=cst_colors,s=30,marker='o')
stacked_axs.set_xticks(data['x_loc'].values,minor=False)
stacked_axs.set_yticks([0,5,10,15,20,25,30,35,40,45,50,55,60],minor=False)
stacked_axs.set_xticklabels(data['subCST'].values,fontsize=14)
stacked_axs.set_yticklabels([0,5,10,15,20,25,30,35,40,45,50,55,60],minor=False,fontsize=14)

stacked_axs.set_ylabel('Odds vaginal pH >4.5',fontsize=16)
stacked_axs.set_xlabel('subCST',fontsize=16)
stacked_axs.vlines(x=data['x_loc'],ymin=data['Low'],ymax=data['High'],linewidth=1,colors=cst_colors,clip_on=False)
stacked_axs.scatter(x=data['x_loc'],y=data['Estimate'],c=cst_colors,s=30,marker='o')
 
stacked_fig.savefig("pH_logodds.pdf")