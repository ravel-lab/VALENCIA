#!/usr/bin/python3

import matplotlib.pyplot as plt
import pandas as pd
import sys
import numpy as np
import matplotlib

matplotlib.rc('font', serif='Helvetica Neue')
matplotlib.rc('axes',linewidth=2)

data = pd.read_csv(sys.argv[1],sep=",")
data_sel = data[data['Number_clust']=='Sel']
data_options = data[data['Number_clust'] !='Sel']

stacked_fig, stacked_axs = plt.subplots(1,2, figsize=(12,6), facecolor='w', edgecolor='k')
stacked_fig.subplots_adjust(left=0.10,right=0.9,bottom=0.1,top=0.9, wspace=0.25,hspace=0.25)

stacked_axs = stacked_axs.ravel()

stacked_axs[0].plot(data_options['Number_clust'],data_options['Silh_score'], linestyle='-', marker='o')
stacked_axs[0].set_xlabel("Number of Clusters",fontsize="14")
stacked_axs[0].set_ylabel("Silhouette score",fontsize="14")
stacked_axs[0].set_ylim(0,0.4)
stacked_axs[0].axhline(y=data_sel['Silh_score'].values)
stacked_axs[0].text(x=12,y=data_sel['Silh_score'].values+0.01,s=r"selected solution$^*$")

stacked_axs[1].plot(data_options['Number_clust'],data_options['Db_score'], linestyle='-', marker='o')
stacked_axs[1].set_xlabel("Number of Clusters",fontsize="14")
stacked_axs[1].set_ylabel("Davies-Bouldin Index",fontsize="14")
stacked_axs[1].set_ylim(0.5,2.5)
stacked_axs[1].axhline(y=data_sel['Db_score'].values)
stacked_axs[1].text(x=12,y=data_sel['Db_score'].values+0.05,s=r"selected solution$^*$")

stacked_fig.savefig("Cluster_number_scores.pdf")
