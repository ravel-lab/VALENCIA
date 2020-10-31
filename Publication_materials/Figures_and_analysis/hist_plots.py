#!/usr/bin/env python3

#importing packages to be used
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import seaborn as sns
import matplotlib.patches as mpatches

data_in = pd.read_csv(sys.argv[1],sep=",")

ref_in = pd.read_csv(sys.argv[2],sep=",")

print(data_in.head())
print(ref_in.head())

stacked_fig, stacked_axs = plt.subplots(1,1, figsize=(8,8), facecolor='w', edgecolor='k')
stacked_fig.subplots_adjust(left=0.1,right=0.9,bottom=0.1,top=0.9, wspace=0.2)

num_bins = 20

n_ref, bins_ref, patches_ref = stacked_axs.hist(ref_in['score'], num_bins,density =True,color='#8b8b8b')
n, bins, patches = stacked_axs.hist(data_in['score'], num_bins,density =True,color='#db6565',alpha=0.5)
stacked_axs.set_xlim(0,1)


stacked_fig.savefig("adol_hist.pdf")