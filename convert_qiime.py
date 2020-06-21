#!/usr/bin/env python3
import pandas as pd
import sys


#reading in the qiime taxon key
taxon_key = pd.read_csv(sys.argv[1],sep=",",index_col=0)

taxon_key.columns = ['k','p','c','o','f','g','s']
taxon_key = taxon_key[taxon_key.columns[::-1]]

#replacing taxon_key with 
#reading in the table of counts
counts_table = pd.read_csv(sys.argv[2],sep=",",index_col=0)

#function that determines the highest level of taxonimc specifity and then formats the condensed name
#only provides species assignments for focal taxa used by valencia
def taxon_condense(row):

    row = row.T

    first_nonan = row.first_valid_index()

    if first_nonan == 's':
        if row.loc[['g']].values[0] in ['Lactobacillus','Prevotella','Gardnerella','Atopobium','Sneathia']:
            taxon_name = "%s_%s" %(row.loc[['g']].values[0],row.loc[['s']].values[0])
        else:
            taxon_name = "g_%s" %(row.loc[['g']].values[0])
    elif first_nonan == 'g':
        taxon_name = "g_%s" %(row.loc[['g']].values[0])
    elif first_nonan == 'f':
        taxon_name = "f_%s" %(row.loc[['f']].values[0])
    elif first_nonan == 'o':
        taxon_name = "o_%s" %(row.loc[['o']].values[0])
    elif first_nonan == 'c':
        taxon_name = "c_%s" %(row.loc[['c']].values[0])
    elif first_nonan == 'p':
        taxon_name = "p_%s" %(row.loc[['p']].values[0])
    elif first_nonan == 'k':
        taxon_name = "k_%s" %(row.loc[['k']].values[0])
    else:
        taxon_name = "None"

    return taxon_name

#applying function to each row of the taxa key file
taxon_key['taxa'] = taxon_key.apply(lambda x : taxon_condense(x), axis=1)

#manual correction of names, these should be checked by looking at the ASV sequences and see how they match to the new name 
taxon_key['taxa'] = taxon_key['taxa'].replace({'g_Gardnerella':'Gardnerella_vaginalis','Lactobacillus_acidophilus/casei/crispatus/gallinarum':'Lactobacillus_crispatus'
                                                ,'Lactobacillus_fornicalis/jensenii':'Lactobacillus_jensenii','g_Escherichia/Shigella':'g_Escherichia.Shigella'
                                                ,'Lactobacillus_gasseri/johnsonii':'Lactobacillus_gasseri'})

#creating a dataframe for merging with just the information from the new condense column
taxon_merge = taxon_key[['taxa']]
#merging the counts table with the taxa table
counts_table_named = pd.merge(left=taxon_merge,right=counts_table,right_index=True,left_index=True,how="inner")
#grouping asvs with the same name and summing
counts_table_named = counts_table_named.groupby('taxa').sum()
#transposing to a table with samples are rows and counts as columns
counts_table_named = counts_table_named.T

#sorting the table by the study wide read count for each taxa
counts_table_named = counts_table_named.reindex(counts_table_named.sum().sort_values(ascending=False).index, axis=1)
#summing the read counts for each sample to be used by valencia in calculation of relative abundance
counts_table_named['read_count'] = counts_table_named.sum(axis=1)
#moving read count to first column
read_count_column = counts_table_named.pop('read_count')
counts_table_named.insert(0,'read_count',read_count_column)

#output two files, one with the new condenses table and one that matches the condensed taxa name back to the original ASV
counts_table_named.to_csv("taxon_table_asv_merged.csv",sep=",",index_label="sampleID")
taxon_key.to_csv("asv_condensed_taxa_names.csv",sep=",")