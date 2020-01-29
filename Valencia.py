#!/usr/bin/python

#importing packages to be used
import pandas as pd
import numpy as np
import sys

#removing the  warning that pandas has for when you do complicated things
pd.options.mode.chained_assignment = None  # default='warn'

#!/usr/bin/python

#ComSTaT (COMmunity Stat Type clAssifying Tool)
###The purpose of this tool is to classify vaginal microbial communities into community state types (CSTs)
###in a standardized and repeatable fashion, clusters are based on the 13,000+ women data
###A total of 13 CSTs are considered: CST I-A, I-B, II, III-A, III-B, IV-A, IV-B, IV-C0, IV-C1, IV-C2, IV-C3, IV-C4, V 

###The basic approach is to define the center of these CSTs, develop error profiles based on in cluster variation, 
###place new samples based on similarity to the centers, estimate confidence based on the error profile

####This script tests new samples based on similarity to already defined cluster centroids

#defining function to determine yue-clayton theta
def yue_distance(row, median):

    #creating a counting variable to index the median list
    taxon_count = 0
    
    #creating lists to iterativly store output
    median_times_obs = []
    median_minus_obs_sq = []
    
    #looping through the row and calculating product and difference squared between row data and median data
    for taxon_abund in row:

        #calculate p * q
        median_times_obs.append(median[taxon_count]*taxon_abund)

        #calculate p-q squared
        median_minus_obs_sq.append((median[taxon_count]-taxon_abund)**2)

        taxon_count += 1
    
    #calculate sum p* q
    product = np.nansum(median_times_obs)   

    #calculate sum p-q squared
    diff_sq = np.nansum(median_minus_obs_sq)
    #calculate yue_med_dist
    yue_med_dist = product / (diff_sq + product)
     
    #return the value of yue distance
    return yue_med_dist 

#### defining fuction to determine the penalized similarity score
def penalized_simil_score(row):

    if row['sim_subCST'] == 'I-A':
        row = row.drop('I-B_sim')
    elif row['sim_subCST'] == 'I-B':
        row = row.drop('I-A_sim')
    elif row['sim_subCST'] == 'III-A':
        row = row.drop('III-B_sim')
    elif row['sim_subCST'] == 'III-B':
        row = row.drop('III-A_sim')
    
    row = row.drop('sim_subCST')

    similarity_scores = list(row)
    similarity_scores.sort()

    score_len = len(similarity_scores)
    penalized_score = similarity_scores[score_len-1] * (similarity_scores[score_len-1]-similarity_scores[score_len-2]) ** (1./2)

    return penalized_score

#list of subCSTs 
CSTs = ['I-A','I-B','II','III-A','III-B','IV-A','IV-B','IV-C0','IV-C1','IV-C2','IV-C3','IV-C4','V']

#reading in the input CST centroids
reference_centroids = pd.read_csv(sys.argv[1],sep=',')

#reading in table of samples to be tested against the centroids
sample_data_OG = pd.read_csv(sys.argv[2],sep=',')

#forcing the sample data to have the same number of columns as the reference 
combined_data = pd.concat([sample_data_OG,reference_centroids], ignore_index=True,sort=False)
sample_data = combined_data[:-13].fillna(0)
sample_data = sample_data.drop(['sub_CST'],axis=1)
print(sample_data.shape)

reference_centroids = combined_data.tail(13).fillna(0)

#reference_centroids = reference_centroids.drop(['sampleID','CST','Sample_number_for_SRA','barcode','Type','Project','total_reads'],axis=1).set_index('sub_CST')
reference_centroids = reference_centroids.drop(["sampleID","read_count"],axis=1).set_index('sub_CST')


#converting all of the read counts to relative abundance data for use in determining community stability
sample_data_rel = sample_data[sample_data.columns[2:]].div(sample_data['read_count'],axis=0)

#sample_data_rel = sample_data[sample_data.columns[8:]].div(100,axis=0)
sample_data_rel = pd.concat([sample_data[sample_data.columns[0:2]],sample_data_rel],axis=1)
print(sample_data_rel)
print(reference_centroids)

#activate this line if already converted to relative abundance
#reference_input_rel = reference_input
#print(sample_data_rel)

for CST in CSTs:

    sample_data_OG['%s_sim' %(CST)] = sample_data_rel.apply(lambda x : yue_distance(x[2:],reference_centroids.loc[CST]), axis=1)

#outputting the acquired data with the new variability measure
#print(sample_data.columns[-13:])

sample_data_OG['sim_subCST'] = sample_data_OG.iloc[:,-13:].idxmax(axis=1)
sample_data_OG['sim_subCST'] = sample_data_OG['sim_subCST'].str.replace('_sim',"")
sample_data_OG['penalized_score'] = sample_data_OG.apply(lambda x : penalized_simil_score(x[-14:]), axis=1)
sample_data_OG['score'] = sample_data_OG.iloc[:,-15:-2].max(axis=1)


sample_data_OG['sim_CST'] = sample_data_OG['sim_subCST'].replace({'I-A':'I','I-B':'I','III-A':'III','III-B':'III','IV-C0':'IV-C','IV-C1':'IV-C','IV-C2':'IV-C','IV-C3':'IV-C','IV-C4':'IV-C'})

sample_data_OG.to_csv("%s_StR_CST.csv" %(sys.argv[2].split(".")[0]),index=None)














