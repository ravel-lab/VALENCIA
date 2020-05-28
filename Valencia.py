#!/usr/local/bin/python3

#importing packages to be used with error handling
try:
    import pandas as pd
except:
    print("Required package pandas not available")
    exit()
try:
    import numpy as np
except:
    print("Required package numpy not available")
    exit()
try:
    import argparse
except:
    print("Required package argparse not available")
    exit()

#setting up expected arguments and help messages
parser = argparse.ArgumentParser(description="VALENCIA is a tool to classify vaginal microbial community into community state types")
#required arguments
required = parser.add_argument_group("Required arguments")
required.add_argument("-ref","--reference", help="Complete path to reference centroids file")
required.add_argument("-i", "--input", help="Path to input CSV file with the first column header 'sampleID', second column total read counts with header 'read_count', and remaining columns containing taxa read counts",required=True)
required.add_argument("-o","--output", help="Output csv file prefix")

#optional arguments
parser.add_argument("-p","--plot", default=None, help="File prefix for boxplots of similarty of each sample to assigned CST, if not provided, no plot will be generated")

#reading arguments into parser
args = parser.parse_args()

#removing the  warning that pandas has for when you do complicated things
pd.options.mode.chained_assignment = None  # default='warn'

#VALENCIA
###The purpose of this tool is to classify vaginal microbial communities into community state types (CSTs)
###in a standardized and repeatable fashion, clusters are based on the 13,000+ women data
###A total of 13 CSTs are considered: CST I-A, I-B, II, III-A, III-B, IV-A, IV-B, IV-C0, IV-C1, IV-C2, IV-C3, IV-C4, V 
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

#### defining fuction to determine the penalized similarity score, not currently in use
def penalized_simil_score(row):

    if row['subCST'] == 'I-A':
        row = row.drop('I-B_sim')
    elif row['subCST'] == 'I-B':
        row = row.drop('I-A_sim')
    elif row['subCST'] == 'III-A':
        row = row.drop('III-B_sim')
    elif row['subCST'] == 'III-B':
        row = row.drop('III-A_sim')
    row = row.drop('subCST')
    similarity_scores = list(row)
    similarity_scores.sort()
    score_len = len(similarity_scores)
    penalized_score = similarity_scores[score_len-1] * (similarity_scores[score_len-1]-similarity_scores[score_len-2]) ** (1./2)

    return penalized_score

#list of subCSTs 
CSTs = ['I-A','I-B','II','III-A','III-B','IV-A','IV-B','IV-C0','IV-C1','IV-C2','IV-C3','IV-C4','V']

#reading in the input CST centroids
try:
    reference_centroids = pd.read_csv(args.reference,sep=',')
except:
    print('Please provide a valid path to the reference centroids using -r')
    exit()

try:
#reading in table of samples to be tested against the centroids
    sample_data_OG = pd.read_csv(args.input,sep=',')
except:
    print('Please provide a valid path to the input file using -i')
    exit()

#checking if first two columns have appropriate headers
if list(sample_data_OG.columns[0:2]) != ['sampleID','read_count']:
    print('Input file expected to be a CSV with first two column headers: sampleID,read_count')
    exit()

#forcing the sample data to have the same number of columns as the reference and in the same order 
combined_data = pd.concat([sample_data_OG,reference_centroids], ignore_index=True,sort=False)
sample_data = combined_data[:-13].fillna(0)
sample_data = sample_data.drop(['sub_CST'],axis=1)
reference_centroids = combined_data.tail(13).fillna(0)
reference_centroids = reference_centroids.drop(["sampleID","read_count"],axis=1).set_index('sub_CST')

#converting all of the read counts to relative abundance data and adding back the first two columns
sample_data_rel = sample_data[sample_data.columns[2:]].div(sample_data['read_count'],axis=0)
sample_data_rel = pd.concat([sample_data[sample_data.columns[0:2]],sample_data_rel],axis=1)

#loop measuring the similarity of each sample to each subCST centroid using yue + clayon theta
for CST in CSTs:

    sample_data_OG['%s_sim' %(CST)] = sample_data_rel.apply(lambda x : yue_distance(x[2:],reference_centroids.loc[CST]), axis=1)

#outputting the acquired data with the new variability measure
#identify for each sample, which subCST was most similar, then correcting name of subCST to remove _sim
sample_data_OG['subCST'] = sample_data_OG.iloc[:,-13:].idxmax(axis=1)
sample_data_OG['subCST'] = sample_data_OG['subCST'].str.replace('_sim',"")

#applying function to calculate a penalized score, not currently in use
#sample_data_OG['penalized_score'] = sample_data_OG.apply(lambda x : penalized_simil_score(x[-14:]), axis=1)

#store the similarity between each sample and its as=ssigned 
sample_data_OG['score'] = sample_data_OG.iloc[:,-14:-1].max(axis=1)

#determine higher order CST assignment based on subCST assignment
sample_data_OG['CST'] = sample_data_OG['subCST'].replace({'I-A':'I','I-B':'I','III-A':'III','III-B':'III','IV-C0':'IV-C','IV-C1':'IV-C','IV-C2':'IV-C','IV-C3':'IV-C','IV-C4':'IV-C'})

#output the assignments in new CSV file
sample_data_OG.to_csv("%s.csv" %(args.output),index=None)

#plotting distributions of similarity scores for each subCST agaist that for the reference dataset
if args.plot != None:

    #loading matplotlib module
    try:    
        import matplotlib.pyplot as plt
    except:
        print("Required package matplotlib not available, either install or do not use --plot option")

    #defining figure size
    similarity_fig, similarity_axs = plt.subplots(1,1, figsize=(10,6), facecolor='w', edgecolor='k')
    similarity_fig.subplots_adjust(left=0.15,right=0.85,bottom=0.15,top=0.9,hspace = 0.4, wspace=0.4)
    
    #creating x axis location variables
    loc=0
    #creating blank list to store subCSTs that are actually present in the data
    CST_labels = list()
    
    #values for the average and stdev of the similarity scores of reference samples to their assigned CSTs
    ref_ave = {'I-A':0.995678,'I-B':0.858443,'II': 0.811734,'III-A':0.972983,'III-B':0.810466,'IV-A':0.718100,'IV-B':0.659592,'IV-C0':0.321432,'IV-C1':0.745991,'IV-C2':0.695203,'IV-C3':0.758735,'IV-C4':0.681502,'V':0.734973}
    ref_std = {'I-A':0.007615,'I-B':0.151048,'II':0.161987,'III-A':0.058698,'III-B':0.133990,'IV-A':0.150829,'IV-B':0.162577,'IV-C0':0.137558,'IV-C1':0.231891,'IV-C2':0.246322,'IV-C3':0.213653,'IV-C4':0.186462,'V':0.140276}

    #establishing the CST color scheme
    CST_color_scheme = {'I-A':'#ff6868','I-B':'#ffd4da','II':'#b4ff68','III-A':'#ffbc6b','III-B':'#e4a67b','IV-A':'#c1adec','IV-B':'#91a8ed',
                    'IV-C0':'#989898','IV-C1':'#ffc0cb','IV-C2':'#a8e5e5','IV-C3':'#9acc9a','IV-C4':'#800080','V':'#ffff71'}

    #plotting a boxplot, ref average and ref standard deviation for each subCST
    for CST in CSTs:
        
        #checking if subCST present in data
        if len(sample_data_OG[sample_data_OG['subCST']==CST].index) > 0:

            boxprops = dict(linewidth=1, color="k")
            medianprops = dict(linewidth=1,color="k")
            
            #plotting box plot of similarity distribution
            similarity_axs.boxplot(x=sample_data_OG[sample_data_OG['subCST'] == CST].score,positions=[loc-0.25],notch=True,widths=[0.25],patch_artist=True,boxprops=boxprops,medianprops=medianprops)
            #plotting the reference average and stdev
            similarity_axs.scatter(x=loc,y=ref_ave[CST],s=15,marker='d',c=CST_color_scheme[CST])
            similarity_axs.plot([loc,loc],[ref_ave[CST]+ref_std[CST],ref_ave[CST]-ref_std[CST]],c=CST_color_scheme[CST])
            CST_labels.append(CST)

            loc+=1
        
        else:
            continue

    #setting plot labels etc
    similarity_axs.set_xticklabels(CST_labels)
    similarity_axs.set_ylim(0,1)
    similarity_axs.set_xlabel("subCST")
    similarity_axs.set_ylabel("Similarity to assigned subCST")
    #saving plot
    similarity_fig.savefig("%s.pdf" %(args.plot))
