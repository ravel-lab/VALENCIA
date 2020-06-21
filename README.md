![](https://github.com/ravel-lab/VALENCIA/blob/master/valencia_logo.png)
# Tutorial for VALENCIA

##### Author: Michael France

*Please contact the [author](mailto:mfrance@som.umaryland.edu?subject=[GitHub]%20VALENCIA) if you have any questions, suggestions, or if would like to contribute data to further develop VALENCIA. Thank you and Enjoy!*

&nbsp;
**VALENCIA** is a nearest-centroid based algorithm for the classification of human vaginal microbial community into state types based on their taxonomic composition. Samples are individually assigned to community state types (CSTs) based on their similarity to a set of thirteen reference centroids. The assignments provided by VALENCIA are reproducible and comparable across studies. Reference centroids were identified and defined using a dataset of >13,000 vaginal microbial taxonomic compositions from >1,900 North American women. We have done fairly extensive work to validate the usage of VALENCIA on vaginal microbial communities (see publication) including those that were derived from sequencing the V1V3, V3V4 and V4 regions. We have also looked at how it performs on samples from adolescent girls and post menopausal women, as well as on reproductive age African women. 

---
**Required dependencies** 

> Python 3.0+ with pandas, numpy, argparse, and matplotlib modules


Taxa names need to match used by VALENCIA for proper CST assignment. We use species level assignments for the *Lactobacillus*, *Gardnerella*, *Prevotella*, *Atopobium* and *Sneathia*. These appear as “Genus_species” format (e.g. “Lactobacillus_crispatus”. All other taxa are summarized to the genus or higher level. These appear as “g_Genus” or “f_Family” (e.g. “g_Bifidobacterium”). It is important that the major vaginal taxa’s names match the reference. This includes any critical in the definition of a CST, which appear below. *Prevotella* often causes problems due variations in naming conventions. If you are having difficulty with this taxa it is okay to combine all the data into "g_Prevotella", just make sure to use the reference centroid with the matching change .

---
**Running VALENCIA** 

VALENCIA has the follow required arguments: 

- -ref, --reference : path to the reference centroids file (provided)

- -i, --input : CSV file of the sample dataset with column 1 named “sampleID” containing unique sample names and column 2 named “read_count” containing each sample’s total read count. The remaining columns are expected to be taxa read counts, with the taxa name as the header.

- -o, --output : File prefix to store the output, all of the added CST information will appear as the final columns of the dataset

And the following optional arguments:
- -p, --plot : File prefix to store diagnostic plot in, not generated unless this argument is provided

The script is run as follows:
 
python3 /path/to/Valencia.py -ref /path/to/CST_profiles_012920.csv -i /path/to/test_dataset.csv -o /path/to/test_out -p /path/to/test
 
 
---
**CST architecture** 

-**CST I** communities are dominated by *L. crispatus* and include subtypes
  - CST **I-A** almost completely *L. crispatus*
  - CST **I-B** less *L. crispatus* but still majority

-**CST II** communities are dominated by *L. gasseri*

-**CST III** communities are dominated by *L. iners* and include subtypes
  - CST **III-A** almost completely *L. iners*
  - CST **III-B** less *L. iners* but still majority

-**CST IV** communities that have a low relative abundance of *Lactobacillus* spp. and include subtypes
  - CST **IV-A** contains a high to moderate relative abundance of BVAB1 and *G. vaginalis*
  - CST **IV-B** contains a high to moderate relative abundance of *G. vaginalis* and *A. vaginae* 
  - CST **IV-C** contains low relative abundances of *G. vaginalis*, BVAB1, and *Lactobacillus* spp. and includes
  
      - CST **IV-C0** relatively even community with *Prevotella* spp.
      - CST **IV-C1** dominated by *Streptococcus* spp.
      - CST **IV-C2** dominated by *Enterococcus* spp.
      - CST **IV-C3** dominated by *Bifidobacterium* spp.
      - CST **IV-C4** dominated by *Staphylococcus* spp.

-**CST V** communities are dominated by *L. jensenii*

---
**qiime2 conversion script**

This script takes ouput from qiime2 and converts it to a format suitable for VALENCIA. The two files expected as input are the ASV taxon names key and the ASV read count table. The expected output is a new table which contains the samples as rows and the condensed taxa as columns. This file can be used as input to VALENCIA. The second file provides a link between the original taxa names and the condensed names.

The script is run as follows:

python3 /path/to/convert_qiime.py /path/to/taxon_key.csv /path/to/asv_count_table.csv 


