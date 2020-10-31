library(vegan)
library(cluster)
set.seed(3)
setwd("~/bioinformatics_analysis/HMP_16S/Valencia/Cluster_testing/")
data_in = read.csv("allsamples_taxononomic_composition.csv",sep=",")

all.1k <- rrarefy(data_in[,7:ncol(data_in)], 1000)
all.1k <- all.1k[, order(colSums(all.1k), decreasing = T)]
all.1k_2 <- cbind(data_in[,1:6], all.1k)
nreads <- rowSums(all.1k)
relab.1k <- cbind(data_in[,1:6], all.1k/nreads)

source("vegdist.R")
dyn.load("vegdist.so")
#bray-curtis
bc1 <- vegdist(relab.1k[,7:ncol(relab.1k)], method ="bray", useShrinkage = TRUE)
hc1 <- hclust(bc1, method="ward.D")

#9 cuts for CSTs
nClrs <- 9
memb1 <- cutree(hc1,k=nClrs) 
memb1 <- as.factor(memb1)
table(memb1)

cluster1k <- cbind(memb1,cluster1k)
colnames(cluster1k)[1] <- "CST_9"

