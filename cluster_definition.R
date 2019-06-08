#Load packages
library("ape")
library("ggtree")
library("dplyr")
library("plyr")
library("reshape")
library("stringr")
library("tidyr")
library("ggplot2")
#library("Biostrings")
library("RColorBrewer")
library("phylobase")
library("data.table")
library("cluster")
#library("seqinr")
library("haplotypes")
# library("reshape2")
# library("igraph")
# library("phangorn")
# library("robustbase")
library("tidyverse")
library("dendextend")
library("e1071")


# Loading data

geo_alignment <- read.fas("K:/TBRU/GENERAL/Users/Gygli_Sebastian/bio/phd/projects/Georgia_genomes/combined/Keystone_Georgia_Ret_and_Pro_MDR_variable_positions.fasta")


# creating distance matrix

dist_geo_combined <- distance(geo_alignment, indels = "missing")

# clustering 

clusters_geo_combined <- agnes(dist_geo_combined, diss = TRUE, method = "average")

# cutting tree at 5 SNPs

clusters_geo_combined_5snps <- as.data.frame(cutree(as.hclust(clusters_geo_combined), h = 5))

colnames(clusters_geo_combined_5snps) <- "cluster_id"

# Reformating agnes output to aggregate the strain identifyers belonging to the same cluster

clusters_geo_combined_5snps <- as.data.frame(aggregate(rn ~ clusters_geo_combined_5snps$cluster_id, data =  setDT(clusters_geo_combined_5snps, keep.rownames = TRUE)[], toString))

#rename columns

colnames(clusters_geo_combined_5snps) <- c("cluster_id", "Gnr")

# Getting cluster list

results <- vector("list", length(clusters_geo_combined_5snps$cluster_id))

for (i in 1:length(clusters_geo_combined_5snps$cluster_id))
{
  results[[i]] <- length(unlist(strsplit(clusters_geo_combined_5snps$Gnr[i], ", ")))
}

results_vector <- as.vector(unlist(results))

final_clusters <- clusters_geo_combined_5snps[which(results_vector >=3),]
clusters_2_members <- clusters_geo_combined_5snps[which(results_vector ==2),]
final_clusters$cluster_id <- as.character(final_clusters$cluster_id)
final_clusters$Gnr <- as.character(final_clusters$Gnr)

list_clusters <- vector("list", length(final_clusters$cluster_id))

for (j in 1:length(final_clusters$cluster_id)) 
{
  list_clusters[[j]] <- unlist(strsplit(final_clusters$Gnr[j], ", "))
}

names(list_clusters) <- c(1:length(final_clusters$cluster_id))

# Make dataframe from list of final clusters
clusters_5snps <- data.frame(
  Taxa = unlist(list_clusters), Cluster_Nr = rep(names(list_clusters), lapply(list_clusters, length)), stringsAsFactors = FALSE)


# write.table(clusters_5snps, file = "K:/TBRU/GENERAL/Users/Gygli_Sebastian/bio/phd/projects/Georgia_genomes/combined/Keystone_Georgia_Ret_and_Pro_MDR_clusters", col.names = TRUE, quote = FALSE,
#             row.names = FALSE)


