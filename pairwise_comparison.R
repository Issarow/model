library(ape)
library(adegenet)
library(ggplot2)
library(reshape2)
#library(seriation)
library(igraph)
library(dplyr)

#MyseqDNA<-fasta2DNAbin("/home/chacha/Downloads/multiple_44fasta/alinged_33stbosh.fasta")

#MyseqDNA_L1<-fasta2DNAbin("/home/chacha/Documents/all178_outputs/alinged_L1.fasta")
MyseqDNA_L4<-fasta2DNAbin("/home/cissarow/Documents/isolates187/clust187/threshold0/alingnedthrezero.fasta")
#MyseqDNA_L4<-fasta2DNAbin("/home/cissarow/Documents/isolates187/cluster10_alingned.fasta")
MyseqDNA_L4

Myseq.distance1<-dist.dna(MyseqDNA_L4, model="N", pairwise.deletion=FALSE)
#Myseq.distance2<-dist.dna(MyseqDNA_L2, model="N", pairwise.deletion=FALSE)
#Myseq.distance4<-dist.dna(MyseqDNA_L4, model="N", pairwise.deletion=FALSE)
Myseq.distance1
mean(Myseq.distance1)
median(Myseq.distance1)
hist(Myseq.distance1, breaks = 2000, xlab = "Pairwise SNP distance", main = "", xlim = c(0,200)) #hist(Myseq.distance2)
#hist(Myseq.distance4)
#boxplot(Myseq.distance1, Myseq.distance2, Myseq.distance4, xlab="Lineages", ylab="Genetic distance")

Myseq.mat<-as.matrix(Myseq.distance1)
min(Myseq.mat)
#subset(Myseq.mat<4)

#Myseq.mat<5
hist(Myseq.mat)
#boxplot(Myseq.mat)
longData<-melt(Myseq.mat)
longData
ggplot(longData, aes(x = Var1, y = Var2)) + 
  geom_raster(aes(fill=value)) + #scale_fill_gradient(low="grey90", high="red") +
  labs(x="Y", y="X") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))

Myseq.mat <- as.data.frame(as.matrix(Myseq.distance1))
Myseq.mat
#table.paint(Myseq.mat, cleg = 0, clabel.row = 0.5, clabel.col = 0.5)
write.table(Myseq.mat, file="/home/cissarow/Documents/isolates187/cluster10matrix.txt", sep="\t")

