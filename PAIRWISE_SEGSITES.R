setwd("/Users/Tash/Dropbox (MMRU)/Backup/Papers/DELETIONS_KIPKORIR/ANALYSIS/PAIRWISE_COMP/FASTA_FILES/DEL")
library(ape)
library(igraph)
library(adegenet)
library(pheatmap)
library(pegas)
library (ggplot2)
library(reshape2)
library(scales)
library(Rmisc)
library(gdata)
library(gridExtra)
library(cowplot)
####Import files
AlldelFasta<-fasta2DNAbin("all.del.mafft.fasta")
AlldelFast.D<-dist.dna(AlldelFasta, model ="N", pairwise.deletion=FALSE)
AlldelFast.D.Mat<-as.matrix(AlldelFast.D)
write.table(AlldelFast.D.Mat, file="AllD.mafft5.txt", sep="\t")




Lin2DNA<-fasta2DNAbin("Lin2.fasta")
Lin2NegDNA<-fasta2DNAbin("Lin2Neg.fasta")
Lin2PosDNA<-fasta2DNAbin("Lin2Pos.fasta")
Lin4DNA<-fasta2DNAbin("Lin4.fasta")
Lin4NegDNA<-fasta2DNAbin("Lin4Neg.fasta")
Lin4PosDNA<-fasta2DNAbin("Lin4Pos.fasta")
Lin2Lin4NegDNA<-fasta2DNAbin("Lin2Lin4Neg.fasta")
Lin2Lin4PosDNA<-fasta2DNAbin("Lin2Lin4Pos.fasta")
NegDNA<-fasta2DNAbin("Neg.fasta")
PosDNA<-fasta2DNAbin("Pos.fasta")
###Distance matrices
AllD<-dist.dna(AllDNA, model="N", pairwise.deletion=TRUE)
Lin2D<-dist.dna(Lin2DNA, model="N", pairwise.deletion=TRUE)
Lin2NegD<-dist.dna(Lin2NegDNA, model="N", pairwise.deletion=TRUE)
Lin2PosD<-dist.dna(Lin2PosDNA, model="N", pairwise.deletion=TRUE)
Lin4D<-dist.dna(Lin4DNA, model="N", pairwise.deletion=TRUE)
Lin4NegD<-dist.dna(Lin4NegDNA, model="N", pairwise.deletion=TRUE)
Lin4PosD<-dist.dna(Lin4PosDNA, model="N", pairwise.deletion=TRUE)
Lin2Lin4NegD<-dist.dna(Lin2Lin4NegDNA, model="N", pairwise.deletion=TRUE)
Lin2Lin4PosD<-dist.dna(Lin2Lin4PosDNA, model="N", pairwise.deletion=TRUE)
NegD<-dist.dna(NegDNA, model="N",pairwise.deletion=TRUE)
PosD<-dist.dna(PosDNA, model="N",pairwise.deletion=TRUE)
AllMat<-as.matrix(AllD)
write.table(AllMat, file="Allmat.txt", sep="\t")
Lin2Mat<-as.matrix(Lin2D)
write.table(Lin2Mat, file="Lin2mat.txt", sep="\t")
Lin2NegMat<-as.matrix(Lin2NegD)
write.table(Lin2NegMat, file="Lin2Negmat.txt", sep="\t")
Lin2PosMat<-as.matrix(Lin2PosD)
write.table(Lin2PosMat, file="Lin2Posmat.txt", sep="\t")
Lin4Mat<-as.matrix(Lin4D)
write.table(Lin4Mat, file="Lin4mat.txt", sep="\t")
Lin4NegMat<-as.matrix(Lin4NegD)
write.table(Lin4NegMat, file="Lin4Negmat.txt", sep="\t")
Lin4PosMat<-as.matrix(Lin4PosD)
write.table(Lin4PosMat, file="Lin4Posmat.txt", sep="\t")
Lin2Lin4NegMat<-as.matrix(Lin2Lin4NegD)
write.table(Lin4NegMat, file="Lin2Lin4Negmat.txt", sep="\t")
Lin2Lin4PosMat<-as.matrix(Lin2Lin4PosD)
write.table(Lin4PosMat, file="Lin2Lin4Posmat.txt", sep="\t")
NegMat<-as.matrix(NegD)
write.table(NegMat, file="Negmat.txt", sep="\t")
PosMat<-as.matrix(PosD)
write.table(PosMat, file="Posmat.txt", sep="\t")

####FOR BOXPLOTS
AllMat<-as.matrix(AllD)
AllMatUTri<-upper.tri(AllMat)
AllMatUTriNo<-AllMat[AllMatUTri]
AllMatUTriNo.df<-as.data.frame(AllMatUTriNo)

Lin2Mat<-as.matrix(Lin2D)
Lin2MatUTri<-upper.tri(Lin2Mat)
Lin2MatUTriNo<-Lin2Mat[Lin2MatUTri]
Lin2MatUTriNo.df<-as.data.frame(Lin2MatUTriNo)

Lin4Mat<-as.matrix(Lin4D)
Lin4MatUTri<-upper.tri(Lin4Mat)
Lin4MatUTriNo<-Lin2Mat[Lin4MatUTri]
Lin4MatUTriNo.df<-as.data.frame(Lin4MatUTriNo)

Lin2NegMat<-as.matrix(Lin2NegD)
Lin2NegMatUTri<-upper.tri(Lin2NegMat)
Lin2NegMatUTriNo<-Lin2NegMat[Lin2NegMatUTri]
Lin2NegMatUTriNo.df<-as.data.frame(Lin2NegMatUTriNo)

Lin2PosMat<-as.matrix(Lin2PosD)
Lin2PosMatUTri<-upper.tri(Lin2PosMat)
Lin2PosMatUTriNo<-Lin2PosMat[Lin2PosMatUTri]
Lin2PosMatUTriNo.df<-as.data.frame(Lin2PosMatUTriNo)

Lin4NegMat<-as.matrix(Lin4NegD)
Lin4NegMatUTri<-upper.tri(Lin4NegMat)
Lin4NegMatUTriNo<-Lin4NegMat[Lin4NegMatUTri]
Lin4NegMatUTriNo.df<-as.data.frame(Lin4NegMatUTriNo)

Lin4PosMat<-as.matrix(Lin4PosD)
Lin4PosMatUTri<-upper.tri(Lin4PosMat)
Lin4PosMatUTriNo<-Lin4PosMat[Lin4PosMatUTri]
Lin4PosMatUTriNo.df<-as.data.frame(Lin4PosMatUTriNo)

NegMat<-as.matrix(NegD)
NegMatUTri<-upper.tri(NegMat)
NegMatUTriNo<-NegMat[NegMatUTri]
NegMatUTriNo.df<-as.data.frame(NegMatUTriNo)

PosMat<-as.matrix(PosD)
PosMatUTri<-upper.tri(PosMat)
PosMatUTriNo<-PosMat[PosMatUTri]
PosMatUTriNo.df<-as.data.frame(PosMatUTriNo)

####For statistical measure of diffs between pairwise genetic diffs
#### THis is the right test. Turned correct (coninuity correction) and paired, off
### Between All HIVpos and HIVneg
wilcox.test(PosMatUTriNo, NegMatUTriNo, correct=TRUE, paired=FALSE)
quantile(PosMatUTriNo)
median(PosMatUTriNo)

quantile(NegMatUTriNo)
median(NegMatUTriNo)
### Lin2Pos and neg
wilcox.test(Lin2PosMatUTriNo, Lin2NegMatUTriNo, correct=FALSE, paired=FALSE)
### Lin4Pos and neg
wilcox.test(Lin4PosMatUTriNo, Lin4NegMatUTriNo, correct=FALSE, paired=FALSE)
### Lin4 (all) and Lin2 (all)
wilcox.test(Lin2MatUTriNo, Lin4MatUTriNo, correct=FALSE, paired=FALSE)
### Lin4 (pos) and Lin2 (Pos)
wilcox.test(Lin4PosMatUTriNo, Lin2PosMatUTriNo, correct=FALSE, paired=FALSE)
### Lin4 (neg) and Lin2 (neg)
wilcox.test(Lin4NegMatUTriNo, Lin2NegMatUTriNo, correct=FALSE, paired=FALSE)


colnames(AllMatUTriNo.df)<-c('All')
colnames(Lin2NegMatUTriNo.df)<-c('Lin2Neg')
colnames(Lin2PosMatUTriNo.df)<-c('Lin2Pos')
colnames(Lin4NegMatUTriNo.df)<-c('Lin4neg')
colnames(Lin4PosMatUTriNo.df)<-c('Lin4Pos')
colnames(NegMatUTriNo.df)<-c('AllNeg')
colnames(PosMatUTriNo.df)<-c('AllPos')

###Name the column so can merge on something
####Will be the number of observations for each matrix - so diffs number
####NEEDS MANUAL CHANGING
colno<-c(1:14365)
AllMatUTriNo.df<-data.frame(AllMatUTriNo.df,colno)
colno<-c(1:120)
Lin2NegMatUTriNo.df<-data.frame(Lin2NegMatUTriNo.df,colno)
colno<-c(1:378)
Lin2PosMatUTriNo.df<-data.frame(Lin2PosMatUTriNo.df,colno)
colno<-c(1:1711)
Lin4NegMatUTriNo.df<-data.frame(Lin4NegMatUTriNo.df,colno)
colno<-c(1:1596)
Lin4PosMatUTriNo.df<-data.frame(Lin4PosMatUTriNo.df,colno)
colno<-c(1:3003)
NegMatUTriNo.df<-data.frame(NegMatUTriNo.df,colno)
colno<-c(1:4095)
PosMatUTriNo.df<-data.frame(PosMatUTriNo.df,colno)

View(PosMatUTriNo.df)
###Need to merge separately - or sequentially. Longest first for each seq
###Find lengthn by dim(mat)
###Check for NAs
Lin2NegPos.Merge.df<-merge(Lin2PosMatUTriNo.df, Lin2NegMatUTriNo.df, by="colno", all=TRUE)
Lin4NegPos.Merge.df<-merge(Lin4NegMatUTriNo.df, Lin4PosMatUTriNo.df, by="colno", all=TRUE)
AllNegPos.merge.df<-merge(PosMatUTriNo.df, NegMatUTriNo.df, by="colno", all=TRUE)
Lin2Lin4.merge.df<-merge(Lin4NegPos.Merge.df,Lin2NegPos.Merge.df,by="colno", all=TRUE)
All.merge.df<-merge(AllNegPos.merge.df,Lin2Lin4.merge.df,by="colno", all=TRUE)

####RESHAPE THE DATA FOR BOXPLOT
All.merge.df.m<-melt(All.merge.df,id.vars='colno', measure.vars=c('AllNeg','AllPos','Lin4Pos','Lin4neg','Lin2Pos','Lin2Neg'))
Lin2Lin4.merge.df.m<-melt(Lin2Lin4.merge.df,id.vars='colno', measure.vars=c('Lin4Pos','Lin4neg','Lin2Pos','Lin2Neg'))
NegPos.merge.df.m<-melt(AllNegPos.merge.df,id.vars='colno', measure.vars=c('AllPos','AllNeg'))
View(NegPos.merge.df.m)

###Reorder before plotting
All.merge.df.m$variable <- factor(All.merge.df.m$variable, levels = c("AllNeg","AllPos","Lin2Neg","Lin2Pos","Lin4neg", "Lin4Pos"))
Lin2Lin4.merge.df.m$variable <- factor(Lin2Lin4.merge.df.m$variable, levels = c("Lin2Neg","Lin2Pos","Lin4neg", "Lin4Pos"))
NegPos.merge.df.m$variable <- factor(NegPos.merge.df.m$variable, levels = c("AllNeg","AllPos"))


NegPos_boxplot<-ggplot() + 
  geom_boxplot(data=NegPos.merge.df.m, aes(x=variable, y=value, colour=variable, fill=variable)) +
  geom_point(data=NegPos.merge.df.m, aes(x=variable, y=value), alpha=0.4) +
  theme_classic() + 
  scale_y_continuous(limits=c(5,2000),
                     breaks=c(0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000),
                     labels=comma, name="Number of SNP differences") + 
  scale_color_manual(values=c("black", "black")) + 
  scale_x_discrete(labels=c("HIV \n uninfected", "HIV \n infected")) +
  scale_fill_manual(values=c("white", "grey")) +
  theme(axis.text.x=element_text(size=16)) +
  theme(axis.title.y=element_text(size=16)) +
  theme(axis.text.y=element_text(size=16)) 
NegPos_boxplot

All_boxplot<-ggplot() + 
  #geom_point(data=All.merge.df.m, aes(x=variable, y=value, colour=variable, fill=variable), alpha=0.5, position="jitter",alpha=0.4) +
  geom_boxplot(data=All.merge.df.m, aes(x=variable, y=value, colour=variable, fill=variable), alpha=0.75) +
  scale_y_continuous(limits=c(0,1905),
                     breaks=c(0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900),
                     labels=comma, name="Number of SNP differences") + 
  scale_colour_manual("variable",values=c("AllPos" = "black", "AllNeg"="black",
                      "Lin2Neg" = "blue", "Lin2Pos" = "dodgerblue",
                      "Lin4neg" = "firebrick", "Lin4Pos" = "orangered"), guide=FALSE) +
  scale_fill_manual("variable",values=c("AllPos" = "white", "AllNeg"="black",
                                          "Lin2Neg" = "blue", "Lin2Pos" = "white",
                                          "Lin4neg" = "firebrick", "Lin4Pos" = "white"), 
                    breaks=c("AllNeg", "AllPos", "Lin2Neg", "Lin2Pos","Lin4neg", "Lin4Pos"),
                    labels=c("All\nHIV uninfected", "All\nHIV infected","Lineage 2\nHIV uninfected",
                             "Lineage 2\nHIV infected",
                             "Lineage 4\nHIV uninfected","Lineage 4\nHIV uninfected")) +
  scale_x_discrete(labels=c("","","","","",""),name="") +
  theme(axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=12)) +
  theme(axis.text.y=element_text(size=12)) +
  theme(axis.ticks.x=element_blank()) +
  theme(legend.title.align=0.5) +
  theme(legend.key.width=unit(1,"line")) +
  theme(legend.key.height=unit(2,"line"))
All_boxplot

ggsave(file="All_boxplot_FINAL_LAND_LEG.pdf", All_boxplot, width =297, height =210, units=c("mm"), dpi=300)
ggsave(file="All_boxplot_FINAL_PORT_LEG.pdf", All_boxplot, width =210, height =297, units=c("mm"), dpi=300)


All.merge.df.m$variable<-factor(All.merge.df.m$variable,levels=c("AllNeg", "AllPos","Lin4neg", "Lin4Pos","Lin2Neg", "Lin2Pos"))

All_boxplot_axis<-ggplot() + 
  #geom_point(data=All.merge.df.m, aes(x=variable, y=value, colour=variable, fill=variable), alpha=0.5, position="jitter",alpha=0.4) +
  geom_boxplot(data=All.merge.df.m, aes(x=variable, y=value, colour=variable, fill=variable), alpha=0.75) +
  scale_y_continuous(limits=c(0,1905),
                     breaks=c(0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900),
                     labels=comma, name="Number of SNP differences") + 
  scale_colour_manual("variable",values=c("AllPos" = "black", "AllNeg"="black",
                                          "Lin2Neg" = "blue", "Lin2Pos" = "dodgerblue",
                                          "Lin4neg" = "firebrick", "Lin4Pos" = "orangered"), guide=FALSE) +
  scale_fill_manual("variable",values=c("AllPos" = "white", "AllNeg"="black",
                                        "Lin2Neg" = "blue", "Lin2Pos" = "white",
                                        "Lin4neg" = "firebrick", "Lin4Pos" = "white"), 
                    breaks=c("AllNeg", "AllPos", "Lin2Neg", "Lin2Pos","Lin4neg", "Lin4Pos"), guide=FALSE) +
  scale_x_discrete(labels=c("All\nHIV uninfected", "All\nHIV infected","Lineage 2\nHIV uninfected","Lineage 2\nHIV infected",
                                        "Lineage 4\nHIV uninfected","Lineage 4\nHIV infected"),name="") +
  theme(axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=12)) +
  theme(axis.text.y=element_text(size=12)) 
All_boxplot_axis

ggsave(file="All_boxplot_FINAL_LAND_AXIS.pdf", All_boxplot_axis, width =297, height =210, units=c("mm"), dpi=300)
ggsave(file="All_boxplot_FINAL_PORT_AXIS.pdf", All_boxplot_axis, width =210, height =297, units=c("mm"), dpi=300)




####SEGSITE
All_SS<-seg.sites(AllDNA)
Neg_SS<-seg.sites(NegDNA)
Pos_SS<-seg.sites(PosDNA)
Lin2Neg_SS<-seg.sites(Lin2NegDNA)
Lin2Pos_SS<-seg.sites(Lin2PosDNA)
Lin4Neg_SS<-seg.sites(Lin4NegDNA)
Lin4Pos_SS<-seg.sites(Lin4PosDNA)

All.df_SS<-as.data.frame(All_SS)
Neg.df_SS<-as.data.frame(Neg_SS)
Pos.df_SS<-as.data.frame(Pos_SS)
Lin2Neg.df_SS<-as.data.frame(Lin2Neg_SS)
Lin2Pos.df_SS<-as.data.frame(Lin2Pos_SS)
Lin4Neg.df_SS<-as.data.frame(Lin4Neg_SS)
Lin4Pos.df_SS<-as.data.frame(Lin4Pos_SS)

All<-ggplot() +  
  geom_line(data=Neg.df_SS, aes(x=Neg_SS), stat="density",
            adjust=0.1, colour="black", size=1) +
  geom_line(data=Pos.df_SS, aes(x=Pos_SS), stat="density",
                           adjust=0.1, colour="black",  size=1,linetype="twodash" ) +
  scale_y_continuous(expand = c(0, 0), name="Density", limits=c(0,7e-05), 
                     breaks=c(0,1e-05,2e-05,3e-05,4e-05,5e-05,6e-05,7e-05)) +
  scale_x_continuous(limits=c(0,30000), name="SNP index",
                     breaks=c(0,5000,10000,15000,20000,25000,30000)) +
  theme(axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=12)) +
  theme(axis.text.y=element_text(size=12)) +
  theme(axis.title.x=element_text(size=12)) 
All

Lin2<-ggplot() +  
  geom_line(data=Lin2Neg.df_SS, aes(x=Lin2Neg_SS), stat="density", 
            adjust=0.1, colour="blue", size=1) +
  geom_line(data=Lin2Pos.df_SS, aes(x=Lin2Pos_SS), stat="density", 
            adjust=0.1, colour="dodgerblue", size=1, linetype="twodash") +
  scale_y_continuous(expand = c(0, 0), name="Density", limits=c(0,7e-05), 
                     breaks=c(0,1e-05,2e-05,3e-05,4e-05,5e-05,6e-05,7e-05)) +
  scale_x_continuous(limits=c(0,30000), name="SNP index",
                     breaks=c(0,5000,10000,15000,20000,25000,30000)) +
  theme(axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=12)) +
  theme(axis.text.y=element_text(size=12)) +
  theme(axis.title.x=element_text(size=12)) 
Lin2

Lin4<-ggplot() +  
  geom_line(data=Lin4Neg.df_SS, aes(x=Lin4Neg_SS), stat="density", 
            adjust=0.1, colour="firebrick", size=1) +
  geom_line(data=Lin4Pos.df_SS, aes(x=Lin4Pos_SS), stat="density", 
            adjust=0.1, colour="orangered", size=1, linetype="twodash") +
  scale_y_continuous(expand = c(0, 0), name="Density", limits=c(0,7e-05), 
                     breaks=c(0,1e-05,2e-05,3e-05,4e-05,5e-05,6e-05,7e-05)) +
  scale_x_continuous(limits=c(0,30000), name="SNP index",
                     breaks=c(0,5000,10000,15000,20000,25000,30000)) +
  theme(axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=12)) +
  theme(axis.text.y=element_text(size=12)) +
  theme(axis.title.x=element_text(size=12)) 
Lin4

All_2_4<-grid.arrange(All, arrangeGrob(Lin2,Lin4, ncol=2), nrow=2)

###Save plots
ggsave(file="DENS_ALL_2_4.pdf", All_2_4, width =210, height =297, units=c("mm"), dpi=300)
ggsave(file="DENS_ALL_2_4..jpeg", All_2_4, width =210, height =297, units=c("mm"), dpi=300)
ggsave(file="DENS_ALL_2_4..tiff", All_2_4, width =210, height =297, units=c("mm"), dpi=300)

All_2<-ggplot() +  
  geom_line(data=Neg.df_SS, aes(x=Neg_SS), stat="density",
            adjust=0.2, colour="black", size=1) +
  geom_line(data=Pos.df_SS, aes(x=Pos_SS), stat="density",
            adjust=0.2, colour="black",  size=1,linetype="twodash" ) +
  scale_y_continuous(expand = c(0, 0), name="Density", limits=c(0,5e-05), 
                     breaks=c(0,1e-05,2e-05,3e-05,4e-05,5e-05)) +
  scale_x_continuous(limits=c(1,30000), name="SNP index",
                     breaks=c(1,4000,8000,12000,
                              16000,20000,24000,28000)) +
  theme(axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=12)) +
  theme(axis.text.y=element_text(size=12)) +
  theme(axis.title.x=element_text(size=12)) 
All_2

ggsave(file="ALL_ONLY.pdf", All_2, width =297, height =210, units=c("mm"), dpi=300)
ggsave(file="ALL_ONLY_LAND.pdf", All_2, width =210, height =297, units=c("mm"), dpi=300)
ggsave(file="ALL_ONLY.jpeg", All_2, width =297, height =210, units=c("mm"), dpi=300)
ggsave(file="ALL_ONLY.tiff", All_2, width =297, height =210, units=c("mm"), dpi=300)

