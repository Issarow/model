library(gplots)
library(csv)
#library(pheatmap)
library(NMF)
library(RColorBrewer)

df <- read.csv("/home/cissarow/Desktop/chambusoData/modfyd.csv", sep="\t", header=T)
head(df)
##df$Sample.ID
hv <- as.matrix(df[2:24])
head(hv)
#rownames(hv) <- df$Sample.ID
annotation = data.frame(HIV = factor(df$HIV.status), Age = factor(df$Age.group), Cancer = factor(df$Tumour.stage))
#HIV = C("navy", "darkgreen")
#names(HIV) = c("negative", "positive")
#head(annotation)
aheatmap(hv, scale = "row",  annRow = annotation) 
#hv2

library(csv)
#library(gplots)
library(ComplexHeatmap)
library(RColorBrewer)
#library(NMF)

df <- read.csv("/home/cissarow/Desktop/chambusoData/modfyd.csv", sep="\t", header=T)
head(df)

rownames(df) <- df$Sample.ID
rownames(df)
df_matrix <- data.matrix(df[5:27])
head(df_matrix)
df2 <- scale(df_matrix )
head(df2)
annot_df <- data.frame(Age = df$Age.group, HIV = df$HIV.status, Tumar = df$Tumour.stage)
col = list( Age = c("Below 30" = "black", "30 to 40" = "violet", "Above 40" = "purple"), HIV = c("Negative" = "blue", "Positive" = "red"), Tumar = c("CIN 1" = "green", "CIN 2" = "skyblue", "CIN 3" = "orange", "Invasive" = "firebrick"))
annotation <- HeatmapAnnotation(annot_df, col = col)
#col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
Heatmap(df2, name = "HPV", top_annotation = annotation, width = unit(12, "cm"), show_row_names = FALSE, col = circlize::colorRamp2(c(0, 1), c("blue", "red")))



