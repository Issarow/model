library(ggtree)
library(ape)
library(ggplot2)
library("phytools")
#col = c("0/0"="grey90", "1/1"="darkblue", "0/1"="red", "0/2"="blue", "1/2"="green" , "2/2"="brown", "0/3"="orange", "1/3"="violet", "2/3"="yellow", "3/3"="lightblue") 

Tree_File = "/home/cissarow/Documents/isolates187/RAxML_187bestTree.TEST"


Tree = read.tree(Tree_File)
Tree <- root(Tree,"G25780",resolve.root= FALSE)     #reroot tree using ancestor

Tree2 <- ladderize(Tree)
p <- ggtree(Tree2, layout = "rectangular", size=.6) #+ 
  #geom_text2(aes(subset=!isTip, label=node), hjust=-.3) 
p <- groupClade(p, c(202, 197, 200))
p <- ggtree(p, aes(color=group)) + scale_color_manual(values=c("black", "blue", "red", "violet")) + 
  #geom_hilight(node=202, fill="blue", alpha=.3) +
  #geom_hilight(node=197, fill="red", alpha=.3) + 
  #geom_hilight(node=200, fill="violet", alpha=.5) +
  geom_treescale(fontsize=4, x=0.1, y=-1, offset=-3)
#p

DR_file = "/home/cissarow/Documents/isolates187/mutations187.txt"
DR_annotation = read.table(DR_file,sep='\t',header=T,stringsAsFactors = F)
#head(DR_annotation)

#pdf("/home/cissarow/Documents/isolates187/mutation187_tree.pdf")
gheatmap(p, DR_annotation, offset=0.003, width = 3.0, font.size = 2.6, hjust = -0.02, colnames=TRUE, colnames_position = "top",colnames_angle = 60) + 
  #scale_fill_manual(values=c("0/0"="gray90", "1/1"="darkblue"), labels= c("No mutations", "Mutations")) +
  scale_fill_manual(values=c("white", "darkblue", "forestgreen", "skyblue", "orange" , "violet", "brown", "purple", "greenyellow", "purple", "black",  "pink"), breaks = c("0/0", "1/1", "0/1", "0/2", "1/2" , "2/2", "0/3", "1/3", "2/3", "1/4", "2/4", "3/3"),  labels = c("", "RIF", "INH", "EMB", "STR", "PZA", "ETO", "FLQ", "AMK", "KM", "Pre-XDR/XDR", "HIV_POS")) +
  
  geom_tiplab(align=TRUE, linesize=.2, size=1.2, color="black") +
  theme(legend.position="none") #+ geom_hilight(node=35, fill="purple") + theme(legend.position=c(0.1,0.7))


#head(DR_annotation)

#pdf("/home/cissarow/Documents/isolates187/mutationXDR187_tree.pdf")
#gheatmap(p, DR_annotation, offset=0.008, width = 0.6, font.size = 2.5, hjust = 0.02, colnames=TRUE, colnames_position = "top",colnames_angle = 60) #+ 
  #scale_fill_manual(values=c("0/0"="gray90", "1/1"="darkblue"), labels= c("No mutations", "Mutations")) +
  #scale_fill_manual(values=c("white", "darkblue", "forestgreen", "blue", "orange" , "magenta", "brown", "skyblue", "greenyellow", "black"), breaks = c("0/0", "1/1", "0/1", "0/2", "1/2" , "2/2", "0/3", "1/3", "2/3", "0/4"),  labels = c("", "RIF", "INH", "EMB", "STR", "PZA", "ETO", "FLQ", "AMK", "Pre-XDR")) +
  
  dev.off()

