library("ape")
library("Biostrings")
library("ggplot2")
library("ggtree")
library("tidyverse")
library("treeio")

#beast trees 
tree = read.beast("concat41.tree")

#matK gene tree 70 taxa
#ggtree(root(tree,13), layout = )+ geom_tiplab(fontface='italic', size = 2)+ xlim(0, 3)+ geom_text2(size = 1.75,aes(label=round(as.numeric(posterior), 1), subset=as.numeric(posterior)> 0.5, x=branch, vjust=-1.1, hjust=1)) + geom_treescale(fontsize=1.75)

#annotated with fossil priors
#ggtree(tree)+ geom_tiplab(fontface='italic', size = 2)+ xlim(0, 35)+ geom_text2(size = 1.75,aes(label=round(as.numeric(posterior), 1), subset=as.numeric(posterior)> 0.5, x=branch, vjust=-1, hjust=2)) + geom_treescale(fontsize=1.75)

tree@data[isTip(tree, tree@data$node),]$length_0.95_HPD <- NA
# + geom_text(aes(label=node))
p <- ggtree(tree)+ geom_tiplab(fontface='italic', size = 2)+ scale_x_continuous(labels = abs)+ theme_tree2() + geom_text2(size = 1.75,aes(label=round(as.numeric(posterior), 1), subset=as.numeric(posterior)> 0.5, x=branch, vjust=-1, hjust=2)) + geom_treescale(fontsize=1.75) + geom_range(range='height_0.95_HPD', color='red', alpha=.6, size=2)
p <-revts(p)
p
                                                                                                                                                                                                                                                  
                                  
p <-flip(p, 45, 55)
p

p<-flip(p, 44, 68)
p


p<- rotate(p,56)
p
