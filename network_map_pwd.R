#Bioconductor has many packages for working with DNA sequences. Install the ShortRead package with

#source("http://bioconductor.org/biocLite.R")
#biocLite("ShortRead")

library(ShortRead)
library(DECIPHER)
library(phangorn)
library(reshape)

# no Pattern, Whole Directory
fasta <- readFasta("/Users/solomon.champion/Desktop/FASTA/loci/ITS/HPDL")

#Then write the object out to fasta file

writeFasta(fasta, "/Users/solomon.champion/Desktop/FASTA/loci/ITS/ITS.fasta")

# specify the path to the FASTA file (in quotes)
fasta <- "/Users/solomon.champion/Desktop/FASTA/loci/ITS/ITS.fasta"

# load the sequences from the file
# change "DNA" to "RNA" or "AA" if necessary
seqs <- readDNAStringSet(fasta)

# look at some of the sequences (optional)
seqs

# nucleotide sequences need to be in the same orientation
# if they are not, then they can be reoriented (optional)
seqs <- OrientNucleotides(seqs)

# perform the alignment
aligned <- AlignSeqs(seqs)

# view the alignment in a browser (optional)
BrowseSeqs(aligned, highlight=0)

# write the alignment to a new FASTA file
writeXStringSet(aligned, file="/Users/solomon.champion/Desktop/FASTA/loci/ITS/ITS_aligned.fasta")



ML <- read.phyDat("/Users/solomon.champion/Desktop/FASTA/loci/ITS/ITS_aligned.fasta", format = "fasta", type = "DNA")

pairwise <- dist.ml(ML, model = "JC69", exclude = "none", bf = NULL, Q = NULL,
                    k = 1L, shape = 1)


treeUPGMA  <- upgma(pairwise)
# treeNJ  <- NJ(dm)

layout(matrix(c(1,2), 2, 1), height=c(1,2))
par(mar = c(0,0,2,0)+ 0.1)
plot(treeUPGMA, main="UPGMAtree")
write.tree(treeUPGMA, "/Users/solomon.champion/Desktop/FASTA/loci/ITS/ITS.tree")

#make pwd into matrix
pairwise_matrix <- as.matrix(pairwise)

write.csv(pairwise_matrix, "pairwise.csv")


##############################################
# Data Dictionary
# variable name = Species + Accession #
# data type = Continous value
# description = Pairwise Distance from Internally Transcribed Spacer Nucleotide Sequences

pairwise <- read.csv("/Users/solomon.champion/Desktop/pwd/pairwise.csv")
dist_m <- as.matrix(dist(pairwise))
dist_mi <- 1/dist_m # one over, as qgraph takes similarity matrices as input
library(qgraph)
jpeg('ITS.jpg', width=1000, height=1000, unit='px')
qgraph(dist_mi, layout='spring', vsize=3)
dev.off()


#################################################
