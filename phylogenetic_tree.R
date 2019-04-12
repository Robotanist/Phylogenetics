#Bioconductor has many packages for working with DNA sequences. Install the ShortRead package with

#source("http://bioconductor.org/biocLite.R")
#biocLite("ShortRead")

library(ShortRead)
library(DECIPHER)
library(phangorn)
library(muscle)

# reads in whole directory
fasta <- readFasta("/Users/solomon.champion/Desktop/FASTA/loci/matK/HPDL/")

#Then write the object out
writeFasta(fasta, "/Users/solomon.champion/Desktop/FASTA/loci/matK/matK.fasta")


fasta <- "/Users/solomon.champion/Desktop/FASTA/loci/matK/matK.fasta"

# load the sequences from the file
# change "DNA" to "RNA" or "AA" if necessary
seqs <- readDNAStringSet(fasta)

# look at some of the sequences (optional)
seqs

# nucleotide sequences need to be in the same orientation
# if they are not, then they can be reoriented (optional)
seqs <- OrientNucleotides(seqs)

# perform alignment with MUSCLE
aligned <- muscle(seqs)

aligned <- as.phyDat(aligned, type = "DNA", levels = NULL, return.index = TRUE)

# perform alignment with DECIPHER package
# aligned <- AlignSeqs(seqs)

# view the alignment in a browser (optional)
# BrowseSeqs(aligned, highlight=0)

# write the alignment to a new FASTA file (optional)
# writeXStringSet(aligned, file="Users/solomon.champion/Desktop/FASTA/loci/matK/matK_aligned.fasta")

# calculates pairwise distance
pairwise <- dist.ml(aligned, model = "F81", exclude = "none", bf = NULL, Q = NULL,
                    k = 1L, shape = 1)

# bootstrapping

treeUPGMA  <- upgma(pairwise)

fit <- pml(treeUPGMA, aligned)

fit <- optim.pml(fit, rearrangements="NNI")
set.seed(123)
bs <- bootstrap.pml(fit, bs=100, optNni=TRUE)
treeBS <- plotBS(fit$tree,bs)

# open tree in FigTree
write.tree(treeBS, "Users/solomon.champion/Desktop/FASTA/loci/matK/matK.tree")
