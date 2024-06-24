
#make tree with 16s data

# Packages and functions ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(phangorn); packageVersion("phangorn")
# library(msa); packageVersion("msa")
library(ape); packageVersion("ape")
library(seqinr); packageVersion("seqinr")
library(DECIPHER); packageVersion("DECIPHER")

setwd("C:/Users/becca/Dropbox/data_working/16S_soils")
print("packages read in")

# Read in phyloseq object from first script output ####
ps3 <- readRDS("./Output/ps_preLULU_clustered_97_lowcountsremoved.RDS")


# simplify ASV names
seqs <- rownames(tax_table(ps3))
names(seqs) <- paste0("ASV_",1:length(seqs)) # This propagates to the tip labels of the tree

# create DNAStringSet object
seqs_StringSet <- DNAStringSet(seqs)

print("starting alignment")

# Multiple sequence alignment  ####
decipher_alignment <- DECIPHER::AlignSeqs(seqs_StringSet, processors = parallel::detectCores() - 1,verbose = TRUE) # DECIPHER method
adj_decipher_alignment <- DECIPHER::AdjustAlignment(decipher_alignment,processors = parallel::detectCores() - 1)


saveRDS(adj_decipher_alignment,"16S_dna_alignment_DECIPHER.RDS")

adj_decipher_alignment <- readRDS("16S_dna_alignment_DECIPHER.RDS")

print("alignment done")

# Convert to various formats
phang.align <- as.phyDat(as.character(adj_decipher_alignment), type = "DNA")
dnabin.align <- as.DNAbin(adj_decipher_alignment)


# distance - maximum likelihood ####
dm <- DistanceMatrix(adj_decipher_alignment,type = "dist",correction = "Jukes-Cantor",
                     processors = parallel::detectCores() - 1,verbose = TRUE)

dm <- readRDS("16S_ML_Distance.RDS")



# dm.TN93 <- dist.dna(dnabin.align,model = "TN93")
# missing <- which(is.na(dm.TN93))
# possible_positions <- round(missing / length(labels(dm.TN93)),0)
# 
# labels(dm.TN93)[possible_positions]

print("distances done")

#save
saveRDS(dm,"16S_ML_Distance.RDS")

# Initial neighbor-joining tree ####
treeNJ <- njs(dm) # Note, tip order != sequence order

# save progress
saveRDS(treeNJ, "16S_treeNJ.RDS")

treeNJ <- readRDS("./Output/16S_treeNJ.RDS")

g_NJ <-  plot_tree(treeNJ,"Neighbor Joining (NJ)")

# Estimate model parameters ####
#fit <-  pml(treeNJ, data=as.phyDat(muscle_alignment))

#fit = pml(treeNJ, data=phang.align)

#save
#saveRDS(fit,"16S_fit_treeNJ_muscle.RDS")


#fit$tree$tip.label <- seqs



treeNJ$tip.label #this is asv 1,2,3 etc good

#make taxa names match the tree tip names (asv1, 2, etc)
taxaname <- names(seqs) <- paste0("ASV_",1:length(seqs)) 
taxa_names(ps3) <- taxaname


# add tree to phyloseq object ####
ps5 <- phyloseq(tax_table(tax_table(ps3)),
                otu_table(otu_table(ps3)),
                sample_data(sample_data(ps3)),
                phy_tree(treeNJ))
 


# Save updated phyloseq object with tree

saveRDS(ps5, "ps_not-cleaned_w_tree.RDS")

print("added to phlyo")

