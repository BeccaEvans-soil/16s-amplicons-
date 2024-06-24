# setup
library(dada2)
library(tidyverse)
library(ShortRead)
FWD <- "GTGCCAGCMGCCGCGGTAA" # Sequence of FWD primer (515F)
REV <- "GGACTACHVGGGTWTCTAAT"  # Sequence of REV primer (806R)


setwd("C:/Users/becca/Dropbox/data_working/16S_soils")

# load ps object
ps <- readRDS("./Output/16S_ps_not-cleaned.RDS")

# get ASVs into a DNAStringSet
seqs <- taxa_names(ps) %>% DNAStringSet()

# try glomming taxa at genus level
genus <- tax_glom(ps,"Genus")
ntaxa(genus)

# export ASVs to fasta
names(seqs) <- paste0("ASV_",seq_along(seqs))
writeFasta(seqs,'ASVs.fasta')

# cluster at 97% using vsearch in terminal
# vsearch --cluster_fast ASVs.fasta --centroids ASV_centroids_97.fasta --id .97 --strand both --clusters ASV_97_clusters.fasta --relabel_self

# read in clustered version 
asvs_97 <- readFasta("./Output/ASV_centroids_97.fasta")


# merge taxa at 97% similarity clustered by vsearch, keeping centroid as reference sequence
clusters_to_merge <- list.files(pattern = "ASV_97_clusters.fasta.*.seqs")

# back up ps object
ps2 <- ps
# in for-loop, read in each cluster's sequences, and progressively merge them together, overwriting
# the ps object with each iteration
for(i in seq_along(clusters_to_merge)){
  ps <- merge_taxa(ps,readLines(clusters_to_merge[i])[readLines(clusters_to_merge[i]) %in% taxa_names(ps)])
}
# Check on what happened for sanity
ntaxa(ps2); ntaxa(ps) # got rid of ~30,000 ASVs by merging 
sum(sample_sums(ps2)); sum(sample_sums(ps)) # kept same number of overall reads (additive merge)

# save new version of ps object
saveRDS(ps,"./ps_preLULU_clustered_97.RDS")

ps3 <- subset_taxa(ps,taxa_sums(ps) >= 100)
ntaxa(ps3)

saveRDS(ps3,"./ps_preLULU_clustered_97_lowcountsremoved.RDS")





