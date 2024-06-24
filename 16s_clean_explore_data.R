# -----------------------------------------------------------------------------#
# Cleaning up the phyloseq object
# Author: Geoffrey Zahn Rebecca Evans
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     phyloseq v 1.40.0
#                     ShortRead v 1.54.0
#                     Biostrings v 2.64.0
#                     adegenet v 2 .1 10
#                     readxl v 1.4.1
#                     janitor v 2.1.0
#                     microbiome v 1.20.0
# -----------------------------------------------------------------------------#

#################################################################################
#                               Main workflow                                   #
#        Remove non-bacteria, chloroplast and mitochondrial sequences           #
#                                                                               #
#################################################################################

# SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(ShortRead); packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")
library(adegenet); packageVersion("adegenet")
library(readxl); packageVersion("readxl")
library(janitor); packageVersion("janitor")
library(microbiome); packageVersion("microbiome")

# source("./R/googlemap_styling.R")
# seed
set.seed(789)


setwd("C:/Users/becca/Dropbox/data_working/16S_soils")


# DATA ####
ps <- readRDS("./Output/ps_not-cleaned_w_tree.RDS") # change to non-phylogeny stuff

# check positive control efficacy...
pos <- ps %>% 
  subset_samples(control == "Positive")
y <- tax_table(pos)[which(taxa_sums(pos) > 0),5:7] %>% as.data.frame()


x <- as.data.frame(otu_table(pos))[,which(taxa_sums(pos) > 0)]
y$abund <- as.numeric(x[1,])

y %>% group_by(Genus) %>% summarize(total=sum(abund)) %>% arrange(desc(total)) %>% mutate(relabund=total/sum(total)) %>% 
  mutate(Genus=factor(Genus,levels=pluck(.,"Genus"))) %>% 
  ggplot(aes(x=Genus,y=relabund)) +
  geom_col() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=60,hjust=1)) +
  labs(y="Relative abundance",title="Bacterial positive control sample")
ggsave("./Output/16S_Bacterial_positive_control_taxonomy.png")

# REMOVE NON-BACTERIA ####
ps <- subset_taxa(ps, Kingdom == "Bacteria")
# tax <- tax_table(ps)

ps <- subset_taxa(ps,Class != "Chloroplast")
ps <- subset_taxa(ps, taxa_sums(ps) > 0)
ps <- subset_samples(ps, sample_sums(ps) > 0)


# Save DNA sequences apart from rownames (from subsetted ps object)
seqs <- taxa_names(ps)
seqs <- DNAStringSet(seqs)
saveRDS(seqs,"./Output/16S_ASV_reference_sequences.RDS")


saveRDS(ps, file = "./Output/16S_clean_phyloseq_object_tree.RDS")
