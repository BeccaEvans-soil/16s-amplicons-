# -----------------------------------------------------------------------------#
# 16S DADA2 pipeline
# Processing raw amplicon reads
# Author: Geoffrey Zahn Becca Evans
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     dada2 v 1.24.0
#                     decontam v 1.16.0
#                     phyloseq v 1.40.0
#                     Biostrings v 2.64.0
#                     patchwork v 1.1.1
#                     readxl v 1.4.1
#                     janitor::clean_names() v 2.1.0
# -----------------------------------------------------------------------------#

#################################################################################
#                               Main workflow                                   #
# Filter and trim, denoise, sample inferrence, chimera and contaminant removal, # 
# taxonomic assignment, combine sequence table and metadata                     #
#################################################################################

# PACKAGES, SCRIPTS, AND SETUP ####

# why each package (put in onboarding document)
library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(patchwork); packageVersion("patchwork")
library(readxl); packageVersion("readxl")

# PARSE FILE PATHS ####

# File parsing - 

path <- "C:/Users/rebecca.c.evans/Dropbox/data_working/16S_soils/new" # CHANGE to the directory containing your adaptor-free demultiplexed fastq files when using your own data
filtpath <- file.path(path, "filtered_new") # Filtered files go into the filtered/ subdirectory
if(!file_test("-d", filtpath)) dir.create(filtpath) # make directory for filtered fqs if not already present


# Here, the pattern is set to match the forward reads tagged as "...pass_1.fastq.gz" in their filenames; 
# "...pass_2.fastq.gz" for reverse reads
list.files(path)
# Your data may differ, using "F" and "R" in the filenames, or something similar..
# Be sure to change that pattern to fit your own files when using your own data
fns <- sort(list.files(file.path(path), full.names = TRUE, pattern = "_R1_001.fastq")) # make pattern match your FWD reads the star means anything else after
rns <- sort(list.files(file.path(path), full.names = TRUE, pattern = "_R2_001.fastq"))



sample.names <- basename(fns) %>% str_remove("_L001_R1_001.fastq")
# sample.names <- basename(fns) %>% str_split("_") %>% map_chr(1)

# visualize a couple of fwd read quality profiles to help select reasonable filtration parameters
# you can select any number of files here...
# as-is, this just shows the Fwd and Rev quality profiles for the 1st and 2nd files
p1 <- plotQualityProfile(fns[1:2]) + ggtitle("Example forward reads")
p2 <- plotQualityProfile(rns[1:2]) + ggtitle("Example reverse reads")

# display and save the plots
p1 / p2
ggsave("./Output/figs/unfiltered_quality_plots.png",dpi=300,height = 6,width = 6)

# FILTER AND TRIM ####

# here, we decide what filenames we want our filtered reads to be called
# in this case, we use the base name of the sample and save the filtered files into the "filtered_new" subdirectory
filts_f <- file.path(path, "filtered_new", paste0(sample.names, "_FWD_filt.fastq.gz"))
filts_r <- file.path(path, "filtered_new", paste0(sample.names, "_REV_filt.fastq.gz"))

# this is the actual qualit control step
# These values are informed by our quality plots
out <- filterAndTrim(fns, filts_f, rns, filts_r, # input and output file names as denoted above
                     maxN=0, # uncalled bases are currently not supported in dada2
                     maxEE=c(2,2), # refers to the maximum expected errors allowed
                     truncQ=2, # special value denoting "end of good quality sequence" (optional)
                     rm.phix=TRUE, # automatically remove PhiX spike-in reads from sequencing center
                     truncLen = c(250,200), # refers to the lengths at which to truncate Fwd and Rev reads, respectively 
                     compress=TRUE, # compress output files with gzip
                     multithread=FALSE) # On Windows set multithread=FALSE

# save the filtration info in case we need it later
saveRDS(out, "./Output/16S_trackreads.RDS")
# out <- readRDS("./output/16S_trackreads.RDS")

# In case some samples may have had zero reads pass QC, reassign filts
filts_f <- sort(list.files(filtpath, full.names = TRUE,pattern = "FWD"))
filts_r <- sort(list.files(filtpath, full.names = TRUE,pattern = "REV"))
sample.names <- basename(filts_f) %>% str_remove("_FWD_filt.fastq.gz")


# sanity check  comparison of before and after filtration
p3 <- plotQualityProfile(rns[1:2]) + ggtitle("Unfiltered")
p4 <- plotQualityProfile(filts_r[1:2])+ ggtitle("Filtered") + coord_cartesian(xlim = c(0,300))
p3 / p4
ggsave("./Output/figs/16S_filtered_quality_comparison.png",dpi=300,height = 6,width = 6)

# LEARN ERROR RATES ####

# learn errors
set.seed(123) # "random" seed for reproducibility
errF <- learnErrors(filts_f, multithread=FALSE, MAX_CONSIST = 20,verbose = 1) # set multithread = FALSE on Windows
errR <- learnErrors(filts_r, multithread=FALSE, MAX_CONSIST = 20,verbose = 1) # set multithread = FALSE on Windows

saveRDS(errF,"./Output/16S_errF.RDS")
saveRDS(errR,"./Output/16S_errR.RDS")
errF <- readRDS("./Output/16S_errF.RDS")
errR <- readRDS("./Output/16S_errR.RDS")

# sanity check for error model
# explain what to look for in the plot! 
plotErrors(errF, nominalQ=FALSE)
ggsave("./Output/figs/16S_error_model.png",dpi=300,height = 6,width = 6)
plotErrors(errR, nominalQ=FALSE)

###########################################################
# IF DATA WON'T FIT IN MEMORY
# Can do derep and dada one sample at a time in a for-loop
names(filts_f) <- sample.names
names(filts_r) <- sample.names
dada_f <- vector("list",length(sample.names))
names(dada_f) <- sample.names
dada_r <- vector("list",length(sample.names))
names(dada_r) <- sample.names
mergers <- vector("list",length(sample.names))
names(mergers) <- sample.names

for(sam in sample.names){
  cat("Processing: ", sam, "\n")
  derep_f <- derepFastq(filts_f[[sam]])
  derep_r <- derepFastq(filts_r[[sam]])
  dada_f[[sam]] <- dada(derep_f,err = errF,multithread=TRUE)
  dada_r[[sam]] <- dada(derep_r,err = errR,multithread=TRUE)
  class(dada_f)
  merger <- mergePairs(dada_f[[sam]],derep_f,dada_r[[sam]],derep_r)
  mergers[[sam]] <- merger
}
##############################################################

# dereplication
 derepF <- derepFastq(filts_f, verbose=TRUE)
 saveRDS(derepF,"./Output/16S_derepF.RDS")
 derepR <- derepFastq(filts_r, verbose=TRUE)
 saveRDS(derepR,"./Output/16S_derepR.RDS")
 
 #derepF <- readRDS("./Output/16S_derepF.RDS")
 #derepR <- readRDS("./Output/16S_derepR.RDS")

# Name the derep-class objects by the sample names
if(identical(map_chr(strsplit(basename(filts_f), "_FWD_filt"), 1), map_chr(strsplit(basename(filts_r), "_REV_filt"), 1))){
  names(derepF) <- sample.names
  names(derepR) <- sample.names
} else {
  stop("Make sure fwd and rev files are in same order!")
}  

 
 ####HPC
 #files needed errF errR derepF, derepR
 
 library(tidyverse); packageVersion("tidyverse")
 library(dada2); packageVersion("dada2")
 library(decontam); packageVersion("decontam")
 library(phyloseq); packageVersion("phyloseq")
 library(Biostrings); packageVersion("Biostrings")
 library(patchwork); packageVersion("patchwork")
 
 print("loaded libraries")
 
 # DADA SAMPLE INFERRENCE #####
 
 derepF <- readRDS("16S_derepF.RDS")
 derepR  <- readRDS("16S_derepR.RDS")
 
 print("loaded RDS")
 
 errF <- readRDS("16S_errF.RDS")
 errR <- readRDS("16S_errR.RDS")
 
 print("loaded err")
 
 dadaFs <- dada(derepF, err=errF, multithread=FALSE, selfConsist = TRUE, pool = "pseudo") # set multithread = FALSE on Windows
 
 dadaRs <- dada(derepR, err=errR, multithread=FALSE, selfConsist = TRUE, pool = "pseudo") # set multithread = FALSE on Windows
 
 saveRDS(dadaFs,"16S_dadaFs.RDS")
 saveRDS(dadaRs,"16S_dadaRs.RDS")
 
 
 #################################
 
 

# SAMPLE INFERRENCE ####
#dadaFs <- dada(derepF, err=errF, multithread=FALSE, selfConsist = TRUE, pool = "pseudo") # set multithread = FALSE on Windows
 #dadaRs <- dada(derepR, err=errR, multithread=TRUE, selfConsist = TRUE, pool = "pseudo") # set multithread = FALSE on Windows
#saveRDS(dadaFs,"Output/16S_dadaFs.RDS")
 #saveRDS(dadaRs,"Output/16S_dadaRs.RDS")

# rm(derepF);rm(derepR)

 ############################
 #start here after HPC 
 #load in files
 
dadaF <- readRDS("16S_dadaFs.RDS")
dadaR <- readRDS("16S_dadaRs.RDS") 
 
# MERGE FWD and REV READS ####
#mergers <- mergePairs(dadaF, filts_f, dadaR, filts_r, verbose=TRUE)

#this one works and fits dada protocol 
#use this 
merg <- mergePairs(dadaF, derepF, dadaR, derepR, verbose=TRUE)

#rename
mergers <- merg


  
# MAKE SEQUENCE TABLE ####
seqtab <- makeSequenceTable(mergers)

# REMOVE CHIMERAS ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
saveRDS(seqtab.nochim,"./Output/16S_seqtab.nochim.RDS")
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)



############START HPC again########
#mergers and chimera take too long, do on HPC
#code for HPC

library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(patchwork); packageVersion("patchwork")

print("loaded libraries")


#read in data
derepF <- readRDS("C:/Users/becca/Dropbox/data_working/16S_soils/new/dada_derep/16S_derepF.RDS")
derepR  <- readRDS("C:/Users/becca/Dropbox/data_working/16S_soils/new/dada_derep/16S_derepR.RDS")

print("read in derep")

dada_F <- readRDS("C:/Users/becca/Dropbox/data_working/16S_soils/new/dada_derep/16S_dadaFs.RDS")
dada_R <- readRDS("C:/Users/becca/Dropbox/data_working/16S_soils/new/dada_derep/16S_dadaRs.RDS")

print("read in dada")


#merge files

mergers <- mergePairs(dada_F, derepF, dada_R, derepR, verbose=TRUE)
saveRDS(mergers, "./Output/mergers.RDS")


#make seqtab

seqtab <- makeSequenceTable(mergers)

print("seq table made")

# REMOVE CHIMERAS ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim,"seqtab.nochim.RDS")

print("done with chimeras")

#################DONE WITH HPC################

##START AGAIN######

#read in seq tab no chim

seqtab.nochim <- readRDS("./Output/16S_seqtab.nochim.clean.RDS")
mergers <- readRDS("./Output/mergers.RDS")

# reassign "out" to remove any missing reads
out <- readRDS("./Output/16S_trackreads.RDS")
out <- out[as.data.frame(out)$reads.out > 0,]

# TRACK READS THROUGH PIPELINE ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dada_F, getN), sapply(dada_R, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names %>% str_split("_") %>% map_chr(1)
track = as.data.frame(track)
track$total.loss.proportion = (track[,1]-track$nonchim)/track[,1]
head(track)

write.csv(track, file = "./Output/16S_read_counts_at_each_step.csv", row.names = TRUE)

# IMPORT METADATA ####

# import and clean
#meta <- read.csv("C:/Users/becca/Dropbox/data_working/16S_soils/new/meta_nwere_2018_2023.csv") %>% 
 # janitor::clean_names() 
#%>%  dplyr::filter(amplicon == "16S") # just bacteria samples for this


# subset to match seq table sample names
#meta <- meta[meta$sample_name %in% (sample.names %>% str_split("_") %>% map_chr(1)), ]
#row.names(seqtab.nochim) <- row.names(seqtab.nochim) %>% str_split("_") %>% map_chr(1)
# reorder metadata to match seqtab
#df <- data.frame(seqtab_rows=row.names(seqtab.nochim),
                # sample_name=row.names(seqtab.nochim))

#df2 <- left_join(meta,df,by="sample_name")
#row.names(df2) <- df2$sample_name
#meta <- df2[row.names(seqtab.nochim),]
#row.names(meta) <- meta$sample_name
#identical(row.names(meta),row.names(seqtab.nochim))



#############
############


# IMPORT METADATA ####
meta<- read.csv("./Output/meta_nwere_2018_2023.csv")

str(meta) #sample id must be character

meta$n.treat <- factor(meta$n_treat)
meta$plot.unq <- factor(meta$plot_unq)
meta$block <- factor(meta$block)
str(meta)

#bind asv with meta

#give names
rownames <- rownames(seqtab.nochim)

print(rownames)

row.names(seqtab.nochim) <- c("E10S","E11S", "E12S", "E13S","E14S","E15S",
                             "E16S","E17S","E18S","E1S","E2S","E3S","E4S","E5S",
                             "E6S","E7S","E8S","E9S","H10S","H11S","H12S",
                             "H1S","H2S","H3S","H4S","H5S","H6S","H7S","H8S",
                             "H9S","N10S","N11S","N12S","N1S","N2S","N3S",
                             "N4S","N5S","N6S","N7S","N8S","N9S",
                             "W10S","W11S","W12S","W13S","W14S","W15S",
                             "W16S","W17S","W18S","W1S","W2S","W3S","W4S",
                             "W5S","W6S","W7S","W8S","W9S")



#back to cleaning up seqtab
                               
# Remove all seqs with fewer than 100 nucleotides (if any) ####
keeper_esvs <- nchar(names(as.data.frame(seqtab.nochim))) > 99
seqtab.nochim <- seqtab.nochim[,keeper_esvs]

# remove newly singleton taxa
seqtab.nochim <- seqtab.nochim[,colSums(seqtab.nochim) > 1]

# save cleaned up seqtab
saveRDS(seqtab.nochim,"./Output/16S_seqtab.nochim.clean.RDS")

##############
copyseqtab.nochim <- seqtab.nochim
Sample_ID <- rownames(copyseqtab.nochim)
data <- cbind(Sample_ID,copyseqtab.nochim)

seqtab.nochim <-data 
saveRDS(seqtab.nochim,"seqtab.nochim.samplenames.RDS")
seqtab_nochim <- readRDS("seqtab.nochim.samplenames.RDS")

#meta with ASV by sample
meta_asv_notclean<- left_join(meta,seqtab.nochim,by="Sample_ID", copy=TRUE)

#####################
#SKIP#
# Find and remove contaminants ####
# Find and remove contaminants ####
contams = isContaminant(seqtab.nochim, neg = (meta$control == "Negative"), normalize = TRUE)
contams[contams$contaminant,]
table(contams$contaminant) # how many taxa are contaminants?
write.csv(contams, file = "./Output/16S_likely_contaminants.csv", row.names = TRUE)

# remove contaminant sequences and control samples from both tables, respectively ####
seqtab.nochim = seqtab.nochim[,(which(contams$contaminant != TRUE))]
seqtab.nochim = seqtab.nochim[(meta$control != "Negative"),]
meta = meta[(meta$control != "Negative"),]
dim(seqtab.nochim)
##################
 getwd()


# ASSIGN TAXONOMY ####

##########################
###HPC

library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(patchwork); packageVersion("patchwork")
library(readxl); packageVersion("readxl")
library(seqinr)

print("libraries loaded")

#read in data
seqtab.nochim <- readRDS("./Output/16S_seqtab.nochim.clean.RDS")

taxa <- read.fasta(file="rdp_train_set_16.fa.gz")

print("seq tab and taxa loaded")


# Use RDP training set for 16S
taxa <- assignTaxonomy(seqtab.nochim, "rdp_train_set_16.fa.gz", multithread=22)

# Save intermediate taxonomy file
saveRDS(taxa, file = "16S_RDP_Taxonomy_from_dada2.RDS")

print("taxonomy done")

# add_species names
taxa <- addSpecies(taxa, "rdp_species_assignment_16.fa.gz")

# Save completed taxonomy file
saveRDS(taxa, file = "16S_RDP_Taxonomy_from_dada2_sp.RDS")

print("saved files")

##############################

#on local computer
#use clean seq tab no chim without sample id column
seqtab.nochim.clean <- readRDS("./Output/16S_seqtab.nochim.clean.RDS")

# Use RDP training set for 16S
taxa <- assignTaxonomy(seqtab.nochim.clean, "C:/Users/becca/Dropbox/data_working/16S_soils/rdp_train_set_16.fa.gz", multithread=22)

# Save intermediate taxonomy file
saveRDS(taxa, file = "./Output/16S_RDP_Taxonomy_from_dada2.RDS")

# add_species names
taxa <- addSpecies(taxa, "./Output/rdp_species_assignment_16.fa.gz")

# Save completed taxonomy file
saveRDS(taxa, file = "./Output/16S_RDP_Taxonomy_from_dada2_sp.RDS")

#############################################################

#read in data from HPC

taxa <- readRDS("./Output/16S_RDP_Taxonomy_from_dada2_sp.RDS")

# inspect taxonomy
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


# Hand off to Phyloseq ####

#need clean seqtab with row names same as meta data "13S" "11S" etc. 

seqtab_a <- readRDS("./Output/16S_seqtab.nochim.clean.RDS")
seqtab.nochim_rownames <- seqtab_a
#okay this has the right row names


#get meta to have the same row names
rownames(meta) <- meta$Sample_ID


otu <- otu_table(seqtab.nochim,taxa_are_rows = FALSE)
tax <- tax_table(taxa)
met <- sample_data(meta)
sample_names(met) <- meta$Sample_ID

ps <- phyloseq(otu,met,tax, treeNJ)
saveRDS(ps,"./Output/16S_ps_not-cleaned.RDS")
ps


