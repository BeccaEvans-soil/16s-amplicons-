# -----------------------------------------------------------------------------#
# 16s SOILS
#
#Author: Geoffrey Zahn Rebecca Evans
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     vegan v 2.6.4
#                     phyloseq v 1.42.0
#                     broom v 1.0.3
# -----------------------------------------------------------------------------#

# SETUP ####

# packages 
library(tidyverse); packageVersion("tidyverse")
library(vegan); packageVersion("vegan")
library(phyloseq); packageVersion("phyloseq")
library(broom); packageVersion("broom")
library(lmerTest); packageVersion("lmerTest")
library(corncob); packageVersion("corncob")
library(microbiome)
library(emmeans)
library(microViz)
library(BAT)
library(BiodiversityR)

#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"), update = FALSE)




dodge <- position_dodge(.3)

setwd("C:/Users/becca/Dropbox/data_working/16S_soils")

# data
ps <- readRDS("./Output/16S_clean_phyloseq_object_tree.RDS")




####ALL SITES  (skip if analyzing primary and secondary seperatly)

#select just 0 and 16 to compare the composition evenly 
ps <- microViz:: ps_filter(ps, n.treat !="8kgN")


######################
#fix order of sites and ntreat

todf <- function(ps) {
  sd <- sample_data(ps)
  return(as(sd,"data.frame"))
}

ps_df<-todf(ps)

ps_df$n.treat<- factor(ps_df$n.treat, levels = c("0kgN","8kgN","16kgN"))
ps_df$site<- factor(ps_df$site, levels = c('East', "West", "Norway", "Harmony"))

str(ps_df)

sample_data(ps) <- as.data.frame(ps_df)

##################

#select pri and sec 
pri_ps <- microViz:: ps_filter(ps, succession=="Primary")
sec_ps <- microViz:: ps_filter(ps, succession=="Secondary")  

meta_pri <- meta(pri_ps)
meta_sec <- meta(sec_ps)


#alpha
ntaxa(pri_ps)
a <- otu_table(pri_ps)

plot_richness(pri_ps)

# ALPHA-DIV ESTIMATES ####
alpha <- estimate_richness(pri_ps) %>% 
  dplyr::select(Observed,Shannon, Simpson)
#calculate evenness
H <- alpha$Shannon
S1 <- alpha$Observed
S <- log(S1)
evenness <- H/S
evenness



#add alpha scores to ps_sample data
as.data.frame(alpha)
alpha$evenness <- evenness
str(alpha)
head(alpha)

#Sample_ID <- sample_data(ps)$Sample_ID
#alpha$Sample_ID <- Sample_ID


sam.new <- data.frame(New_var = sample(alpha))


# Turn into `sample_data` 
sam.new <- sample_data(sam.new)
head(sam.new)

# Merge with original phyloseq object
sec_ps <- merge_phyloseq(sec_ps, sam.new)

write.csv(sec_ps@sam_data, "sec_alpha_16s.csv")



meta_numeric <- read.csv("sec_alpha_16s.csv")


###########double check in vegan
####use vegan to calculate alpha

# convert the sample_data() within a phyloseq object to a vegan compatible data object
pssd2veg <- function(pri_ps) {
  sd <- sample_data(pri_ps)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a phyloseq object to a vegan compatible data object
psotu2veg <- function(pri_ps) {
  OTU <- otu_table(pri_ps)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

###get OTU table 
OTU1 = as(otu_table(pri_ps), "matrix")

str(OTU1)

data_richness <- estimateR(OTU1) # calculate richness and Chao1 using vegan package

t(data_richness)

data_shannon <- diversity(OTU1, index = "shannon")   # calculate Shannon index using vegan package


data_evenness <-  data_shannon/ log(specnumber(OTU1))        # calculate evenness index using vegan package


data_alphadiv <- cbind(data_shannon, data_evenness) # combine all indices in one data table


##########################


##########################
#########################pri 

###MODELS###

#########LMER PCR

#select the parameters of interest
subset_pca_numeric_env <- meta_numeric %>% dplyr::select(X21_d13c, X21_mg.c.gsoil, X21_d15n, X21_mg.n.gsoil, percent.moist2019,pH2019)

#prcomp is in base stats package
#https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/prcomp
pca_model <- prcomp(subset_pca_numeric_env , scale=TRUE)



par(mfrow=c(1,1))
biplot(pca_model)
summary(pca_model)
vectors <- pca_model$rotation #main PCA
#write.csv(vectors, "vectors.csv")


#summary
pca_model$x

# making variables from each of the PCs.
meta_numeric$e1 <- pca_model$x[,1]  
meta_numeric$e2 <- pca_model$x[,2] 
meta_numeric$e3 <- pca_model$x[,3] 
meta_numeric$e4 <- pca_model$x[,4] 
meta_numeric$e5 <- pca_model$x[,5] 
meta_numeric$e6 <- pca_model$x[,6] 


#Model 
#primary observed full 
names(meta_numeric)
str(meta_numeric)
mod <- lmerTest::lmer(New_var.Observed~weevil*n.treat+
                        e1+e2+e3+e4+e5+e6+(1|block_), data=meta_numeric)


summary(mod)
anova(mod)

step(mod)

mod <- lmerTest::lmer(New_var.Observed~e3+(1|block_), data=meta_numeric)
summary(mod)
anova(mod)

#post hoc
emmeans(mod, specs = pairwise ~ weevil:succession, type = "response")


###SHANNON


names(meta_numeric)
mod <- lmerTest::lmer(New_var.Shannon~weevil+n.treat+e1+e2+e3+e4+e5+e6+(1|block_), data=meta_numeric)
summary(mod)
anova(mod)

step(mod)





###evenness 
mod <- lmerTest::lmer(New_var.evenness~weevil+n.treat+e1+e2+e3+e4+e5+e6+(1|block_), data=meta_numeric)
summary(mod)
anova(mod)

step(mod)

anova(mod)







######################MODELS PRI#####################

#read in ps again to have the 8 in pri
ps <- readRDS("./Output/16S_clean_phyloseq_object_tree.RDS")

#alpha

######################
#fix order of sites and ntreat

todf <- function(ps) {
  sd <- sample_data(ps)
  return(as(sd,"data.frame"))
}

ps_df<-todf(ps)

ps_df$n.treat<- factor(ps_df$n.treat, levels = c("0kgN","8kgN","16kgN"))
ps_df$site<- factor(ps_df$site, levels = c('East', "West", "Norway", "Harmony"))

str(ps_df)

sample_data(ps) <- as.data.frame(ps_df)

##################

#alpha
ntaxa(ps)
otu_table(ps)[1:5, 1:5]

# ALPHA-DIV ESTIMATES ####
alpha <- estimate_richness(ps) %>% 
  dplyr::select(Observed,Shannon, Simpson)
#calculate evenness
H <- alpha$Shannon
S1 <- alpha$Observed
S <- log(S1)
evenness <- H/S
evenness



#add alpha scores to ps_sample data
as.data.frame(alpha)
alpha$evenness <- evenness
str(alpha)
head(alpha)

#Sample_ID <- sample_data(ps)$Sample_ID
#alpha$Sample_ID <- Sample_ID


sam.new <- data.frame(New_var = sample(alpha))


# Turn into `sample_data` 
sam.new <- sample_data(sam.new)
head(sam.new)

# Merge with original phyloseq object
ps <- merge_phyloseq(ps, sam.new)

#check
head(sample_data(ps))

str(sample_data(ps))

######################MODELS PRI #####################
################
#####################MODELS
meta_numeric <- read.csv("meta_16s_numeric.csv")
meta_numeric <- meta_numeric %>% filter(succession== "Primary")


#select the parameters of interest
subset_pca_numeric_env <- meta_numeric %>% dplyr::select(X21_d13c, X21_mg.c.gsoil, X21_d15n, X21_mg.n.gsoil, percent.moist2019,pH2019)

#prcomp is in base stats package
#https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/prcomp
pca_model <- prcomp(subset_pca_numeric_env , scale=TRUE)



par(mfrow=c(1,1))
biplot(pca_model)
summary(pca_model)
vectors <- pca_model$rotation #main PCA
#write.csv(vectors, "vectors.csv")


#summary
pca_model$x

# making variables from each of the PCs.
meta_numeric$e1 <- pca_model$x[,1]  
meta_numeric$e2 <- pca_model$x[,2] 
meta_numeric$e3 <- pca_model$x[,3] 
meta_numeric$e4 <- pca_model$x[,4] 
meta_numeric$e5 <- pca_model$x[,5] 
meta_numeric$e6 <- pca_model$x[,6] 



#Model 
#primary observed full 
names(meta_numeric)
mod <- lmerTest::lmer(New_var.Observed~weevil*n.treat+
                        e1+e2+e3+e4+e5+e6+(1|block), data=meta_numeric)

summary(mod)
anova(mod)

step(mod)

plot(meta_numeric$New_var.Observed~meta_numeric$numberwillow)

##NONE



###SHANNON


names(meta_numeric)
mod <- lmerTest::lmer(New_var.Shannon~weevil*n.treat+e1+e2+e3+e4+e5+e6+(1|block_), data=meta_numeric)
summary(mod)
anova(mod)

step(mod)

mod <- lmerTest::lmer(New_var.Shannon~weevil*n.treat+(1|block_), data=meta_numeric)
summary(mod)
anova(mod)


###evenness 
mod <- lmerTest::lmer(New_var.evenness~weevil*n.treat+e1+e2+e3+e4+e5+e6+(1|block_), data=meta_numeric)
summary(mod)
anova(mod)

step(mod)
mod <- lmerTest::lmer(New_var.Shannon~weevil+e3+n.treat+weevil:n.treat+(1|block_), data=meta_numeric)
summary(mod)
anova(mod)

###############################
####Secondary models

ps <- readRDS("./Output/16S_clean_phyloseq_object_tree.RDS")

#alpha

######################
#fix order of sites and ntreat

todf <- function(ps) {
  sd <- sample_data(ps)
  return(as(sd,"data.frame"))
}

ps_df<-todf(ps)

ps_df$n.treat<- factor(ps_df$n.treat, levels = c("0kgN","8kgN","16kgN"))
ps_df$site<- factor(ps_df$site, levels = c('East', "West", "Norway", "Harmony"))

str(ps_df)

sample_data(ps) <- as.data.frame(ps_df)

meta_numeric <- read.csv("meta_16s_numeric.csv")
meta_numeric <- meta_numeric %>% filter(succession== "Secondary")


sec_ps <- microViz:: ps_filter(ps, succession=="Secondary")  



##################

#alpha
ntaxa(sec_ps)
otu_table(sec_ps)[1:5, 1:5]

# ALPHA-DIV ESTIMATES ####
alpha <- estimate_richness(sec_ps) %>% 
  dplyr::select(Observed,Shannon, Simpson)
#calculate evenness
H <- alpha$Shannon
S1 <- alpha$Observed
S <- log(S1)
evenness <- H/S
evenness



#add alpha scores to ps_sample data
as.data.frame(alpha)
alpha$evenness <- evenness
str(alpha)
head(alpha)

#Sample_ID <- sample_data(ps)$Sample_ID
#alpha$Sample_ID <- Sample_ID


sam.new <- data.frame(New_var = sample(alpha))


# Turn into `sample_data` 
sam.new <- sample_data(sam.new)
head(sam.new)

# Merge with original phyloseq object
sec_ps <- merge_phyloseq(sec_ps, sam.new)


######################MODELS sec #####################
################
#####################MODELS


#select the parameters of interest
subset_pca_numeric_env <- meta_numeric %>% dplyr::select(X21_d13c, X21_mg.c.gsoil, X21_d15n, X21_mg.n.gsoil, percent.moist2019,pH2019)

#prcomp is in base stats package
#https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/prcomp
pca_model <- prcomp(subset_pca_numeric_env , scale=TRUE)



par(mfrow=c(1,1))
biplot(pca_model)
summary(pca_model)
vectors <- pca_model$rotation #main PCA
#write.csv(vectors, "vectors.csv")


#summary
pca_model$x

# making variables from each of the PCs.
meta_numeric$e1 <- pca_model$x[,1]  
meta_numeric$e2 <- pca_model$x[,2] 
meta_numeric$e3 <- pca_model$x[,3] 
meta_numeric$e4 <- pca_model$x[,4] 
meta_numeric$e5 <- pca_model$x[,5] 
meta_numeric$e6 <- pca_model$x[,6] 



#Model 
#primary observed full 
names(meta_numeric)
mod <- lmerTest::lmer(New_var.Observed~weevil*n.treat+
                        e1+e2+e3+e4+e5+e6+(1|block), data=meta_numeric)

summary(mod)
anova(mod)

step(mod)

mod <- lmerTest::lmer(New_var.Observed~weevil+(1|block), data=meta_numeric)
summary(mod)
anova(mod)




###SHANNON


names(meta_numeric)
mod <- lmerTest::lmer(New_var.Shannon~weevil*n.treat+e1+e2+e3+e4+e5+e6+(1|block_), data=meta_numeric)
summary(mod)
anova(mod)

step(mod)

mod <- lmerTest::lmer(New_var.Shannon~weevil+(1|block_), data=meta_numeric)
summary(mod)
anova(mod)


###evenness 
mod <- lmerTest::lmer(New_var.evenness~weevil*n.treat+e1+e2+e3+e4+e5+e6+(1|block_), data=meta_numeric)
summary(mod)
anova(mod)

step(mod)
mod <- lmerTest::lmer(New_var.evenness~weevil+(1|block_), data=meta_numeric)
summary(mod)
anova(mod)



###FIGS
####################
#compare primary vs secondary figures 


names(meta)
std.error <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
all.summary <- meta_numeric %>% group_by(n.treat, weevil,succession)%>%
  summarise_at(c("New_var.Observed",       "New_var.Simpson"      ,  "New_var.Shannon"    ,    "New_var.evenness"), 
               funs(mean, std.error))

names(all.summary)

#plot pri vs sec, ntreat, weevil  observed
obs_fig_weevil_n_<- ggplot(data=all.summary, aes(x=n.treat,y=New_var.Observed_mean, group=weevil))+
  geom_point(aes(color=weevil),position=dodge, size=4)+
  facet_grid(~succession)+
  geom_errorbar(aes(ymin=New_var.Observed_mean- New_var.Observed_std.error,
                    ymax=New_var.Observed_mean+ New_var.Observed_std.error, 
                    width=0.4), position=dodge)+
  scale_color_manual(values=c("grey", "black"), name="", breaks=c("NW","W"), labels=c("No Weevil", "Weevil"))+
  xlab(expression(Nitrogen~addition~(kg~N~ha^{-1}~yr^{-1})))+
  ylab("Bacteria Richness")+
  theme(axis.title.y = element_text(size=13))+
  theme(axis.text.y=element_text(size=11))+
  theme(axis.title.x=element_text(size=13))+
  theme(axis.text.x=element_text(size=11))+
  theme(strip.text.x = element_text(size = 13))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave("16s_soil_observed_by_weevil_n_treat_succession.png",height = 3.5,width = 5)


#plot pri vs sec, ntreat, weevil  shannon
shannon_fig_weevil_n_<- ggplot(data=all.summary, aes(x=n.treat,y=New_var.Shannon_mean, group=weevil))+
  geom_point(aes(color=weevil),position=dodge, size=4)+
  facet_grid(~succession)+
  geom_errorbar(aes(ymin=New_var.Shannon_mean- New_var.Shannon_std.error,
                    ymax=New_var.Shannon_mean+ New_var.Shannon_std.error, 
                    width=0.4), position=dodge)+
  scale_color_manual(values=c("grey", "black"), name="", breaks=c("NW","W"), labels=c("No Weevil", "Weevil"))+
  xlab(expression(Nitrogen~addition~(kg~N~ha^{-1}~yr^{-1})))+
  ylab("Bacteria Shannon Diversity (H')")+
  theme(axis.title.y = element_text(size=13))+
  theme(axis.text.y=element_text(size=11))+
  theme(axis.title.x=element_text(size=13))+
  theme(axis.text.x=element_text(size=11))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(strip.text.x = element_blank())
  

ggsave("16s_soil_shannon_by_weevil_n_treat_succession.png",height = 3.5,width = 5)


#plot pri vs sec, ntreat, weevil  eveness
even_fig_weevil_n_<- ggplot(data=all.summary, aes(x=n.treat,y=New_var.evenness_mean, group=weevil))+
  geom_point(aes(color=weevil),position=dodge, size=4)+
  facet_grid(~succession)+
  geom_errorbar(aes(ymin=New_var.evenness_mean- New_var.evenness_std.error,
                    ymax=New_var.evenness_mean+ New_var.evenness_std.error, 
                    width=0.4), position=dodge)+
  scale_color_manual(values=c("grey", "black"), name="", breaks=c("NW","W"), labels=c("No Weevil", "Weevil"))+
  xlab(expression(Nitrogen~addition~(kg~N~ha^{-1}~yr^{-1})))+
  ylab("Bacteria Evenness")+
  theme(axis.title.y = element_text(size=13))+
  theme(axis.text.y=element_text(size=11))+
  theme(axis.title.x=element_text(size=13))+
  theme(axis.text.x=element_text(size=11))+
  theme(strip.text.x = element_blank())


ggsave("16s_soil_even_by_weevil_n_treat_succession.png",height = 3.5,width = 5)

############################
#################################
###move to vegan
veganotu = function(ps) {
  require("vegan")
  OTU = otu_table(ps)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}

bac_vegan <- veganotu(ps)
write.csv(bac_vegan,"otu.csv")
otu <- read.csv("otu.csv")

data(dune)
str(dune)
dune$area = sample(1:100, 20) 
mod.w = specaccum(dune, "random", w = area)

mod.orig = specaccum(dune, "random")

mod.w = specaccum(otu, "random", w = otu$numberwillow)

sac <- specaccum(bac_vegan)
sac[[3]]

data(BCI)

dune$area = sample(1:100, 20) 
pool <- specpool(bac_vegan)

acc <- specaccum(dune, method = "rarefaction",permutations=100)

data <- data.frame(Sites=acc$sites, Richness=acc$richness, SD=acc$sd)

#####
library(BAT)
hill(bac_vegan, raref=1,q=0, runs=100)
hill(bac_vegan, q=1, runs=100)


