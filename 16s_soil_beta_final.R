###BETA 16s soil 

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
library(xtable)
library(AICcPermanova)
library(microViz)
library(pairwiseAdonis)

#BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"), update = FALSE)

dodge <- position_dodge(.3)

setwd("C:/Users/becca/Dropbox/data_working/16S_soils")

# data
ps <- readRDS("./Output/16S_clean_phyloseq_object_tree.RDS")


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


####ALL SITES 

#select just 0 and 16 to compare the composition evenly 
ps <- microViz:: ps_filter(ps, n.treat !="8kgN")




# transform raw counts to relative abundance ####
all_ps_ra <- transform_sample_counts(ps, fun = function(x){x/sum(x)})




# Beta-diversity ####


##################
###PRI and SEC
#################

# pull out components for stats
asv_all <- otu_table(all_ps_ra) %>% as("matrix") %>% as.data.frame()
meta_all <- read.csv("meta_16s_numeric.csv")
meta_all <- meta_all %>% filter(n.treat != "8kgN")


dist_all <- phyloseq::distance(all_ps_ra, method="bray")
ord_all <- ordinate(all_ps_ra, method="NMDS", distance= dist_all)

##get pca for env var

subset_pca_numeric <- meta_all %>% select(X21_d13c, X21_mg.c.gsoil, X21_d15n, X21_mg.n.gsoil, percent.moist2019,pH2019)

#prcomp is in base stats package
#https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/prcomp
pca_model <- prcomp(subset_pca_numeric , scale=TRUE)



par(mfrow=c(1,1))
biplot(pca_model)
summary(pca_model)
vectors <- pca_model$rotation #main PCA


#summary
pca_model$x
#to see that x is an array of all the data values for PCs 1-8. PC1 is pca_pri_model$x[,1] 


# making variables from each of the PCs.
meta_all$e1 <- pca_model$x[,1]  #do this for each of the 8 PCs; it might be easier for you to give the dataframe a shorter name
meta_all$e2 <- pca_model$x[,2] 
meta_all$e3 <- pca_model$x[,3] 
meta_all$e4 <- pca_model$x[,4] 
meta_all$e5 <- pca_model$x[,5] 
meta_all$e6 <- pca_model$x[,6] 




#run permanova with all interactions
names(meta_all)

permanova_mm_allinteract <- vegan::adonis2(as.dist(dist_all)~meta_all$n.treat*meta_all$weevil*meta_all$succession+ meta_all$e1+
                                             meta_all$e2+meta_all$e3+meta_all$e4+meta_all$e5+meta_all$e6,
                                           strata=meta_all$block_)
permanova_mm_allinteract
AICc_permanova2(permanova_mm_allinteract)


permanova_mm_allinteract_simple <- vegan::adonis2(as.dist(dist_all)~meta_all$n.treat+meta_all$weevil+meta_all$succession+ 
                                                    meta_all$weevil:meta_all$succession,
                                           strata=meta_all$block_)
permanova_mm_allinteract_simple
AICc_permanova2(permanova_mm_allinteract_simple)


write.csv(as.data.frame(permanova_mm_allinteract_simple), "soil16s_perm.csv")

##NMDS

nmds <- metaMDS(asv_all, distance="bray")
data.scores <- as.data.frame(scores(nmds)$sites)
data.scores$weevil <-meta_all$weevil 
data.scores$n.treat <-meta_all$n.treat
data.scores$block_ <-meta_all$block_
data.scores$succession<-meta_all$succession


str(data.scores)
data.scores$weevil <- as.factor(data.scores$weevil)
data.scores$n.treat <- as.factor(data.scores$n.treat)
data.scores$succession <- as.factor(data.scores$succession)
str(data.scores)


todf <- function(ps) {
  sd <- sample_data(ps)
  return(as(sd,"data.frame"))
}

ps_df<-todf(ps)
data.scores$n.treat<- factor(data.scores$n.treat, levels = c("0kgN","8kgN","16kgN"))
ps_df$site<- factor(ps_df$site, levels = c('East', "West", "Norway", "Harmony"))

str(ps_df)

sample_data(ps) <- as.data.frame(ps_df)

#n weevil and weevil and succession
n_weevil_succession_nmds = ggplot(data.scores, aes(x = NMDS1, y = NMDS2, )) + 
  geom_point(size=4, aes(color=weevil))  +
  stat_ellipse(geom = "polygon",
               aes(fill = succession), 
               alpha = 0.25)+
  scale_color_manual(values=c("grey", "black"), name="", breaks=c("NW","W"), labels=c("No Weevil", "Weevil"))+
  #ylim(c(-1.5,1.5))+
  #xlim(c(-1.5,1.5))+
  theme(axis.title.y = element_text(size=13))+
  theme(axis.text.y=element_text(size=11))+
  theme(axis.title.x=element_text(size=13))+
  theme(axis.text.x=element_text(size=11))+
  theme(strip.text.x = element_text(size = 13))

####DCA

a <- decorana(dist_all, iweigh = 1)
summary(a)

dca <- scores(a, choices = c(1,2))

data.scores <- as.data.frame(dca)
data.scores$weevil <-meta_all$weevil 
data.scores$n.treat <-meta_all$n.treat
data.scores$block_ <-meta_all$block_
data.scores$succession<-meta_all$succession


str(data.scores)
data.scores$weevil <- as.factor(data.scores$weevil)
data.scores$n.treat <- as.factor(data.scores$n.treat)
data.scores$succession <- as.factor(data.scores$succession)
str(data.scores)

data.scores$n.treat<- factor(data.scores$n.treat, levels = c("0kgN","8kgN","16kgN"))



#plot
n_weevil_succession_nmds = ggplot(data.scores, aes(x = DCA1, y = DCA2, )) + 
  geom_point(size=4, aes(shape = n.treat, color=weevil))  +
  stat_ellipse(geom = "polygon",
               aes(fill = succession), 
               alpha = 0.25)+
  scale_color_manual(values=c("grey", "black"), name="", breaks=c("NW","W"), labels=c("No Weevil", "Weevil"))+
  #ylim(c(-1.5,1.5))+
  #xlim(c(-1.5,1.5))+
  theme(axis.title.y = element_text(size=13))+
  theme(axis.text.y=element_text(size=11))+
  theme(axis.title.x=element_text(size=13))+
  theme(axis.text.x=element_text(size=11))+
  theme(strip.text.x = element_text(size = 13))

###SEPERATE PRI AND SEC 
#seperate
pri_ps <- microViz:: ps_filter(ps, succession=="Primary")
sec_ps <- microViz:: ps_filter(ps, succession=="Secondary") 

meta_pri <- meta(pri_ps)
meta_sec <- meta(sec_ps)

#pri
pri_ps_ra <- transform_sample_counts(pri_ps, fun = function(x){x/sum(x)})


dist_pri <- phyloseq::distance(pri_ps_ra, method="bray")
ord_pri <- ordinate(pri_ps_ra, method="NMDS", distance= dist_pri)


# pull out components for stats
asv_pri <- otu_table(pri_ps_ra) %>% as("matrix") %>% as.data.frame()
meta_pri <- sample_data(pri_ps) %>% as.data.frame()

write.csv(meta_pri,"metapri.csv")
meta_pri <- read.csv("metapri.csv")


#sec
sec_ps_ra <- transform_sample_counts(sec_ps, fun = function(x){x/sum(x)})


dist_sec <- phyloseq::distance(sec_ps_ra, method="bray")
ord_sec <- ordinate(sec_ps_ra, method="NMDS", distance= dist_sec)


# pull out components for stats
asv_sec <- otu_table(sec_ps_ra) %>% as("matrix") %>% as.data.frame()
meta_sec <- sample_data(sec_ps) %>% as.data.frame()

write.csv(meta_sec,"metasec.csv")
meta_sec <- read.csv("metasec.csv")


###JUST PRI



##get pca for env var

subset_pca_numeric <- meta_pri %>% select(X21_d13c, X21_mg.c.gsoil, X21_d15n, X21_mg.n.gsoil, percent.moist2019,pH2019)

#prcomp is in base stats package
#https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/prcomp
pca_model <- prcomp(subset_pca_numeric , scale=TRUE)



par(mfrow=c(1,1))
biplot(pca_model)
summary(pca_model)
vectors <- pca_model$rotation #main PCA


#summary
pca_model$x
#to see that x is an array of all the data values for PCs 1-8. PC1 is pca_pri_model$x[,1] 


# making variables from each of the PCs.
meta_pri$e1 <- pca_model$x[,1]  #do this for each of the 8 PCs; it might be easier for you to give the dataframe a shorter name
meta_pri$e2 <- pca_model$x[,2] 
meta_pri$e3 <- pca_model$x[,3] 
meta_pri$e4 <- pca_model$x[,4] 
meta_pri$e5 <- pca_model$x[,5] 
meta_pri$e6 <- pca_model$x[,6] 




#run permanova with all interactions


permanova_mm_priinteract <- vegan::adonis2(as.dist(dist_pri)~meta_pri$n.treat*meta_pri$weevil+ meta_pri$e1+
                                             meta_pri$e2+meta_pri$e3+meta_pri$e4+meta_pri$e5+meta_pri$e6,
                                           strata=meta_pri$block_)
permanova_mm_priinteract
AICc_permanova2(permanova_mm_priinteract)




#simple
permanova_mm_priinteract_simple <- vegan::adonis2(as.dist(dist_pri)~meta_pri$n.treat+meta_pri$weevil+
                                             meta_pri$e2+meta_pri$e3,
                                           strata=meta_pri$block_)


permanova_mm_priinteract_simple
AICc_permanova2(permanova_mm_priinteract_simple)
write.csv(as.data.frame(permanova_mm_priinteract_simple), "soil16s_perm.csv")



#pc4= moist and d13c
moist_model <- lmerTest::lmer(percent.moist2019~n.treat+weevil+site+(1|block_), data=meta_pri)
summary(moist_model)

d15n_model <- lmerTest::lmer(X21_d15n~n.treat+weevil+site+(1|block_), data=meta_pri)
summary(d15n_model)


##means


std.error <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

moist_summary <- meta_pri %>% dplyr::group_by(site)%>%
  summarise_at("percent.moist2019", funs(mean, std.error))


moist <- ggplot(moist_summary, aes(x= site, y = mean)) + 
  geom_point(size=4)  +
  geom_errorbar(aes(ymin=mean- std.error,
                    ymax=mean+ std.error, 
                    width=0.4), position=dodge)+
  #ylim(c(-1.5,1.5))+
  #xlim(c(-1.5,1.5))+
  theme(axis.title.y = element_text(size=13))+
  theme(axis.text.y=element_text(size=11))+
  theme(axis.title.x=element_text(size=13))+
  theme(axis.text.x=element_text(size=11))+
  ylab("Mean soil moisture (%)")+
  theme(strip.text.x = element_text(size = 13))


##NMDS

nmds <- metaMDS(asv_pri, distance="bray")
data.scores <- as.data.frame(scores(nmds)$sites)
data.scores$weevil <-meta_pri$weevil 
data.scores$n.treat <-meta_pri$n.treat
data.scores$block_ <-meta_pri$block_
data.scores$site <-meta_pri$site


str(data.scores)
data.scores$weevil <- as.factor(data.scores$weevil)
data.scores$n.treat <- as.factor(data.scores$n.treat)
data.scores$site <- as.factor(data.scores$site)
str(data.scores)


todf <- function(pri_ps) {
  sd <- sample_data(pri_ps)
  return(as(sd,"data.frame"))
}

pri_ps_df<-todf(pri_ps)
data.scores$n.treat<- factor(data.scores$n.treat, levels = c("0kgN","8kgN","16kgN"))
pri_ps_df$site<- factor(pri_ps_df$site, levels = c('East', "West"))


sample_data(pri_ps) <- as.data.frame(pri_ps_df)

#weevil in pri
weevil_succession_nmds = ggplot(data.scores, aes(x = NMDS1, y = NMDS2, )) + 
  geom_point(size=4, aes(shape = weevil, color=site))  +
  stat_ellipse(geom = "polygon",
               aes(fill = site), 
               alpha = 0.25)+
  scale_color_manual(values=c("grey", "black"), name="", breaks=c("NW","W"), labels=c("No Weevil", "Weevil"))+
  #ylim(c(-1.5,1.5))+
  #xlim(c(-1.5,1.5))+
  theme(axis.title.y = element_text(size=13))+
  theme(axis.text.y=element_text(size=11))+
  theme(axis.title.x=element_text(size=13))+
  theme(axis.text.x=element_text(size=11))+
  theme(strip.text.x = element_text(size = 13))



##DCA

pri <- decorana(asv_pri, iweigh = 1)
summary(pri)

dca <- scores(pri, choices = c(1,2))

data.scores_pri <- as.data.frame(dca)
data.scores_pri$weevil <-meta_all$weevil 
data.scores_pri$n.treat <-meta_all$n.treat
data.scores_pri$block_ <-meta_all$block_




data.scores_pri$weevil <- as.factor(data.scores$weevil)
data.scores_pri$n.treat <- as.factor(data.scores$n.treat)
str(data.scores_pri)

data.scores_pri$n.treat<- factor(data.scores_pri$n.treat, levels = c("0kgN","8kgN","16kgN"))

##PLOT

n_weevil_pri_= ggplot(data.scores_pri, aes(x = DCA1, y = DCA2 )) + 
  geom_point(size=4, aes(color=weevil, shape=n.treat))  +
  stat_ellipse(geom = "polygon",
               aes(fill = n.treat), 
               alpha = 0.25)+
  scale_color_manual(values=c("grey", "black"), name="", breaks=c("NW","W"), labels=c("No Weevil", "Weevil"))+
  #ylim(c(-1.5,1.5))+
  #xlim(c(-1.5,1.5))+
  theme(axis.title.y = element_text(size=13))+
  theme(axis.text.y=element_text(size=11))+
  theme(axis.title.x=element_text(size=13))+
  theme(axis.text.x=element_text(size=11))+
  theme(strip.text.x = element_text(size = 13))


n_weevil_pri_= ggplot(data.scores_pri, aes(x = DCA1, y = DCA2 )) + 
  geom_point(size=4, aes(color=weevil, shape=n.treat))  +
  stat_ellipse(geom = "polygon",
               aes(fill = weevil), 
               alpha = 0.25)+
  scale_color_manual(values=c("grey", "black"), name="", breaks=c("NW","W"), labels=c("No Weevil", "Weevil"))+
  #ylim(c(-1.5,1.5))+
  #xlim(c(-1.5,1.5))+
  theme(axis.title.y = element_text(size=13))+
  theme(axis.text.y=element_text(size=11))+
  theme(axis.title.x=element_text(size=13))+
  theme(axis.text.x=element_text(size=11))+
  theme(strip.text.x = element_text(size = 13))



################ ###SECONDARY 
##get pca for env var

subset_pca_numeric <- meta_sec %>% select(X21_d13c, X21_mg.c.gsoil, X21_d15n, X21_mg.n.gsoil, percent.moist2019,pH2019)
#prcomp is in base stats package
#https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/prcomp
pca_model<- prcomp(subset_pca_numeric , scale=TRUE)



par(mfrow=c(1,1))
biplot(pca_model)
summary(pca_model)
vectors <- pca_model$rotation #main PCA


#summary
pca_model$x
#to see that x is an array of all the data values for PCs 1-8. PC1 is pca_pri_model$x[,1] 


# making variables from each of the PCs.
meta_sec$e1 <- pca_model$x[,1]  #do this for each of the 8 PCs; it might be easier for you to give the dataframe a shorter name
meta_sec$e2 <- pca_model$x[,2] 
meta_sec$e3 <- pca_model$x[,3] 
meta_sec$e4 <- pca_model$x[,4] 
meta_sec$e5 <- pca_model$x[,5] 
meta_sec$e6 <- pca_model$x[,6] 

sec_ps_ra <- transform_sample_counts(sec_ps, fun = function(x){x/sum(x)})
# pull out components for stats
asv_sec <- otu_table(sec_ps_ra) %>% as("matrix") %>% as.data.frame()
#meta_sec <- sample_data(sec_ps) %>% as.data.frame()

#PERMANOVA 

#run permanova with all interactions


permanova_mm_priinteract <- vegan::adonis2(as.dist(dist_sec)~meta_sec$n.treat*meta_sec$weevil+ meta_sec$e1+
                                             meta_sec$e2+meta_sec$e3+meta_sec$e4+meta_sec$e5+meta_sec$e6,
                                           strata=meta_sec$block_)
permanova_mm_priinteract
AICc_permanova2(permanova_mm_priinteract)




#simple
permanova_mm_priinteract_simple <- vegan::adonis2(as.dist(dist_sec)~meta_sec$n.treat+meta_sec$weevil+ meta_sec$e4,
                                                  strata=meta_sec$block_)


permanova_mm_priinteract_simple
AICc_permanova2(permanova_mm_priinteract_simple)
write.csv(as.data.frame(permanova_mm_priinteract_simple), "soil16s_perm.csv")

#pc3=ph

ph_model <- lmerTest::lmer(pH2019~n.treat+weevil+site+(1|block_), data=meta_pri)
summary(ph_model)
str(meta_pri)

#pc4= moist and d15n
moist_model <- lmerTest::lmer(percent.moist2019~n.treat+weevil+site+(1|block_), data=meta_pri)
summary(moist_model)

#pc4= moist and d13c
moist_model <- lmerTest::lmer(percent.moist2019~n.treat+weevil+site+(1|block_), data=meta_sec)
summary(moist_model)

d15n_model <- lmerTest::lmer(X21_d15n~n.treat+weevil+site+(1|block_), data=meta_sec)
summary(d15n_model)



##NMDS

nmds <- metaMDS(asv_sec, distance="bray")
data.scores <- as.data.frame(scores(nmds)$sites)
data.scores$weevil <-meta_sec$weevil 
data.scores$n.treat <-meta_sec$n.treat
data.scores$block_ <-meta_sec$block_
data.scores$site <-meta_sec$site


str(data.scores)
data.scores$weevil <- as.factor(data.scores$weevil)
data.scores$n.treat <- as.factor(data.scores$n.treat)
data.scores$site <- as.factor(data.scores$site)
str(data.scores)


todf <- function(sec_ps) {
  sd <- sample_data(sec_ps)
  return(as(sd,"data.frame"))
}

sec_ps_df<-todf(sec_ps)
data.scores$n.treat<- factor(data.scores$n.treat, levels = c("0kgN","16kgN"))
sec_ps_df$site<- factor(sec_ps_df$site, levels = c('Norway', "Harmony"))


sample_data(sec_ps) <- as.data.frame(sec_ps_df)

#weevil in sec
n_treat = ggplot(data.scores, aes(x = NMDS1, y = NMDS2, )) + 
  geom_point(size=4, aes(shape = n.treat))  +
  stat_ellipse(geom = "polygon",
               aes(fill = n.treat), 
               alpha = 0.25)+
  #scale_color_manual(values=c("grey", "black"), name="", breaks=c("NW","W"), labels=c("No Weevil", "Weevil"))+
  #ylim(c(-1.5,1.5))+
  #xlim(c(-1.5,1.5))+
  theme(axis.title.y = element_text(size=13))+
  theme(axis.text.y=element_text(size=11))+
  theme(axis.title.x=element_text(size=13))+
  theme(axis.text.x=element_text(size=11))+
  theme(strip.text.x = element_text(size = 13))


##DCA

sec <- decorana(asv_sec, iweigh = 1)
summary(sec)

dca_sec <- scores(sec, choices = c(1,2))

data.scores_sec <- as.data.frame(dca_sec)
data.scores_sec $weevil <-meta_sec$weevil 
data.scores_sec $n.treat <-meta_sec$n.treat
data.scores_sec $block_ <-meta_sec$block_
data.scores_sec $succession<-meta_sec$succession



data.scores_sec$weevil <- as.factor(data.scores_sec$weevil)
data.scores_sec$n.treat <- as.factor(data.scores_sec$n.treat)
data.scores_sec$succession <- as.factor(data.scores_sec$succession)
str(data.scores_sec)

data.scores_sec$n.treat<- factor(data.scores_sec$n.treat, levels = c("0kgN","16kgN"))



n_weevil_sec= ggplot(data.scores_sec, aes(x = DCA1, y = DCA2, )) + 
  geom_point(size=4, aes(color=weevil, shape=n.treat))  +
  stat_ellipse(geom = "polygon",
               aes(fill = n.treat), 
               alpha = 0.25)+
  scale_color_manual(values=c("grey", "black"), name="", breaks=c("NW","W"), labels=c("No Weevil", "Weevil"))+
  #ylim(c(-1.5,1.5))+
  #xlim(c(-1.5,1.5))+
  theme(axis.title.y = element_text(size=13))+
  theme(axis.text.y=element_text(size=11))+
  theme(axis.title.x=element_text(size=13))+
  theme(axis.text.x=element_text(size=11))+
  theme(strip.text.x = element_text(size = 13))


################
########CORNCOB
####################
###################################
##################pri vs sec

# merge taxa at genus level
# data
ps <- readRDS("./Output/16S_clean_phyloseq_object_tree.RDS")


#read in data and get rid of 8
ps <- microViz:: ps_filter(ps, n.treat !="8kgN")

#cleann up
tax_table(ps)[, colnames(tax_table(ps))] <- gsub(tax_table(ps)[, colnames(tax_table(ps))],   
                                                 pattern = "[a-z]__", replacement = "")
phyloseq::tax_table(ps)[1:7,1:7]




# merge taxa at genus level
ps_genus <- ps %>% 
  tax_glom(taxrank = "Genus")

# Clean up ASV names to show taxonomy
ASV_names <- otu_table(ps) %>% colnames()
#ASV_taxa <- otu_to_taxonomy(ASV_names,ps,level = c("Phylum","Class","Order","Family","Genus"))
ASV_taxa <- otu_to_taxonomy(ASV_names,ps,level = c("Family","Genus"))

genus_names <- otu_table(ps_genus) %>% colnames()


#genus_taxa <- otu_to_taxonomy(genus_names,ps_genus,level = c("Phylum","Class","Order","Family","Genus"))

genus_taxa <- otu_to_taxonomy(genus_names,ps_genus,level = c("Family","Genus"))

# CORNCOB DIFFABUND #### succession
ps_genus@sam_data$succession %>% unique

# set levels so "primary" is intercept
ps_genus@sam_data$succession<- 
  ps_genus@sam_data$succession %>% factor(levels = c( "Primary"  , "Secondary"))


# use raw count data for corncob
da_analysis_succession <- differentialTest(formula = ~ succession, #abundance
                                           phi.formula = ~ 1, #dispersion
                                           formula_null = ~ 1, #mean
                                           phi.formula_null = ~ 1,
                                           test = "Wald", boot = FALSE,
                                           data = ps_genus,
                                           fdr_cutoff = 0.05,
                                           full_output = TRUE)

da_analysis_succession$significant_models[1:5]
da_analysis_succession$significant_taxa



#test how significant 
#We can examine a subset of the p-values of our tests using:


da_analysis_succession[1:5]


####################

#asv to taxa
otu_to_taxonomy(OTU = da_analysis_succession$significant_taxa, data = ps_genus)


p6 <- 
  plot(da_analysis_succession, level = c("Family", "Genus")) +
  theme(legend.position = 'none') 



sig_taxa <- da_analysis_succession$significant_taxa %>% otu_to_taxonomy(data=ps_genus)


pri_weevil_vs_noweevil_sigtaxa_16s_soil<- as.data.frame(sig_taxa)
write.csv(pri_weevil_vs_noweevil_sigtaxa_16s_soil, "pri_weevil_vs_noweevil_sigtaxa_16s_soils.csv")


#get rel abundance of these taxa
tax_table(ps) <- cbind(tax_table(ps), 
                           rownames(tax_table(ps)))

colnames(tax_table(ps)) <- 
  c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTUID")

tax_table(ps)
relabun.ps <- transform_sample_counts(ps,function(x) x / sum(x))
ps_genus <- tax_glom(relabun.ps, taxrank = "Genus", NArm = FALSE)


gen_all_taxa<- ps_genus@tax_table
write.csv(gen_all_taxa ,"gen_all_taxa.csv")



ps_genusP <- subset_taxa(ps_genus, Genus %in% c(
                                                'Streptomyces',
                                                'Bradyrhizobium',
                                                'Segetibacter',
                                                'Gemmatimonas',
                                                'Sphingomonas',
                                                'Actinomycetospora',
                                                'Angustibacter',
                                                'Chitinophaga',
                                                'Cupriavidus',
                                                'Methylobacterium',
                                                'Tumebacillus',
                                                'Aciditerrimonas',
                                                'Jatrophihabitans',
                                                'Rhodococcus',
                                                'Hyphomicrobium',
                                                'Parachlamydia',
                                                'Aquicella',
                                                'Parafilimonas',
                                                'Psychroglaciecola',
                                                'Pseudomonas',
                                                'Flavisolibacter',
                                                'Ferruginibacter',
                                                'Hymenobacter',
                                                'Cystobacter',
                                                'Kribbella',
                                                'Acidisoma',
                                                'Actinoplanes',
                                                'Conexibacter',
                                                'Nitrospira',
                                                'Flavobacterium',
                                                'Chryseobacterium',
                                                'Methylocella',
                                                'Spirosoma',
                                                'Nitrolancea',
                                                'Ramlibacter',
                                                'Rhizobacter',
                                                'Romboutsia',
                                                'Phycicoccus',
                                                'Paenibacillus',
                                                'Legionella',
                                                'Roseiarcus',
                                                'Clostridium',
                                                'Neochlamydia',
                                                'Pseudolabrys',
                                                'Thermoflexus',
                                                'Actinotalea',
                                                'Labilithrix',
                                                'Roseococcus',
                                                'Thermasporomyces',
                                                'Fibrella',
                                                'Indibacter',
                                                'Flexithrix',
                                                'Coxiella',
                                                'Sporocytophaga',
                                                'Bosea',
                                                'Elstera',
                                                'Dongia',
                                                'Clostridium',
                                                'Amaricoccus',
                                                'Kitasatospora',
                                                'Deinococcus',
                                                'Tepidisphaera',
                                                'Phaselicystis',
                                                'Stella',
                                                'Rhodomicrobium',
                                                'Sporichthya',
                                                'Oligoflexus',
                                                'Rubellimicrobium',
                                                'Psychrosinus',
                                                'Altererythrobacter',
                                                'Luteolibacter',
                                                'Geobacter',
                                                'Turicibacter',
                                                'Thermomarinilinea'))


genus.df <- psmelt(ps_genusP)
head(genus.df)


MySummary <- genus.df %>%
  group_by(Genus) %>%
  summarize(mean_abund = mean(Abundance, na.rm=TRUE)) 
head(MySummary)

write.csv(MySummary, "16s_prisec_relabund.csv")

#this file
bacrelabundprisec <- read.csv("16s_prisec_relabund.csv")






############PRIMARY 

# CORNCOB DIFFABUND #### weevil primary 
# merge taxa at genus level


pri_ps <- microViz:: ps_filter(ps, succession=="Primary")



pri_ps_genus <- pri_ps %>% 
  tax_glom(taxrank = "Genus")

pri_ps_genus@sam_data$weevil %>% unique

# set levels so "weevil" is intercept
pri_ps_genus@sam_data$weevil <- 
  pri_ps_genus@sam_data$weevil %>% factor(levels = c("W","NW"))


# use raw count data for corncob
da_analysis_weevil <- differentialTest(formula = ~ weevil, #abundance
                                       phi.formula = ~ 1, #dispersion
                                       formula_null = ~ 1, #mean
                                       phi.formula_null = ~ 1,
                                       test = "Wald", boot = FALSE,
                                       data = pri_ps_genus,
                                       fdr_cutoff = 0.05,
                                       full_output = TRUE)


p4 <- 
  plot(da_analysis_weevil, level = c("Family", "Genus")) +
  theme(legend.position = 'none') 

#saveRDS(p4,"./Output/pri_16S_diffabund_weevil.RDS")



sig_taxa <- da_analysis_weevil$significant_taxa %>% otu_to_taxonomy(data=pri_ps_genus)


#get taxonomy of sig taxa
pri_weevil_vs_noweevil_sig_taxa <- otu_to_taxonomy(OTU=da_analysis_weevil$significant_taxa, data=pri_ps_genus)

pri_weevil_vs_noweevil_sig_taxa_16s_soil <- as.data.frame((pri_weevil_vs_noweevil_sig_taxa))
write.csv(pri_weevil_vs_noweevil_sig_taxa_16s_soil, "pri_weevil_vs_noweevil_sig_taxa_16s_soil.csv")



#get rel abundance of these taxa





tax_table(pri_ps_genus) <- cbind(tax_table(pri_ps_genus), 
                       rownames(tax_table(pri_ps_genus)))

#colnames(tax_table(pri_ps_genus)) <- 
  c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTUID")

tax_table(pri_ps)
relabun.ps <- transform_sample_counts(pri_ps,function(x) x / sum(x))
pri_ps_genus <- tax_glom(relabun.ps, taxrank = "Genus", NArm = FALSE)


gen_pri_taxa<- ps_genus@tax_table
write.csv(gen_pri_taxa ,"gen_pri_taxa.csv")



pri_ps_genusP <- subset_taxa(pri_ps_genus, Genus %in% c('Rhizobacter','Roseiarcus'))


genus.df <- psmelt(pri_ps_genusP)
head(genus.df)


MySummary <- genus.df %>%
  group_by(Genus) %>%
  summarize(mean_abund = mean(Abundance, na.rm=TRUE)) 
head(MySummary)

write.csv(MySummary, "16s_pri_relabund.csv")

#this file
bacrelabundpri <- read.csv("16s_pri_relabund.csv")

p <- ggplot(bacrelabundpri,aes(x=Genus, y=mean_abund))+
 geom_bar(stat="identity")+ coord_flip() +labs(y='Relative Abundance',
                                                          x='Bacteria Genus')+

theme(axis.text.x=element_text(size=10))+
  theme(axis.title.y = element_text(size=13))+
  theme(axis.text.y=element_text(size=10))+
  theme(axis.title.x=element_text(size=13))+
  theme(axis.text.x=element_text(size=11))
p



# CORNCOB DIFFABUND #### secondary N treat

sec_ps <- microViz:: ps_filter(ps, succession=="Secondary") 

#sec_ps<- transform_sample_counts(sec_ps, fun = function(x){x/sum(x)})


# merge taxa at genus level
sec_ps_genus <- sec_ps %>% 
  tax_glom(taxrank = "Genus")


sec_ps_genus@sam_data$n.treat %>% unique

# set levels of n so "0" is intercept
sec_ps_genus@sam_data$n.treat <- 
  sec_ps_genus@sam_data$n.treat %>% factor(levels = c("0kgN", "16kgN"))


# use raw count data for corncob
da_analysis_n.treat <- differentialTest(formula = ~ n.treat, #abundance
                                        phi.formula = ~ 1, #dispersion
                                        formula_null = ~ 1, #mean
                                        phi.formula_null = ~ 1,
                                        test = "Wald", boot = FALSE,
                                        data = sec_ps_genus,
                                        fdr_cutoff = 0.05,
                                        full_output = TRUE)


p5 <- 
  plot(da_analysis_n.treat, level = c("Family", "Genus")) +
  theme(legend.position = 'none') 

#saveRDS(p5,"./Output/sec_16S_diffabund_ntreat.RDS")



#get taxonomy of sig taxa
sig_taxa_sec <- da_analysis_n.treat$significant_taxa %>% otu_to_taxonomy(data=sec_ps_genus)


sig_taxa_sec_df<- as.data.frame((sig_taxa_sec))
write.csv(sig_taxa_sec_df, "sig_taxa_sec.csv")


#get rel abundance of these taxa
tax_table(sec_ps_genus) <- cbind(tax_table(sec_ps_genus), 
                                 rownames(tax_table(sec_ps_genus)))

#colnames(tax_table(sec_ps_genus)) <- 
  #c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTUID")

tax_table(sec_ps)
relabun.ps <- transform_sample_counts(sec_ps,function(x) x / sum(x))
sec_ps_genus <- tax_glom(relabun.ps, taxrank = "Genus", NArm = FALSE)


gen_sec_taxa<- sec_ps_genus@tax_table
write.csv(gen_sec_taxa ,"gen_sec_taxa.csv")



sec_ps_genusP <- subset_taxa(sec_ps_genus, Genus %in% c('Aquicella','Vampirovibrio', "Magnetospirillum"))


genus.df <- psmelt(sec_ps_genusP )
head(genus.df)


MySummary <- genus.df %>%
  group_by(Genus) %>%
  summarize(mean_abund = mean(Abundance, na.rm=TRUE)) 
head(MySummary)

write.csv(MySummary, "16s_sec_relabund.csv")

#this file
bacrelabundsec <- read.csv("16s_sec_relabund.csv")

p <- ggplot(bacrelabundsec,aes(x=reorder(Genus, mean_abund),y=mean_abund))+
  geom_bar(stat="identity")+ coord_flip() +labs(y='Relative Abundance',
                                                x='Bacteria Genus')
theme(axis.text.x=element_text(size=10))+
  theme(axis.title.y = element_text(size=13))+
  theme(axis.text.y=element_text(size=10))+
  theme(axis.title.x=element_text(size=13))+
  theme(axis.text.x=element_text(size=11))
p


###Look at taxa

#top 20 pri and sec

gen_all <- 
  ps %>% 
  tax_glom("Genus", NArm = FALSE) %>% 
  transform_sample_counts(function(x){x/sum(x)})

gen_all %>% taxa_sums()/100%>% unname() # look at sums for families

write.csv(gen_all@tax_table, "taxa_all_16s.csv")

# top 20 genera
toptwenty_gen_all <-
  taxa_sums(gen_all) %>%
  sort(decreasing = TRUE) %>%
  head(20) %>%
  names() # returns the OTU/ASV names


# subset to those top twenty taxa phyloseq
gen_all_top_20 <-
  gen_all%>%
  subset_taxa(taxa_names(gen_all) %in% toptwenty_gen_all)

gen_all_top_20
# this gives you a phyloseq object with the twenty most abundant genera
# that you can play with as you see fit


#get the csv of top 20  into csv
gen_all_20_16ssoil <- gen_all_top_20@otu_table
write.csv(gen_all_20_16ssoil ,"gen_all_20_16ssoil.csv")


gen_all_20_16ssoil_taxa <- gen_all_top_20@tax_table
write.csv(gen_all_20_16ssoil_taxa ,"gen_all_20_16ssoil_taxa.csv")



#########PRI


gen_pri <- 
  pri_ps %>% 
  tax_glom("Genus", NArm = FALSE) %>% 
  transform_sample_counts(function(x){x/sum(x)})

gen_pri %>% taxa_sums()/100%>% unname() # look at sums for families


# top 20 genera
toptwenty_gen_pri <-
  taxa_sums(gen_pri) %>%
  sort(decreasing = TRUE) %>%
  head(20) %>%
  names() # returns the OTU/ASV names


# subset to those top twenty taxa phyloseq
gen_pri_top_20 <-
  gen_pri%>%
  subset_taxa(taxa_names(gen_pri) %in% toptwenty_gen_pri)

gen_pri_top_20
# this gives you a phyloseq object with the twenty most abundant genera
# that you can play with as you see fit


#get the csv of top 20  into csv
gen_pri_20_16ssoil <- gen_pri_top_20@otu_table
write.csv(gen_pri_20_16ssoil ,"gen_pri_20_16ssoil.csv")


gen_pri_20_16ssoil_taxa <- gen_pri_top_20@tax_table
write.csv(gen_pri_20_16ssoil_taxa ,"gen_pri_20_16ssoil_taxa.csv")


#### look at taxa

ex3 <- subset_taxa(pri_ps, Genus=="Roseiarcus")
#ex4<-subset_samples(ex3, succession=="Primary")
df<-as.data.frame(ex3@tax_table)


ex4 <- subset_taxa(pri_ps, Genus=="Rhizobacter")
#ex4<-subset_samples(ex3, succession=="Primary")
df<-as.data.frame(ex4@tax_table)

#####SEC

gen_sec <- 
  sec_ps %>% 
  tax_glom("Genus", NArm = FALSE) %>% 
  transform_sample_counts(function(x){x/sum(x)})

gen_sec %>% taxa_sums()/100%>% unname() # look at sums for families


# top 20 genera
toptwenty_gen_sec <-
  taxa_sums(gen_sec) %>%
  sort(decreasing = TRUE) %>%
  head(20) %>%
  names() # returns the OTU/ASV names


# subset to those top twenty taxa phyloseq
gen_sec_top_20 <-
  gen_sec%>%
  subset_taxa(taxa_names(gen_sec) %in% toptwenty_gen_sec)


gen_sec_taxa <- gen_sec@tax_table
write.csv(gen_sec_taxa  ,"gen_sec_taxa.csv")

gen_sec_top_20
# this gives you a phyloseq object with the twenty most abundant genera
# that you can play with as you see fit


#get the csv of top 20  into csv
gen_sec_20_16ssoil <- gen_sec_top_20@otu_table
write.csv(gen_sec_20_16ssoil ,"gen_sec_20_16ssoil.csv")


gen_sec_20_16ssoil_taxa <- gen_sec_top_20@tax_table
write.csv(gen_sec_20_16ssoil_taxa ,"gen_sec_20_16ssoil_taxa.csv")



#### look at taxa

ex3 <- subset_taxa(sec_ps, Genus=="Aquicella")
#ex4<-subset_samples(ex3, succession=="Primary")
df<-as.data.frame(ex3@tax_table)
df

ex4 <- subset_taxa(sec_ps, Genus=="Vampirovibrio")
#ex4<-subset_samples(ex3, succession=="Primary")
df<-as.data.frame(ex4@tax_table)
df

ex5 <- subset_taxa(sec_ps, Genus=="Magnetospirillum")
#ex4<-subset_samples(ex3, succession=="Primary")
df<-as.data.frame(ex4@tax_table)
df
