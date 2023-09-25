## Dovrefjell 2001-2022

# set WD
#setwd("D:/Users/thvneste/Documents/GLORIA/Data_doverfjell")
setwd("E:/Users/thvneste/Backup/D/Documents/GLORIA/Data_doverfjell")

# packages
library(lme4)
library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(mgcv)
#options(repos = c(
  #rtrees = 'https://daijiang.r-universe.dev',
  #CRAN = 'https://cloud.r-project.org'))
#install.packages("rtrees")
library(rtrees)
library(ape)
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("ggtree", force = TRUE)
library(ggtree)
library(tidyr)
library(tidyverse)
library(vegan)
library(picante)
library(FD)
library(phytools)
library(Rphylopars)
library(car)
library(pez)
library(ggbeeswarm)
library(egg)
library(effects)
library(AICcmodavg)
library(MuMIn)
library(hillR)
library(fossil)
library(iNEXT)
library(adiv)
library(betapart)
library(reshape2)
library(BAT)

## summit sections ######################################################################################################################################
sites <- read.csv("SAS_veg.csv", header = T)
sites <- sites %>% remove_rownames %>% column_to_rownames(var="id")
species <- colnames(sites)
species <- replace(species,species=="Arctostaphylos_uva.ursi","Arctostaphylos_uva-ursi")
species <- replace(species,species=="Vaccinium_vitis.idaea","Vaccinium_vitis-idaea")
colnames(sites) <- species

# taxonomic diversity
TR <- specnumber(sites)
TD <- diversity(sites, "shannon")
TE <- TD/log(TR)
Tax <- data.frame(TR,TD,TE)

Tax1 <- hill_taxa(sites)
Tax2 <- hill_taxa(sites, q = 2)
Tax3 <- hill_taxa(sites, q = 3)

TD2 <- diversity(sites, "simpson")
TD3 <- diversity(sites, "inv")
TD4 <- fisher.alpha(sites)
TD.extra <- data.frame(TD,TD2,TD3,TD4)

TD.full <- speciesdiv(sites, method = "full", tol = 1e-8)
TD.full <- data.frame(TD.full)
TD.full$id <- rownames(TD.full)
TD.full <- TD.full %>% separate(id, c("year", "summit", "aspect", "plot" , "elevation"))
write.csv(TD.full,"TD_SAS.csv")
TD.full$year <- as.numeric(TD.full$year)
TD.full$elevation <- as.numeric(TD.full$elevation)

# phylogentic diversity

sp_list_df(sp_list = species, taxon = "plant")

tree = get_tree(sp_list = species,
                taxon = "plant",
                scenario = "at_basal_node",
                show_grafted = TRUE)
tree$tip.label <- sub('[*]','',tree$tip.label)
ggtree(tree, branch.length='none', layout='circular') + geom_tiplab(size = 3, hjust = -.05)

sites3 <- sites[rowSums(sites)>2,]
PR <- psr(samp = sites3, tree = tree)
PD <- psv(samp = sites3, tree = tree)
PE <- pse(samp = sites3, tree = tree)
Phy <- data.frame(PR$PSR,PD$PSVs,PE$PSEs)

Phy1 <- hill_phylo(sites3, tree)
Phy2 <- hill_phylo(sites3, tree, q = 2)
Phy3 <- hill_phylo(sites3, tree, q = 3)

# functional diversity
traits <- read.csv("SAS_traits.csv", header = T)
traits <- traits %>% remove_rownames %>% column_to_rownames(var="species")   

# Imputation
# Decomposing phylogenetic distance matrix derived from tree into a set of orthogonal vectors
x <- PVR::PVRdecomp(tree)
# Extract the PVRs
pvrs <- x@Eigen$vectors

# Combine traits and PVRs
traits.pvrs <- cbind(traits, pvrs)

# Imputation using missForest (note that this function have other arguments, see details in Stekhoven and Buehlmann 2012)
RF.imp <- missForest::missForest(traits.pvrs, maxiter = 15, ntree = 100, variablewise = FALSE)

# OOBerror
RF.imp$OOBerror

# Imputed dataset
traits.imp <- RF.imp$ximp[, seq_len(ncol(traits)), drop = FALSE]
traits.imp

sites2 <- sites[rowSums(sites)>0,]
Fun1 <- hill_func(sites2, traits.imp)
Fun2 <- hill_func(sites2, traits.imp, q = 2)
Fun3 <- hill_func(sites2, traits.imp, q = 3)
Fun1 <- Fun1[1,]

sites2 <- sites[rowSums(sites)>0,]
FDi <- dbFD(scale(traits.imp), sites2, corr = "cailliez")
Fct <- data.frame(FDi$FRic,FDi$RaoQ,FDi$FEve)
cwm <- FDi$CWM
cwm$id <- rownames(sites2)
cwm <- cwm %>% separate(id, c("year", "summit", "aspect", "plot" , "elevation"))
x <- traits.imp$plant_height
names(x)<-rownames(traits.imp)
phylosig(tree, x, method="lambda",test=T) # phylogentic signal
phylosig(tree, x, method="K",test=T)
x <- traits.imp$seed_mass
names(x)<-rownames(traits.imp)
phylosig(tree, x, method="lambda",test=T) # phylogentic signal
phylosig(tree, x, method="K",test=T)
x <- traits.imp$sla
names(x)<-rownames(traits.imp)
phylosig(tree, x, method="lambda",test=T) # phylogentic signal
phylosig(tree, x, method="K",test=T)

Tax$id <- rownames(Tax)
Fct$id <- rownames(sites2)
Phy$id <- rownames(PR)
indices <- left_join(Tax,Fct, by = "id")
indices <- left_join(indices,Phy, by = "id")
colnames(indices) <- c("Taxric","Taxdiv","Taxeve", "id", "Funric","Fundiv","Funeve","Phyric","Phydiv","Phyeve")

indices <- indices %>% separate(id, c("year", "summit", "aspect", "plot" , "elevation"))
indices$year <- as.numeric(indices$year)
indices$elevation <- as.numeric(indices$elevation)
write.csv(indices,"indices_SAS.csv")

Tax1 <- data.frame(Tax1)
Tax1$id <- rownames(Tax1)
Fun1 <- data.frame(Fun1)
Fun1$id <- rownames(Fun1)
Phy1 <- data.frame(Phy1)
Phy1$id <- rownames(Phy1)
indices_hill <- left_join(Tax1,Fun1, by = "id")
indices_hill <- left_join(indices_hill,Phy1, by = "id")
indices_hill <- indices_hill %>% separate(id, c("year", "summit", "aspect", "plot" , "elevation"))

# thermophilization
sites <- read.csv("SAS_veg.csv", header = T)
sites <- sites %>% remove_rownames %>% column_to_rownames(var="id")
S <- read.csv("SAS_S.csv", header = T)
S <- S %>% remove_rownames %>% column_to_rownames(var="species")
cwm.S <- as.data.frame(BAT::cwm(sites, S, abund = T, na.rm = T))
cwm.S$id <- rownames(cwm.S)
cwm.S <- cwm.S %>% separate(id, c("year", "summit", "aspect", "plot" , "elevation"))
cwm.S$year <- as.numeric(cwm.S$year)
cwm.S$plot <- as.numeric(cwm.S$plot)
cwm.S$elevation <- as.numeric(cwm.S$elevation)

indices <- read.csv("indices_SAS.csv", header = T)
indices$year <- as.numeric(indices$year)
indices$plot <- as.numeric(indices$plot)
indices$elevation <- as.numeric(indices$elevation)
indices <- left_join(indices,cwm.S, by = c("year", "summit", "aspect", "plot" , "elevation"))
indices <- indices %>% arrange(summit, aspect, plot, year)

diff <- indices %>% select(summit, aspect, plot, year, Taxric, Taxdiv, Funric, Fundiv, Phyric, Phydiv, S) %>% filter(year==2001 | year==2022) %>% group_by(summit, aspect, plot) %>%
        mutate(dTaxric = (Taxric - lag(Taxric))/(year -lag(year)), dTaxdiv = (Taxdiv - lag(Taxdiv))/(year -lag(year)),
               dFunric = (Funric - lag(Funric))/(year -lag(year)), dFundiv = (Fundiv - lag(Fundiv))/(year -lag(year)),
               dPhyric = (Phyric - lag(Phyric))/(year -lag(year)), dPhydiv = (Phydiv - lag(Phydiv))/(year -lag(year)),
               dS = (S - lag(S))/(year - lag(year)))

# general plots
ggplot(indices, aes(x = summit, y = Taxric, fill = year)) + 
  geom_violin(trim = FALSE, position=position_dodge(width = 0.9), alpha=0.5) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1),
               geom="pointrange", color="black",
               shape = 18, size = 0.75,
               position = position_dodge(width = 0.9)) +
  theme_bw() +
  labs(x = "Elevation", y = "Taxonomic richness", color = "Year")

ggplot(indices, aes(x = elevation, y = Fundiv, color = year)) + geom_boxplot()
ggplot(indices, aes(x = elevation, y = Phydiv, color = year)) + geom_boxplot()

ggplot(indices, aes(x = year, y = Funric, color = elevation)) + geom_point()

# models
indices <- read.csv("indices_SAS.csv", header = T)
indices$year <- as.numeric(indices$year)
indices$elevation <- as.numeric(indices$elevation)
m_TR <- glmer(Taxric ~ scale(year) + scale(elevation) + scale(year)*scale(elevation) + (1|summit/aspect), family=poisson(), data=indices, na.action = na.omit)
fixef(m_TR)
Anova(m_TR)
r.squaredGLMM(m_TR)
m_TD <- lmer(Taxdiv ~ scale(year) + scale(elevation) + scale(year)*scale(elevation) + (1|summit/aspect), data=indices, na.action = na.omit)
fixef(m_TD)
Anova(m_TD)
r.squaredGLMM(m_TD)

m_FR <- lmer(Funric ~ scale(year) + scale(elevation) + scale(year)*scale(elevation) + (1|summit/aspect), data=indices, na.action = na.omit)
fixef(m_FR)
Anova(m_FR)
r.squaredGLMM(m_FR)
m_FD <- lmer(Fundiv ~ scale(year) + scale(elevation) + scale(year)*scale(elevation) + (1|summit/aspect), data=indices, na.action = na.omit)
fixef(m_FD)
Anova(m_FD)
r.squaredGLMM(m_FD)
m_FE <- lmer(Funeve ~ scale(year) + scale(elevation) + scale(year)*scale(elevation) + (1|summit/aspect), data=indices, na.action = na.omit)
fixef(m_FE)
Anova(m_FE)
r.squaredGLMM(m_FE)

m_PR <- lmer(Phyric ~ scale(year) + scale(elevation) + scale(year)*scale(elevation) + (1|summit/aspect), data=indices, na.action = na.omit)
fixef(m_PR)
Anova(m_PR)
r.squaredGLMM(m_PR)
m_PD <- lmer(Phydiv ~ scale(year) + scale(elevation) + scale(year)*scale(elevation) + (1|summit/aspect), data=indices, na.action = na.omit)
fixef(m_PD)
Anova(m_PD)
r.squaredGLMM(m_PD)
m_PE <- lmer(Phyeve ~ scale(year) + scale(elevation) + scale(year)*scale(elevation) + (1|summit/aspect), data=indices, na.action = na.omit)
fixef(m_PE)
Anova(m_PE)
r.squaredGLMM(m_PE)

m_TD <- lmer(GiniSimpson ~ scale(year) + scale(elevation) + scale(year)*scale(elevation) + (1|summit/aspect), data=TD.full, na.action = na.omit)
fixef(m_TD)
Anova(m_TD)
r.squaredGLMM(m_TD)
m_TD <- lmer(Margalef ~ scale(year) + scale(elevation) + scale(year)*scale(elevation) + (1|summit/aspect), data=TD.full, na.action = na.omit)
fixef(m_TD)
Anova(m_TD)
r.squaredGLMM(m_TD)
m_TD <- lmer(Menhinick ~ scale(year) + scale(elevation) + scale(year)*scale(elevation) + (1|summit/aspect), data=TD.full, na.action = na.omit)
fixef(m_TD)
Anova(m_TD)
r.squaredGLMM(m_TD)

## permanent 1-m² plots ################################################################################################################################
species <- read.csv("species.csv", header = T)
names <- as.vector(read.csv("species.csv", header = T))
sites <- read.csv("vp2022.csv", header = T)
sites <- subset(sites, select=c(id,Antdio:Vacvit))
sites <- sites %>% remove_rownames %>% column_to_rownames(var="id")
colnames(sites) <- names$species

# taxonomic diversity
TR <- specnumber(sites)
TD <- diversity(sites, "shannon")
TE <- TD/log(TR)
Tax <- data.frame(TR,TD,TE)

TD.full <- speciesdiv(sites, method = "full", tol = 1e-8)
TD.full <- data.frame(TD.full)
TD.full$id <- rownames(TD.full)
TD.full <- TD.full %>% separate(id, c("year", "summit", "aspect", "plot" , "elevation"))
write.csv(TD.full,"TD_plots.csv")
TD.full$year <- as.numeric(TD.full$year)
TD.full$elevation <- as.numeric(TD.full$elevation)

# pylogentic diversity
list <- sp_list_df(sp_list = species, taxon = "plant")

tree = get_tree(sp_list = species,
                     taxon = "plant",
                     scenario = "at_basal_node",
                     show_grafted = TRUE)
tree$tip.label <- sub('[*]','',tree$tip.label)
ggtree(tree, branch.length='none', layout='circular') + geom_tiplab(size = 3, hjust = -.05)

PR <- psr(samp = sites, tree = tree)
PD <- psv(samp = sites, tree = tree)
PE <- pse(samp = sites, tree = tree)
Phy <- data.frame(PR$PSR,PD$PSVs,PE$PSEs)

# functional diversity
traits <- read.csv("traits2.csv", header = T)
traits <- traits %>% remove_rownames %>% column_to_rownames(var="species")
rownames(traits) <- names$species
sites2 <- sites[rowSums(sites)>0,]
FDi <- dbFD(scale(traits), sites2, corr = "cailliez")
Fct <- data.frame(FDi$FRic,FDi$RaoQ,FDi$FEve)

x <- traits$plant_height
names(x)<-rownames(traits)
phylosig(tree, x, method="lambda",test=T) # phylogenetic signal
phylosig(tree, x, method="K",test=T)
x <- traits$seed_mass
names(x)<-rownames(traits)
phylosig(tree, x, method="lambda",test=T) # phylogenetic signal
phylosig(tree, x, method="K",test=T)
x <- traits$sla
names(x)<-rownames(traits)
phylosig(tree, x, method="lambda",test=T) # phylogenetic signal
phylosig(tree, x, method="K",test=T)

Tax$id <- rownames(Tax)
Fct$id <- rownames(sites2)
Phy$id <- rownames(PR)
indices <- left_join(Tax,Fct, by = "id")
indices <- left_join(indices,Phy, by = "id")
colnames(indices) <- c("Taxric","Taxdiv","Taxeve", "id", "Funric","Fundiv","Funeve","Phyric","Phydiv","Phyeve")

indices <- indices %>% separate(id, c("year", "summit", "aspect", "plot" , "elevation"))
indices$year <- as.numeric(indices$year)
indices$elevation <- as.numeric(indices$elevation)
write.csv(indices,"indices_plots.csv")

# models
indices <- read.csv("indices_plots.csv", header = T)
indices$year <- as.numeric(indices$year)
indices$elevation <- as.numeric(indices$elevation)
m_TR <- glmer(Taxric ~ scale(year) + scale(elevation) + scale(year)*scale(elevation) + (1|summit/aspect), family=poisson(), data=indices, na.action = na.omit)
fixef(m_TR)
Anova(m_TR)
r.squaredGLMM(m_TR)
m_TD <- lmer(Taxdiv ~ scale(year) + scale(elevation) + scale(year)*scale(elevation) + (1|summit/aspect), data=indices, na.action = na.omit)
fixef(m_TD)
Anova(m_TD)
r.squaredGLMM(m_TD)
m_TE <- lmer(Taxeve ~ scale(year) + scale(elevation) + scale(year)*scale(elevation) + (1|summit/aspect), data=indices, na.action = na.omit)
fixef(m_TE)
Anova(m_TE)
r.squaredGLMM(m_TE)

m_FR <- lmer(Funric ~ scale(year) + scale(elevation) + scale(year)*scale(elevation) + (1|summit/aspect), data=indices, na.action = na.omit)
fixef(m_FR)
Anova(m_FR)
r.squaredGLMM(m_FR)
m_FD <- lmer(Fundiv ~ scale(year) + scale(elevation) + scale(year)*scale(elevation) + (1|summit/aspect), data=indices, na.action = na.omit)
fixef(m_FD)
Anova(m_FD)
r.squaredGLMM(m_FD)
m_FE <- lmer(Funeve ~ scale(year) + scale(elevation) + scale(year)*scale(elevation) + (1|summit/aspect), data=indices, na.action = na.omit)
fixef(m_FE)
Anova(m_FE)
r.squaredGLMM(m_FE)

m_PR <- lmer(Phyric ~ scale(year) + scale(elevation) + scale(year)*scale(elevation) + (1|summit/aspect), data=indices, na.action = na.omit)
fixef(m_PR)
Anova(m_PR)
r.squaredGLMM(m_PR)
m_PD <- lmer(Phydiv ~ scale(year) + scale(elevation) + scale(year)*scale(elevation) + (1|summit/aspect), data=indices, na.action = na.omit)
fixef(m_PD)
Anova(m_PD)
r.squaredGLMM(m_PD)
m_PE <- lmer(Phyeve ~ scale(year) + scale(elevation) + scale(year)*scale(elevation) + (1|summit/aspect), data=indices, na.action = na.omit)
fixef(m_PE)
Anova(m_PE)
r.squaredGLMM(m_PE)

## plots
indices <- read.csv("indices_SAS.csv", header = T)
m <- glmer(Taxric ~ elevation + (1|summit/aspect), family = poisson(), data = indices)
newdat <- expand.grid(elevation = seq(min(indices$elevation), max(indices$elevation), by = .5))
newdat$pred <- predictSE(m, newdata = newdat, re.form=~0, type = "response")$fit
newdat$se <- predictSE(m, newdata = newdat, re.form=~0, type = "response")$se.fit

p1 <- ggplot() + 
  geom_point(data=indices, aes(elevation, Taxric), alpha = 0.3, size = 2, color = "#fde725") + 
  geom_line(data=newdat, aes(x=elevation, y=pred), color="#fde725") +
  geom_ribbon(data=newdat, aes(x=elevation, ymin=pred-se, ymax=pred+se), alpha= 0.3, fill="#fde725") +
  labs(x="Elevation (m)", y="Taxonomic richness") +
  ggtitle("(a)") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0), axis.title.x=element_blank())


m <- glmer(Taxric ~ year + (1|summit/aspect), family = poisson(), data = indices)
newdat <- expand.grid(year = seq(min(indices$year), max(indices$year), by = .1))
newdat$pred <- predictSE(m, newdata = newdat, re.form=~0, type = "response")$fit
newdat$se <- predictSE(m, newdata = newdat, re.form=~0, type = "response")$se.fit

p2 <- ggplot() + 
  geom_point(data=indices, aes(year, Taxric), alpha = 0.3, size = 2, color = "#fde725") + 
  labs(x="Year", y="Taxonomic richness") +
  ggtitle("(a)") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0), axis.title.x=element_blank())


m <- lmer(Taxdiv ~ elevation + (1|summit/aspect), data = indices)
newdat <- expand.grid(elevation = seq(min(indices$elevation), max(indices$elevation), by = .5))
newdat$pred <- predictSE(m, newdata = newdat, re.form=~0, type = "response")$fit
newdat$se <- predictSE(m, newdata = newdat, re.form=~0, type = "response")$se.fit

p3 <- ggplot() + 
  geom_point(data=indices, aes(elevation, Taxdiv), alpha = 0.3, size = 2, color = "#fde725") + 
  geom_line(data=newdat, aes(x=elevation, y=pred), color="#fde725") +
  geom_ribbon(data=newdat, aes(x=elevation, ymin=pred-se, ymax=pred+se), alpha= 0.3, fill="#fde725") +
  labs(x="Elevation (m)", y="Taxonomic diversity") +
  ggtitle("(d)") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0), axis.title.x=element_blank())

m <- lmer(Taxdiv ~ year + (1|summit/aspect), data = indices)
newdat <- expand.grid(year = seq(min(indices$year), max(indices$year), by = .1))
newdat$pred <- predictSE(m, newdata = newdat, re.form=~0, type = "response")$fit
newdat$se <- predictSE(m, newdata = newdat, re.form=~0, type = "response")$se.fit

p4 <- ggplot() + 
  geom_point(data=indices, aes(year, Taxdiv), alpha = 0.3, size = 2, color = "#fde725") + 
  geom_line(data=newdat, aes(x=year, y=pred), color="#fde725") +
  geom_ribbon(data=newdat, aes(x=year, ymin=pred-se, ymax=pred+se), alpha= 0.3, fill="#fde725") +
  labs(x="Year", y="Taxonomic diversity") +
  ggtitle("(d)") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0), axis.title.x=element_blank())

m <- lmer(Funric ~ elevation + (1|summit/aspect), data = indices)
newdat <- expand.grid(elevation = seq(min(indices$elevation), max(indices$elevation), by = .5))
newdat$pred <- predictSE(m, newdata = newdat, re.form=~0, type = "response")$fit
newdat$se <- predictSE(m, newdata = newdat, re.form=~0, type = "response")$se.fit

p5 <- ggplot() + 
  geom_point(data=indices, aes(elevation, Funric), alpha = 0.3, size = 2, color = "#21918c") + 
  geom_line(data=newdat, aes(x=elevation, y=pred), color="#21918c") +
  geom_ribbon(data=newdat, aes(x=elevation, ymin=pred-se, ymax=pred+se), alpha= 0.3, fill="#21918c") +
  labs(x="Elevation (m)", y="Functional richness") +
  ggtitle("(b)") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0), axis.title.x=element_blank())

indices1 <- subset(indices, elevation==1161)
m <- lmer(Funric ~ year + (1|aspect), data = indices1)
newdat <- expand.grid(year = seq(min(indices$year), max(indices$year), by = .1))
newdat$pred1 <- predictSE(m, newdata = newdat, re.form=~0, type = "response")$fit
newdat$se1 <- predictSE(m, newdata = newdat, re.form=~0, type = "response")$se.fit
indices2 <- subset(indices, elevation==1418)
m <- lmer(Funric ~ year + (1|aspect), data = indices2)
newdat2 <- expand.grid(year = seq(min(indices$year), max(indices$year), by = .1))
newdat$pred2 <- predictSE(m, newdata = newdat2, re.form=~0, type = "response")$fit
newdat$se2 <- predictSE(m, newdata = newdat2, re.form=~0, type = "response")$se.fit
indices3 <- subset(indices, elevation==1651)
m <- lmer(Funric ~ year + (1|aspect), data = indices3)
newdat3 <- expand.grid(year = seq(min(indices$year), max(indices$year), by = .1))
newdat$pred3 <- predictSE(m, newdata = newdat3, re.form=~0, type = "response")$fit
newdat$se3 <- predictSE(m, newdata = newdat3, re.form=~0, type = "response")$se.fit
indices4 <- subset(indices, elevation==1845)
m <- lmer(Funric ~ year + (1|aspect), data = indices4)
newdat4 <- expand.grid(year = seq(min(indices$year), max(indices$year), by = .1))
newdat$pred4 <- predictSE(m, newdata = newdat4, re.form=~0, type = "response")$fit
newdat$se4 <- predictSE(m, newdata = newdat4, re.form=~0, type = "response")$se.fit

colors <- c("1161 m" = "#90d743", "1418 m" = "#35b779", "1651 m" = "#21918c", "1845 m" = "#31688e")
p6 <- ggplot() + 
  geom_point(data=indices1, aes(year, Funric, color = "1161 m"), alpha = 0.3, size = 2) + 
  geom_line(data=newdat, aes(x=year, y=pred1, color="1161 m")) +
  geom_ribbon(data=newdat, aes(x=year, ymin=pred1-se1, ymax=pred1+se1), alpha= 0.3, fill="#90d743") +
  geom_point(data=indices2, aes(year, Funric, color = "1418 m"), alpha = 0.3, size = 2) + 
  geom_line(data=newdat, aes(x=year, y=pred2, color="1418 m")) +
  geom_ribbon(data=newdat, aes(x=year, ymin=pred2-se2, ymax=pred2+se2), alpha= 0.3, fill="#35b779") +
  geom_point(data=indices3, aes(year, Funric, color = "1651 m"), alpha = 0.3, size = 2) + 
  geom_line(data=newdat, aes(x=year, y=pred3, color="1651 m")) +
  geom_ribbon(data=newdat, aes(x=year, ymin=pred3-se3, ymax=pred3+se3), alpha= 0.3, fill="#21918c") +
  geom_point(data=indices4, aes(year, Funric, color = "1845 m"), alpha = 0.3, size = 2) + 
  geom_line(data=newdat, aes(x=year, y=pred4, color="1845 m")) +
  geom_ribbon(data=newdat, aes(x=year, ymin=pred4-se4, ymax=pred4+se4), alpha= 0.3, fill="#31688e") +
  labs(x="Year", y="Functional richness", color = "Elevation (m)") +
  ggtitle("(b)") +
  scale_color_manual(values = colors) +
  guides(color = guide_legend(override.aes = list(size = 1))) +
  theme_bw() +
  theme(legend.position = c(0, 1), legend.justification = c(0, 1), legend.title = element_blank(), legend.box.margin = margin(0, 0, 0, 0), legend.text=element_text(size=7),
        legend.background = element_rect(fill="white",
                                         size=0.5,
                                         linetype="solid",
                                         colour ="gray"),
        legend.key.size = unit(0.2, "cm"), legend.spacing.y = unit(0, "pt"), legend.spacing.x = unit(0, "pt"), plot.title=element_text(hjust=0), axis.title.x=element_blank())

m <- lmer(Fundiv ~ elevation + (1|summit/aspect), data = indices)
newdat <- expand.grid(elevation = seq(min(indices$elevation), max(indices$elevation), by = .5))
newdat$pred <- predictSE(m, newdata = newdat, re.form=~0, type = "response")$fit
newdat$se <- predictSE(m, newdata = newdat, re.form=~0, type = "response")$se.fit

p7 <- ggplot() + 
  geom_point(data=indices, aes(elevation, Fundiv), alpha = 0.3, size = 2, color = "#21918c") + 
  geom_line(data=newdat, aes(x=elevation, y=pred), color="#21918c") +
  geom_ribbon(data=newdat, aes(x=elevation, ymin=pred-se, ymax=pred+se), alpha= 0.3, fill="#21918c") +
  labs(x="Elevation (m)", y="Functional diversity") +
  ggtitle("(e)") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0))

m <- lmer(Fundiv ~ year + (1|summit/aspect), data = indices)
newdat <- expand.grid(year = seq(min(indices$year), max(indices$year), by = .1))
newdat$pred <- predictSE(m, newdata = newdat, re.form=~0, type = "response")$fit
newdat$se <- predictSE(m, newdata = newdat, re.form=~0, type = "response")$se.fit

p8 <- ggplot() + 
  geom_point(data=indices, aes(year, Fundiv), alpha = 0.3, size = 2, color = "#21918c") + 

  labs(x="Year", y="Functional diversity") +
  ggtitle("(e)") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0))

m <- lmer(Phyric ~ elevation + (1|summit/aspect), data = indices)
newdat <- expand.grid(elevation = seq(min(indices$elevation), max(indices$elevation), by = .5))
newdat$pred <- predictSE(m, newdata = newdat, re.form=~0, type = "response")$fit
newdat$se <- predictSE(m, newdata = newdat, re.form=~0, type = "response")$se.fit

p9 <- ggplot() + 
  geom_point(data=indices, aes(elevation, Phyric), alpha = 0.3, size = 2, color = "#440154") + 
  geom_line(data=newdat, aes(x=elevation, y=pred), color="#440154") +
  geom_ribbon(data=newdat, aes(x=elevation, ymin=pred-se, ymax=pred+se), alpha= 0.3, fill="#440154") +
  labs(x="Elevation (m)", y="Phylogenetic richness") +
  ggtitle("(c)") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0), axis.title.x=element_blank())

m <- lmer(Phyric ~ year + (1|summit/aspect), data = indices)
newdat <- expand.grid(year = seq(min(indices$year), max(indices$year), by = .1))
newdat$pred <- predictSE(m, newdata = newdat, re.form=~0, type = "response")$fit
newdat$se <- predictSE(m, newdata = newdat, re.form=~0, type = "response")$se.fit

p10 <- ggplot() + 
  geom_point(data=indices, aes(year, Phyric), alpha = 0.3, size = 2, color = "#440154") + 
  geom_line(data=newdat, aes(x=year, y=pred), color="#440154") +
  geom_ribbon(data=newdat, aes(x=year, ymin=pred-se, ymax=pred+se), alpha= 0.3, fill="#440154") +
  labs(x="Year", y="Phylogenetic richness") +
  ggtitle("(c)") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0), axis.title.x=element_blank())

m <- lmer(Phydiv ~ elevation + (1|summit/aspect), data = indices)
newdat <- expand.grid(elevation = seq(min(indices$elevation), max(indices$elevation), by = .5))
newdat$pred <- predictSE(m, newdata = newdat, re.form=~0, type = "response")$fit
newdat$se <- predictSE(m, newdata = newdat, re.form=~0, type = "response")$se.fit

p11 <- ggplot() + 
  geom_point(data=indices, aes(elevation, Phydiv), alpha = 0.3, size = 2, color = "#440154") + 
  geom_line(data=newdat, aes(x=elevation, y=pred), color="#440154") +
  geom_ribbon(data=newdat, aes(x=elevation, ymin=pred-se, ymax=pred+se), alpha= 0.3, fill="#440154") +
  labs(x="Elevation (m)", y="Phylogenetic diversity") +
  ggtitle("(f)") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0), axis.title.x=element_blank())

m <- lmer(Phydiv ~ year + (1|summit/aspect), data = indices)
newdat <- expand.grid(year = seq(min(indices$year), max(indices$year), by = .1))
newdat$pred <- predictSE(m, newdata = newdat, re.form=~0, type = "response")$fit
newdat$se <- predictSE(m, newdata = newdat, re.form=~0, type = "response")$se.fit

p12 <- ggplot() + 
  geom_point(data=indices, aes(year, Phydiv), alpha = 0.3, size = 2, color = "#440154") + 
  labs(x="Year", y="Phylogenetic diversity") +
  ggtitle("(f)") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0), axis.title.x=element_blank())

p <- ggarrange(p1,p5,p9,p3,p7,p11, ncol = 3, nrow = 2)
png(file="elevation_SAS.png",width=16,height=10,units="cm",res=300, pointsize=12)
p
dev.off()

p <- ggarrange(p2,p6,p10,p4,p8,p12, ncol = 3, nrow = 2)
png(file="year_SAS.png",width=16,height=10,units="cm",res=300, pointsize=12)
p
dev.off()

# SAS
indices <- read.csv("indices_SAS.csv", header = T)
indices$year <- as.factor(indices$year)
indices$elevation <- as.factor(indices$elevation)
p1<-ggplot(indices, aes(elevation, Taxric, col = year)) + # y: iq
  geom_point(alpha = 0.5, position = position_dodge(0.6)) +
  stat_summary(fun = mean, geom = 'point', size = 3, position = position_dodge(0.6)) + 
  stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', 
               width = 0, size = 1, position = position_dodge(0.6)) +
  ggtitle("(a)") +
  theme_bw() +
  theme(legend.position = "none", plot.title=element_text(hjust=0), axis.title.x=element_blank(), axis.text.x=element_blank()) +
  labs(x="Elevation",y="Taxonomic richness",color="Year")

p2<-ggplot(indices, aes(elevation, Taxdiv, col = year)) + # y: iq
  geom_point(alpha = 0.5, position = position_dodge(0.6)) +
  stat_summary(fun = mean, geom = 'point', size = 3, position = position_dodge(0.6)) + 
  stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', 
               width = 0, size = 1, position = position_dodge(0.6)) +
  ggtitle("(b)") +
  theme_bw() +
  theme(legend.position = c(1,1), legend.justification = c(1,1), legend.title=element_blank(), legend.text=element_text(size=7), legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"), plot.title=element_text(hjust=0), axis.title.x=element_blank(), axis.text.x=element_blank()) +
  labs(x="Elevation",y="Taxonomic diversity",color="Year")

p3<-ggplot(indices, aes(elevation, Funric, col = year)) + # y: iq
  geom_point(alpha = 0.5, position = position_dodge(0.6)) +
  stat_summary(fun = mean, geom = 'point', size = 3, position = position_dodge(0.6)) + 
  stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', 
               width = 0, size = 1, position = position_dodge(0.6)) +
  ggtitle("(c)") +
  theme_bw() +
  theme(legend.position = "none", plot.title=element_text(hjust=0), axis.title.x=element_blank(), axis.text.x=element_blank()) +
  labs(x="Elevation",y="Functional richness",color="Year")

p4<-ggplot(indices, aes(elevation, Fundiv, col = year)) + # y: iq
  geom_point(alpha = 0.5, position = position_dodge(0.6)) +
  stat_summary(fun = mean, geom = 'point', size = 3, position = position_dodge(0.6)) + 
  stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', 
               width = 0, size = 1, position = position_dodge(0.6)) +
  ggtitle("(d)") +
  theme_bw() +
  theme(legend.position = "none", plot.title=element_text(hjust=0), axis.title.x=element_blank(), axis.text.x=element_blank()) +
  labs(x="Elevation",y="Functional diversity",color="Year")

p5<-ggplot(indices, aes(elevation, Phyric, col = year)) + # y: iq
  geom_point(alpha = 0.5, position = position_dodge(0.6)) +
  stat_summary(fun = mean, geom = 'point', size = 3, position = position_dodge(0.6)) + 
  stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', 
               width = 0, size = 1, position = position_dodge(0.6)) +
  ggtitle("(e)") +
  theme_bw() +
  theme(legend.position = "none", plot.title=element_text(hjust=0)) +
  labs(x="Elevation (m)",y="Phylogenetic richness",color="Year")

p6<-ggplot(indices, aes(elevation, Phydiv, col = year)) + # y: iq
  geom_point(alpha = 0.5, position = position_dodge(0.6)) +
  stat_summary(fun = mean, geom = 'point', size = 3, position = position_dodge(0.6)) + 
  stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', 
               width = 0, size = 1, position = position_dodge(0.6)) +
  ggtitle("(f)") +
  theme_bw() +
  theme(legend.position = "none", plot.title=element_text(hjust=0)) +
  labs(x="Elevation (m)",y="Phylogenetic diversity",color="Year")

p <- ggarrange(p1,p2,p3,p4,p5,p6, ncol = 2, nrow = 3)
png(file="indices_SAS.png",width=16,height=23,units="cm",res=300, pointsize=12)
p
dev.off()

# plots
indices <- read.csv("indices_plots.csv", header = T)
indices$year <- as.factor(indices$year)
indices$elevation <- as.factor(indices$elevation)
p1<-ggplot(indices, aes(elevation, Taxric, col = year)) + # y: iq
  geom_point(alpha = 0.5, position = position_dodge(0.6)) +
  stat_summary(fun = mean, geom = 'point', size = 3, position = position_dodge(0.6)) + 
  stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', 
               width = 0, size = 1, position = position_dodge(0.6)) +
  ggtitle("(a)") +
  theme_bw() +
  theme(legend.position = "none", plot.title=element_text(hjust=0), axis.title.x=element_blank(), axis.text.x=element_blank()) +
  labs(x="Elevation",y="Taxonomic richness",color="Year")

p2<-ggplot(indices, aes(elevation, Taxdiv, col = year)) + # y: iq
  geom_point(alpha = 0.5, position = position_dodge(0.6)) +
  stat_summary(fun = mean, geom = 'point', size = 3, position = position_dodge(0.6)) + 
  stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', 
               width = 0, size = 1, position = position_dodge(0.6)) +
  ggtitle("(b)") +
  theme_bw() +
  theme(legend.position = c(1,1), legend.justification = c(1,1), legend.title=element_blank(), legend.text=element_text(size=7), legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"), plot.title=element_text(hjust=0), axis.title.x=element_blank(), axis.text.x=element_blank()) +
  labs(x="Elevation",y="Taxonomic diversity",color="Year")

p3<-ggplot(indices, aes(elevation, Funric, col = year)) + # y: iq
  geom_point(alpha = 0.5, position = position_dodge(0.6)) +
  stat_summary(fun = mean, geom = 'point', size = 3, position = position_dodge(0.6)) + 
  stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', 
               width = 0, size = 1, position = position_dodge(0.6)) +
  ggtitle("(c)") +
  theme_bw() +
  theme(legend.position = "none", plot.title=element_text(hjust=0), axis.title.x=element_blank(), axis.text.x=element_blank()) +
  labs(x="Elevation",y="Functional richness",color="Year")

p4<-ggplot(indices, aes(elevation, Fundiv, col = year)) + # y: iq
  geom_point(alpha = 0.5, position = position_dodge(0.6)) +
  stat_summary(fun = mean, geom = 'point', size = 3, position = position_dodge(0.6)) + 
  stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', 
               width = 0, size = 1, position = position_dodge(0.6)) +
  ggtitle("(d)") +
  theme_bw() +
  theme(legend.position = "none", plot.title=element_text(hjust=0), axis.title.x=element_blank(), axis.text.x=element_blank()) +
  labs(x="Elevation",y="Functional diversity",color="Year")

p5<-ggplot(indices, aes(elevation, Phyric, col = year)) + # y: iq
  geom_point(alpha = 0.5, position = position_dodge(0.6)) +
  stat_summary(fun = mean, geom = 'point', size = 3, position = position_dodge(0.6)) + 
  stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', 
               width = 0, size = 1, position = position_dodge(0.6)) +
  ggtitle("(e)") +
  theme_bw() +
  theme(legend.position = "none", plot.title=element_text(hjust=0)) +
  labs(x="Elevation (m)",y="Phylogenetic richness",color="Year")

p6<-ggplot(indices, aes(elevation, Phydiv, col = year)) + # y: iq
  geom_point(alpha = 0.5, position = position_dodge(0.6)) +
  stat_summary(fun = mean, geom = 'point', size = 3, position = position_dodge(0.6)) + 
  stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', 
               width = 0, size = 1, position = position_dodge(0.6)) +
  ggtitle("(f)") +
  theme_bw() +
  theme(legend.position = "none", plot.title=element_text(hjust=0)) +
  labs(x="Elevation (m)",y="Phylogenetic diversity",color="Year")

p <- ggarrange(p1,p2,p3,p4,p5,p6, ncol = 2, nrow = 3)
png(file="indices_plots.png",width=16,height=23,units="cm",res=300, pointsize=12)
p
dev.off()

## correltation plots
indices <- read.csv("indices_SAS.csv", header = T)
indices <- na.omit(data.frame(indices$Taxric,indices$Taxdiv,indices$Funric,indices$Fundiv,indices$Phyric,indices$Phydiv))
colnames(indices) <- c("Taxric","Taxdiv","Funric","Fundiv","Phyric","Phydiv")
M <-cor(indices)
corrplot(M, method ="number")

indices <- read.csv("indices_plot.csv", header = T)
indices <- na.omit(data.frame(indices$Taxric,indices$Taxdiv,indices$Taxeve,indices$Funric,indices$Fundiv,indices$Funeve,indices$Phyric,indices$Phydiv,indices$Phyeve))
colnames(indices) <- c("Taxric","Taxdiv","Taxeve","Funric","Fundiv","Funeve","Phyric","Phydiv","Phyeve")
M <-cor(indices)
corrplot(M, method ="number")

# diversity indices plots
TD.full <- read.csv("TD_SAS.csv", header = T)
TD.full$year <- as.factor(TD.full$year)
TD.full$elevation <- as.factor(TD.full$elevation)
p1<-ggplot(TD.full, aes(elevation, GiniSimpson, col = year)) + # y: iq
  geom_point(alpha = 0.5, position = position_dodge(0.6)) +
  stat_summary(fun = mean, geom = 'point', size = 3, position = position_dodge(0.6)) + 
  stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', 
               width = 0, size = 1, position = position_dodge(0.6)) +
  ggtitle("(a)") +
  theme_bw() +
  theme(legend.position = "none", plot.title=element_text(hjust=0), axis.title.x=element_blank(), axis.text.x=element_blank()) +
  labs(x="Elevation (m)",y="Gini-Simpson index",color="Year")

p2<-ggplot(TD.full, aes(elevation, Shannon, col = year)) + # y: iq
  geom_point(alpha = 0.5, position = position_dodge(0.6)) +
  stat_summary(fun = mean, geom = 'point', size = 3, position = position_dodge(0.6)) + 
  stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', 
               width = 0, size = 1, position = position_dodge(0.6)) +
  ggtitle("(b)") +
  theme_bw() +
  theme(legend.position = c(1,1), legend.justification = c(1,1), legend.title=element_blank(), legend.text=element_text(size=7), legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"), plot.title=element_text(hjust=0), axis.title.x=element_blank(), axis.text.x=element_blank()) +
  labs(x="Elevation (m)",y="Shannon index",color="Year")

p3<-ggplot(TD.full, aes(elevation, Margalef, col = year)) + # y: iq
  geom_point(alpha = 0.5, position = position_dodge(0.6)) +
  stat_summary(fun = mean, geom = 'point', size = 3, position = position_dodge(0.6)) + 
  stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', 
               width = 0, size = 1, position = position_dodge(0.6)) +
  ggtitle("(c)") +
  theme_bw() +
  theme(legend.position = "none", plot.title=element_text(hjust=0)) +
  labs(x="Elevation (m)",y="Margalef's index",color="Year")

p4<-ggplot(TD.full, aes(elevation, Menhinick, col = year)) + # y: iq
  geom_point(alpha = 0.5, position = position_dodge(0.6)) +
  stat_summary(fun = mean, geom = 'point', size = 3, position = position_dodge(0.6)) + 
  stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', 
               width = 0, size = 1, position = position_dodge(0.6)) +
  ggtitle("(d)") +
  theme_bw() +
  theme(legend.position = "none", plot.title=element_text(hjust=0)) +
  labs(x="Elevation (m)",y="Menhinick's index",color="Year")

p <- ggarrange(p1,p2,p3,p4, ncol = 2, nrow = 2)
png(file="diversity_SAS.png",width=16,height=16,units="cm",res=300, pointsize=12)
p
dev.off()

TD.full <- read.csv("TD_plots.csv", header = T)
TD.full$year <- as.factor(TD.full$year)
TD.full$elevation <- as.factor(TD.full$elevation)
p1<-ggplot(TD.full, aes(elevation, GiniSimpson, col = year)) + # y: iq
  geom_point(alpha = 0.5, position = position_dodge(0.6)) +
  stat_summary(fun = mean, geom = 'point', size = 3, position = position_dodge(0.6)) + 
  stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', 
               width = 0, size = 1, position = position_dodge(0.6)) +
  ggtitle("(a)") +
  theme_bw() +
  theme(legend.position = "none", plot.title=element_text(hjust=0), axis.title.x=element_blank(), axis.text.x=element_blank()) +
  labs(x="Elevation (m)",y="Gini-Simpson index",color="Year")

p2<-ggplot(TD.full, aes(elevation, Shannon, col = year)) + # y: iq
  geom_point(alpha = 0.5, position = position_dodge(0.6)) +
  stat_summary(fun = mean, geom = 'point', size = 3, position = position_dodge(0.6)) + 
  stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', 
               width = 0, size = 1, position = position_dodge(0.6)) +
  ggtitle("(b)") +
  theme_bw() +
  theme(legend.position = c(1,1), legend.justification = c(1,1), legend.title=element_blank(), legend.text=element_text(size=7), legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"), plot.title=element_text(hjust=0), axis.title.x=element_blank(), axis.text.x=element_blank()) +
  labs(x="Elevation (m)",y="Shannon index",color="Year")

p3<-ggplot(TD.full, aes(elevation, Margalef, col = year)) + # y: iq
  geom_point(alpha = 0.5, position = position_dodge(0.6)) +
  stat_summary(fun = mean, geom = 'point', size = 3, position = position_dodge(0.6)) + 
  stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', 
               width = 0, size = 1, position = position_dodge(0.6)) +
  ggtitle("(c)") +
  theme_bw() +
  theme(legend.position = "none", plot.title=element_text(hjust=0)) +
  labs(x="Elevation (m)",y="Margalef's index",color="Year")

p4<-ggplot(TD.full, aes(elevation, Menhinick, col = year)) + # y: iq
  geom_point(alpha = 0.5, position = position_dodge(0.6)) +
  stat_summary(fun = mean, geom = 'point', size = 3, position = position_dodge(0.6)) + 
  stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', 
               width = 0, size = 1, position = position_dodge(0.6)) +
  ggtitle("(d)") +
  theme_bw() +
  theme(legend.position = "none", plot.title=element_text(hjust=0)) +
  labs(x="Elevation (m)",y="Menhinick's index",color="Year")

p <- ggarrange(p1,p2,p3,p4, ncol = 2, nrow = 2)
png(file="diversity_plots.png",width=16,height=16,units="cm",res=300, pointsize=12)
p
dev.off()

# CWM plots
p1<-ggplot(cwm, aes(elevation, plant_height, col = year)) + # y: iq
  geom_point(alpha = 0.5, position = position_dodge(0.6)) +
  stat_summary(fun = mean, geom = 'point', size = 3, position = position_dodge(0.6)) + 
  stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', 
               width = 0, size = 1, position = position_dodge(0.6)) +
  ggtitle("(a)") +
  theme_bw() +
  theme(legend.position = "none", plot.title=element_text(hjust=0), axis.title.x=element_blank(), axis.text.x=element_blank()) +
  labs(x="Elevation (m)",y="CWM plant height (m)",color="Year")

p2<-ggplot(cwm, aes(elevation, seed_mass, col = year)) + # y: iq
  geom_point(alpha = 0.5, position = position_dodge(0.6)) +
  stat_summary(fun = mean, geom = 'point', size = 3, position = position_dodge(0.6)) + 
  stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', 
               width = 0, size = 1, position = position_dodge(0.6)) +
  ggtitle("(b)") +
  theme_bw() +
  theme(legend.position = "none", plot.title=element_text(hjust=0)) +
  labs(x="Elevation (m)",y="CWM seed mass (mg)",color="Year")

p3<-ggplot(cwm, aes(elevation, sla, col = year)) + # y: iq
  geom_point(alpha = 0.5, position = position_dodge(0.6)) +
  stat_summary(fun = mean, geom = 'point', size = 3, position = position_dodge(0.6)) + 
  stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', 
               width = 0, size = 1, position = position_dodge(0.6)) +
  ggtitle("(c)") +
  theme_bw() +
  theme(legend.position = c(1,1), legend.justification = c(1,1), legend.title=element_blank(), legend.text=element_text(size=7), legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"), plot.title=element_text(hjust=0), axis.title.x=element_blank(), axis.text.x=element_blank()) +
  labs(x="Elevation (m)",y="CWM SLA (mm²/mg)",color="Year")

p <- ggarrange(p1,p2,p3, ncol = 3, nrow = 1)
png(file="CWM_SAS.png",width=25,height=20,units="cm",res=300, pointsize=12)
p
dev.off()

# phylogenetic diversity as SES MPD
test <- ses.mpd(sites2, dist(traits.imp), null.model = "taxa.labels", abundance.weighted = FALSE, runs = 100)
id<-rownames(sites2)
test<-data.frame(id,test$mpd.obs.z)
test <- test %>% separate(id, c("year", "summit", "aspect", "plot" , "elevation"))
test$year <- as.numeric(test$year)
test$elevation <- as.numeric(test$elevation)
m_PR <- lmer(test.mpd.obs.z ~ scale(year) + scale(elevation) + scale(year)*scale(elevation) + (1|summit/aspect), data=test, na.action = na.omit)
fixef(m_PR)
Anova(m_PR)
r.squaredGLMM(m_PR)

test <- ses.mpd(sites3, cophenetic(tree), null.model = "taxa.labels", abundance.weighted = FALSE, runs = 100)
id<-rownames(sites3)
test<-data.frame(id,test$mpd.obs.z)
test <- test %>% separate(id, c("year", "summit", "aspect", "plot" , "elevation"))
test$year <- as.numeric(test$year)
test$elevation <- as.numeric(test$elevation)
m_PR <- lmer(test.mpd.obs.z ~ scale(year) + scale(elevation) + scale(year)*scale(elevation) + (1|summit/aspect), data=test, na.action = na.omit)
fixef(m_PR)
Anova(m_PR)
r.squaredGLMM(m_PR)

## species list
sites <- read.csv("SAS_veg.csv", header = T)
sites <- melt(sites, id.vars=c("id"))
sites <- sites %>% separate(id, c("year", "summit", "aspect", "plot" , "elevation"))
sites <- dcast(sites, variable ~ year, value.var="value", fun.aggregate = max)
write.csv(sites,"species_table.csv")