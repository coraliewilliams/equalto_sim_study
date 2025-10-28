
#### Meta-analysis 2023 ####

setwd("~/PhD work/meta-analysis/Code")



# Load packages

library(metafor)
library(tidyverse)
library(pacman)
library(orchaRd)
library(rotl)
library(ape)
library(ggplot2)
library(ggtext)
library(patchwork)

# Load data

ef <- read.csv("effect_sizes.csv") 

options(na.action = "na.pass")

## calculate absolute ef

yi_ab <- abs(ef$yi)
ef_ab <- ef
ef_ab$yi <- yi_ab  #converts effect sizes to absolute values


# construct phylogentic tree matrix for use as random factor

species <- unique(ef$species) # list of unique species in meta-analysis
species <- as.character(species) # change to character object
taxa <- tnrs_match_names(species)
tree <- rotl::tol_induced_subtree(taxa$ott_id)
tree$tip.label <-
  strip_ott_ids(tree$tip.label, remove_underscores = TRUE) # change ids to the names from dataset

# calculate correlations between all species = cor

tree2 <- compute.brlen(tree)
cor <- vcv(tree2, cor = T)

# add phylo random factor

ef$phylo <- ef$species
ef_ab$phylo <- ef_ab$species

# Subgroups for each trait category   

# directional 

phenolics <- subset(ef, trait_category == "phenolics")
LDMC <- subset(ef, trait_category == "LDMC")
C <- subset(ef, trait_category == "C")
N <- subset(ef, trait_category == "N")
terpenes <- subset(ef, trait_category == "terpenes")
thickness <- subset(ef, trait_category == "thickness")
toughness <- subset(ef, trait_category == "toughness")
SLA <- subset(ef, trait_category == "SLA")

# absolute

phenolics_ab <- subset(ef_ab, trait_category == "phenolics")
LDMC_ab <- subset(ef_ab, trait_category == "LDMC")
C_ab <- subset(ef_ab, trait_category == "C")
N_ab <- subset(ef_ab, trait_category == "N")
terpenes_ab <- subset(ef_ab, trait_category == "terpenes")
thickness_ab <- subset(ef_ab, trait_category == "thickness")
toughness_ab <- subset(ef_ab, trait_category == "toughness")
SLA_ab <- subset(ef_ab, trait_category == "SLA")

##############################################################################
    
#                     M E T A - A N A L Y S E S                              #

##############################################################################


# Meta-analyses for 8 leaf traits
# random factors = ACC (study), experiment (study site), species, phylogeny  R = list(species = cor) (phylogeny) and individual effect size (ef)


#--------------------------------#
       ### thickness ###
#--------------------------------#


# thickness

thick_ma<- rma.mv(yi,vi,data=ef, subset=(trait_category=="thickness"), method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                 R = list(phylo = cor), control = list(optimizer = "optim"))

summary(thick_ma)


# thickness absolute

thick_ab_ma<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="thickness"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"))

summary(thick_ab_ma)


#--------------------------------#
       ### toughness ###
#--------------------------------#

# toughness

tough_ma<- rma.mv(yi,vi,data=ef, subset=(trait_category=="toughness"), method= "REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"))

summary(tough_ma)



# toughness absolute

tough_ab_ma<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="toughness"), method="REML",
                  random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                  R = list(phylo = cor), control = list(optimizer = "optim"))

summary(tough_ab_ma)


#--------------------------------#
        ### terpenes ###
#--------------------------------#

# terpenes 

terp_ma<- rma.mv(yi,vi,data=ef, subset=(trait_category=="terpenes"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"))

summary(terp_ma)

# terpenes absolute

terp_ab_ma<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="terpenes"),method="REML",
                  random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                  R = list(phylo = cor), control = list(optimizer = "optim"))

summary(terp_ab_ma)


#--------------------------------#
           ### LDMC ###
#--------------------------------#

# LDMC

LDMC_ma<- rma.mv(yi,vi,data=ef, subset=(trait_category=="LDMC"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"))

summary(LDMC_ma)


# LDMC absolute

LDMC_ab_ma<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="LDMC"),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                 R = list(phylo = cor), control = list(optimizer = "optim"))

summary(LDMC_ab_ma)


#--------------------------------#
           ### SLA ###
#--------------------------------#


# SLA

SLA_ma<- rma.mv(yi,vi,data=ef, subset=(trait_category=="SLA"),method="REML",
             random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"))

summary(SLA_ma)


# SLA absolute

SLA_ab_ma<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="SLA"),method="REML",
                random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                 R = list(phylo = cor), control = list(optimizer = "optim"))

summary(SLA_ab_ma)


#--------------------------------#
           ### C ###
#--------------------------------#

# C

C_ma<- rma.mv(yi,vi,data=ef, subset=(trait_category=="C"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"))

summary(C_ma)


# C absolute

C_ab_ma<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="C"),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                 R = list(phylo = cor), control = list(optimizer = "optim"))

summary(C_ab_ma)



#--------------------------------#
            ### N ###
#--------------------------------#

# N

N_ma<- rma.mv(yi,vi,data=ef, subset=(trait_category=="N"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"))

summary(N_ma)



# N absolute

N_ab_ma<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="N"),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                 R = list(phylo = cor), control = list(optimizer = "optim"))

summary(N_ab_ma)


#--------------------------------#
        ### phenolics ###
#--------------------------------#

# phenolics

phenolics_ma<- rma.mv(yi,vi,data=ef, subset=(trait_category=="phenolics"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"))


summary(phenolics_ma)

# phenolics absolute

phenolics_ab_ma<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="phenolics"),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                 R = list(phylo = cor), control = list(optimizer = "optim"))

summary(phenolics_ab_ma)


# phenolics by type

phenolics$defence_type_specific <- as.factor(phenolics$defence_type_specific) 

phenolics$defence_type_specific <- factor(phenolics$defence_type_specific, levels = c("total phenolics", "hydrolysable tannins", "condensed tannins", "flavonoids", "lignins"))

phenolics_type <- rma.mv(yi,vi,data=phenolics, method="REML",
                         random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                         R = list(phylo = cor), control = list(optimizer = "optim"),
                         mods = ~ defence_type_specific-1)


phenolics_type



# -------------------------#
### prediction intervals ###
# -------------------------#

predict(tough_ma)
predict(thick_ma)
predict(LDMC_ma)
predict(SLA_ma)
predict(terp_ma)
predict(phenolics_ma)
predict(N_ma)
predict(C_ma)


predict(tough_ab_ma)
predict(thick_ab_ma)
predict(LDMC_ab_ma)
predict(SLA_ab_ma)
predict(terp_ab_ma)
predict(phenolics_ab_ma)
predict(N_ab_ma)
predict(C_ab_ma)

#-------------------#
### orchard plots ###
#-------------------#

# directional #

toughnesso <- orchard_plot(tough_ma, group = "ACC", data = toughness, xlab = "Standardised mean difference", transfm = "none", angle = 45, k =FALSE, legend.pos = "none") +
ylim(-5, 5) +
annotate(geom="text", x= 1.5, y= 3.5,  label = "k(N) = 20(3)", size = 6) +
annotate(geom="text", x= 1.5, y= -3.5,  label = "Toughness", size = 6) +
  scale_fill_manual(values = "#56B4E9") + scale_colour_manual(values = "#56B4E9") +
   theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
               axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank())

toughnesso 


thicknesso <- orchard_plot(thick_ma, group = "ACC", data = thickness, xlab = "Standardised mean difference", transfm = "none", angle = 45, k =FALSE, legend.pos = "none") +
ylim(-5, 5) +
annotate(geom="text", x= 1.5, y= 3.5,  label = "k(N) = 20(3)", size = 6) +
annotate(geom="text", x= 1.5, y= -3.5,  label = "Thickness", size = 6) +
  scale_fill_manual(values = "#E69F00") + scale_colour_manual(values = "#E69F00") +
   theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
               axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank())

thicknesso 

LDMCo <- orchard_plot(LDMC_ma, group = "ACC", data = LDMC, xlab = "Standardised mean difference", transfm = "none", angle = 45, k =FALSE, legend.pos = "none") +
ylim(-5, 5) +
annotate(geom="text", x= 1.5, y= 3.5,  label = "k(N) = 119(9)", size = 6) +
annotate(geom="text", x= 1.5, y= -3.5,  label = "LDMC", size = 6) +
  scale_fill_manual(values = "#F0E442") + scale_colour_manual(values = "#F0E442") +
   theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank())
LDMCo 


SLAo <- orchard_plot(SLA_ma, group = "ACC", data = SLA, xlab = "Standardised mean difference", transfm = "none", angle = 45, k =FALSE, legend.pos = "none") +
ylim(-5, 5) +
annotate(geom="text", x= 1.5, y= 3.5,  label = "k(N) = 251(17)", size = 6) +
annotate(geom="text", x= 1.5, y= -3.5,  label = "SLA", size = 6) +
  scale_fill_manual(values = "#0072B2") + scale_colour_manual(values = "#0072B2") +
   theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
         axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank())
SLAo


terpeneso <- orchard_plot(terp_ma, group = "ACC", data = terpenes, xlab = "Standardised mean difference", transfm = "none", angle = 45, k =FALSE, legend.pos = "none") +
ylim(-5, 5) +
annotate(geom="text", x= 1.5, y= 3.5,  label = "k(N) = 24(6)", size = 6) +
annotate(geom="text", x= 1.5, y= -3.5,  label = "Terpenoids", size = 6) +
  scale_fill_manual(values = "#CC79A7") + scale_colour_manual(values = "#CC79A7") +
   theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
           axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank())

terpeneso 

phenolicso <- orchard_plot(phenolics_ma, group = "ACC", data = phenolics, xlab = "Standardised mean difference", transfm = "none", angle = 45, k =FALSE, legend.pos = "none") +
ylim(-5, 5) +
annotate(geom="text", x= 1.5, y= 3.5,  label = "k(N) = 228(13)", size = 6) +
annotate(geom="text", x= 1.5, y= -3.5,  label = "Phenolics", size = 6) +
  scale_fill_manual(values = "#D55E00" ) + scale_colour_manual(values = "#D55E00") +
   theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank())

phenolicso 


No <- orchard_plot(N_ma, group = "ACC", data = N, xlab = "Standardised mean difference", transfm = "none", angle = 45, k =FALSE, legend.pos = "none") +
ylim(-5, 5) +
annotate(geom="text", x= 1.5, y= 3.5,  label = "k(N) = 206(27)", size = 6) +
annotate(geom="text", x= 1.5, y= -3.5,  label = "Nitrogen", size = 6) +
 scale_fill_manual(values = "#009E73" ) + scale_colour_manual(values = "#009E73") +
   theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
          axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank())

No 

Co <- orchard_plot(C_ma, group = "ACC", data = C, xlab = "Sandardised mean difference", transfm = "none", angle = 45, k =FALSE, legend.pos = "none") +
ylim(-5, 5) +
annotate(geom="text", x= 1.5, y= 3.5,  label = "k(N) = 139(11)", size = 6) +
annotate(geom="text", x= 1.5, y= -3.5,  label = "Carbon", size = 6) +
  scale_fill_manual(values = "#999999" ) + scale_colour_manual(values = "#999999" ) +
   theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_blank())
Co 


thicknesso / toughnesso / 
LDMCo / SLAo /
terpeneso  / phenolicso  /
No / Co 


# absolute #


thicknesso_ab <- orchard_plot(thick_ab_ma, group = "ACC", data = thickness_ab, xlab = "Standardised mean difference", transfm = "none", angle = 45, k =FALSE, legend.pos = "none") +
ylim(-1, 5)  +
annotate(geom="text", x= 1.5, y= 3.5,  label = "k(N) = 20(3)", size = 6) +
annotate(geom="text", x= 1.5, y= -0.6,  label = "Thickness", size = 6) +
  scale_fill_manual(values = "#E69F00") + scale_colour_manual(values = "#E69F00") +
   theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
          axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank())

thicknesso_ab 

toughnesso_ab <- orchard_plot(tough_ab_ma, group = "ACC", data = toughness_ab, xlab = "Standardised mean difference", transfm = "none", angle = 45, k =FALSE, legend.pos = "none") +
ylim(-1, 5)  +
annotate(geom="text", x= 1.5, y= 3.5,  label = "k(N) = 20(3)", size = 6) +
annotate(geom="text", x= 1.5, y= -0.6,  label = "Toughness", size = 6) +
  scale_fill_manual(values = "#56B4E9") + scale_colour_manual(values = "#56B4E9") +
   theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
          axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank())

toughnesso_ab 

LDMCo_ab <- orchard_plot(LDMC_ab_ma, group = "ACC", data = LDMC_ab, xlab = "Standardised mean difference", transfm = "none", angle = 45, k =FALSE, legend.pos = "none") +
ylim(-1, 5)  +
annotate(geom="text", x= 1.5, y= 3.5,  label = "k(N) = 119(9)", size = 6) +
annotate(geom="text", x= 1.5, y= -0.6,  label = "LDMC", size = 6) +
  scale_fill_manual(values = "#F0E442") + scale_colour_manual(values = "#F0E442") +
   theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
          axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank())

LDMCo_ab 


SLAo_ab <- orchard_plot(SLA_ab_ma, group = "ACC", data = SLA_ab, xlab = "Standardised mean difference", transfm = "none", angle = 45, k =FALSE, legend.pos = "none") +
ylim(-1, 5)  +
annotate(geom="text", x= 1.5, y= 3.5,  label = "k(N) = 251(17)", size = 6) +
annotate(geom="text", x= 1.5, y= -0.6,  label = "SLA", size = 6) +
  scale_fill_manual(values = "#0072B2") + scale_colour_manual(values = "#0072B2") +
   theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
          axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank())

SLAo_ab 

terpeneso_ab <- orchard_plot(terp_ab_ma, group = "ACC", data = terpenes_ab, xlab = "Standardised mean difference", transfm = "none", angle = 45, k =FALSE, legend.pos = "none") +
ylim(-1, 5)  +
annotate(geom="text", x= 1.5, y= 3.5,  label = "k(N) = 24(6)", size = 6) +
annotate(geom="text", x= 1.5, y= -0.6,  label = "Terpenoids", size = 6) +
  scale_fill_manual(values = "#CC79A7") + scale_colour_manual(values = "#CC79A7") +
   theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
           axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank())

terpeneso_ab

phenolicso_ab <- orchard_plot(phenolics_ab_ma, group = "ACC", data = phenolics_ab, xlab = "Standardised mean difference", transfm = "none", angle = 45, k =FALSE, legend.pos = "none") +
ylim(-1, 5)  +
annotate(geom="text", x= 1.5, y= 3.5,  label = "k(N) = 228(13)", size = 6) +
annotate(geom="text", x= 1.5, y= -0.6,  label = "Phenolics", size = 6) +
  scale_fill_manual(values = "#D55E00" ) + scale_colour_manual(values = "#D55E00") +
   theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
          axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank())


phenolicso_ab 

No_ab <- orchard_plot(N_ab_ma, group = "ACC", data = N_ab, xlab = "Standardised mean difference", transfm = "none", angle = 45, k =FALSE, legend.pos = "none") +
ylim(-1, 5)  +
annotate(geom="text", x= 1.5, y= 3.5,  label = "k(N) = 206(27)", size = 6) +
annotate(geom="text", x= 1.5, y= -0.6,  label = "Nitrogen", size = 6) +
 scale_fill_manual(values = "#009E73" ) + scale_colour_manual(values = "#009E73") +
   theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
          axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank())

No_ab 

Co_ab <- orchard_plot(C_ab_ma, group = "ACC", data = C_ab, xlab = "Absolute standardised mean difference", transfm = "none", angle = 45, k =FALSE, legend.pos = "none") +
ylim(-1, 5)  +
annotate(geom="text", x= 1.5, y= 3.5,  label = "k(N) = 139(11)", size = 6) +
annotate(geom="text", x= 1.5, y= -0.6,  label = "Carbon", size = 6) +
  scale_fill_manual(values = "#999999" ) + scale_colour_manual(values = "#999999" ) +
   theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_blank())
Co_ab 



thicknesso_ab / toughnesso_ab / 
LDMCo_ab / SLAo_ab /
terpeneso_ab  / phenolicso_ab  /
No_ab / Co_ab 


##############################################################################

#                      M E T A - R E G R E S S I O N S                      #

##############################################################################


# meta regressions only performed for SLA, LDMC, C, N and phenolics
# terpenes, thickness and toughness had too few effect sizes to be considered
# SR, PD, density = continuous moderators
# experiment type, plant age, N-fixing neighbours = categorical moderators

# Meta regression performed separately for both directional and absolute effect sizes

### C A T E G O R I C A L  M O D E R A T O R S ###


#--------------------------------#
###    N-fixing neighbours     ###
#--------------------------------#

options(na.action = "na.omit")

# SLA - directional

SLA_N <- rma.mv (yi,vi,data=ef, subset=(trait_category=="SLA"), method="REML",
                 mods = ~ Nfixing,
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                 R = list(phylo = cor), control = list(optimizer = "optim"))

summary(SLA_N)

# SLA - absolute

SLA_N_ab <- rma.mv (yi,vi,data=ef_ab, subset=(trait_category=="SLA"), method="REML",
                 mods = ~ Nfixing-1,
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                 R = list(phylo = cor), control = list(optimizer = "optim"))

summary(SLA_N_ab)

# Phenolics - directional

phenolics_N <- rma.mv (yi,vi,data=ef, subset=(trait_category=="phenolics"),method="REML",
                 mods = ~ Nfixing-1,
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                 R = list(phylo = cor), control = list(optimizer = "optim"))

summary(phenolics_N)

# Phenolics - absolute

phenolics_N_ab <- rma.mv (yi,vi,data=ef_ab, subset=(trait_category=="phenolics"),method="REML",
                 mods = ~ Nfixing-1,
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                 R = list(phylo = cor), control = list(optimizer = "optim"))

summary(phenolics_N_ab)

# LDMC - directional

LDMC_N <- rma.mv (yi,vi,data=ef, subset=(trait_category=="LDMC"),method="REML",
                       mods = ~ Nfixing-1,
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                       R = list(phylo = cor), control = list(optimizer = "optim"))

summary(LDMC_N)

# LDMC - absolute

LDMC_N_ab <- rma.mv (yi,vi,data=ef_ab, subset=(trait_category=="LDMC"),method="REML",
                       mods = ~ Nfixing-1,
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                       R = list(phylo = cor), control = list(optimizer = "optim"))

summary(LDMC_N_ab)


# Nitrogen - directional (with phylogentic matrix argument)

N_N <- rma.mv (yi,vi,data=ef, subset=(trait_category=="N"),method="REML",
                       mods = ~ Nfixing-1,
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                       R = list(phylo = cor), control = list(optimizer = "optim"))

summary(N_N)


# Nitrogen - directional (without phylogenetic argument)  
  
N_N2 <- rma.mv (yi,vi,data=ef, subset=(trait_category=="N"),method="REML",
                       mods = ~ Nfixing-1,
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                        control = list(optimizer = "optim"))

summary(N_N2)

# Nitrogen - absolute

N_N_ab <- rma.mv (yi,vi,data=ef_ab, subset=(trait_category=="N"),method="REML",
                       mods = ~ Nfixing-1,
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                       R = list(phylo = cor), control = list(optimizer = "optim"))

summary(N_N_ab)

# carbon - directional

C_N <- rma.mv (yi,vi,data=ef, subset=(trait_category=="C"),method="REML",
                       mods = ~ Nfixing-1,
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                       R = list(phylo = cor), control = list(optimizer = "optim"))

summary(C_N)

# carbon - absolute

C_N_ab <- rma.mv (yi,vi,data=ef_ab, subset=(trait_category=="C"),method="REML",
                       mods = ~ Nfixing-1,
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                       R = list(phylo = cor), control = list(optimizer = "optim"))

summary(C_N_ab)

#--------------------------------#
###         Study type         ###
#--------------------------------#

# SLA - directional

SLA_ST<- rma.mv (yi,vi,data=ef, subset=(trait_category=="SLA"),method="REML",
                 mods = ~ study_type-1,
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                 R = list(phylo = cor), control = list(optimizer = "optim"))

summary(SLA_ST)

# SLA - absolute

SLA_ST_ab<- rma.mv (yi,vi,data=ef_ab, subset=(trait_category=="SLA"),method="REML",
                 mods = ~ study_type-1,
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                 R = list(phylo = cor), control = list(optimizer = "optim"))

summary(SLA_ST_ab)

# phenolics - directional

phenolics_ST <- rma.mv (yi,vi,data=ef, subset=(trait_category=="phenolics"),method="REML",
                       mods = ~ study_type-1,
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                       R = list(phylo = cor), control = list(optimizer = "optim"))

summary(phenolics_ST)

# phenolics - absolute

phenolics_ST_ab <- rma.mv (yi,vi,data=ef_ab, subset=(trait_category=="phenolics"),method="REML",
                       mods = ~ study_type-1,
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                       R = list(phylo = cor), control = list(optimizer = "optim"))

summary(phenolics_ST_ab)

# LDMC - directional

LDMC_ST <- rma.mv (yi,vi,data=ef, subset=(trait_category=="LDMC"),method="REML",
                  mods = ~ study_type-1,
                  random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                  R = list(phylo = cor), control = list(optimizer = "optim"))

summary(LDMC_ST)

# LDMC - absolute

LDMC_ST_ab <- rma.mv (yi,vi,data=ef_ab, subset=(trait_category=="LDMC"),method="REML",
                  mods = ~ study_type-1,
                  random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                  R = list(phylo = cor), control = list(optimizer = "optim"))

summary(LDMC_ST_ab)


# Nitrogen - directional

N_ST <- rma.mv (yi,vi,data=ef, subset=(trait_category=="N"),method="REML",
               mods = ~ study_type,
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"))

summary(N_ST)

# Nitrogen - absolute

N_ST_ab <- rma.mv (yi,vi,data=ef_ab, subset=(trait_category=="N"),method="REML",
               mods = ~ study_type-1,
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"))

summary(N_ST_ab)

# Carbon - directional

C_ST <- rma.mv (yi,vi,data=ef, subset=(trait_category=="C"),method="REML",
               mods = ~ study_type-1,
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"))

summary(C_ST)

# Carbon - absolute

C_ST_ab <- rma.mv (yi,vi,data=ef_ab, subset=(trait_category=="C"),method="REML",
               mods = ~ study_type-1,
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"))

summary(C_ST_ab)


#--------------------------------#
###      plant lifestage       ###
#--------------------------------#

# SLA - directional

SLA_LS <- rma.mv (yi,vi,data=ef, subset=(trait_category=="SLA"),method="REML",
                 mods = ~ plant_lifestage_2-1,
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                 R = list(phylo = cor), control = list(optimizer = "optim"))

summary(SLA_LS)

# SLA - absolute

SLA_LS_ab <- rma.mv (yi,vi,data=ef_ab, subset=(trait_category=="SLA"),method="REML",
                 mods = ~ plant_lifestage_2-1,
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                 R = list(phylo = cor), control = list(optimizer = "optim"))

summary(SLA_LS_ab)


# Phenolics - directional

phenolics_LS <- rma.mv (yi,vi,data=ef, subset=(trait_category=="phenolics"),method="REML",
                       mods = ~ plant_lifestage_2-1,
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                       R = list(phylo = cor), control = list(optimizer = "optim"))

summary(phenolics_LS)

# Phenolics - absolute

phenolics_LS_ab <- rma.mv (yi,vi,data=ef_ab, subset=(trait_category=="phenolics"),method="REML",
                       mods = ~ plant_lifestage_2-1,
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                       R = list(phylo = cor), control = list(optimizer = "optim"))

summary(phenolics_LS_ab)

# LDMC - directional

LDMC_LS <- rma.mv (yi,vi,data=ef, subset=(trait_category=="LDMC"),method="REML",
                  mods = ~ plant_lifestage_2-1,
                  random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                  R = list(phylo = cor), control = list(optimizer = "optim"))

summary(LDMC_LS)

# LDMC - absolute

LDMC_LS_ab <- rma.mv (yi,vi,data=ef_ab, subset=(trait_category=="LDMC"),method="REML",
                  mods = ~ plant_lifestage_2-1,
                  random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                  R = list(phylo = cor), control = list(optimizer = "optim"))

summary(LDMC_LS_ab)

# Nitrogen - directional

N_LS <- rma.mv (yi,vi,data=ef, subset=(trait_category=="N"),method="REML",
               mods = ~ plant_lifestage_2-1,
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"))

summary(N_LS)

# Nitrogen - absolute

N_LS_ab <- rma.mv (yi,vi,data=ef_ab, subset=(trait_category=="N"),method="REML",
               mods = ~ plant_lifestage_2-1,
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"))

summary(N_LS_ab)


# Carbon - directional

C_LS <- rma.mv (yi,vi,data=ef, subset=(trait_category=="C"),method="REML",
               mods = ~ plant_lifestage_2-1,
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"))

summary(C_LS)

# Carbon - absolute

C_LS_ab <- rma.mv (yi,vi,data=ef_ab, subset=(trait_category=="C"),method="REML",
               mods = ~ plant_lifestage_2-1,
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"))

summary(C_LS_ab)


### C O N T I N I O U S   M O D E R A T O R S ###

#--------------------------------#
###          density           ###
#--------------------------------#

# SLA - directional

# excluded experiments with extreme density values

SLA_density <- rma.mv(yi,vi,data=ef, subset=(trait_category=="SLA" & experiment != "IDENT Montreal"),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                 R = list(phylo = cor), mods = ~ density, control = list(optimizer = "optim"))

summary(SLA_density)

SLA_density_r <- regplot(SLA_density, mod="density", pi=TRUE, refline=0, 
                    xlab="trees per m2", ylab="SMD", bg= "#009e73", main = "SLA", ylim =c(-5,5))

# SLA - absolute

SLA_density_ab <- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="SLA" & experiment != "IDENT Montreal"),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                 R = list(phylo = cor), mods = ~ density, control = list(optimizer = "optim"))

summary(SLA_density_ab)

SLA_density_r_ab <- regplot(SLA_density_ab, mod="density", pi=TRUE, refline=0, 
                    xlab="trees per m2", ylab="SMD", bg= "#009e73", main = "SLA", ylim =c(-5,5))


# phenolics - directional 

phenolics_density <- rma.mv(yi,vi,data=ef, subset=(trait_category=="phenolics"),method="REML",
                      random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                      R = list(phylo = cor), mods = ~ density,  control = list(optimizer = "optim"))

summary(phenolics_density)

phenolics_density_r <- regplot(phenolics_density, mod="density", pi=TRUE, refline=0, 
                         xlab="trees per m2", ylab="SMD", bg= "#009e73", main = "phenolics", ylim =c(-5,5))

# phenolics - absolute 

phenolics_density_ab <- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="phenolics"),method="REML",
                      random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                      R = list(phylo = cor), mods = ~ density,  control = list(optimizer = "optim"))

summary(phenolics_density_ab)

phenolics_density_r_ab <- regplot(phenolics_density_ab, mod="density", pi=TRUE, refline=0, 
                         xlab="trees per m2", ylab="SMD", bg= "#009e73", main = "phenolics", ylim =c(-5,5))

# LDMC - directional

LDMC_density <- rma.mv(yi,vi,data=ef, subset=(trait_category=="LDMC"),method="REML",
                            random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                            R = list(phylo = cor), mods = ~ density, control = list(optimizer = "optim"))
summary(LDMC_density)

LDMC_density_r <- regplot(LDMC_density, mod="density", pi=TRUE, refline=0, 
                               xlab="trees per m2", ylab="SMD", bg= "#009e73", main = "LDMC", ylim =c(-5,5))

# LDMC - absolute

LDMC_density_ab <- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="LDMC"),method="REML",
                            random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                            R = list(phylo = cor), mods = ~ density, control = list(optimizer = "optim"))
summary(LDMC_density_ab)

LDMC_density_r_ab <- regplot(LDMC_density_ab, mod="density", pi=TRUE, refline=0, 
                               xlab="trees per m2", ylab="SMD", bg= "#009e73", main = "LDMC", ylim =c(-5,5))


# Nitrogen - directional

# excluded experiments with extreme density values

N_density <- rma.mv(yi,vi,data=ef, subset=(trait_category=="N" & experiment != "IDENT Cloquet"),method="REML",
                            random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                            R = list(phylo = cor), mods = ~ density, control = list(optimizer = "optim"))

summary(N_density)

N_density_r <- regplot(N_density, mod="density", pi=TRUE, refline=0, 
                               xlab="trees per m2", ylab="SMD", bg= "#009e73", main = "N", ylim =c(-5,5))

# Nitrogen - absolute

N_density_ab <- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="N" & experiment != "IDENT Cloquet"),method="REML",
                            random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                            R = list(phylo = cor), mods = ~ density, control = list(optimizer = "optim"))

summary(N_density_ab)

N_density_r_ab <- regplot(N_density_ab, mod="density", pi=TRUE, refline=0, 
                               xlab="trees per m2", ylab="SMD", bg= "#009e73", main = "N", ylim =c(-5,5))


# Carbon - directional

# excluded experiments with extreme density values

C_density <- rma.mv(yi,vi,data=ef, subset=(trait_category=="C" & experiment != "IDENT Cloquet"),method="REML",
                            random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                            R = list(phylo = cor), mods = ~ density, control = list(optimizer = "optim"))
summary(C_density)

C_density_r <- regplot(C_density, mod="density", pi=TRUE, refline=0, 
                               xlab="trees per m2", ylab="SMD", bg= "#009e73", main = "C", ylim =c(-5,5))

# Carbon - absolute
  
C_density_ab <- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="C" & experiment != "IDENT Cloquet"),method="REML",
                            random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                            R = list(phylo = cor), mods = ~ density, control = list(optimizer = "optim"))
summary(C_density_ab)

C_density_r_ab <- regplot(C_density_ab, mod="density", pi=TRUE, refline=0, 
                               xlab="trees per m2", ylab="SMD", bg= "#009e73", main = "C", ylim =c(-5,5))


#--------------------------------#
###      species richness      ###
#--------------------------------#

# SLA - directional

SLA_SR <- rma.mv(yi,vi,data=ef, subset=(trait_category=="SLA"),method="REML",
                   random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                   R = list(phylo = cor), mods = ~ SR, control = list(optimizer = "optim"))

summary(SLA_SR)

SLA_SR_r <- regplot(SLA_SR, mod="SR", pi=TRUE, refline=0, labsize=0., 
                    xlab="Species richness", ylab="Standardised Mean Difference", bg= "#009e73", main = "SLA", ylim =c(-5,5)) 

# SLA - absolute

SLA_SR_ab <- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="SLA"),method="REML",
                   random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                   R = list(phylo = cor), mods = ~ SR, control = list(optimizer = "optim"))

summary(SLA_SR_ab)

SLA_SR_r_ab <- regplot(SLA_SR, mod="SR", pi=TRUE, refline=0, labsize=0., 
                    xlab="Species richness", ylab="Standardised Mean Difference", bg= "#009e73", main = "SLA", ylim =c(-5,5)) 



# phenolics - directional

phenolics_SR <- rma.mv(yi,vi,data=ef, subset=(trait_category=="phenolics"),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                 R = list(phylo = cor), mods = ~ SR, control = list(optimizer = "optim"))

summary(phenolics_SR)

phenolics_SR_r <- regplot(phenolics_SR, mod="SR", pi=TRUE, refline=0, 
                    xlab="Species richness", ylab="SMD", bg= "#009e73", main = "phenolics", ylim =c(-5,5))

# phenolics - absolute

phenolics_SR_ab <- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="phenolics"),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                 R = list(phylo = cor), mods = ~ SR, control = list(optimizer = "optim"))

summary(phenolics_SR_ab)

phenolics_SR_r_ab <- regplot(phenolics_SR_ab, mod="SR", pi=TRUE, refline=0, 
                    xlab="Species richness", ylab="SMD", bg= "#009e73", main = "phenolics", ylim =c(-5,5))


# LDMC - directional

LDMC_SR <- rma.mv(yi,vi,data=ef, subset=(trait_category=="LDMC"),method="REML",
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                       R = list(phylo = cor), mods = ~ SR, control = list(optimizer = "optim"))

summary(LDMC_SR)

LDMC_SR_r <- regplot(LDMC_SR, mod="SR", pi=TRUE, refline=0, 
                          xlab="Species richness", ylab="SMD", bg= "#009e73", main = "LDMC", ylim =c(-5,5))

# LDMC - absolute

LDMC_SR_ab <- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="LDMC"),method="REML",
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                       R = list(phylo = cor), mods = ~ SR, control = list(optimizer = "optim"))

summary(LDMC_SR_ab)

LDMC_SR_r_ab <- regplot(LDMC_SR_ab, mod="SR", pi=TRUE, refline=0, 
                          xlab="Species richness", ylab="SMD", bg= "#009e73", main = "LDMC", ylim =c(-5,5))


# Nitrogen - directional

N_SR <- rma.mv(yi,vi,data=ef, subset=(trait_category=="N"),method="REML",
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                       R = list(phylo = cor), mods = ~ SR, control = list(optimizer = "optim"))

summary(N_SR)

N_SR_r <- regplot(N_SR, mod="SR", pi=TRUE, refline=0, 
                          xlab="Species richness", ylab="SMD", bg= "#009e73", main = "N", ylim =c(-5,5))

# Nitrogen - absolute

N_SR_ab <- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="N"),method="REML",
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                       R = list(phylo = cor), mods = ~ SR, control = list(optimizer = "optim"))

summary(N_SR_ab)

N_SR_r_ab <- regplot(N_SR_ab, mod="SR", pi=TRUE, refline=0, 
                          xlab="Species richness", ylab="SMD", bg= "#009e73", main = "N", ylim =c(-5,5))

# Carbon - directional

C_SR <- rma.mv(yi,vi,data=ef, subset=(trait_category=="C"),method="REML",
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                       R = list(phylo = cor), mods = ~ SR, control = list(optimizer = "optim"))

summary(C_SR)

C_SR_r <- regplot(C_SR, mod="SR", pi=TRUE, refline=0, 
                          xlab="Species richness", ylab="SMD", bg= "#009e73", main = "C", ylim =c(-5,5))
# Carbon - absolute

C_SR_ab <- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="C"),method="REML",
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                       R = list(phylo = cor), mods = ~ SR, control = list(optimizer = "optim"))

summary(C_SR_ab)

C_SR_r_ab <- regplot(C_SR_ab, mod="SR", pi=TRUE, refline=0, 
                          xlab="Species richness", ylab="SMD", bg= "#009e73", main = "C", ylim =c(-5,5))

#--------------------------------#
###   phylogentic diversity    ###
#--------------------------------#

# SLA - directional

SLA_PD2 <- rma.mv(yi,vi,data=ef, subset=(trait_category=="SLA"),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                 R = list(phylo = cor), mods = ~ PD2, control = list(optimizer = "optim"))

summary(SLA_PD2)

SLA_PD2_r <- regplot(SLA_PD2, mod="PD2", pi=TRUE, refline=0, 
                    xlab="Phylogenetic diversity", ylab="SMD", bg= "#009e73", main = "SLA", ylim =c(-5,5))
# SLA - absolute

SLA_PD2_ab <- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="SLA"),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                 R = list(phylo = cor), mods = ~ PD2, control = list(optimizer = "optim"))

summary(SLA_PD2_ab)

SLA_PD2_r_ab <- regplot(SLA_PD2_ab, mod="PD2", pi=TRUE, refline=0, 
                    xlab="Phylogenetic diversity", ylab="SMD", bg= "#009e73", main = "SLA", ylim =c(-5,5))

# phenolics - directional

phenolics_PD2 <- rma.mv(yi,vi,data=ef, subset=(trait_category=="phenolics"),method="REML",
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                       R = list(phylo = cor), mods = ~ PD2, control = list(optimizer = "optim"))

summary(phenolics_PD2)

phenolics_PD2_r <- regplot(phenolics_PD2, mod="PD2", pi=TRUE, refline=0, 
                          xlab="Phylogenetic diversity", ylab="SMD", bg= "#009e73", main = "phenolics", ylim =c(-5,5))

# phenolics - absolute

phenolics_PD2_ab <- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="phenolics"),method="REML",
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                       R = list(phylo = cor), mods = ~ PD2, control = list(optimizer = "optim"))

summary(phenolics_PD2_ab)

phenolics_PD2_r_ab <- regplot(phenolics_PD2_ab, mod="PD2", pi=TRUE, refline=0, 
                          xlab="Phylogenetic diversity", ylab="SMD", bg= "#009e73", main = "phenolics", ylim =c(-5,5))


# LDMC - directional

LDMC_PD2 <- rma.mv(yi,vi,data=ef, subset=(trait_category=="LDMC"),method="REML",
                  random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                  R = list(phylo = cor), mods = ~ PD2, control = list(optimizer = "optim"))

summary(LDMC_PD2)

LDMC_PD2_r <- regplot(LDMC_PD2, mod="PD2", pi=TRUE, refline=0, 
                     xlab="Phylogenetic diversity", ylab="SMD", bg= "#009e73", main = "LDMC", ylim =c(-5,5))

# LDMC - absolute

LDMC_PD2_ab <- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="LDMC"),method="REML",
                  random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                  R = list(phylo = cor), mods = ~ PD2, control = list(optimizer = "optim"))

summary(LDMC_PD2_ab)

LDMC_PD2_r_ab <- regplot(LDMC_PD2_ab, mod="PD2", pi=TRUE, refline=0, 
                     xlab="Phylogenetic diversity", ylab="SMD", bg= "#009e73", main = "LDMC", ylim =c(-5,5))


# Nitrogen - directional

N_PD2 <- rma.mv(yi,vi,data=ef, subset=(trait_category=="N" ) ,method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), mods = ~ PD2, control = list(optimizer = "optim"))

summary(N_PD2)

N_PD2_r <- regplot(N_PD2, mod="PD2", pi=TRUE, refline=0, 
                  xlab="Phylogenetic diversity", ylab="SMD", bg= "#009e73", main = "N", ylim =c(-5,5))

# Nitrogen - absolute

N_PD2_ab <- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="N" ) ,method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), mods = ~ PD2, control = list(optimizer = "optim"))

summary(N_PD2_ab)

N_PD2_r_ab <- regplot(N_PD2_ab, mod="PD2", pi=TRUE, refline=0, 
                  xlab="Phylogenetic diversity", ylab="SMD", bg= "#009e73", main = "N", ylim =c(-5,5))

# Carbon - directional

C_PD2 <- rma.mv(yi,vi,data=ef, subset=(trait_category=="C"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), mods = ~ PD2, control = list(optimizer = "optim"))

summary(C_PD2)

C_PD2_r <- regplot(C_PD2, mod="PD2", pi=TRUE, refline=0, 
                  xlab="Phylogenetic diversity", ylab="SMD", bg= "#009e73", main = "C", ylim =c(-5,5))

# Carbon - absolute

C_PD2_ab <- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="C"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), mods = ~ PD2, control = list(optimizer = "optim"))

summary(C_PD2_ab)

C_PD2_r_ab <- regplot(C_PD2_ab, mod="PD2", pi=TRUE, refline=0, 
                  xlab="Phylogenetic diversity", ylab="SMD", bg= "#009e73", main = "C", ylim =c(-5,5))


##############################################################################

#                             SENSITIVITY ANALYSIS                           #
 
##############################################################################


# What species are over-represented ?
  
# C

 # 15.2 % Betula pendula
 # 9.4 % Fagus sylvetcia
 # 5.1 % Quercus robur
 # All other species < 5 %

# LDMC

 # 10.1 % Betula pendula
 # 8.4 % Fagus sylvetica
 # 5 % Quercus robur
 # All other species < 5 % 

# N

 # 10.1 % Betula pendula
 # 6.2 % Fagus sylvetica
 # All other species < 5 % 


# phenolics

 # 67.1 % Betula pendula
 # 7.1 % Fagus slyvetica

# SLA

 # 8.7 % Fagus sylvetica
 # 7.1 % Betula pendula


# Terpenes 

 # 24 % Plantago lancelota
 # 24 % Pinus halpensis
 # 24 % Rosmarinus officentalis
 # 12 % Cistus albidus


# Thickness

 # 80 % Betula Pendula
 # 15 % Quercus robur
 # 5 % Pinus pinaster 

# Toughness

 # 80 % Betula pendula
 # 15 % Quercus robur
 # 5 % Shorea leprosula



##### Re-analysis without Betula pendula #####

# Bpn = Betula pendula no

#--------------------------#
###  C - Betula Pendula  ###
#--------------------------#

## meta-analysis ##

# directional

C_Bpn<- rma.mv(yi,vi,data=ef, subset=(trait_category=="C" & species != "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"))

summary(C_Bpn)

# absolute

C_Bpn_ab <- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="C" & species != "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"))

summary(C_Bpn_ab)

## meta-regressions ##


# Phylogenetic diversity

# directional

C_Bpn_PD2<- rma.mv(yi,vi,data=ef, subset=(trait_category=="C" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ PD2)

summary(C_Bpn_PD2)

# absolute 

C_Bpn_PD2_ab<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="C" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ PD2)

summary(C_Bpn_PD2_ab)


# species richness

# directional

C_Bpn_SR<- rma.mv(yi,vi,data=ef, subset=(trait_category=="C" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ SR)

summary(C_Bpn_SR)

# absolute

C_Bpn_SR_ab<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="C" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ SR)

summary(C_Bpn_SR_ab)


# density

# directional

C_Bpn_density<- rma.mv(yi,vi,data=ef, subset=(trait_category=="C" & species != "Betula pendula" & experiment != "IDENT Cloquet"),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ density)

summary(C_Bpn_density)

# absolute  

C_Bpn_density_ab<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="C" & species != "Betula pendula" & experiment != "IDENT Cloquet"),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ density)

summary(C_Bpn_density_ab)

# study type

# directional

C_Bpn_ST<- rma.mv(yi,vi,data=ef, subset=(trait_category=="C" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ study_type-1)

summary(C_Bpn_ST)

# absolute

C_Bpn_ST_ab<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="C" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ study_type-1)

summary(C_Bpn_ST_ab)

# Nfixing

# directional

C_Bpn_N<- rma.mv(yi,vi,data=ef, subset=(trait_category=="C" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ Nfixing-1)

summary(C_Bpn_N)

# absolute 

C_Bpn_N_ab<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="C" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ Nfixing-1)

summary(C_Bpn_N_ab)


# ontogeny

# directional

C_Bpn_LS<- rma.mv(yi,vi,data=ef, subset=(trait_category=="C" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ plant_lifestage_2-1)

summary(C_Bpn_LS)

# absolute

C_Bpn_LS_ab<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="C" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ plant_lifestage_2-1)

summary(C_Bpn_LS_ab)


#-----------------------------#
### LDMC - Betula pendula ###
#-----------------------------#

## meta-analysis ##

# directional

LDMC_Bpn<- rma.mv(yi,vi,data=ef, subset=(trait_category=="LDMC" & species != "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"))

summary(LDMC_Bpn)

# absolute

LDMC_Bpn_ab <- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="LDMC" & species != "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"))

summary(LDMC_Bpn_ab)

## meta-regressions ##


# Phylogenetic diversity

# directional

LDMC_Bpn_PD2<- rma.mv(yi,vi,data=ef, subset=(trait_category=="LDMC" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ PD2)

summary(LDMC_Bpn_PD2)

# absolute 

LDMC_Bpn_PD2_ab<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="LDMC" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ PD2)

summary(LDMC_Bpn_PD2_ab)


# species richness

# directional

LDMC_Bpn_SR<- rma.mv(yi,vi,data=ef, subset=(trait_category=="LDMC" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ SR)

summary(LDMC_Bpn_SR)

# absolute

LDMC_Bpn_SR_ab<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="LDMC" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ SR)

summary(LDMC_Bpn_SR_ab)


# density

# directional

LDMC_Bpn_density<- rma.mv(yi,vi,data=ef, subset=(trait_category=="LDMC" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ density)

summary(LDMC_Bpn_density)

# absolute  

LDMC_Bpn_density_ab<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="LDMC" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ density)

summary(LDMC_Bpn_density_ab)

# study type

# directional

LDMC_Bpn_ST<- rma.mv(yi,vi,data=ef, subset=(trait_category=="LDMC" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ study_type-1)

summary(LDMC_Bpn_ST)

# absolute

LDMC_Bpn_ST_ab<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="LDMC" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ study_type-1)

summary(LDMC_Bpn_ST_ab)

# Nfixing

# directional

LDMC_Bpn_N<- rma.mv(yi,vi,data=ef, subset=(trait_category=="LDMC" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ Nfixing-1)

summary(LDMC_Bpn_N)

# absolute 

LDMC_Bpn_N_ab<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="LDMC" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ Nfixing-1)

summary(LDMC_Bpn_N_ab)


# ontogeny

# directional

LDMC_Bpn_LS<- rma.mv(yi,vi,data=ef, subset=(trait_category=="LDMC" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ plant_lifestage_2-1)

summary(LDMC_Bpn_LS)

# absolute

LDMC_Bpn_LS_ab<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="LDMC" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ plant_lifestage_2-1)

summary(LDMC_Bpn_LS_ab)

#--------------------------#
### N - Betula Pendula ###
#--------------------------#

## meta-analysis ##

# directional

N_Bpn<- rma.mv(yi,vi,data=ef, subset=(trait_category=="N" & species != "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species , ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"))

summary(N_Bpn)

# absolute

N_Bpn_ab <- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="N" & species != "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"))

summary(N_Bpn_ab)

## meta-regressions ##


# Phylogenetic diversity

# directional

N_Bpn_PD2<- rma.mv(yi,vi,data=ef, subset=(trait_category=="N" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ PD2)

summary(N_Bpn_PD2)

# absolute 

N_Bpn_PD2_ab<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="N" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ PD2)

summary(N_Bpn_PD2_ab)


# species richness

# directional

N_Bpn_SR<- rma.mv(yi,vi,data=ef, subset=(trait_category=="N" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ SR)

summary(N_Bpn_SR)

# absolute

N_Bpn_SR_ab<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="N" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ SR)

summary(N_Bpn_SR_ab)


# density

# directional

N_Bpn_density<- rma.mv(yi,vi,data=ef, subset=(trait_category=="N" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ density)

summary(N_Bpn_density)

# absolute  

N_Bpn_density_ab<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="N" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ density)

summary(N_Bpn_density_ab)

# study type

# directional

N_Bpn_ST<- rma.mv(yi,vi,data=ef, subset=(trait_category=="N" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ study_type-1)

summary(N_Bpn_ST)

orchard_plot(N_Bpn_ST,  mod = "study_type", group = "ACC", data = ef, xlab = "Standardised mean difference", transfm = "none", angle = 45)


# absolute

N_Bpn_ST_ab<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="N" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ study_type-1)

summary(N_Bpn_ST_ab)

# Nfixing

# directional

N_Bpn_N<- rma.mv(yi,vi,data=ef, subset=(trait_category=="N" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ Nfixing-1)

summary(N_Bpn_N)

# absolute 

N_Bpn_N_ab<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="N" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ Nfixing-1)

summary(N_Bpn_N_ab)


# ontogeny

# directional

N_Bpn_LS<- rma.mv(yi,vi,data=ef, subset=(trait_category=="N" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ plant_lifestage_2-1)

summary(N_Bpn_LS)

# absolute

N_Bpn_LS_ab<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="N" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ plant_lifestage_2-1)

summary(N_Bpn_LS_ab)

#----------------------------------#
###  phenolics - Betula pendula  ###
#----------------------------------#

## meta-analysis ##

# directional

phenolics_Bpn<- rma.mv(yi,vi,data=ef, subset=(trait_category=="phenolics" & species != "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"))

summary(phenolics_Bpn)

# absolute

phenolics_Bpn_ab <- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="phenolics" & species != "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"))

summary(phenolics_Bpn_ab)

## meta-regressions ##


# Phylogenetic diversity
  
# directional

phenolics_Bpn_PD2<- rma.mv(yi,vi,data=ef, subset=(trait_category=="phenolics" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ PD2)

summary(phenolics_Bpn_PD2)

# absolute 

phenolics_Bpn_PD2_ab<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="phenolics" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ PD2)

summary(phenolics_Bpn_PD2_ab)


# species richness

# directional

phenolics_Bpn_SR<- rma.mv(yi,vi,data=ef, subset=(trait_category=="phenolics" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ SR)

summary(phenolics_Bpn_SR)

# absolute

phenolics_Bpn_SR_ab<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="phenolics" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ SR)

summary(phenolics_Bpn_SR_ab)


# density

# directional

phenolics_Bpn_density<- rma.mv(yi,vi,data=ef, subset=(trait_category=="phenolics" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ density)

summary(phenolics_Bpn_density)

# absolute  

phenolics_Bpn_density_ab<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="phenolics" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ density)

summary(phenolics_Bpn_density_ab)

# study type

# directional

phenolics_Bpn_ST<- rma.mv(yi,vi,data=ef, subset=(trait_category=="phenolics" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ study_type-1)

summary(phenolics_Bpn_ST)

# absolute

phenolics_Bpn_ST_ab<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="phenolics" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ study_type-1)

summary(phenolics_Bpn_ST_ab)

orchard_plot(phenolics_Bpn_ST_ab,  mod = "study_type", group = "ACC", data = ef, xlab = "Standardised mean difference", transfm = "none", angle = 45)


# Nfixing

# directional

phenolics_Bpn_N<- rma.mv(yi,vi,data=ef, subset=(trait_category=="phenolics" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ Nfixing-1)

summary(phenolics_Bpn_N)

# absolute 

phenolics_Bpn_N_ab<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="phenolics" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ Nfixing-1)

summary(phenolics_Bpn_N_ab)


# ontogeny

# directional

phenolics_Bpn_LS<- rma.mv(yi,vi,data=ef, subset=(trait_category=="phenolics" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ plant_lifestage_2-1)

summary(phenolics_Bpn_LS)

# absolute

phenolics_Bpn_LS_ab<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="phenolics" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
               R = list(phylo = cor), control = list(optimizer = "optim"),  mods = ~ plant_lifestage_2-1)

summary(phenolics_Bpn_LS_ab)

#----------------------------------#
### Thickness - Betula pendula ###
#----------------------------------#

# directional

thick_Bpn<- rma.mv(yi,vi,data=ef, subset=(trait_category=="thickness" & species != "Betula pendula"),method="REML",
            random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
             R = list(phylo = cor),  control = list(optimizer = "optim"))

summary(thick_Bpn)

# absolute

thick_Bpn_ab<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="thickness" & species != "Betula pendula"),method="REML",
            random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
             R = list(phylo = cor),  control = list(optimizer = "optim"))

summary(thick_Bpn_ab)

#----------------------------------#
### Toughness - Betula pendula ###
#----------------------------------#

# directional

tough_Bpn<- rma.mv(yi,vi,data=ef, subset=(trait_category=="toughness" & species != "Betula pendula"),method="REML",
            random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
             R = list(phylo = cor),  control = list(optimizer = "optim"))

summary(tough_Bpn)

# absolute

tough_Bpn_ab<- rma.mv(yi,vi,data=ef_ab, subset=(trait_category=="toughness" & species != "Betula pendula"),method="REML",
            random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
             R = list(phylo = cor),  control = list(optimizer = "optim"))

summary(tough_Bpn_ab)


#-----------------------#
###    funnel plots   ###
#-----------------------#


f_thick <- funnel(thick_ma, xlim = c(-5, 5), ylim= c(2, 0), main = "Thickness")

f_tough <- funnel(tough_ma, xlim = c(-5, 5), ylim= c(2, 0), main = "Toughness")

f_LDMC <- funnel(LDMC_ma, xlim = c(-5, 5), ylim= c(2, 0), main = "LDMC")

f_SLA <- funnel(SLA_ma, xlim = c(-5, 5), ylim= c(2, 0), main = "SLA")

f_terp <- funnel(terp_ma, xlim = c(-5, 5), ylim= c(2, 0), main = "Terpenoids")

f_phenolics <- funnel(phenolics_ma, xlim = c(-5, 5), ylim= c(2, 0), main = "Phenolics")

f_N <- funnel(N_ma, xlim = c(-5, 5), ylim= c(2, 0), main = "Nitrogen")

f_C <- funnel(C_ma, xlim = c(-5, 5), ylim= c(2, 0), main = "Carbon")



###### Publication bias tests #####


# 1) sample size bias 
#run meta-regressions with standard error as moderator 

ef$sei <- sqrt(ef$vi) # calculate sampling error

# SLA

# removed regression data as this had very large variance

SLA_pb1 <- rma.mv(yi,vi,data=ef, subset=(trait_category=="SLA" & ACC != "18" & ACC != "57" & EF != "ef877"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"),
             mods= sei)

summary(SLA_pb1)

SLA_pb1_r <- regplot(SLA_pb1, pi=TRUE, refline=0, 
                    xlab="sampling error", ylab="SMD", bg= "#009e73", main = "SLA")

trimfill(SLA_pb1)

# phenolics

phenolics_pb1 <- rma.mv(yi,vi,data=ef, subset=(trait_category=="phenolics"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"),
             mods = ~ sei)


summary(phenolics_pb1)


phenolics_pb1_r <- regplot(phenolics_pb1, pi=TRUE, refline=0, 
                    xlab="sampling error", ylab="SMD", bg= "#009e73", main = "Phenolics")


# LDMC

LDMC_pb1 <- rma.mv(yi,vi,data=ef, subset=(trait_category=="LDMC"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"),
             mods = ~ sei)


summary(LDMC_pb1)


LDMC_pb1_r <- regplot(LDMC_pb1, pi=TRUE, refline=0, 
                    xlab="sampling error", ylab="SMD", bg= "#009e73", main = "LDMC")


# N

N_pb1 <- rma.mv(yi,vi,data=ef, subset=(trait_category=="N" & ACC != c("18", "57") & EF != c("ef1215", "ef901",	"ef1235")), method="REML", 
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"),
             mods = ~ sei)

#here i have remvoved studies 18 and 57 which had massive variences 

N_pb1 <- rma.mv(yi,vi,data=ef, subset=(trait_category=="N" & sei < 2), method="REML", 
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"),
             mods = ~ sei)

# removed outlier here with argument  # sei < 2 #

summary(N_pb1)

N_pb1_r <- regplot(N_pb1, mod ="sei", pi=TRUE, refline=0, 
    xlab="sampling error", ylab="SMD", bg= "#009e73", main = "Nitrogen")
  
  
# C
  

C_pb1 <- rma.mv(yi,vi,data=ef, subset=(trait_category=="C" & sei < 2),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"),
             mods = ~ sei)


summary(C_pb1)


C_pb1_r <- regplot(C_pb1, pi=TRUE, refline=0, 
                    xlab="sampling error", ylab="SMD", bg= "#009e73", main = "Carbon")


# terpenes


terpenes_pb1 <- rma.mv(yi,vi,data=ef, subset=(trait_category=="terpenes"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"),
             mods = ~ sei)


summary(terpenes_pb1)


terpenes_pb1_r <- regplot(terpenes_pb1, pi=TRUE, refline=0, 
                    xlab="sampling error", ylab="SMD", bg= "#009e73", main = "terpenes")


# toughness



toughness_pb1 <- rma.mv(yi,vi,data=ef, subset=(trait_category=="toughness"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"),
             mods = ~ sei)


summary(toughness_pb1)


toughness_pb1_r <- regplot(toughness_pb1, pi=TRUE, refline=0, 
                    xlab="sampling error", ylab="SMD", bg= "#009e73", main = "toughness")


# thickness



thickness_pb1 <- rma.mv(yi,vi,data=ef, subset=(trait_category=="thickness"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"),
             mods = ~ sei)


summary(thickness_pb1)


thickness_pb1_r <- regplot(thickness_pb1, pi=TRUE, refline=0, 
                    xlab="sampling error", ylab="SMD", bg= "#009e73", main = "thickness")

## no convergance 



# 2) publication year bias

# run meta-regressions with year as moderator 


# SLA


SLA_pb2 <- rma.mv(yi,vi,data=ef, subset=(trait_category=="SLA"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"),
             mods= ~ year)

summary(SLA_pb2)

SLA_pb2_r <- regplot(SLA_pb2, pi=TRUE, refline=0, 
                    xlab="year", ylab="SMD", bg= "#E69F00", main = "SLA")


# phenolics

phenolics_pb2 <- rma.mv(yi,vi,data=ef, subset=(trait_category=="phenolics"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"),
             mods = ~ year)


summary(phenolics_pb2)


phenolics_pb2_r <- regplot(phenolics_pb2, pi=TRUE, refline=0, 
                    xlab="year", ylab="SMD", bg= "#E69F00", main = "Phenolics")


# LDMC

LDMC_pb2 <- rma.mv(yi,vi,data=ef, subset=(trait_category=="LDMC"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"),
             mods = ~ year)


summary(LDMC_pb2)


LDMC_pb2_r <- regplot(LDMC_pb2, pi=TRUE, refline=0, 
                    xlab="year", ylab="SMD", bg= "#E69F00", main = "LDMC")


# Nitrogen


N_pb2 <- rma.mv(yi,vi,data=ef, subset=(trait_category=="N"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"),
             mods = ~ year)


summary(N_pb2)

N_pb2_r <- regplot(N_pb2, pi=TRUE, refline=0, 
                    xlab="year", ylab="SMD", bg= "#E69F00", main = "Nitrogen")
  
  
# Carbon
  

C_pb2 <- rma.mv(yi,vi,data=ef, subset=(trait_category=="C"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"),
             mods = ~ year)


summary(C_pb2)


C_pb2_r <- regplot(C_pb2, pi=TRUE, refline=0, 
                    xlab="year", ylab="SMD", bg= "#E69F00", main = "Carbon")


# terpenes

terpenes_pb2 <- rma.mv(yi,vi,data=ef, subset=(trait_category=="terpenes"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"),
             mods = ~ year)


summary(terpenes_pb2)


terpenes_pb2_r <- regplot(terpenes_pb2, pi=TRUE, refline=0, 
                    xlab="year", ylab="SMD", bg= "#E69F00", main = "terpenes")


# toughness


oughness_pb2 <- rma.mv(yi,vi,data=ef, subset=(trait_category=="toughness"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"),
             mods = ~ year)


summary(toughness_pb2)


toughness_pb2_r <- regplot(toughness_pb2, pi=TRUE, refline=0, 
                    xlab="year", ylab="SMD", bg= "#E69F00", main = "toughness")


# thickness

thickness_pb2 <- rma.mv(yi,vi,data=ef, subset=(trait_category=="thickness"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
              R = list(phylo = cor), control = list(optimizer = "optim"),
             mods = ~ year)


summary(thickness_pb2)


thickness_pb2_r <- regplot(thickness_pb2, pi=TRUE, refline=0, 
                    xlab="year", ylab="SMD", bg= "#E69F00", main = "thickness")


### Additional orchard plots for figures ###

# Fig 2.

phenolics_type <- rma.mv(yi,vi,data=phenolics, method="REML",
                         random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                         R = list(phylo = cor), control = list(optimizer = "optim"),
                         mods = ~ defence_type_specific-1)

phenolics_type


colours2 = c("#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00")

orchard_plot(phenolics_type,  mod = "defence_type_specific", group = "EF", data = phenolics, xlab = "Standardised mean difference", transfm = "none", angle = 0, k =FALSE)+
  scale_fill_manual(values = colours2) + 
  scale_colour_manual(values = colours2) +
    theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))


# Fig. 3

SLA_SR <- rma.mv(yi,vi,data=ef, subset=(trait_category=="SLA"),method="REML",
                   random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phtthttp://127.0.0.1:11457/graphics/plot_zoom_png?width=1030&height=616p://127.0.0.1:11457/graphics/plot_zoom_png?width=1030&height=616hylo, ~1 | EF),
                   R = list(phylo = cor), mods = ~ SR, control = list(optimizer = "optim"))

SLA_SR_r <- regplot(SLA_SR, mod="SR", pi=TRUE, refline=0, labsize=0., 
                    xlab="Species richness", ylab="Standardised Mean Difference", bg= "#009e73", main = "SLA", ylim =c(-5,5)) 


# Fig. 4


phenolics_N <- rma.mv (yi,vi,data=ef, subset=(trait_category=="phenolics"),method="REML",
                 mods = ~ Nfixing-1,
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                 R = list(phylo = cor), control = list(optimizer = "optim"))

colours2 = c("#009E73", "#F0E442",)

phenolics_N_o <- orchard_plot(phenolics_N,  mod = "Nfixing", group = "EF", data = phenolics, xlab = "Standardised mean difference", transfm = "none", angle = 0, k =FALSE)+
  ylim(-5, 5) +
  scale_fill_manual(values = colours2) + 
  scale_colour_manual(values = colours2) +
  annotate(geom="text", x= 2.5, y= -3.5,  label = "Phenolics", size = 8) +
    theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))

phenolics_N_o 


N_N <- rma.mv (yi,vi,data=ef, subset=(trait_category=="N" & study_type =="observational"),method="REML",
                 mods = ~ Nfixing-1,
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                 R = list(phylo = cor), control = list(optimizer = "optim"))

colours2 = c("#009E73", "#F0E442")

N_N_o <-  orchard_plot(N_N,  mod = "Nfixing", group = "EF", data = N, xlab = "Standardised mean difference", transfm = "none", angle = 0, k =FALSE)+
  ylim(-5, 5) +
  scale_fill_manual(values = colours2) + 
  scale_colour_manual(values = colours2) +
    annotate(geom="text", x= 2.5, y= -3.5,  label = "Nitrogen", size = 8) +
    theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_blank())
N_N_o 

phenolics_N_o + N_N_o 

# Fig. S5

N_N <- rma.mv (yi,vi,data=ef, subset=(trait_category=="N"),method="REML",
                       mods = ~ Nfixing-1,
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                       R = list(phylo = cor), control = list(optimizer = "optim"))

colours2 = c("#009E73", "#F0E442")

orchard_plot(N_N,  mod = "Nfixing", group = "EF", data = N, xlab = "Standardised mean difference", transfm = "none", angle = 0, k =FALSE)+
  scale_fill_manual(values = colours2) + 
  scale_colour_manual(values = colours2) +
    theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))


# Nitrogen (without phylogenetic argument)  
  
N_N2 <- rma.mv (yi,vi,data=ef, subset=(trait_category=="N"),method="REML",
                       mods = ~ Nfixing-1,
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | phylo, ~1 | EF),
                       control = list(optimizer = "optim"))

summary(N_N2)

orchard_plot(N_N2,  mod = "Nfixing", group = "EF", data = N, xlab = "Standardised mean difference", transfm = "none", angle = 0, k =FALSE)+
  scale_fill_manual(values = colours2) + 
  scale_colour_manual(values = colours2) +
    theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))



