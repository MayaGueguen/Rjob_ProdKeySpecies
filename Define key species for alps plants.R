################################################################################################################
# R script to define species that are potential key contributor to a given ecosystem service.                  #
#                                                                                                              #
# Step 1: The goal is to model a given ecosystem service (Y) according to Environmental (E) and                #
# Socio-Economic conditions (SE) and species richness (R); to check its relevance according to its             #
# explanatory power and to save its Akaike Information Criterion (AIC_M0) as a reference for the next step     #
#                                                                                                              #
# Step2: The goal is to identify species key for the studied ecosystem service (Y) adding each candidate       #
# species (presence-absence) as an additional explanatory variable to M0 to compute model M1 and its           #
# associated AIC (AIC_M1). Finally, a species is declared as a key potential contributor to the ecosystem      #
# service if ΔAIC (AIC_M0-AIC_M1) > 4 and if its partial effect is positive (positive coefficient in the model)#
#                                                                                                              #
# Here we apply this framework to define the "key plant species" for total plant biomass                         # 
# Code by Eva Maire (eva.maire@umontpellier.fr). Last update on 2018-06-26                                     #
# This code provides results of analysis used in Maire, E., Villéger, S., Graham, N., Hoey, A., Cinner, J.,    #
# Ferse, S., Aliaume, C, Booth, D., Feary, D., Kulbicki, M., Sandin, S., Vigliola, L. & Mouillot, D. (2018).   # 
# Community-wide scan identifies plant species associated to coral reef services globally. Proceedings  B.      #
#                                                                                                              #
################################################################################################################

rm(list=ls())

#loading required libraries
require(corrplot)
require(MuMIn)
require(lme4)
require(gam)
require(optimx)
library(foreach)
library(foreign)
library(ade4)
library(ggplot2)
library(ggthemes)
library(reshape2)

# setting working directory
my_path <- "/home/gueguen/Documents/_SCRIPTS/2018_09_NPP_Functio/2_Alps_version/" #  <= "Biomass" folder
setwd(my_path)

###################################################################################################################################
### LOADING AND PREPARING DATA
###################################################################################################################################

load("Data_for_Wilfried.RData")
dim(commplots) ## sites x species
dim(env.data) ## sites x env
dim(ss.trts) ## species x traits

## Add slope values and Chamois classification
suppData = read.dbf("PlotsAlps_slope_classifChamois.dbf")
env.data = unique(merge(env.data, suppData[, c("longL93fic","latL93fich","slope")], by = c("longL93fic","latL93fich")))
env.data = unique(merge(env.data, suppData[, c("longL93fic","latL93fich","chamois")], by = c("longL93fic","latL93fich")))

# table(env.data$Forests, env.data$lnus)
# table(env.data$Forest_trts, env.data$lnus)
# table(env.data$Forest_both, env.data$lnus)
# 
# table(env.data$lnus, env.data$chamois)

## --------------------------------------------------------------------------------------------------------------------------------
## (i) Loading complete dataset of environnemental drivers, productivity, etc
nrow(env.data) ## 11172 plots are available
rownames(env.data) = env.data$numchrono

## Compute annual mean precipitation
env.data$prec_mean = rowMeans(env.data[, c(paste0("prec_0", 1:9), "prec_10", "prec_11", "prec_12")])

## Keep only variables of interest
env.data = env.data[, c("Prec",
                        "Forest_both",
                        "Grassland_both",
                        "Mountain_belts" ,
                        "lnus",
                        "Elevation",
                        "slope",
                        "temp_mean",
                        "prec_mean",
                        "aridity",
                        "evi_mn",
                        "X_savi_mn",
                        "msavi_mn",
                        "ndmi_mn",
                        "ndvi_mn",
                        "NPP_Mean")]

## Remove NAs
env.data = na.exclude(env.data) ## 10948

## All the quantitative covariables are standardised. 
for (i in 1:5) env.data[, i] = as.factor(env.data[, i])
for (i in 6:ncol(env.data)) env.data[, i] = as.vector(scale(env.data[, i]))
summary(env.data)
boxplot(env.data, las = 2)

## Check for correlation between variables
corrplot(cor(env.data[, 6:ncol(env.data)]), method = "number"
         , type = "lower", diag = F, tl.cex = 0.8
         , order = "FPC")


## --------------------------------------------------------------------------------------------------------------------------------
## (ii) Loading the matrix which describes presence/absence of species for each site
## Transform into 0/1
presence_plant_species = commplots
presence_plant_species[which(presence_plant_species > 0)] = 1

## We excluded species present on less than 1% of the sites
limit <- round(0.01 * nrow(env.data), 0) ## 109
plant_sp <- colnames(presence_plant_species)

occurrence <- apply(presence_plant_species, 2, sum)
length(which(occurrence >= limit)) # we retained 444 plant species as candidates
required_threshold <- which(occurrence >= limit)

## Keeping only candidates species and their presence/absence 
candidates_species <- plant_sp[required_threshold]
presence_candidates <- presence_plant_species[,required_threshold]
dim(presence_candidates)

## --------------------------------------------------------------------------------------------------------------------------------
## Verifying that the 2 datasets have the same rows
env_data = env.data[which(rownames(env.data) %in% rownames(presence_candidates)), ]
presence_candidates = presence_candidates[which(rownames(presence_candidates) %in% rownames(env.data)), ]
length(is.element(rownames(env_data),rownames(presence_candidates))) # 10948 OK

## Reorder lines and compute Species Richness
env_data = env_data[order(rownames(env_data)), ]
presence_candidates = presence_candidates[rownames(env_data), ]
env_data$SR = rowSums(presence_candidates, na.rm = T)


# dim(presence_candidates) ## sites x species
# dim(env_data) ## sites x env
# dim(ss.trts) ## species x traits

## --------------------------------------------------------------------------------------------------------------------------------
# (iii) Loading functional traits of species
Functional_traits = ss.trts
summary(Functional_traits)
boxplot(Functional_traits, las = 2)

Functional_traits$SEEDM = log(Functional_traits$SEEDM)
Functional_traits$PL_VEG_H = log(Functional_traits$PL_VEG_H)
Functional_traits$PL_REPR_H = log(Functional_traits$PL_REPR_H)
Functional_traits$LDMC = log(Functional_traits$LDMC)
Functional_traits$SLA = log(Functional_traits$SLA)
boxplot(Functional_traits, las = 2)

dim(Functional_traits)
dim(na.exclude(Functional_traits))
dim(na.exclude(Functional_traits[, c("LDMC","SLA","PL_VEG_H","PL_REPR_H","SEEDM","DISP_VITTOZ","WOODY")])) ## rm GROWTHF


###################################################################################################################################
#DEFINING KEY SPECIES
###################################################################################################################################

# par(mfrow=c(2,3))
# plot(temp_mean ~ lnus, env_data)
# plot(prec_mean ~ lnus, env_data)
# plot(slope ~ lnus, env_data)
# plot(temp_mean ~ Mountain_belts, env_data)
# plot(prec_mean ~ Mountain_belts, env_data)
# plot(slope ~ Mountain_belts, env_data)

vari_prod = c("ndvi_mn", "X_savi_mn", "msavi_mn", "evi_mn", "ndmi_mn", "NPP_Mean")
vari_hab = c("Mountain_belts", "lnus")
vari_env = c("temp_mean", "prec_mean", "aridity")
vari_topo = "slope"
vari_SP = c("SR")
vari_all = c(vari_hab, vari_env, vari_topo, vari_SP)

comb_vari = unlist(sapply(1:length(vari_all), function(x)
  apply(combn(vari_all, x), 2, function(y)
    paste0(y, collapse = "."))))
comb_vari = comb_vari[-grep("prec_mean.aridity", comb_vari)]
comb_vari = comb_vari[-grep("temp_mean.aridity", comb_vari)]

comb_mod = expand.grid(prod = vari_prod, expl = comb_vari, stringsAsFactors = F)

## --------------------------------------------------------------------------------------------------------------------------------
## STEP 1: Computing initial plant biomass model (M0)
COMPAR = foreach(vp = comb_mod$prod, ve = comb_mod$expl, .combine = "rbind") %do%
{
  cat("\n Calculate model M0 for ", vp, " with variables ", ve, "...\n")

  M0 = as.formula(paste0(vp, " ~ ", gsub("[.]", " + ", ve)))
  mod = "glm"
  type_mod = strsplit(mod, "_")[[1]][1]
  supp_arg = ifelse(type_mod == "lmer"
                    , paste0(", control = lmerControl(optimizer = 'optimx', calc.derivs = FALSE,"
                             , "optCtrl = list(method = 'nlminb', starttests = FALSE, kkt = FALSE))")
                    , "")
  eval(parse(text = paste0("M0 = ", type_mod, "(M0, data = env_data", supp_arg, ")")))
  
  res = data.frame(PROD = vp, VAR = ve, TYPE_MOD = mod, AIC = AIC(M0), R2 = r.squaredGLMM(M0)[2])
  return(res)
}
COMPAR_melt = melt(COMPAR, id.vars = c("PROD", "VAR", "TYPE_MOD"))
COMPAR_melt$NB_VAR = sapply(COMPAR_melt$VAR, function(x) length(strsplit(as.character(x), "[.]")[[1]]))
COMPAR_melt$GROUP1 = "Mountain_belts"
COMPAR_melt$GROUP1[grep("^lnus", COMPAR_melt$VAR)] = "lnus"
COMPAR_melt$GROUP1[grep("Mountain_belts.lnus", COMPAR_melt$VAR)] = "Mountain_belts.lnus"
COMPAR_melt$GROUP2 = "no_slope"
COMPAR_melt$GROUP2[grep("slope", COMPAR_melt$VAR)] = "slope"
COMPAR_melt$GROUP3 = "no_SR"
COMPAR_melt$GROUP3[grep("SR", COMPAR_melt$VAR)] = "SR"
head(COMPAR_melt)

TAB_plot = COMPAR_melt[which(COMPAR_melt$variable == "AIC"),]
ggplot(TAB_plot, aes(x = PROD, y = value, group = VAR, color = GROUP3)) + #interaction(GROUP1,GROUP2))) +
  geom_line() +
  facet_grid(~ NB_VAR, scales = "free_y") +
  scale_color_discrete("") +
  theme_fivethirtyeight() +
  theme(axis.text.x = element_text(angle = 90))

ggplot(TAB_plot, aes(x = VAR, y = value)) +
  geom_boxplot() +
  theme_fivethirtyeight() +
  theme(axis.text.x = element_text(angle = 90))


# ## --------------------------------------------------------------------------------------------------------------------------------
# # STEP 1: Computing initial plant biomass model (M0)
# COMPAR = foreach(vari = vari_prod, .combine = "rbind") %do%
# {
#   M0 = as.formula(paste0(vari, " ~ Mountain_belts + lnus + slope + temp_mean + prec_mean"))
#   M0_glm = as.formula(paste0(vari, " ~ Mountain_belts + lnus + slope + temp_mean + prec_mean"))
#   M0_gam = as.formula(paste0(vari, " ~ Mountain_belts + lnus + s(slope,4) + s(temp_mean,4) + s(prec_mean,4)"))
#   M0_lmer_1 = as.formula(paste0(vari, " ~ Mountain_belts + lnus + slope + temp_mean + prec_mean + (1 | Mountain_belts) + (1 | lnus)"))
#   M0_lmer_2 = as.formula(paste0(vari, " ~ Mountain_belts + lnus + slope + temp_mean + prec_mean + (1 | Mountain_belts) + (temp_mean | lnus)"))
#   M0_lmer_3 = as.formula(paste0(vari, " ~ Mountain_belts + lnus + slope + temp_mean + prec_mean + (1 | Mountain_belts) + (temp_mean + prec_mean | lnus)"))
#   M0_lmer_4 = as.formula(paste0(vari, " ~ Mountain_belts + lnus + slope + temp_mean + prec_mean + (temp_mean | Mountain_belts) + (temp_mean | lnus)"))
#   M0_lmer_5 = as.formula(paste0(vari, " ~ Mountain_belts + lnus + slope + temp_mean + prec_mean + (temp_mean | Mountain_belts) + (temp_mean + prec_mean | lnus)"))
#   
#   mods = c("glm", "gam", paste0("lmer_", 1:5))
#   for (mod in mods)
#   {
#     type_mod = strsplit(mod, "_")[[1]][1]
#     eval(parse(text = paste0("M0_", mod, " = ", type_mod, "(M0_", mod, ", data = env_data"
#                              , ifelse(type_mod == "lmer", ", control = lmerControl(optimizer = 'optimx', calc.derivs = FALSE, optCtrl = list(method = 'nlminb', starttests = FALSE, kkt = FALSE))", "")
#                              , ")")))
#   }
#   
#   eval(parse(text = paste0("res = data.frame(VARI = vari"
#                            , paste0(", AIC_", mods, " = AIC(M0_", mods, "), R2_"
#                                     , mods, " = r.squaredGLMM(M0_", mods, ")[2]"
#                                     # , mods, " = M0_", mods, "$converged"
#                                     , collapse = "")
#                            , ")" )))
#   return(res)
# }
# COMPAR_melt = melt(COMPAR)
# COMPAR_melt = COMPAR_melt[grep("AIC", COMPAR_melt$variable), ]
# COMPAR_melt$VARI = factor(COMPAR_melt$VARI, c("ndvi_mn", "msavi_mn","X_savi_mn",
#                                               "evi_mn", "ndmi_mn", "NPP_Mean"))
# ggplot(COMPAR_melt, aes(x = VARI, y = value, group = variable, color = variable)) +
#   geom_line(lwd = 1) +
#   labs(x = "", y = "R squared") +
#   scale_color_discrete("", labels = c("AIC_glm" = "GLM : MountainBelts + LU + Slope + Temp + Precip"
#                                       , "AIC_gam" = "GAM : MountainBelts + LU + s(Slope,4) + s(Temp,4) + s(Precip,4)"
#                                       , "AIC_lmer_1" = "LMER : GLM + (1 | Mountain_belts) + (1 | LU)"
#                                       , "AIC_lmer_2" = "LMER : GLM + (1 | Mountain_belts) + (Temp | LU)"
#                                       , "AIC_lmer_3" = "LMER : GLM + (1 | Mountain_belts) + (Temp + Precip | LU)"
#                                       , "AIC_lmer_4" = "LMER : GLM + (Temp | Mountain_belts) + (Temp | LU)"
#                                       , "AIC_lmer_5" = "LMER : GLM + (Temp | Mountain_belts) + (Temp + Precip | LU)")) +
#   theme_fivethirtyeight()

## For a specific variable (i.e. NPP)
# anova(M0_lmer_1, M0_lmer_2) ## 2
# anova(M0_lmer_2, M0_lmer_3) ## 3
# anova(M0_lmer_2, M0_lmer_4) ## 4
# anova(M0_lmer_3, M0_lmer_5) ## 5
# anova(M0_lmer_4, M0_lmer_5) ## 5



###################################################################################################################################

COMPAR = foreach(vari = c("evi_mn",
                          "X_savi_mn",
                          "msavi_mn",
                          "ndmi_mn",
                          "ndvi_mn",
                          "NPP_Mean"), .combine = "rbind") %do%
                          {
                            M0_glm_1 = as.formula(paste0(vari, " ~ SR"))
                            M0_glm_2 = as.formula(paste0(vari, " ~ SR + Mountain_belts"))
                            M0_glm_3 = as.formula(paste0(vari, " ~ SR + lnus"))
                            M0_glm_4 = as.formula(paste0(vari, " ~ SR + Mountain_belts + lnus"))
                            M0_gam = as.formula(paste0(vari, " ~ s(SR,4)"))
                            M0_lmer_1 = as.formula(paste0(vari, " ~ SR + (1 | Mountain_belts)"))
                            M0_lmer_2 = as.formula(paste0(vari, " ~ SR + (1 | lnus)"))
                            M0_lmer_3 = as.formula(paste0(vari, " ~ SR + (1 | Mountain_belts) + (1 | lnus)"))
                            
                            
                            mods = c(paste0("glm_", 1:4), "gam", paste0("lmer_", 1:3))
                            for (mod in mods)
                            {
                              type_mod = strsplit(mod, "_")[[1]][1]
                              eval(parse(text = paste0("M0_", mod, " = ", type_mod, "(M0_", mod, ", data = env_data"
                                                       , ifelse(type_mod == "lmer", ", control = lmerControl(optimizer = 'optimx', calc.derivs = FALSE, optCtrl = list(method = 'nlminb', starttests = FALSE, kkt = FALSE))", "")
                                                       , ")")))
                            }
                            
                            eval(parse(text = paste0("res = data.frame(VARI = vari"
                                                     , paste0(", AIC_", mods, " = AIC(M0_", mods, "), R2_"
                                                              , mods, " = r.squaredGLMM(M0_", mods, ")[2]"
                                                              , collapse = "")
                                                     , ")" )))
                            return(res)
                          }

COMPAR_melt = melt(COMPAR)
COMPAR_melt = COMPAR_melt[grep("AIC", COMPAR_melt$variable), ]
COMPAR_melt$VARI = factor(COMPAR_melt$VARI, c("ndvi_mn", "msavi_mn","X_savi_mn",
                                              "evi_mn", "ndmi_mn", "NPP_Mean"))
ggplot(COMPAR_melt, aes(x = VARI, y = value, group = variable, color = variable)) +
  geom_line(lwd = 1) +
  labs(x = "", y = "AIC") +
  scale_color_discrete("", labels = c("AIC_glm_1" = "GLM : SR"
                                      , "AIC_glm_2" = "GLM : SR + MountainBelts"
                                      , "AIC_glm_3" = "GLM : SR + LandCover"
                                      , "AIC_glm_4" = "GLM : SR + MountainBelts + LandCover"
                                      , "AIC_gam" = "GAM : SR")) +
  theme_fivethirtyeight()

###################################################################################################################################


## STEP 1: Computing initial plant biomass model (M0)
mod_form = "glm(evi_mn ~ SR + Mountain_belts + lnus, data = env_data)"
mod_form = "glm(NPP_Mean ~ SR + Mountain_belts + lnus + slope + aridity, data = env_data)"
mod_form = paste0("lmer(evi_mn ~ SR + Mountain_belts + lnus + slope + temp_mean + prec_mean +"
                  , "(1 | Mountain_belts) + (1 | lnus)"
                  , ", data = env_data"
                  , ", control = lmerControl(optimizer = 'optimx', calc.derivs = FALSE"
                  , ", optCtrl = list(method = 'nlminb', starttests = FALSE, kkt = FALSE)))")

eval(parse(text = paste0("M0 = ", mod_form)))

## Checking for robustness of M0
AIC(M0)
r.squaredGLMM(M0)[2] ## conditional (variance explained by both fixed and random factors)

#Storing AIC of M0 that will be the threshold to define which species improve the prediction of plant biomass
AIC_M0 <- AIC(M0)

## STEP 2: TESTING THE EFFECT OF EACH plant SPECIES INDIVIDUALLY
##  For each candidate species, an alternative model (M1) is obtained by adding the presence/absence of this plant species to M0
##  AIC of M1 and the effect (positive or negative) of the candidate species are saved.

## Results for all species will be stored in the Test_candidates matrix
Test_candidates <- matrix(NA,length(candidates_species),3)
colnames(Test_candidates) <- c("occurrence","AIC_M1","coeff_sp")
rownames(Test_candidates) <- candidates_species

cat(" \n ==> CANDIDATE ")
for (k in 1:length(candidates_species)) {
  cat(" ", k)
  
  # Defining presence/absence of the species k
  candidate <- as.factor(presence_candidates[,k])
  Test_candidates[k,"occurrence"] <- length(which(candidate==1))
  
  data_with_candidate <- as.data.frame(cbind(env_data,candidate)) 
  
  ## Adding presence/absence of species k to M0 to obtain M1
  mod_form_M1 = sub("env_data", "data_with_candidate", mod_form)
  mod_form_M1 = sub("SR", "candidate + SR", mod_form_M1)
  eval(parse(text = paste0("M1 = ", mod_form_M1)))
  
  Test_candidates[k,"AIC_M1"] <- AIC(M1)
  if (length(grep("lmer", mod_form)) > 0)
  {
    Test_candidates[k,"coeff_sp"] <- unique(unlist(coef(M1)$lnus[grep("candidate", names(coef(M1)$lnus))]))
  } else
  {
    Test_candidates[k,"coeff_sp"] <- unique(coef(M1)[grep("candidate", names(coef(M1)))])
  }
} # end of k
Test_candidates <- as.data.frame(Test_candidates)
# save(Test_candidates, file = "Test_candidates_glm.RData")



###################################################################################################################################
#Determining species that improve the prediction of plant biomass as ΔAIC (AIC_M0 - AIC_M1) > 4 with a positive effect
Delta_AIC <- AIC_M0 - Test_candidates$AIC_M1
better_AIC <- ifelse(Delta_AIC >= 4, 1, 0)
positive_coef <- ifelse(Test_candidates$coeff_sp > 0, 1, 0)

key_species <- which(better_AIC == 1 & positive_coef == 1)

#Extracting matrix that combines "key species" and their performances in the model
Summary_key_species <- cbind(Test_candidates,Delta_AIC)[key_species,]
nrow(Summary_key_species) 
#26 plant species are significantly and positively related to plant biomass



###################################################################################################################################
# END OF SCRIPT
###################################################################################################################################
###################################################################################################################################

func_tab = na.exclude(Functional_traits[, c("SLA","PL_VEG_H","SEEDM")])
dim(func_tab)
TT = dudi.hillsmith(func_tab, scannf = F, nf = 4)
# s.corcircle(TT$c1)
# s.label(TT$li)

ii = c(1,3)
xy_0 = TT$li[, ii]
xy_1 = TT$li[which(rownames(TT$li) %in% rownames(Summary_key_species)), ii]
chu_0 = chull(xy_0)
chu_1 = chull(xy_1)
plot(xy_0)
abline(h = 0, v = 0)
lines(xy_0[c(chu_0,chu_0[1]),])
points(xy_1, col = "red")
lines(xy_1[c(chu_1,chu_1[1]),], col = "red")

