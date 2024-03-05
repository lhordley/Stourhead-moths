##########################
#### user: Lisbeth Hordley
#### date: September 2021
#### info: Stourhead moths: run models testing effect of habitat and surrounding landscape

library(MuMIn)
library(DHARMa)
library(lme4)
library(StatisticalModels)
library(glmmTMB)
library(dplyr)
library(ggeffects)
library(performance)

rm(list = ls())
options(scipen=999)

moth_hab_data <- read.csv("Data/Stourhead_moths_final.csv", header=TRUE)
land_cover <- read.csv("../Data/Land_cover_stourhead_2015.csv", header=TRUE)

# change orientation of land cover data
colnames(land_cover) <- c("site","scale", "prop_arable","prop_broadleaf","prop_conifer","prop_grassland")
land_cover <- reshape(land_cover, idvar = "site", timevar = "scale", direction = "wide")

# merge in land cover data
moth_hab_data <- merge(moth_hab_data, land_cover, by.x="gridref", by.y="site", all.x=TRUE)

# Make sure variables are in right format
moth_hab_data[moth_hab_data == "n/a"] <- NA
moth_hab_data$basal_area <- as.numeric(moth_hab_data$basal_area)
moth_hab_data$canopy_openess <- as.numeric(moth_hab_data$canopy_openess)

summary(moth_hab_data) # basal area 1 NA, canopy openness 11 NAs 
moth_hab_data <- na.omit(moth_hab_data) 
length(unique(moth_hab_data$plot)) # 82 plots

# Check correlations between explanatory variables
round(cor(moth_hab_data[,c("dbh_average","basal_area","per_broadleaf_canopy","complexity_score","canopy_openess",
                           "wind_speed",  "minimum_temperature","cloud_cover_.", "dis_to_edge", "dis_broadleaf", "altitude")]),3)
# no correlations above r>0.7 for continuous variables (only moon cycle isn't included - but this is removed due to high VIF values)


########## SURROUNDING LANDSCAPE MODELS
# This is to determine which parameters (land type) will be included in the habitat variable models
# Following Fuentes-Montemayor et al. 2012 paper methods

# Total richness

tot_rich_arable1 <- glmer(richness ~ prop_arable.500 + (1|subcompartment) + (1|year),
                          family=poisson(), data=moth_hab_data)
summary(tot_rich_arable1)
r.squaredGLMM(tot_rich_arable1) # 0.2485051   
tot_rich_arable2 <- glmer(richness ~ prop_arable.1500 + (1|subcompartment) + (1|year),
                          family=poisson(), data=moth_hab_data)
summary(tot_rich_arable2) 
r.squaredGLMM(tot_rich_arable2) # 0.002046054   
tot_rich_arable3 <- glmer(richness ~ prop_arable.3000 + (1|subcompartment) + (1|year),
                          family=poisson(), data=moth_hab_data)
summary(tot_rich_arable3)    
r.squaredGLMM(tot_rich_arable3) # 0.01033699  
#### ARABLE 500M #### 

tot_rich_broad1 <- glmer(richness ~ prop_broadleaf.500 + (1|subcompartment) + (1|year),
                         family=poisson(), data=moth_hab_data)
summary(tot_rich_broad1) 
r.squaredGLMM(tot_rich_broad1) # 0.02999598    
tot_rich_broad2 <- glmer(richness ~ prop_broadleaf.1500 + (1|subcompartment) + (1|year),
                         family=poisson(), data=moth_hab_data)
summary(tot_rich_broad2)
r.squaredGLMM(tot_rich_broad2) # 0.01131397    
tot_rich_broad3 <- glmer(richness ~ prop_broadleaf.3000 + (1|subcompartment) + (1|year),
                         family=poisson(), data=moth_hab_data)
summary(tot_rich_broad3) 
r.squaredGLMM(tot_rich_broad3) # 0.02333280   
#### BROADLEAF 500M #### 


tot_rich_conifer1 <- glmer(richness ~ prop_conifer.500 + (1|subcompartment) + (1|year),
                           family=poisson(), data=moth_hab_data)
summary(tot_rich_conifer1)   
r.squaredGLMM(tot_rich_conifer1) # 0.0001103585     
tot_rich_conifer2 <- glmer(richness ~ prop_conifer.1500 + (1|subcompartment) + (1|year),
                           family=poisson(), data=moth_hab_data)
summary(tot_rich_conifer2)       
r.squaredGLMM(tot_rich_conifer2) # 0.0001536378     
tot_rich_conifer3 <- glmer(richness ~ prop_conifer.3000 + (1|subcompartment) + (1|year),
                           family=poisson(), data=moth_hab_data)
summary(tot_rich_conifer3)         
r.squaredGLMM(tot_rich_conifer3) # 0.02628789    
#### CONIFER 3000M #### 

tot_rich_grass1 <- glmer(richness ~ prop_grassland.500 + (1|subcompartment) + (1|year),
                         family=poisson(), data=moth_hab_data)
summary(tot_rich_grass1)        
r.squaredGLMM(tot_rich_grass1) # 0.02611387      
tot_rich_grass2 <- glmer(richness ~ prop_grassland.1500 + (1|subcompartment) + (1|year),
                         family=poisson(), data=moth_hab_data)
summary(tot_rich_grass2)                   
r.squaredGLMM(tot_rich_grass2) # 0.0146671      
tot_rich_grass3 <- glmer(richness ~ prop_grassland.3000 + (1|subcompartment) + (1|year),
                         family=poisson(), data=moth_hab_data)
summary(tot_rich_grass3)          
r.squaredGLMM(tot_rich_grass3) # 0.002249187     
#### GRASSLAND 500M #### 

# save results
r2 <- c(r.squaredGLMM(tot_rich_arable1)[1,1], r.squaredGLMM(tot_rich_arable2)[1,1], r.squaredGLMM(tot_rich_arable3)[1,1],
        r.squaredGLMM(tot_rich_broad1)[1,1], r.squaredGLMM(tot_rich_broad2)[1,1], r.squaredGLMM(tot_rich_broad3)[1,1],
        r.squaredGLMM(tot_rich_conifer1)[1,1], r.squaredGLMM(tot_rich_conifer2)[1,1], r.squaredGLMM(tot_rich_conifer3)[1,1],
        r.squaredGLMM(tot_rich_grass1)[1,1], r.squaredGLMM(tot_rich_grass2)[1,1], r.squaredGLMM(tot_rich_grass3)[1,1])
land_cover <- c("Arable.500", "Arable.1500", "Arable.3000", "Broadleaf.500", "Broadleaf.1500", "Broadleaf.3000",
                "Conifer.500", "Conifer.1500", "Conifer.3000", "Grassland.500", "Grassland.1500", "Grassland.3000")
tot_rich_land <- data.frame(Response=rep("Total moth richness"), Land_cover=land_cover, R2=r2)
write.csv(tot_rich_land, file="Results/Total_richness_land_cover.csv", row.names=FALSE)

## Now run models with habitat variables only 
tot_rich_hab <- glmer(richness ~ scale(dbh_average) + scale(basal_area) +
                        scale(per_broadleaf_canopy) + scale(complexity_score) + 
                        scale(canopy_openess) + scale(minimum_temperature) + 
                        scale(cloud_cover_.) + scale(dis_to_edge) + 
                        scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ (1|year) + 
                        (1|subcompartment), data = moth_hab_data, family = poisson(),
                      glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                      na.action = "na.fail")
summary(tot_rich_hab)
check_collinearity(tot_rich_hab) # low correlations after removing moon cycle
testDispersion(tot_rich_hab) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = tot_rich_hab, plot = F)
plot(simulationOutput) ## no assumptions violated!

final_mod2 <- glmer(richness ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + 
                      scale(wind_speed) + scale(dis_broadleaf) + scale(altitude) + 
                      scale(prop_arable.500) + (1|year) + (1|subcompartment), 
                    data = moth_hab_data, family = poisson(),
                    glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                    na.action = "na.fail")
summary(final_mod2)
check_collinearity(final_mod2) # low correlations after removing moon cycle
testDispersion(final_mod2) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod2, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod3 <- glmer(richness ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + 
                      scale(wind_speed) + scale(dis_broadleaf) + scale(altitude) + 
                      scale(prop_broadleaf.500) + (1|year) + (1|subcompartment), 
                    data = moth_hab_data, family = poisson(),
                    glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                    na.action = "na.fail")
summary(final_mod3)
check_collinearity(final_mod3) # low correlations after removing moon cycle
testDispersion(final_mod3) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod3, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod4 <- glmer(richness ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + 
                      scale(wind_speed) + scale(dis_broadleaf) + scale(altitude) + 
                      scale(prop_conifer.3000) + (1|year) + (1|subcompartment), 
                    data = moth_hab_data, family = poisson(),
                    glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                    na.action = "na.fail")
summary(final_mod4)
check_collinearity(final_mod4) # low correlations after removing moon cycle
testDispersion(final_mod4) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod4, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod5 <- glmer(richness ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + 
                      scale(wind_speed) + scale(dis_broadleaf) + scale(altitude) + 
                      scale(prop_grassland.500) + (1|year) + (1|subcompartment), 
                    data = moth_hab_data, family = poisson(),
                    glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                    na.action = "na.fail")
summary(final_mod5)
check_collinearity(final_mod5) # low correlations after removing moon cycle
testDispersion(final_mod5) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod5, plot = F)
plot(simulationOutput) ## no assumptions violated

AIC(tot_rich_hab, final_mod2, final_mod3, final_mod4, final_mod5)
# No land cover variables reduce AIC by >2 - use habitat model
summary(tot_rich_hab)
# complexity score (positive) and canopy openness (negative)

## save all AIC values for supplementary table
AIC_values <- data.frame(model=c("habitat", "arable 500m", "broadleaf 500m", "conifer 3000m", "grassland 500m"),
                 AIC=AIC(tot_rich_hab, final_mod2, final_mod3, final_mod4, final_mod5)$AIC, stringsAsFactors=FALSE) 
write.csv(AIC_values, file="Results/Total_richness_AIC.csv", row.names=FALSE)

## save lowest AIC model output for supplementary table
tot_rich_hab_sum <- as.data.frame(coef(summary(tot_rich_hab))) 
tot_rich_hab_sum$parameters <- row.names(tot_rich_hab_sum)
row.names(tot_rich_hab_sum) <- 1:nrow(tot_rich_hab_sum)
tot_rich_hab_sum$model <- "habitat"
tot_rich_hab_sum$AIC <- AIC(tot_rich_hab)

write.csv(tot_rich_hab_sum, file="Results/Total_richness_habitat.csv", row.names=FALSE)

###################################################################################################################################

# Broadleaf richness

broad_rich_arable1 <- glmer(broadleaf_rich ~ prop_arable.500 + (1|subcompartment) + (1|year),
                            family=poisson(), data=moth_hab_data)
summary(broad_rich_arable1) 
r.squaredGLMM(broad_rich_arable1) # 0.1152524    
broad_rich_arable2 <- glmer(broadleaf_rich ~ prop_arable.1500 + (1|subcompartment) + (1|year),
                            family=poisson(), data=moth_hab_data)
summary(broad_rich_arable2) 
r.squaredGLMM(broad_rich_arable2) # 0.007854100    
broad_rich_arable3 <- glmer(broadleaf_rich ~ prop_arable.3000 + (1|subcompartment) + (1|year),
                            family=poisson(), data=moth_hab_data)
summary(broad_rich_arable3)        
r.squaredGLMM(broad_rich_arable3) # 0.03390211   
#### ARABLE 500M #### (3000m previously)

broad_rich_broad1 <- glmer(broadleaf_rich ~ prop_broadleaf.500 + (1|subcompartment) + (1|year),
                           family=poisson(), data=moth_hab_data)
summary(broad_rich_broad1)    
r.squaredGLMM(broad_rich_broad1) # 0.05081882     
broad_rich_broad2 <- glmer(broadleaf_rich ~ prop_broadleaf.1500 + (1|subcompartment) + (1|year),
                           family=poisson(), data=moth_hab_data)
summary(broad_rich_broad2)   
r.squaredGLMM(broad_rich_broad2) # 0.02466290     
broad_rich_broad3 <- glmer(broadleaf_rich ~ prop_broadleaf.3000 + (1|subcompartment) + (1|year),
                           family=poisson(), data=moth_hab_data)
summary(broad_rich_broad3)      
r.squaredGLMM(broad_rich_broad3) # 0.03161780    
#### BROADLEAF 500M #### (3000m previously)


broad_rich_conifer1 <- glmer(broadleaf_rich ~ prop_conifer.500 + (1|subcompartment) + (1|year),
                             family=poisson(), data=moth_hab_data)
summary(broad_rich_conifer1)           
r.squaredGLMM(broad_rich_conifer1) # 0.01844945      
broad_rich_conifer2 <- glmer(broadleaf_rich ~ prop_conifer.1500 + (1|subcompartment) + (1|year),
                             family=poisson(), data=moth_hab_data)
summary(broad_rich_conifer2)               
r.squaredGLMM(broad_rich_conifer2) # 0.0001670618      
broad_rich_conifer3 <- glmer(broadleaf_rich ~ prop_conifer.3000 + (1|subcompartment) + (1|year),
                             family=poisson(), data=moth_hab_data)
summary(broad_rich_conifer3)                
r.squaredGLMM(broad_rich_conifer3) # 0.00003653343     
#### CONIFER 500M #### (1500m previously)


broad_rich_grass1 <- glmer(broadleaf_rich ~ prop_grassland.500 + (1|subcompartment) + (1|year),
                           family=poisson(), data=moth_hab_data)
summary(broad_rich_grass1)            
r.squaredGLMM(broad_rich_grass1) # 0.003853630       
broad_rich_grass2 <- glmer(broadleaf_rich ~ prop_grassland.1500 + (1|subcompartment) + (1|year),
                           family=poisson(), data=moth_hab_data)
summary(broad_rich_grass2)                    
r.squaredGLMM(broad_rich_grass2) # 0.008908418       
broad_rich_grass3 <- glmer(broadleaf_rich ~ prop_grassland.3000 + (1|subcompartment) + (1|year),
                           family=poisson(), data=moth_hab_data)
summary(broad_rich_grass3)        
r.squaredGLMM(broad_rich_grass3) # 0.04922648      
#### GRASSLAND 3000M #### 

# save results
r2 <- c(r.squaredGLMM(broad_rich_arable1)[1,1], r.squaredGLMM(broad_rich_arable2)[1,1], r.squaredGLMM(broad_rich_arable3)[1,1],
        r.squaredGLMM(broad_rich_broad1)[1,1], r.squaredGLMM(broad_rich_broad2)[1,1], r.squaredGLMM(broad_rich_broad3)[1,1],
        r.squaredGLMM(broad_rich_conifer1)[1,1], r.squaredGLMM(broad_rich_conifer2)[1,1], r.squaredGLMM(broad_rich_conifer3)[1,1],
        r.squaredGLMM(broad_rich_grass1)[1,1], r.squaredGLMM(broad_rich_grass2)[1,1], r.squaredGLMM(broad_rich_grass3)[1,1])
land_cover <- c("Arable.500", "Arable.1500", "Arable.3000", "Broadleaf.500", "Broadleaf.1500", "Broadleaf.3000",
                "Conifer.500", "Conifer.1500", "Conifer.3000", "Grassland.500", "Grassland.1500", "Grassland.3000")
broad_rich_land <- data.frame(Response=rep("Broadleaf moth richness"), Land_cover=land_cover, R2=r2)
write.csv(broad_rich_land, file="Results/Moths/Broadleaf_richness_land_cover.csv", row.names=FALSE)

## Now run models with habitat variables only 
broad_rich_hab <- glmer(broadleaf_rich ~ scale(dbh_average) + scale(basal_area) +
                          scale(per_broadleaf_canopy) + scale(complexity_score) + 
                          scale(canopy_openess) + scale(minimum_temperature) + 
                          scale(cloud_cover_.) + scale(dis_to_edge) + 
                          scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ (1|year) + 
                          (1|subcompartment), data = moth_hab_data, family = poisson(),
                        glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                        na.action = "na.fail")
summary(broad_rich_hab)
check_collinearity(broad_rich_hab) # low correlations after removing moon cycle
testDispersion(broad_rich_hab) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = broad_rich_hab, plot = F)
plot(simulationOutput) ## no assumptions violated!

final_mod2 <- glmer(broadleaf_rich ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + scale(wind_speed) + 
                      scale(dis_broadleaf) + scale(altitude)+ scale(prop_arable.500) + (1|year) + 
                      (1|subcompartment), data = moth_hab_data, family = poisson(),
                    glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                    na.action = "na.fail")
check_collinearity(final_mod2) # low correlations after removing moon cycle
testDispersion(final_mod2) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod2, plot = F)
plot(simulationOutput) ## no assumptions violated!

final_mod3 <- glmer(broadleaf_rich ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + scale(wind_speed) + 
                      scale(dis_broadleaf) + scale(altitude)+ scale(prop_broadleaf.500) + (1|year) + 
                      (1|subcompartment), data = moth_hab_data, family = poisson(),
                    glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                    na.action = "na.fail")
check_collinearity(final_mod3) # low correlations after removing moon cycle
testDispersion(final_mod3) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod3, plot = F)
plot(simulationOutput) ## no assumptions violated!

final_mod4 <- glmer(broadleaf_rich ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + scale(wind_speed) + 
                      scale(dis_broadleaf) + scale(altitude)+ scale(prop_conifer.500) + (1|year) + 
                      (1|subcompartment), data = moth_hab_data, family = poisson(),
                    glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                    na.action = "na.fail")
check_collinearity(final_mod4) # low correlations after removing moon cycle
testDispersion(final_mod4) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod4, plot = F)
plot(simulationOutput) ## no assumptions violated!

final_mod5 <- glmer(broadleaf_rich ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + scale(wind_speed) + 
                      scale(dis_broadleaf) + scale(altitude)+ scale(prop_grassland.3000) + (1|year) + 
                      (1|subcompartment), data = moth_hab_data, family = poisson(),
                    glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                    na.action = "na.fail")
check_collinearity(final_mod5) # low correlations after removing moon cycle
testDispersion(final_mod5) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod5, plot = F)
plot(simulationOutput) ## no assumptions violated!

AIC(broad_rich_hab, final_mod2, final_mod3, final_mod4, final_mod5)
# conifer model reduces AIC by ~3 
summary(final_mod4)
# complexity score (positive)

## save all AIC values for supplementary table
AIC_values <- data.frame(model=c("habitat", "arable 500m", "broadleaf 500m", "conifer 500m", "grassland 3000m"),
                         AIC=AIC(broad_rich_hab, final_mod2, final_mod3, final_mod4, final_mod5)$AIC, stringsAsFactors=FALSE) 
write.csv(AIC_values, file="Results/Broadleaf_richness_AIC.csv", row.names=FALSE)

## save lowest AIC model output for supplementary table
final_mod4_sum <- as.data.frame(coef(summary(final_mod4))) 
final_mod4_sum$parameters <- row.names(final_mod4_sum)
row.names(final_mod4_sum) <- 1:nrow(final_mod4_sum)
final_mod4_sum$model <- "conifer 500m"
final_mod4_sum$AIC <- AIC(final_mod4)

write.csv(final_mod4_sum, file="Results/Broadleaf_richness_habitat.csv", row.names=FALSE)

###################################################################################################################################

# Broadleaf ONLY richness

broad_rich_arable1 <- glmer(broadleaf_only_rich ~ prop_arable.500 + (1|subcompartment) + (1|year),
                            family=poisson(), data=moth_hab_data)
summary(broad_rich_arable1) 
r.squaredGLMM(broad_rich_arable1) # 0.07253447      
broad_rich_arable2 <- glmer(broadleaf_only_rich ~ prop_arable.1500 + (1|subcompartment) + (1|year),
                            family=poisson(), data=moth_hab_data)
summary(broad_rich_arable2) 
r.squaredGLMM(broad_rich_arable2) # 0.0006206497     
broad_rich_arable3 <- glmer(broadleaf_only_rich ~ prop_arable.3000 + (1|subcompartment) + (1|year),
                            family=poisson(), data=moth_hab_data)
summary(broad_rich_arable3)        
r.squaredGLMM(broad_rich_arable3) # 0.02585461    
#### ARABLE 500M ####

broad_rich_broad1 <- glmer(broadleaf_only_rich ~ prop_broadleaf.500 + (1|subcompartment) + (1|year),
                           family=poisson(), data=moth_hab_data)
summary(broad_rich_broad1)    
r.squaredGLMM(broad_rich_broad1) # 0.0007029386      
broad_rich_broad2 <- glmer(broadleaf_only_rich ~ prop_broadleaf.1500 + (1|subcompartment) + (1|year),
                           family=poisson(), data=moth_hab_data)
summary(broad_rich_broad2)   
r.squaredGLMM(broad_rich_broad2) # 0.001825777      
broad_rich_broad3 <- glmer(broadleaf_only_rich ~ prop_broadleaf.3000 + (1|subcompartment) + (1|year),
                           family=poisson(), data=moth_hab_data)
summary(broad_rich_broad3)      
r.squaredGLMM(broad_rich_broad3) # 0.004081388     
#### BROADLEAF 3000m #### 


broad_rich_conifer1 <- glmer(broadleaf_only_rich ~ prop_conifer.500 + (1|subcompartment) + (1|year),
                             family=poisson(), data=moth_hab_data)
summary(broad_rich_conifer1)           
r.squaredGLMM(broad_rich_conifer1) # 0.001596111       
broad_rich_conifer2 <- glmer(broadleaf_only_rich ~ prop_conifer.1500 + (1|subcompartment) + (1|year),
                             family=poisson(), data=moth_hab_data)
summary(broad_rich_conifer2)               
r.squaredGLMM(broad_rich_conifer2) # 0.002539227       
broad_rich_conifer3 <- glmer(broadleaf_only_rich ~ prop_conifer.3000 + (1|subcompartment) + (1|year),
                             family=poisson(), data=moth_hab_data)
summary(broad_rich_conifer3)                
r.squaredGLMM(broad_rich_conifer3) # 0.01234100      
#### CONIFER 3000M #### 


broad_rich_grass1 <- glmer(broadleaf_only_rich ~ prop_grassland.500 + (1|subcompartment) + (1|year),
                           family=poisson(), data=moth_hab_data)
summary(broad_rich_grass1)            
r.squaredGLMM(broad_rich_grass1) # 0.008869970        
broad_rich_grass2 <- glmer(broadleaf_only_rich ~ prop_grassland.1500 + (1|subcompartment) + (1|year),
                           family=poisson(), data=moth_hab_data)
summary(broad_rich_grass2)                    
r.squaredGLMM(broad_rich_grass2) # 0.01697158        
broad_rich_grass3 <- glmer(broadleaf_only_rich ~ prop_grassland.3000 + (1|subcompartment) + (1|year),
                           family=poisson(), data=moth_hab_data)
summary(broad_rich_grass3)        
r.squaredGLMM(broad_rich_grass3) # 0.03007713       
#### GRASSLAND 3000M #### 

# save results
r2 <- c(r.squaredGLMM(broad_rich_arable1)[1,1], r.squaredGLMM(broad_rich_arable2)[1,1], r.squaredGLMM(broad_rich_arable3)[1,1],
        r.squaredGLMM(broad_rich_broad1)[1,1], r.squaredGLMM(broad_rich_broad2)[1,1], r.squaredGLMM(broad_rich_broad3)[1,1],
        r.squaredGLMM(broad_rich_conifer1)[1,1], r.squaredGLMM(broad_rich_conifer2)[1,1], r.squaredGLMM(broad_rich_conifer3)[1,1],
        r.squaredGLMM(broad_rich_grass1)[1,1], r.squaredGLMM(broad_rich_grass2)[1,1], r.squaredGLMM(broad_rich_grass3)[1,1])
land_cover <- c("Arable.500", "Arable.1500", "Arable.3000", "Broadleaf.500", "Broadleaf.1500", "Broadleaf.3000",
                "Conifer.500", "Conifer.1500", "Conifer.3000", "Grassland.500", "Grassland.1500", "Grassland.3000")
broad_rich_land <- data.frame(Response=rep("Broadleaf only moth richness"), Land_cover=land_cover, R2=r2)
write.csv(broad_rich_land, file="Results/Moths/Broadleaf_only_richness_land_cover.csv", row.names=FALSE)

## Now run models with habitat variables only 
broad_rich_hab <- glmer(broadleaf_only_rich ~ scale(dbh_average) + scale(basal_area) +
                          scale(per_broadleaf_canopy) + scale(complexity_score) + 
                          scale(canopy_openess) + scale(minimum_temperature) + 
                          scale(cloud_cover_.) + scale(dis_to_edge) + 
                          scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ (1|year) + 
                          (1|subcompartment), data = moth_hab_data, family = poisson(),
                        glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                        na.action = "na.fail")
summary(broad_rich_hab)
check_collinearity(broad_rich_hab) # low correlations after removing moon cycle
testDispersion(broad_rich_hab) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = broad_rich_hab, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod2 <- glmer(broadleaf_only_rich ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + scale(wind_speed) + 
                      scale(dis_broadleaf) + scale(altitude)+ scale(prop_arable.500) + (1|year) + 
                      (1|subcompartment), data = moth_hab_data, family = poisson(),
                    glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                    na.action = "na.fail")
check_collinearity(final_mod2) # low correlations after removing moon cycle
testDispersion(final_mod2) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod2, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod3 <- glmer(broadleaf_only_rich ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + scale(wind_speed) + 
                      scale(dis_broadleaf) + scale(altitude)+ scale(prop_broadleaf.3000) + (1|year) + 
                      (1|subcompartment), data = moth_hab_data, family = poisson(),
                    glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                    na.action = "na.fail")
check_collinearity(final_mod3) # low correlations after removing moon cycle
testDispersion(final_mod3) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod3, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod4 <- glmer(broadleaf_only_rich ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + scale(wind_speed) + 
                      scale(dis_broadleaf) + scale(altitude)+ scale(prop_conifer.3000) + (1|year) + 
                      (1|subcompartment), data = moth_hab_data, family = poisson(),
                    glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                    na.action = "na.fail")
summary(final_mod4)
check_collinearity(final_mod4) # low correlations after removing moon cycle
testDispersion(final_mod4) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod4, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod5 <- glmer(broadleaf_only_rich ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + scale(wind_speed) + 
                      scale(dis_broadleaf) + scale(altitude)+ scale(prop_grassland.3000) + (1|year) + 
                      (1|subcompartment), data = moth_hab_data, family = poisson(),
                    na.action = "na.fail")
check_collinearity(final_mod5) # low correlations after removing moon cycle
testDispersion(final_mod5) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod5, plot = F)
plot(simulationOutput) ## no assumptions violated

AIC(broad_rich_hab, final_mod2, final_mod3, final_mod4, final_mod5)
# grassland model reduces AIC by ~5
summary(final_mod5)
# no significant habitat variables

## save all AIC values for supplementary table
AIC_values <- data.frame(model=c("habitat", "arable 500m", "broadleaf 3000m", "conifer 3000m", "grassland 3000m"),
                         AIC=AIC(broad_rich_hab, final_mod2, final_mod3, final_mod4, final_mod5)$AIC, stringsAsFactors=FALSE) 
write.csv(AIC_values, file="Results/Broadleaf_only_richness_AIC.csv", row.names=FALSE)

## save lowest AIC model summary
final_mod5_sum <- as.data.frame(coef(summary(final_mod5))) 
final_mod5_sum$parameters <- row.names(final_mod5_sum)
row.names(final_mod5_sum) <- 1:nrow(final_mod5_sum)
final_mod5_sum$model <- "grassland 3000m"
final_mod5_sum$AIC <- AIC(final_mod5)

write.csv(final_mod5_sum, file="Results/Broadleaf_only_richness_habitat.csv", row.names=FALSE)

###################################################################################################################################

# Conifer richness

conifer_rich_arable1 <- glmer(conifer_rich ~ prop_arable.500 + (1|subcompartment) + (1|year),
                              family=poisson(), data=moth_hab_data)
summary(conifer_rich_arable1)   
r.squaredGLMM(conifer_rich_arable1) # 0.08442663     
conifer_rich_arable2 <- glmer(conifer_rich ~ prop_arable.1500 + (1|subcompartment) + (1|year),
                              family=poisson(), data=moth_hab_data)
summary(conifer_rich_arable2)  
r.squaredGLMM(conifer_rich_arable2) # 0.010843373     
conifer_rich_arable3 <- glmer(conifer_rich ~ prop_arable.3000 + (1|subcompartment) + (1|year),
                              family=poisson(), data=moth_hab_data)
summary(conifer_rich_arable3)              
r.squaredGLMM(conifer_rich_arable3) # 0.02470212    
#### ARABLE 500M #### (3000m previously)


conifer_rich_broad1 <- glmer(conifer_rich ~ prop_broadleaf.500 + (1|subcompartment) + (1|year),
                             family=poisson(), data=moth_hab_data)
summary(conifer_rich_broad1)      
r.squaredGLMM(conifer_rich_broad1) # 0.02136088      
conifer_rich_broad2 <- glmer(conifer_rich ~ prop_broadleaf.1500 + (1|subcompartment) + (1|year),
                             family=poisson(), data=moth_hab_data)
summary(conifer_rich_broad2)        
r.squaredGLMM(conifer_rich_broad2) # 0.000004656428      
conifer_rich_broad3 <- glmer(conifer_rich ~ prop_broadleaf.3000 + (1|subcompartment) + (1|year),
                             family=poisson(), data=moth_hab_data)
summary(conifer_rich_broad3)     
r.squaredGLMM(conifer_rich_broad3) # 0.00002412235     
#### BROADLEAF 500M #### 


conifer_rich_conifer1 <- glmer(conifer_rich ~ prop_conifer.500 + (1|subcompartment) + (1|year),
                               family=poisson(), data=moth_hab_data)
summary(conifer_rich_conifer1)                   
r.squaredGLMM(conifer_rich_conifer1) # 0.01287681       
conifer_rich_conifer2 <- glmer(conifer_rich ~ prop_conifer.1500 + (1|subcompartment) + (1|year),
                               family=poisson(), data=moth_hab_data)
summary(conifer_rich_conifer2)               
r.squaredGLMM(conifer_rich_conifer2) # 0.0001407752       
conifer_rich_conifer3 <- glmer(conifer_rich ~ prop_conifer.3000 + (1|subcompartment) + (1|year),
                               family=poisson(), data=moth_hab_data)
summary(conifer_rich_conifer3)                     
r.squaredGLMM(conifer_rich_conifer3) # 0.009829104      
#### CONIFER 500M #### (3000m previously)


conifer_rich_grass1 <- glmer(conifer_rich ~ prop_grassland.500 + (1|subcompartment) + (1|year),
                             family=poisson(), data=moth_hab_data)
summary(conifer_rich_grass1)                   
r.squaredGLMM(conifer_rich_grass1) # 0.01324291        
conifer_rich_grass2 <- glmer(conifer_rich ~ prop_grassland.1500 + (1|subcompartment) + (1|year),
                             family=poisson(), data=moth_hab_data)
summary(conifer_rich_grass2)                    
r.squaredGLMM(conifer_rich_grass2) # 0.009563803        
conifer_rich_grass3 <- glmer(conifer_rich ~ prop_grassland.3000 + (1|subcompartment) + (1|year),
                             family=poisson(), data=moth_hab_data)
summary(conifer_rich_grass3)        
r.squaredGLMM(conifer_rich_grass3) # 0.02215809       
#### GRASSLAND 3000M #### 

# save results
r2 <- c(r.squaredGLMM(conifer_rich_arable1)[1,1], r.squaredGLMM(conifer_rich_arable2)[1,1], r.squaredGLMM(conifer_rich_arable3)[1,1],
        r.squaredGLMM(conifer_rich_broad1)[1,1], r.squaredGLMM(conifer_rich_broad2)[1,1], r.squaredGLMM(conifer_rich_broad3)[1,1],
        r.squaredGLMM(conifer_rich_conifer1)[1,1], r.squaredGLMM(conifer_rich_conifer2)[1,1], r.squaredGLMM(conifer_rich_conifer3)[1,1],
        r.squaredGLMM(conifer_rich_grass1)[1,1], r.squaredGLMM(conifer_rich_grass2)[1,1], r.squaredGLMM(conifer_rich_grass3)[1,1])
land_cover <- c("Arable.500", "Arable.1500", "Arable.3000", "Broadleaf.500", "Broadleaf.1500", "Broadleaf.3000",
                "Conifer.500", "Conifer.1500", "Conifer.3000", "Grassland.500", "Grassland.1500", "Grassland.3000")
conifer_rich_land <- data.frame(Response=rep("Conifer moth richness"), Land_cover=land_cover, R2=r2)
write.csv(conifer_rich_land, file="Results/Moths/Conifer_richness_land_cover.csv", row.names=FALSE)

## Now run models with habitat variables only 
conifer_rich_hab <- glmmTMB(conifer_rich ~ scale(dbh_average) + scale(basal_area) +
                              scale(per_broadleaf_canopy) + scale(complexity_score) + 
                              scale(canopy_openess) + scale(minimum_temperature) + 
                              scale(cloud_cover_.) + scale(dis_to_edge) + 
                              scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ (1|year) + 
                              (1|subcompartment), data = moth_hab_data, family = genpois(),
                            na.action = "na.fail")
summary(conifer_rich_hab)
check_collinearity(conifer_rich_hab) # low correlations after removing moon cycle
testDispersion(conifer_rich_hab) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = conifer_rich_hab, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod2 <- glmmTMB(conifer_rich ~ scale(dbh_average) + scale(basal_area) +
                        scale(per_broadleaf_canopy) + scale(complexity_score) + 
                        scale(canopy_openess) + scale(minimum_temperature) + 
                        scale(cloud_cover_.) + scale(dis_to_edge) + 
                        scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_arable.500) + 
                        (1|year) + (1|subcompartment), data = moth_hab_data, family = genpois(),
                      na.action = "na.fail")
check_collinearity(final_mod2) # low correlations after removing moon cycle
testDispersion(final_mod2) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod2, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod3 <- glmmTMB(conifer_rich ~ scale(dbh_average) + scale(basal_area) +
                        scale(per_broadleaf_canopy) + scale(complexity_score) + 
                        scale(canopy_openess) + scale(minimum_temperature) + 
                        scale(cloud_cover_.) + scale(dis_to_edge) + 
                        scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_broadleaf.500) + 
                        (1|year) + (1|subcompartment), data = moth_hab_data, family = genpois(),
                      na.action = "na.fail")
check_collinearity(final_mod3) # low correlations after removing moon cycle
testDispersion(final_mod3) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod3, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod4 <- glmmTMB(conifer_rich ~scale(dbh_average) + scale(basal_area) +
                        scale(per_broadleaf_canopy) + scale(complexity_score) + 
                        scale(canopy_openess) + scale(minimum_temperature) + 
                        scale(cloud_cover_.) + scale(dis_to_edge) + 
                        scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_conifer.500) + 
                        (1|year) + (1|subcompartment), data = moth_hab_data, family = genpois(),
                      na.action = "na.fail")
check_collinearity(final_mod4) # low correlations after removing moon cycle
testDispersion(final_mod4) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod4, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod5 <- glmmTMB(conifer_rich ~ scale(dbh_average) + scale(basal_area) +
                        scale(per_broadleaf_canopy) + scale(complexity_score) + 
                        scale(canopy_openess) + scale(minimum_temperature) + 
                        scale(cloud_cover_.) + scale(dis_to_edge) + 
                        scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_grassland.3000) + 
                        (1|year) + (1|subcompartment), data = moth_hab_data, family = genpois(),
                      na.action = "na.fail")
check_collinearity(final_mod5) # low correlations after removing moon cycle
testDispersion(final_mod5) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod5, plot = F)
plot(simulationOutput) ## no assumptions violated

AIC(conifer_rich_hab, final_mod2, final_mod3, final_mod4, final_mod5)
# conifer model reduces AIC by ~5
summary(final_mod4)
# dbh (negative) and basal area (positive)

## save all AIC values for supplementary table
AIC_values <- data.frame(model=c("habitat", "arable 500m", "broadleaf 500m", "conifer 500m", "grassland 3000m"),
                         AIC=AIC(conifer_rich_hab, final_mod2, final_mod3, final_mod4, final_mod5)$AIC, stringsAsFactors=FALSE) 
write.csv(AIC_values, file="Results/Conifer_richness_AIC.csv", row.names=FALSE)

## save lowest AIC model summary
final_mod4_sum <- as.data.frame(summary(final_mod4)$coefficients$cond) 
final_mod4_sum$parameters <- row.names(final_mod4_sum)
row.names(final_mod4_sum) <- 1:nrow(final_mod4_sum)
final_mod4_sum$model <- "conifer"
final_mod4_sum$AIC <- AIC(final_mod4)

write.csv(final_mod4_sum, file="Results/Conifer_richness_habitat.csv", row.names=FALSE)


###################################################################################################################################

# Conifer ONLY richness

conifer_rich_arable1 <- glmer(conifer_only_rich ~ prop_arable.500 + (1|subcompartment) + (1|year),
                              family=poisson(), data=moth_hab_data)
summary(conifer_rich_arable1)   
r.squaredGLMM(conifer_rich_arable1) # 0.08138858      
conifer_rich_arable2 <- glmer(conifer_only_rich ~ prop_arable.1500 + (1|subcompartment) + (1|year),
                              family=poisson(), data=moth_hab_data)
summary(conifer_rich_arable2)  
r.squaredGLMM(conifer_rich_arable2) # 0.02027563      
conifer_rich_arable3 <- glmer(conifer_only_rich ~ prop_arable.3000 + (1|subcompartment) + (1|year),
                              family=poisson(), data=moth_hab_data)
summary(conifer_rich_arable3)              
r.squaredGLMM(conifer_rich_arable3) # 0.002673330     
#### ARABLE 500M #### 


conifer_rich_broad1 <- glmer(conifer_only_rich ~ prop_broadleaf.500 + (1|subcompartment) + (1|year),
                             family=poisson(), data=moth_hab_data)
summary(conifer_rich_broad1)      
r.squaredGLMM(conifer_rich_broad1) # 0.00009950146       
conifer_rich_broad2 <- glmer(conifer_only_rich ~ prop_broadleaf.1500 + (1|subcompartment) + (1|year),
                             family=poisson(), data=moth_hab_data)
summary(conifer_rich_broad2)        
r.squaredGLMM(conifer_rich_broad2) # 0.06173057       
conifer_rich_broad3 <- glmer(conifer_only_rich ~ prop_broadleaf.3000 + (1|subcompartment) + (1|year),
                             family=poisson(), data=moth_hab_data)
summary(conifer_rich_broad3)     
r.squaredGLMM(conifer_rich_broad3) # 0.06779949      
#### BROADLEAF 3000M #### 


conifer_rich_conifer1 <- glmer(conifer_only_rich ~ prop_conifer.500 + (1|subcompartment) + (1|year),
                               family=poisson(), data=moth_hab_data)
summary(conifer_rich_conifer1)                   
r.squaredGLMM(conifer_rich_conifer1) # 0.000002869582        
conifer_rich_conifer2 <- glmer(conifer_only_rich ~ prop_conifer.1500 + (1|subcompartment) + (1|year),
                               family=poisson(), data=moth_hab_data)
summary(conifer_rich_conifer2)               
r.squaredGLMM(conifer_rich_conifer2) # 0.002371121        
conifer_rich_conifer3 <- glmer(conifer_only_rich ~ prop_conifer.3000 + (1|subcompartment) + (1|year),
                               family=poisson(), data=moth_hab_data)
summary(conifer_rich_conifer3)                     
r.squaredGLMM(conifer_rich_conifer3) # 0.04413466       
#### CONIFER 3000m #### 


conifer_rich_grass1 <- glmer(conifer_only_rich ~ prop_grassland.500 + (1|subcompartment) + (1|year),
                             family=poisson(), data=moth_hab_data)
summary(conifer_rich_grass1)                   
r.squaredGLMM(conifer_rich_grass1) # 0.01710834         
conifer_rich_grass2 <- glmer(conifer_only_rich ~ prop_grassland.1500 + (1|subcompartment) + (1|year),
                             family=poisson(), data=moth_hab_data)
summary(conifer_rich_grass2)                    
r.squaredGLMM(conifer_rich_grass2) # 0.0010246101         
conifer_rich_grass3 <- glmer(conifer_only_rich ~ prop_grassland.3000 + (1|subcompartment) + (1|year),
                             family=poisson(), data=moth_hab_data)
summary(conifer_rich_grass3)        
r.squaredGLMM(conifer_rich_grass3) # 0.012603988        
#### GRASSLAND 500M #### 

# save results
r2 <- c(r.squaredGLMM(conifer_rich_arable1)[1,1], r.squaredGLMM(conifer_rich_arable2)[1,1], r.squaredGLMM(conifer_rich_arable3)[1,1],
        r.squaredGLMM(conifer_rich_broad1)[1,1], r.squaredGLMM(conifer_rich_broad2)[1,1], r.squaredGLMM(conifer_rich_broad3)[1,1],
        r.squaredGLMM(conifer_rich_conifer1)[1,1], r.squaredGLMM(conifer_rich_conifer2)[1,1], r.squaredGLMM(conifer_rich_conifer3)[1,1],
        r.squaredGLMM(conifer_rich_grass1)[1,1], r.squaredGLMM(conifer_rich_grass2)[1,1], r.squaredGLMM(conifer_rich_grass3)[1,1])
land_cover <- c("Arable.500", "Arable.1500", "Arable.3000", "Broadleaf.500", "Broadleaf.1500", "Broadleaf.3000",
                "Conifer.500", "Conifer.1500", "Conifer.3000", "Grassland.500", "Grassland.1500", "Grassland.3000")
conifer_rich_land <- data.frame(Response=rep("Conifer only moth richness"), Land_cover=land_cover, R2=r2)
write.csv(conifer_rich_land, file="Results/Moths/Conifer_only_richness_land_cover.csv", row.names=FALSE)

## Now run models with habitat variables only 
conifer_rich_hab <- glmmTMB(conifer_only_rich ~ scale(dbh_average) + scale(basal_area) +
                              scale(per_broadleaf_canopy) + scale(complexity_score) + 
                              scale(canopy_openess) + scale(minimum_temperature) + 
                              scale(cloud_cover_.) + scale(dis_to_edge) + 
                              scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ (1|year) + 
                              (1|subcompartment), data = moth_hab_data, family = compois(),
                            na.action = "na.fail")
check_collinearity(conifer_rich_hab) # low correlations after removing moon cycle
testDispersion(conifer_rich_hab) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = conifer_rich_hab, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod2 <- glmmTMB(conifer_only_rich ~ scale(dbh_average) + scale(basal_area) +
                        scale(per_broadleaf_canopy) + scale(complexity_score) + 
                        scale(canopy_openess) + scale(minimum_temperature) + 
                        scale(cloud_cover_.) + scale(dis_to_edge) + 
                        scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_arable.500) + 
                        (1|year) + (1|subcompartment), data = moth_hab_data, family = compois(),
                      na.action = "na.fail")
check_collinearity(final_mod2) # low correlations after removing moon cycle
testDispersion(final_mod2) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod2, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod3 <- glmmTMB(conifer_only_rich ~ scale(dbh_average) + scale(basal_area) +
                        scale(per_broadleaf_canopy) + scale(complexity_score) + 
                        scale(canopy_openess) + scale(minimum_temperature) + 
                        scale(cloud_cover_.) + scale(dis_to_edge) + 
                        scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_broadleaf.3000) + 
                        (1|year) + (1|subcompartment), data = moth_hab_data, family = compois(),
                      na.action = "na.fail")
check_collinearity(final_mod3) # low correlations after removing moon cycle
testDispersion(final_mod3) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod3, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod4 <- glmmTMB(conifer_only_rich ~scale(dbh_average) + scale(basal_area) +
                        scale(per_broadleaf_canopy) + scale(complexity_score) + 
                        scale(canopy_openess) + scale(minimum_temperature) + 
                        scale(cloud_cover_.) + scale(dis_to_edge) + 
                        scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_conifer.3000) + 
                        (1|year) + (1|subcompartment), data = moth_hab_data, family = compois(),
                      na.action = "na.fail")
check_collinearity(final_mod4) # low correlations after removing moon cycle
testDispersion(final_mod4) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod4, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod5 <- glmmTMB(conifer_only_rich ~ scale(dbh_average) + scale(basal_area) +
                        scale(per_broadleaf_canopy) + scale(complexity_score) + 
                        scale(canopy_openess) + scale(minimum_temperature) + 
                        scale(cloud_cover_.) + scale(dis_to_edge) + 
                        scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_grassland.500) + 
                        (1|year) + (1|subcompartment), data = moth_hab_data, family = compois(),
                      na.action = "na.fail")
summary(final_mod5)
check_collinearity(final_mod5) # low correlations after removing moon cycle
testDispersion(final_mod5) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod5, plot = F)
plot(simulationOutput) ## no assumptions violated

AIC(conifer_rich_hab, final_mod2, final_mod3, final_mod4, final_mod5)
# arable and grassland models are equivalent (both ~4 lower than habitat)
summary(final_mod5)

## save all AIC values for supplementary table
AIC_values <- data.frame(model=c("habitat", "arable 500m", "broadleaf 3000m", "conifer 3000m", "grassland 500m"),
                         AIC=AIC(conifer_rich_hab, final_mod2, final_mod3, final_mod4, final_mod5)$AIC, stringsAsFactors=FALSE) 
write.csv(AIC_values, file="Results/Conifer_only_richness_AIC.csv", row.names=FALSE)

## save lowest AIC model summary
final_mod2_sum <- as.data.frame(summary(final_mod2)$coefficients$cond)
final_mod2_sum$parameters <- row.names(final_mod2_sum)
row.names(final_mod2_sum) <- 1:nrow(final_mod2_sum)
final_mod2_sum$model <- "arable 500m"
final_mod2_sum$AIC <- AIC(final_mod2)

final_mod5_sum <- as.data.frame(summary(final_mod5)$coefficients$cond)
final_mod5_sum$parameters <- row.names(final_mod5_sum)
row.names(final_mod5_sum) <- 1:nrow(final_mod5_sum)
final_mod5_sum$model <- "grassland 500m"
final_mod5_sum$AIC <- AIC(final_mod5)

conifer_rich_hab_mods <- rbind(final_mod2_sum, final_mod5_sum)
write.csv(conifer_rich_hab_mods, file="Results/Conifer_only_richness_habitat.csv", row.names=FALSE)

###################################################################################################################################

# Woodland richness

wood_rich_arable1 <- glmer.nb(woodland_rich ~ prop_arable.500 + (1|subcompartment) + (1|year),
                              data=moth_hab_data)
summary(wood_rich_arable1)   
r.squaredGLMM(wood_rich_arable1) # 0.1183356       
wood_rich_arable2 <- glmer.nb(woodland_rich ~ prop_arable.1500 + (1|subcompartment) + (1|year),
                              data=moth_hab_data)
summary(wood_rich_arable2)      
r.squaredGLMM(wood_rich_arable2) # 0.00003728517       
wood_rich_arable3 <- glmer.nb(woodland_rich ~ prop_arable.3000 + (1|subcompartment) + (1|year),
                              data=moth_hab_data)
summary(wood_rich_arable3)                          
r.squaredGLMM(wood_rich_arable3) # 0.001446814      
#### ARABLE 500M #### 


wood_rich_broad1 <- glmer.nb(woodland_rich ~ prop_broadleaf.500 + (1|subcompartment) + (1|year),
                             data=moth_hab_data)
summary(wood_rich_broad1)      
r.squaredGLMM(wood_rich_broad1) # 0.006769927        
wood_rich_broad2 <- glmer.nb(woodland_rich ~ prop_broadleaf.1500 + (1|subcompartment) + (1|year),
                             data=moth_hab_data)
summary(wood_rich_broad2)         
r.squaredGLMM(wood_rich_broad2) # 0.002365153        
wood_rich_broad3 <- glmer.nb(woodland_rich ~ prop_broadleaf.3000 + (1|subcompartment) + (1|year),
                             data=moth_hab_data)
summary(wood_rich_broad3)                
r.squaredGLMM(wood_rich_broad3) # 0.003968475       
#### BROADLEAF 500M #### (1500m previously)

wood_rich_conifer1 <- glmer.nb(woodland_rich ~ prop_conifer.500 + (1|subcompartment) + (1|year),
                               data=moth_hab_data)
summary(wood_rich_conifer1)                            
r.squaredGLMM(wood_rich_conifer1) # 0.0004031764         
wood_rich_conifer2 <- glmer.nb(woodland_rich ~ prop_conifer.1500 + (1|subcompartment) + (1|year),
                               data=moth_hab_data)
summary(wood_rich_conifer2)                       
r.squaredGLMM(wood_rich_conifer2) # 0.001064857         
wood_rich_conifer3 <- glmer.nb(woodland_rich ~ prop_conifer.3000 + (1|subcompartment) + (1|year),
                               data=moth_hab_data)
summary(wood_rich_conifer3)                               
r.squaredGLMM(wood_rich_conifer3) # 0.01641014        
#### CONIFER 3000M #### 


wood_rich_grass1 <- glmer.nb(woodland_rich ~ prop_grassland.500 + (1|subcompartment) + (1|year),
                             data=moth_hab_data)
summary(wood_rich_grass1)                          
r.squaredGLMM(wood_rich_grass1) # 0.01655224          
wood_rich_grass2 <- glmer.nb(woodland_rich ~ prop_grassland.1500 + (1|subcompartment) + (1|year),
                             data=moth_hab_data)
summary(wood_rich_grass2)                          
r.squaredGLMM(wood_rich_grass2) # 0.000000001038204          
wood_rich_grass3 <- glmer.nb(woodland_rich ~ prop_grassland.3000 + (1|subcompartment) + (1|year),
                             data=moth_hab_data)
summary(wood_rich_grass3)              
r.squaredGLMM(wood_rich_grass3) # 0.0005774505         
#### GRASSLAND 500M #### (3000m previously)

# save results
r2 <- c(r.squaredGLMM(wood_rich_arable1)[1,1], r.squaredGLMM(wood_rich_arable2)[1,1], r.squaredGLMM(wood_rich_arable3)[1,1],
        r.squaredGLMM(wood_rich_broad1)[1,1], r.squaredGLMM(wood_rich_broad2)[1,1], r.squaredGLMM(wood_rich_broad3)[1,1],
        r.squaredGLMM(wood_rich_conifer1)[1,1], r.squaredGLMM(wood_rich_conifer2)[1,1], r.squaredGLMM(wood_rich_conifer3)[1,1],
        r.squaredGLMM(wood_rich_grass1)[1,1], r.squaredGLMM(wood_rich_grass2)[1,1], r.squaredGLMM(wood_rich_grass3)[1,1])
land_cover <- c("Arable.500", "Arable.1500", "Arable.3000", "Broadleaf.500", "Broadleaf.1500", "Broadleaf.3000",
                "Conifer.500", "Conifer.1500", "Conifer.3000", "Grassland.500", "Grassland.1500", "Grassland.3000")
wood_rich_land <- data.frame(Response=rep("Woodland moth richness"), Land_cover=land_cover, R2=r2)
write.csv(wood_rich_land, file="Results/Moths/Woodland_richness_land_cover.csv", row.names=FALSE)

## Now run models with habitat variables only 
wood_rich_hab <- glmer(woodland_rich ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ (1|year) + 
                         (1|subcompartment), data = moth_hab_data, family = poisson(),
                       glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                       na.action = "na.fail")
check_collinearity(wood_rich_hab) # low correlations after removing moon cycle
testDispersion(wood_rich_hab) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = wood_rich_hab, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod2 <- glmer(woodland_rich ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + 
                      scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_arable.500) + 
                      (1|year) + (1|subcompartment), data = moth_hab_data, family = poisson())
check_collinearity(final_mod2) # low correlations after removing moon cycle
testDispersion(final_mod2) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod2, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod3 <- glmer(woodland_rich ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + 
                      scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_broadleaf.500) + 
                      (1|year) + (1|subcompartment), data = moth_hab_data, family = poisson())
check_collinearity(final_mod3) # low correlations after removing moon cycle
testDispersion(final_mod3) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod3, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod4 <- glmer(woodland_rich ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + 
                      scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_conifer.3000) + 
                      (1|year) + (1|subcompartment), data = moth_hab_data, family = poisson())
check_collinearity(final_mod4) # low correlations after removing moon cycle
testDispersion(final_mod4) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod4, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod5 <- glmer(woodland_rich ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + 
                      scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_grassland.500) + 
                      (1|year) + (1|subcompartment), data = moth_hab_data, family = poisson())
check_collinearity(final_mod5) # low correlations after removing moon cycle
testDispersion(final_mod5) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod5, plot = F)
plot(simulationOutput) ## no assumptions violated

AIC(wood_rich_hab, final_mod2, final_mod3, final_mod4, final_mod5)
# nothing lower than 2 AIC values
summary(wood_rich_hab)
# complexity score (positive) and canopy openness (negative)

## save all AIC values for supplementary table
AIC_values <- data.frame(model=c("habitat", "arable 500m", "broadleaf 500m", "conifer 3000m", "grassland 500m"),
                         AIC=AIC(wood_rich_hab, final_mod2, final_mod3, final_mod4, final_mod5)$AIC, stringsAsFactors=FALSE) 
write.csv(AIC_values, file="Results/Woodland_richness_AIC.csv", row.names=FALSE)

## save lowest AIC model summary
wood_rich_hab_sum <- as.data.frame(coef(summary(wood_rich_hab))) 
wood_rich_hab_sum$parameters <- row.names(wood_rich_hab_sum)
row.names(wood_rich_hab_sum) <- 1:nrow(wood_rich_hab_sum)
wood_rich_hab_sum$model <- "habitat"
wood_rich_hab_sum$AIC <- AIC(wood_rich_hab)
write.csv(wood_rich_hab_sum, file="Results/Woodland_richness_habitat.csv", row.names=FALSE)

###################################################################################################################################

# Woodland ONLY richness

wood_rich_arable1 <- glmer.nb(woodland_only_rich ~ prop_arable.500 + (1|subcompartment) + (1|year),
                              data=moth_hab_data)
summary(wood_rich_arable1)   
r.squaredGLMM(wood_rich_arable1) # 0.1486495        
wood_rich_arable2 <- glmer.nb(woodland_only_rich ~ prop_arable.1500 + (1|subcompartment) + (1|year),
                              data=moth_hab_data)
summary(wood_rich_arable2)      
r.squaredGLMM(wood_rich_arable2) # 0.008224477        
wood_rich_arable3 <- glmer.nb(woodland_only_rich ~ prop_arable.3000 + (1|subcompartment) + (1|year),
                              data=moth_hab_data)
summary(wood_rich_arable3)                          
r.squaredGLMM(wood_rich_arable3) # 0.003122006       
#### ARABLE 500M #### 


wood_rich_broad1 <- glmer.nb(woodland_only_rich ~ prop_broadleaf.500 + (1|subcompartment) + (1|year),
                             data=moth_hab_data)
summary(wood_rich_broad1)      
r.squaredGLMM(wood_rich_broad1) # 0.02741532         
wood_rich_broad2 <- glmer.nb(woodland_only_rich ~ prop_broadleaf.1500 + (1|subcompartment) + (1|year),
                             data=moth_hab_data)
summary(wood_rich_broad2)         
r.squaredGLMM(wood_rich_broad2) # 0.001618020         
wood_rich_broad3 <- glmer.nb(woodland_only_rich ~ prop_broadleaf.3000 + (1|subcompartment) + (1|year),
                             data=moth_hab_data)
summary(wood_rich_broad3)                
r.squaredGLMM(wood_rich_broad3) # 0.0008166710        
#### BROADLEAF 500M #### 

wood_rich_conifer1 <- glmer.nb(woodland_only_rich ~ prop_conifer.500 + (1|subcompartment) + (1|year),
                               data=moth_hab_data)
summary(wood_rich_conifer1)                            
r.squaredGLMM(wood_rich_conifer1) # 0.009858200          
wood_rich_conifer2 <- glmer.nb(woodland_only_rich ~ prop_conifer.1500 + (1|subcompartment) + (1|year),
                               data=moth_hab_data)
summary(wood_rich_conifer2)                       
r.squaredGLMM(wood_rich_conifer2) # 0.0001875802          
wood_rich_conifer3 <- glmer.nb(woodland_only_rich ~ prop_conifer.3000 + (1|subcompartment) + (1|year),
                               data=moth_hab_data)
summary(wood_rich_conifer3)                               
r.squaredGLMM(wood_rich_conifer3) # 0.008355841         
#### CONIFER 500M #### 


wood_rich_grass1 <- glmer.nb(woodland_only_rich ~ prop_grassland.500 + (1|subcompartment) + (1|year),
                             data=moth_hab_data)
summary(wood_rich_grass1)                          
r.squaredGLMM(wood_rich_grass1) # 0.01442785            
wood_rich_grass2 <- glmer.nb(woodland_only_rich ~ prop_grassland.1500 + (1|subcompartment) + (1|year),
                             data=moth_hab_data)
summary(wood_rich_grass2)                          
r.squaredGLMM(wood_rich_grass2) # 0.0005004298            
wood_rich_grass3 <- glmer.nb(woodland_only_rich ~ prop_grassland.3000 + (1|subcompartment) + (1|year),
                             data=moth_hab_data)
summary(wood_rich_grass3)              
r.squaredGLMM(wood_rich_grass3) # 0.003863443           
#### GRASSLAND 500M #### 

# save results
r2 <- c(r.squaredGLMM(wood_rich_arable1)[1,1], r.squaredGLMM(wood_rich_arable2)[1,1], r.squaredGLMM(wood_rich_arable3)[1,1],
        r.squaredGLMM(wood_rich_broad1)[1,1], r.squaredGLMM(wood_rich_broad2)[1,1], r.squaredGLMM(wood_rich_broad3)[1,1],
        r.squaredGLMM(wood_rich_conifer1)[1,1], r.squaredGLMM(wood_rich_conifer2)[1,1], r.squaredGLMM(wood_rich_conifer3)[1,1],
        r.squaredGLMM(wood_rich_grass1)[1,1], r.squaredGLMM(wood_rich_grass2)[1,1], r.squaredGLMM(wood_rich_grass3)[1,1])
land_cover <- c("Arable.500", "Arable.1500", "Arable.3000", "Broadleaf.500", "Broadleaf.1500", "Broadleaf.3000",
                "Conifer.500", "Conifer.1500", "Conifer.3000", "Grassland.500", "Grassland.1500", "Grassland.3000")
wood_rich_land <- data.frame(Response=rep("Woodland only moth richness"), Land_cover=land_cover, R2=r2)
write.csv(wood_rich_land, file="Results/Moths/Woodland_only_richness_land_cover.csv", row.names=FALSE)

## Now run models with habitat variables only 
wood_rich_hab <- glmer(woodland_only_rich ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ (1|year) + 
                         (1|subcompartment), data = moth_hab_data, family = poisson(),
                       glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                       na.action = "na.fail")
check_collinearity(wood_rich_hab) # low correlations after removing moon cycle
testDispersion(wood_rich_hab) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = wood_rich_hab, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod2 <- glmer(woodland_only_rich ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + 
                      scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_arable.500) + 
                      (1|year) + (1|subcompartment), data = moth_hab_data, family = poisson())
check_collinearity(final_mod2) # low correlations after removing moon cycle
testDispersion(final_mod2) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod2, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod3 <- glmer(woodland_only_rich ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + 
                      scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_broadleaf.500) + 
                      (1|year) + (1|subcompartment), data = moth_hab_data, family = poisson())
check_collinearity(final_mod3) # low correlations after removing moon cycle
testDispersion(final_mod3) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod3, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod4 <- glmer(woodland_only_rich ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + 
                      scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_conifer.500) + 
                      (1|year) + (1|subcompartment), data = moth_hab_data, family = poisson())
check_collinearity(final_mod4) # low correlations after removing moon cycle
testDispersion(final_mod4) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod4, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod5 <- glmer(woodland_only_rich ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + 
                      scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_grassland.500) + 
                      (1|year) + (1|subcompartment), data = moth_hab_data, family = poisson())
check_collinearity(final_mod5) # low correlations after removing moon cycle
testDispersion(final_mod5) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod5, plot = F)
plot(simulationOutput) ## no assumptions violated

AIC(wood_rich_hab, final_mod2, final_mod3, final_mod4, final_mod5)
# conifer and grassland both lower than habitat
summary(final_mod4) # no significant structure variables
summary(final_mod5) # no significant structure variables

## save all AIC values for supplementary table
AIC_values <- data.frame(model=c("habitat", "arable 500m", "broadleaf 500m", "conifer 500m", "grassland 500m"),
                         AIC=AIC(wood_rich_hab, final_mod2, final_mod3, final_mod4, final_mod5)$AIC, stringsAsFactors=FALSE) 
write.csv(AIC_values, file="Results/Woodland_only_richness_AIC.csv", row.names=FALSE)

## save lowest AIC model summary
final_mod4_sum <- as.data.frame(coef(summary(final_mod4))) #selecting full model coefficient averages
final_mod4_sum$parameters <- row.names(final_mod4_sum)
row.names(final_mod4_sum) <- 1:nrow(final_mod4_sum)
final_mod4_sum$model <- "conifer"
final_mod4_sum$AIC <- AIC(final_mod4)

final_mod5_sum <- as.data.frame(coef(summary(final_mod5))) #selecting full model coefficient averages
final_mod5_sum$parameters <- row.names(final_mod5_sum)
row.names(final_mod5_sum) <- 1:nrow(final_mod5_sum)
final_mod5_sum$model <- "grassland"
final_mod5_sum$AIC <- AIC(final_mod5)

wood_rich_hab_mods <- rbind(final_mod4_sum, final_mod5_sum)
write.csv(wood_rich_hab_mods, file="Results/Woodland_only_richness_habitat.csv", row.names=FALSE)

###################################################################################################################################

# Grassland richness

grass_rich_arable1 <- glmer(grass_rich ~ prop_arable.500 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(grass_rich_arable1) 
r.squaredGLMM(grass_rich_arable1) # 0.02023435         
grass_rich_arable2 <- glmer(grass_rich ~ prop_arable.1500 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(grass_rich_arable2)              
r.squaredGLMM(grass_rich_arable2) # 0.03061222        
grass_rich_arable3 <- glmer(grass_rich ~ prop_arable.3000 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(grass_rich_arable3)                          
r.squaredGLMM(grass_rich_arable3) # 0.0009083745       
#### ARABLE 1500M #### 


grass_rich_broad1 <- glmer(grass_rich ~ prop_broadleaf.500 + (1|subcompartment) + (1|year),
                           data=moth_hab_data, family = poisson())
summary(grass_rich_broad1) 
r.squaredGLMM(grass_rich_broad1) # 0.03974579         
grass_rich_broad2 <- glmer(grass_rich ~ prop_broadleaf.1500 + (1|subcompartment) + (1|year),
                           data=moth_hab_data, family = poisson())
summary(grass_rich_broad2)       
r.squaredGLMM(grass_rich_broad2) # 0.08362410         
grass_rich_broad3 <- glmer(grass_rich ~ prop_broadleaf.3000 + (1|subcompartment) + (1|year),
                           data=moth_hab_data, family = poisson())
summary(grass_rich_broad3)               
r.squaredGLMM(grass_rich_broad3) # 0.08055298        
#### BROADLEAF 1500M #### 


grass_rich_conifer1 <- glmer(grass_rich ~ prop_conifer.500 + (1|subcompartment) + (1|year),
                             data=moth_hab_data, family = poisson())
summary(grass_rich_conifer1)                           
r.squaredGLMM(grass_rich_conifer1) # 0.010003346          
grass_rich_conifer2 <- glmer(grass_rich ~ prop_conifer.1500 + (1|subcompartment) + (1|year),
                             data=moth_hab_data, family = poisson())
summary(grass_rich_conifer2)                      
r.squaredGLMM(grass_rich_conifer2) # 0.0004367106          
grass_rich_conifer3 <- glmer(grass_rich ~ prop_conifer.3000 + (1|subcompartment) + (1|year),
                             data=moth_hab_data, family = poisson())
summary(grass_rich_conifer3)                                    
r.squaredGLMM(grass_rich_conifer3) # 0.04190886         
#### CONIFER 3000M #### (1500m previously) 


grass_rich_grass1 <- glmer(grass_rich ~ prop_grassland.500 + (1|subcompartment) + (1|year),
                           data=moth_hab_data, family = poisson())
summary(grass_rich_grass1)                             
r.squaredGLMM(grass_rich_grass1) # 0.02038663           
grass_rich_grass2 <- glmer(grass_rich ~ prop_grassland.1500 + (1|subcompartment) + (1|year),
                           data=moth_hab_data, family = poisson())
summary(grass_rich_grass2)                              
r.squaredGLMM(grass_rich_grass2) # 0.007680437           
grass_rich_grass3 <- glmer(grass_rich ~ prop_grassland.3000 + (1|subcompartment) + (1|year),
                           data=moth_hab_data, family = poisson())
summary(grass_rich_grass3)                   
r.squaredGLMM(grass_rich_grass3) # 0.01022614          
#### GRASSLAND 500M #### (3000m previously)

# save results
r2 <- c(r.squaredGLMM(grass_rich_arable1)[1,1], r.squaredGLMM(grass_rich_arable2)[1,1], r.squaredGLMM(grass_rich_arable3)[1,1],
        r.squaredGLMM(grass_rich_broad1)[1,1], r.squaredGLMM(grass_rich_broad2)[1,1], r.squaredGLMM(grass_rich_broad3)[1,1],
        r.squaredGLMM(grass_rich_conifer1)[1,1], r.squaredGLMM(grass_rich_conifer2)[1,1], r.squaredGLMM(grass_rich_conifer3)[1,1],
        r.squaredGLMM(grass_rich_grass1)[1,1], r.squaredGLMM(grass_rich_grass2)[1,1], r.squaredGLMM(grass_rich_grass3)[1,1])
land_cover <- c("Arable.500", "Arable.1500", "Arable.3000", "Broadleaf.500", "Broadleaf.1500", "Broadleaf.3000",
                "Conifer.500", "Conifer.1500", "Conifer.3000", "Grassland.500", "Grassland.1500", "Grassland.3000")
grass_rich_land <- data.frame(Response=rep("Grassland moth richness"), Land_cover=land_cover, R2=r2)
write.csv(grass_rich_land, file="Results/Moths/Grassland_richness_land_cover.csv", row.names=FALSE)

## Now run models with habitat variables only 
grass_rich_hab <- glmmTMB(grass_rich ~ scale(dbh_average) + scale(basal_area) +
                            scale(per_broadleaf_canopy) + scale(complexity_score) + 
                            scale(canopy_openess) + scale(minimum_temperature) + 
                            scale(cloud_cover_.) + scale(dis_to_edge) + 
                            scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ (1|year) + 
                            (1|subcompartment), data = moth_hab_data, family = genpois(link="log"),
                          na.action = "na.fail")
check_collinearity(grass_rich_hab) # low correlations after removing moon cycle
testDispersion(grass_rich_hab) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = grass_rich_hab, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod2 <- glmmTMB(grass_rich ~ scale(dbh_average) + scale(basal_area) +
                        scale(per_broadleaf_canopy) + scale(complexity_score) + 
                        scale(canopy_openess) + scale(minimum_temperature) + 
                        scale(cloud_cover_.) + scale(dis_to_edge) + 
                        scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_arable.1500) + 
                        (1|year) + (1|subcompartment), data = moth_hab_data, family = genpois(link="log"),
                      na.action = "na.fail")
check_collinearity(final_mod2) # low correlations after removing moon cycle
testDispersion(final_mod2) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod2, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod3 <- glmmTMB(grass_rich ~ scale(dbh_average) + scale(basal_area) +
                        scale(per_broadleaf_canopy) + scale(complexity_score) + 
                        scale(canopy_openess) + scale(minimum_temperature) + 
                        scale(cloud_cover_.) + scale(dis_to_edge) + 
                        scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_broadleaf.1500) + 
                        (1|year) + (1|subcompartment), data = moth_hab_data, family = genpois(link="log"),
                      na.action = "na.fail")
check_collinearity(final_mod3) # low correlations after removing moon cycle
testDispersion(final_mod3) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod3, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod4 <- glmmTMB(grass_rich ~ scale(dbh_average) + scale(basal_area) +
                        scale(per_broadleaf_canopy) + scale(complexity_score) + 
                        scale(canopy_openess) + scale(minimum_temperature) + 
                        scale(cloud_cover_.) + scale(dis_to_edge) + 
                        scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_conifer.3000) + 
                        (1|year) + (1|subcompartment), data = moth_hab_data, family = genpois(link="log"),
                      na.action = "na.fail")
check_collinearity(final_mod4) # low correlations after removing moon cycle
testDispersion(final_mod4) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod4, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod5 <- glmmTMB(grass_rich ~ scale(dbh_average) + scale(basal_area) +
                        scale(per_broadleaf_canopy) + scale(complexity_score) + 
                        scale(canopy_openess) + scale(minimum_temperature) + 
                        scale(cloud_cover_.) + scale(dis_to_edge) + 
                        scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_grassland.500) + 
                        (1|year) + (1|subcompartment), data = moth_hab_data, family = genpois(link="log"),
                      na.action = "na.fail")
check_collinearity(final_mod5) # low correlations after removing moon cycle
testDispersion(final_mod5) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod5, plot = F)
plot(simulationOutput) ## no assumptions violated

AIC(grass_rich_hab, final_mod2, final_mod3, final_mod4, final_mod5)
# habitat model 
summary(grass_rich_hab) # no significant structure variables

## save all AIC values for supplementary table
AIC_values <- data.frame(model=c("habitat", "arable 1500m", "broadleaf 1500m", "conifer 3000m", "grassland 500m"),
                         AIC=AIC(grass_rich_hab, final_mod2, final_mod3, final_mod4, final_mod5)$AIC, stringsAsFactors=FALSE) 
write.csv(AIC_values, file="Results/Grassland_richness_AIC.csv", row.names=FALSE)

## save lowest AIC model summary
grass_rich_hab_sum <- as.data.frame(summary(grass_rich_hab)$coefficients$cond) 
grass_rich_hab_sum$parameters <- row.names(grass_rich_hab_sum)
row.names(grass_rich_hab_sum) <- 1:nrow(grass_rich_hab_sum)
grass_rich_hab_sum$model <- "habitat"
grass_rich_hab_sum$AIC <- AIC(grass_rich_hab)

write.csv(grass_rich_hab_sum, file="Results/Grassland_richness_habitat.csv", row.names=FALSE)


###################################################################################################################################

# Grassland ONLY richness

grass_rich_arable1 <- glmer(grass_only_rich ~ prop_arable.500 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(grass_rich_arable1) 
r.squaredGLMM(grass_rich_arable1) # 0.008841925          
grass_rich_arable2 <- glmer(grass_only_rich ~ prop_arable.1500 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(grass_rich_arable2)              
r.squaredGLMM(grass_rich_arable2) # 0.002523597         
grass_rich_arable3 <- glmer(grass_only_rich ~ prop_arable.3000 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(grass_rich_arable3)                          
r.squaredGLMM(grass_rich_arable3) # 0.01182436        
#### ARABLE 3000M #### 


grass_rich_broad1 <- glmer(grass_only_rich ~ prop_broadleaf.500 + (1|subcompartment) + (1|year),
                           data=moth_hab_data, family = poisson())
summary(grass_rich_broad1) 
r.squaredGLMM(grass_rich_broad1) # 0.02759025          
grass_rich_broad2 <- glmer(grass_only_rich ~ prop_broadleaf.1500 + (1|subcompartment) + (1|year),
                           data=moth_hab_data, family = poisson())
summary(grass_rich_broad2)       
r.squaredGLMM(grass_rich_broad2) # 0.03677629          
grass_rich_broad3 <- glmer(grass_only_rich ~ prop_broadleaf.3000 + (1|subcompartment) + (1|year),
                           data=moth_hab_data, family = poisson())
summary(grass_rich_broad3)               
r.squaredGLMM(grass_rich_broad3) # 0.03548865         
#### BROADLEAF 1500M #### 


grass_rich_conifer1 <- glmer(grass_only_rich ~ prop_conifer.500 + (1|subcompartment) + (1|year),
                             data=moth_hab_data, family = poisson())
summary(grass_rich_conifer1)                           
r.squaredGLMM(grass_rich_conifer1) # 0.011861686           
grass_rich_conifer2 <- glmer(grass_only_rich ~ prop_conifer.1500 + (1|subcompartment) + (1|year),
                             data=moth_hab_data, family = poisson())
summary(grass_rich_conifer2)                      
r.squaredGLMM(grass_rich_conifer2) # 0.00006236446           
grass_rich_conifer3 <- glmer(grass_only_rich ~ prop_conifer.3000 + (1|subcompartment) + (1|year),
                             data=moth_hab_data, family = poisson())
summary(grass_rich_conifer3)                                    
r.squaredGLMM(grass_rich_conifer3) # 0.002167641          
#### CONIFER 500M #### 


grass_rich_grass1 <- glmer(grass_only_rich ~ prop_grassland.500 + (1|subcompartment) + (1|year),
                           data=moth_hab_data, family = poisson())
summary(grass_rich_grass1)                             
r.squaredGLMM(grass_rich_grass1) # 0.00007288708            
grass_rich_grass2 <- glmer(grass_only_rich ~ prop_grassland.1500 + (1|subcompartment) + (1|year),
                           data=moth_hab_data, family = poisson())
summary(grass_rich_grass2)                              
r.squaredGLMM(grass_rich_grass2) # 0.009194530            
grass_rich_grass3 <- glmer(grass_only_rich ~ prop_grassland.3000 + (1|subcompartment) + (1|year),
                           data=moth_hab_data, family = poisson())
summary(grass_rich_grass3)                   
r.squaredGLMM(grass_rich_grass3) # 0.01915639           
#### GRASSLAND 3000M #### 

# save results
r2 <- c(r.squaredGLMM(grass_rich_arable1)[1,1], r.squaredGLMM(grass_rich_arable2)[1,1], r.squaredGLMM(grass_rich_arable3)[1,1],
        r.squaredGLMM(grass_rich_broad1)[1,1], r.squaredGLMM(grass_rich_broad2)[1,1], r.squaredGLMM(grass_rich_broad3)[1,1],
        r.squaredGLMM(grass_rich_conifer1)[1,1], r.squaredGLMM(grass_rich_conifer2)[1,1], r.squaredGLMM(grass_rich_conifer3)[1,1],
        r.squaredGLMM(grass_rich_grass1)[1,1], r.squaredGLMM(grass_rich_grass2)[1,1], r.squaredGLMM(grass_rich_grass3)[1,1])
land_cover <- c("Arable.500", "Arable.1500", "Arable.3000", "Broadleaf.500", "Broadleaf.1500", "Broadleaf.3000",
                "Conifer.500", "Conifer.1500", "Conifer.3000", "Grassland.500", "Grassland.1500", "Grassland.3000")
grass_rich_land <- data.frame(Response=rep("Grassland only moth richness"), Land_cover=land_cover, R2=r2)
write.csv(grass_rich_land, file="Results/Moths/Grassland_only_richness_land_cover.csv", row.names=FALSE)

## Now run models with habitat variables only 
grass_rich_hab <- glmmTMB(grass_only_rich ~ scale(dbh_average) + scale(basal_area) +
                            scale(per_broadleaf_canopy) + scale(complexity_score) + 
                            scale(canopy_openess) + scale(minimum_temperature) + 
                            scale(cloud_cover_.) + scale(dis_to_edge) + 
                            scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ (1|year) + 
                            (1|subcompartment), data = moth_hab_data, family = genpois(link="log"),
                          na.action = "na.fail")
check_collinearity(grass_rich_hab) # low correlations after removing moon cycle
testDispersion(grass_rich_hab) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = grass_rich_hab, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod2 <- glmmTMB(grass_only_rich ~ scale(dbh_average) + scale(basal_area) +
                        scale(per_broadleaf_canopy) + scale(complexity_score) + 
                        scale(canopy_openess) + scale(minimum_temperature) + 
                        scale(cloud_cover_.) + scale(dis_to_edge) + 
                        scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_arable.3000) + 
                        (1|year) + (1|subcompartment), data = moth_hab_data, family = genpois(link="log"),
                      na.action = "na.fail")
check_collinearity(final_mod2) # low correlations after removing moon cycle
testDispersion(final_mod2) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod2, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod3 <- glmmTMB(grass_only_rich ~ scale(dbh_average) + scale(basal_area) +
                        scale(per_broadleaf_canopy) + scale(complexity_score) + 
                        scale(canopy_openess) + scale(minimum_temperature) + 
                        scale(cloud_cover_.) + scale(dis_to_edge) + 
                        scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_broadleaf.1500) + 
                        (1|year) + (1|subcompartment), data = moth_hab_data, family = genpois(link="log"),
                      na.action = "na.fail")
check_collinearity(final_mod3) # low correlations after removing moon cycle
testDispersion(final_mod3) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod3, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod4 <- glmmTMB(grass_only_rich ~ scale(dbh_average) + scale(basal_area) +
                        scale(per_broadleaf_canopy) + scale(complexity_score) + 
                        scale(canopy_openess) + scale(minimum_temperature) + 
                        scale(cloud_cover_.) + scale(dis_to_edge) + 
                        scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_conifer.500) + 
                        (1|year) + (1|subcompartment), data = moth_hab_data, family = genpois(link="log"),
                      na.action = "na.fail")
check_collinearity(final_mod4) # low correlations after removing moon cycle
testDispersion(final_mod4) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod4, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod5 <- glmmTMB(grass_only_rich ~ scale(dbh_average) + scale(basal_area) +
                        scale(per_broadleaf_canopy) + scale(complexity_score) + 
                        scale(canopy_openess) + scale(minimum_temperature) + 
                        scale(cloud_cover_.) + scale(dis_to_edge) + 
                        scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_grassland.3000) + 
                        (1|year) + (1|subcompartment), data = moth_hab_data, family = genpois(link="log"),
                      na.action = "na.fail")
check_collinearity(final_mod5) # low correlations after removing moon cycle
testDispersion(final_mod5) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod5, plot = F)
plot(simulationOutput) ## no assumptions violated

AIC(grass_rich_hab, final_mod2, final_mod3, final_mod4, final_mod5)
# habitat model 
summary(grass_rich_hab) # basal area (positive)

## save all AIC values for supplementary table
AIC_values <- data.frame(model=c("habitat", "arable 3000m", "broadleaf 1500m", "conifer 500m", "grassland 3000m"),
                         AIC=AIC(grass_rich_hab, final_mod2, final_mod3, final_mod4, final_mod5)$AIC, stringsAsFactors=FALSE) 
write.csv(AIC_values, file="Results/Grassland_only_richness_AIC.csv", row.names=FALSE)

## save lowest AIC model summary
grass_rich_hab_sum <- as.data.frame(summary(grass_rich_hab)$coefficients$cond) #selecting full model coefficient averages
grass_rich_hab_sum$parameters <- row.names(grass_rich_hab_sum)
row.names(grass_rich_hab_sum) <- 1:nrow(grass_rich_hab_sum)
grass_rich_hab_sum$model <- "habitat"
grass_rich_hab_sum$AIC <- AIC(grass_rich_hab)
write.csv(grass_rich_hab_sum, file="Results/Grassland_only_richness_habitat.csv", row.names=FALSE)

###################################################################################################################################

# Total abundance

tot_abund_arable1 <- glmer(tot_abund ~ prop_arable.500 + (1|subcompartment) + (1|year),
                           data=moth_hab_data, family = poisson())
summary(tot_abund_arable1)       
r.squaredGLMM(tot_abund_arable1) # 0.2600335           
tot_abund_arable2 <- glmer(tot_abund ~ prop_arable.1500 + (1|subcompartment) + (1|year),
                           data=moth_hab_data, family = poisson())
summary(tot_abund_arable2)                     
r.squaredGLMM(tot_abund_arable2) # 0.07852662         
tot_abund_arable3 <- glmer(tot_abund ~ prop_arable.3000 + (1|subcompartment) + (1|year),
                           data=moth_hab_data, family = poisson())
summary(tot_abund_arable3)                       
r.squaredGLMM(tot_abund_arable3) # 0.01715200        
#### ARABLE 500M #### (1500m previously)

tot_abund_broad1 <- glmer(tot_abund ~ prop_broadleaf.500 + (1|subcompartment) + (1|year),
                          data=moth_hab_data, family = poisson())
summary(tot_abund_broad1)     
r.squaredGLMM(tot_abund_broad1) # 0.1689073          
tot_abund_broad2 <- glmer(tot_abund ~ prop_broadleaf.1500 + (1|subcompartment) + (1|year),
                          data=moth_hab_data, family = poisson())
summary(tot_abund_broad2)             
r.squaredGLMM(tot_abund_broad2) # 0.07988541          
tot_abund_broad3 <- glmer(tot_abund ~ prop_broadleaf.3000 + (1|subcompartment) + (1|year),
                          data=moth_hab_data, family = poisson())
summary(tot_abund_broad3)                  
r.squaredGLMM(tot_abund_broad3) # 0.1314437         
#### BROADLEAF 500M #### (3000m previously)

tot_abund_conifer1 <- glmer(tot_abund ~ prop_conifer.500 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(tot_abund_conifer1)                              
r.squaredGLMM(tot_abund_conifer1) # 0.0002792449           
tot_abund_conifer2 <- glmer(tot_abund ~ prop_conifer.1500 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(tot_abund_conifer2)                       
r.squaredGLMM(tot_abund_conifer2) # 0.005722575           
tot_abund_conifer3 <- glmer(tot_abund ~ prop_conifer.3000 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(tot_abund_conifer3)                                       
r.squaredGLMM(tot_abund_conifer3) # 0.03932071          
#### CONIFER 3000M #### (500m previously)

tot_abund_grass1 <- glmer(tot_abund ~ prop_grassland.500 + (1|subcompartment) + (1|year),
                          data=moth_hab_data, family = poisson())
summary(tot_abund_grass1)                               
r.squaredGLMM(tot_abund_grass1) # 0.02008794            
tot_abund_grass2 <- glmer(tot_abund ~ prop_grassland.1500 + (1|subcompartment) + (1|year),
                          data=moth_hab_data, family = poisson())
summary(tot_abund_grass2)                                 
r.squaredGLMM(tot_abund_grass2) # 0.1259373            
tot_abund_grass3 <- glmer(tot_abund ~ prop_grassland.3000 + (1|subcompartment) + (1|year),
                          data=moth_hab_data, family = poisson())
summary(tot_abund_grass3)              
r.squaredGLMM(tot_abund_grass3) # 0.004131162           
#### GRASSLAND 1500M #### 

# save results
r2 <- c(r.squaredGLMM(tot_abund_arable1)[1,1], r.squaredGLMM(tot_abund_arable2)[1,1], r.squaredGLMM(tot_abund_arable3)[1,1],
        r.squaredGLMM(tot_abund_broad1)[1,1], r.squaredGLMM(tot_abund_broad2)[1,1], r.squaredGLMM(tot_abund_broad3)[1,1],
        r.squaredGLMM(tot_abund_conifer1)[1,1], r.squaredGLMM(tot_abund_conifer2)[1,1], r.squaredGLMM(tot_abund_conifer3)[1,1],
        r.squaredGLMM(tot_abund_grass1)[1,1], r.squaredGLMM(tot_abund_grass2)[1,1], r.squaredGLMM(tot_abund_grass3)[1,1])
land_cover <- c("Arable.500", "Arable.1500", "Arable.3000", "Broadleaf.500", "Broadleaf.1500", "Broadleaf.3000",
                "Conifer.500", "Conifer.1500", "Conifer.3000", "Grassland.500", "Grassland.1500", "Grassland.3000")
tot_abund_land <- data.frame(Response=rep("Total moth abundance"), Land_cover=land_cover, R2=r2)
write.csv(tot_abund_land, file="Results/Moths/Total_abund_land_cover.csv", row.names=FALSE)


## Now run models with habitat variables only 
tot_abund_hab <- glmer.nb(tot_abund ~ scale(dbh_average) + scale(basal_area) +
                            scale(per_broadleaf_canopy) + scale(complexity_score) + 
                            scale(canopy_openess) + scale(minimum_temperature) + 
                            scale(cloud_cover_.) + scale(dis_to_edge) + 
                            scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ (1|year) + 
                            (1|subcompartment), data = moth_hab_data)
check_collinearity(tot_abund_hab) # low correlations after removing moon cycle
testDispersion(tot_abund_hab) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = tot_abund_hab, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod2 <- glmer.nb(tot_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_arable.500) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data)
check_collinearity(final_mod2) # low correlations after removing moon cycle
testDispersion(final_mod2) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod2, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod3 <- glmer.nb(tot_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_broadleaf.500) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data)
check_collinearity(final_mod3) # low correlations after removing moon cycle
testDispersion(final_mod3) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod3, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod4 <- glmer.nb(tot_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_conifer.3000) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data)
check_collinearity(final_mod4) # low correlations after removing moon cycle
testDispersion(final_mod4) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod4, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod5 <- glmer.nb(tot_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_grassland.1500) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data)
check_collinearity(final_mod5) # low correlations after removing moon cycle
testDispersion(final_mod5) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod5, plot = F)
plot(simulationOutput) ## no assumptions violated

AIC(tot_abund_hab, final_mod2, final_mod3, final_mod4, final_mod5)
# habitat model 
summary(tot_abund_hab) # % broadleaf (positive) and openness (negative)

## save all AIC values for supplementary table
AIC_values <- data.frame(model=c("habitat", "arable 500m", "broadleaf 500m", "conifer 3000m", "grassland 1500m"),
                         AIC=AIC(tot_abund_hab, final_mod2, final_mod3, final_mod4, final_mod5)$AIC, stringsAsFactors=FALSE) 
write.csv(AIC_values, file="Results/Total_abund_AIC.csv", row.names=FALSE)

## save lowest AIC model summary
tot_abund_hab_sum <- as.data.frame(coef(summary(tot_abund_hab))) #selecting full model coefficient averages
tot_abund_hab_sum$parameters <- row.names(tot_abund_hab_sum)
row.names(tot_abund_hab_sum) <- 1:nrow(tot_abund_hab_sum)
tot_abund_hab_sum$model <- "habitat"
tot_abund_hab_sum$AIC <- AIC(tot_abund_hab)
write.csv(tot_abund_hab_sum, file="Results/Total_abund_habitat.csv", row.names=FALSE)

###################################################################################################################################

# Broadleaf abundance

broad_abund_arable1 <- glmer(broadleaf_abund ~ prop_arable.500 + (1|subcompartment) + (1|year),
                             data=moth_hab_data, family = poisson())
summary(broad_abund_arable1)       
r.squaredGLMM(broad_abund_arable1) # 0.1313597            
broad_abund_arable2 <- glmer(broadleaf_abund ~ prop_arable.1500 + (1|subcompartment) + (1|year),
                             data=moth_hab_data, family = poisson())
summary(broad_abund_arable2)               
r.squaredGLMM(broad_abund_arable2) # 0.01237000          
broad_abund_arable3 <- glmer(broadleaf_abund ~ prop_arable.3000 + (1|subcompartment) + (1|year),
                             data=moth_hab_data, family = poisson())
summary(broad_abund_arable3)                    
r.squaredGLMM(broad_abund_arable3) # 0.06248194         
#### ARABLE 500M #### (3000m previously)


broad_abund_broad1 <- glmer(broadleaf_abund ~ prop_broadleaf.500 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(broad_abund_broad1)          
r.squaredGLMM(broad_abund_broad1) # 0.1184017           
broad_abund_broad2 <- glmer(broadleaf_abund ~ prop_broadleaf.1500 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(broad_abund_broad2)            
r.squaredGLMM(broad_abund_broad2) # 0.06963682           
broad_abund_broad3 <- glmer(broadleaf_abund ~ prop_broadleaf.3000 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(broad_abund_broad3)                  
r.squaredGLMM(broad_abund_broad3) # 0.08739877          
#### BROADLEAF 500M #### (3000m previously)

broad_abund_conifer1 <- glmer(broadleaf_abund ~ prop_conifer.500 + (1|subcompartment) + (1|year),
                              data=moth_hab_data, family = poisson())
summary(broad_abund_conifer1)                                    
r.squaredGLMM(broad_abund_conifer1) # 0.008029166            
broad_abund_conifer2 <- glmer(broadleaf_abund ~ prop_conifer.1500 + (1|subcompartment) + (1|year),
                              data=moth_hab_data, family = poisson())
summary(broad_abund_conifer2)                   
r.squaredGLMM(broad_abund_conifer2) # 0.01437321            
broad_abund_conifer3 <- glmer(broadleaf_abund ~ prop_conifer.3000 + (1|subcompartment) + (1|year),
                              data=moth_hab_data, family = poisson())
summary(broad_abund_conifer3)                                               
r.squaredGLMM(broad_abund_conifer3) # 0.00002670306           
#### CONIFER 1500M #### (500m previously)


broad_abund_grass1 <- glmer(broadleaf_abund ~ prop_grassland.500 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(broad_abund_grass1)                               
r.squaredGLMM(broad_abund_grass1) # 0.01769671             
broad_abund_grass2 <- glmer(broadleaf_abund ~ prop_grassland.1500 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(broad_abund_grass2)                                      
r.squaredGLMM(broad_abund_grass2) # 0.01209138             
broad_abund_grass3 <- glmer(broadleaf_abund ~ prop_grassland.3000 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(broad_abund_grass3)           
r.squaredGLMM(broad_abund_grass3) # 0.1191595            
#### GRASSLAND 3000M #### 

# save results
r2 <- c(r.squaredGLMM(broad_abund_arable1)[1,1], r.squaredGLMM(broad_abund_arable2)[1,1], r.squaredGLMM(broad_abund_arable3)[1,1],
        r.squaredGLMM(broad_abund_broad1)[1,1], r.squaredGLMM(broad_abund_broad2)[1,1], r.squaredGLMM(broad_abund_broad3)[1,1],
        r.squaredGLMM(broad_abund_conifer1)[1,1], r.squaredGLMM(broad_abund_conifer2)[1,1], r.squaredGLMM(broad_abund_conifer3)[1,1],
        r.squaredGLMM(broad_abund_grass1)[1,1], r.squaredGLMM(broad_abund_grass2)[1,1], r.squaredGLMM(broad_abund_grass3)[1,1])
land_cover <- c("Arable.500", "Arable.1500", "Arable.3000", "Broadleaf.500", "Broadleaf.1500", "Broadleaf.3000",
                "Conifer.500", "Conifer.1500", "Conifer.3000", "Grassland.500", "Grassland.1500", "Grassland.3000")
broad_abund_land <- data.frame(Response=rep("Broadleaf moth abundance"), Land_cover=land_cover, R2=r2)
write.csv(broad_abund_land, file="Results/Moths/Broadleaf_abund_land_cover.csv", row.names=FALSE)


## Now run models with habitat variables only 
broad_abund_hab <- glmer.nb(broadleaf_abund ~ scale(dbh_average) + scale(basal_area) +
                              scale(per_broadleaf_canopy) + scale(complexity_score) + 
                              scale(canopy_openess) + scale(minimum_temperature) + 
                              scale(cloud_cover_.) + scale(dis_to_edge) + 
                              scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ (1|year) + 
                              (1|subcompartment), data = moth_hab_data)
check_collinearity(broad_abund_hab) # low correlations after removing moon cycle
testDispersion(broad_abund_hab) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = broad_abund_hab, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod2 <- glmer.nb(broadleaf_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_arable.500) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data)
check_collinearity(final_mod2) # low correlations after removing moon cycle
testDispersion(final_mod2) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod2, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod3 <- glmer.nb(broadleaf_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_broadleaf.500) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data)
check_collinearity(final_mod3) # low correlations after removing moon cycle
testDispersion(final_mod3) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod3, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod4 <- glmer.nb(broadleaf_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_conifer.1500) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data)
check_collinearity(final_mod4) # low correlations after removing moon cycle
testDispersion(final_mod4) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod4, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod5 <- glmer.nb(broadleaf_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_grassland.3000) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data)
check_collinearity(final_mod5) # low correlations after removing moon cycle
testDispersion(final_mod5) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod5, plot = F)
plot(simulationOutput) ## no assumptions violated

AIC(broad_abund_hab, final_mod2, final_mod3, final_mod4, final_mod5)
# habitat model 
summary(broad_abund_hab) # % broadleaf (positive) and complexity score (positive)

## save all AIC values for supplementary table
AIC_values <- data.frame(model=c("habitat", "arable 500m", "broadleaf 500m", "conifer 1500m", "grassland 3000m"),
                         AIC=AIC(broad_abund_hab, final_mod2, final_mod3, final_mod4, final_mod5)$AIC, stringsAsFactors=FALSE) 
write.csv(AIC_values, file="Results/Broadleaf_abund_AIC.csv", row.names=FALSE)

## save lowest AIC model summary
broad_abund_hab_sum <- as.data.frame(coef(summary(broad_abund_hab))) #selecting full model coefficient averages
broad_abund_hab_sum$parameters <- row.names(broad_abund_hab_sum)
row.names(broad_abund_hab_sum) <- 1:nrow(broad_abund_hab_sum)
broad_abund_hab_sum$model <- "habitat"
broad_abund_hab_sum$AIC <- AIC(broad_abund_hab)
write.csv(broad_abund_hab_sum, file="Results/Broadleaf_abund_habitat.csv", row.names=FALSE)


###################################################################################################################################

# Broadleaf ONLY abundance

broad_abund_arable1 <- glmer(broadleaf_only_abund ~ prop_arable.500 + (1|subcompartment) + (1|year),
                             data=moth_hab_data, family = poisson())
summary(broad_abund_arable1)       
r.squaredGLMM(broad_abund_arable1) # 0.05062779             
broad_abund_arable2 <- glmer(broadleaf_only_abund ~ prop_arable.1500 + (1|subcompartment) + (1|year),
                             data=moth_hab_data, family = poisson())
summary(broad_abund_arable2)               
r.squaredGLMM(broad_abund_arable2) # 0.003363585           
broad_abund_arable3 <- glmer(broadleaf_only_abund ~ prop_arable.3000 + (1|subcompartment) + (1|year),
                             data=moth_hab_data, family = poisson())
summary(broad_abund_arable3)                    
r.squaredGLMM(broad_abund_arable3) # 0.03856640          
#### ARABLE 500M #### 


broad_abund_broad1 <- glmer(broadleaf_only_abund ~ prop_broadleaf.500 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(broad_abund_broad1)          
r.squaredGLMM(broad_abund_broad1) # 0.0007595397            
broad_abund_broad2 <- glmer(broadleaf_only_abund ~ prop_broadleaf.1500 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(broad_abund_broad2)            
r.squaredGLMM(broad_abund_broad2) # 0.002027030            
broad_abund_broad3 <- glmer(broadleaf_only_abund ~ prop_broadleaf.3000 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(broad_abund_broad3)                  
r.squaredGLMM(broad_abund_broad3) # 0.000007670596           
#### BROADLEAF 1500M #### 

broad_abund_conifer1 <- glmer(broadleaf_only_abund ~ prop_conifer.500 + (1|subcompartment) + (1|year),
                              data=moth_hab_data, family = poisson())
summary(broad_abund_conifer1)                                    
r.squaredGLMM(broad_abund_conifer1) # 0.001079006             
broad_abund_conifer2 <- glmer(broadleaf_only_abund ~ prop_conifer.1500 + (1|subcompartment) + (1|year),
                              data=moth_hab_data, family = poisson())
summary(broad_abund_conifer2)                   
r.squaredGLMM(broad_abund_conifer2) # 0.00001487165             
broad_abund_conifer3 <- glmer(broadleaf_only_abund ~ prop_conifer.3000 + (1|subcompartment) + (1|year),
                              data=moth_hab_data, family = poisson())
summary(broad_abund_conifer3)                                               
r.squaredGLMM(broad_abund_conifer3) # 0.02053437            
#### CONIFER 3000M #### 


broad_abund_grass1 <- glmer(broadleaf_only_abund ~ prop_grassland.500 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(broad_abund_grass1)                               
r.squaredGLMM(broad_abund_grass1) # 0.01712169              
broad_abund_grass2 <- glmer(broadleaf_only_abund ~ prop_grassland.1500 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(broad_abund_grass2)                                      
r.squaredGLMM(broad_abund_grass2) # 0.002246130              
broad_abund_grass3 <- glmer(broadleaf_only_abund ~ prop_grassland.3000 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(broad_abund_grass3)           
r.squaredGLMM(broad_abund_grass3) # 0.03370195             
#### GRASSLAND 3000M #### 

# save results
r2 <- c(r.squaredGLMM(broad_abund_arable1)[1,1], r.squaredGLMM(broad_abund_arable2)[1,1], r.squaredGLMM(broad_abund_arable3)[1,1],
        r.squaredGLMM(broad_abund_broad1)[1,1], r.squaredGLMM(broad_abund_broad2)[1,1], r.squaredGLMM(broad_abund_broad3)[1,1],
        r.squaredGLMM(broad_abund_conifer1)[1,1], r.squaredGLMM(broad_abund_conifer2)[1,1], r.squaredGLMM(broad_abund_conifer3)[1,1],
        r.squaredGLMM(broad_abund_grass1)[1,1], r.squaredGLMM(broad_abund_grass2)[1,1], r.squaredGLMM(broad_abund_grass3)[1,1])
land_cover <- c("Arable.500", "Arable.1500", "Arable.3000", "Broadleaf.500", "Broadleaf.1500", "Broadleaf.3000",
                "Conifer.500", "Conifer.1500", "Conifer.3000", "Grassland.500", "Grassland.1500", "Grassland.3000")
broad_abund_land <- data.frame(Response=rep("Broadleaf only moth abundance"), Land_cover=land_cover, R2=r2)
write.csv(broad_abund_land, file="Results/Moths/Broadleaf_only_abund_land_cover.csv", row.names=FALSE)


## Now run models with habitat variables only 
broad_abund_hab <- glmer.nb(broadleaf_only_abund ~ scale(dbh_average) + scale(basal_area) +
                              scale(per_broadleaf_canopy) + scale(complexity_score) + 
                              scale(canopy_openess) + scale(minimum_temperature) + 
                              scale(cloud_cover_.) + scale(dis_to_edge) + 
                              scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ (1|year) + 
                              (1|subcompartment), data = moth_hab_data)
check_collinearity(broad_abund_hab) # low correlations after removing moon cycle
testDispersion(broad_abund_hab) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = broad_abund_hab, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod2 <- glmer.nb(broadleaf_only_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_arable.500) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data)
check_collinearity(final_mod2) # low correlations after removing moon cycle
testDispersion(final_mod2) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod2, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod3 <- glmer.nb(broadleaf_only_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_broadleaf.1500) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data)
check_collinearity(final_mod3) # low correlations after removing moon cycle
testDispersion(final_mod3) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod3, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod4 <- glmer.nb(broadleaf_only_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_conifer.3000) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data)
check_collinearity(final_mod4) # low correlations after removing moon cycle
testDispersion(final_mod4) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod4, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod5 <- glmer.nb(broadleaf_only_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_grassland.3000) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data)
check_collinearity(final_mod5) # low correlations after removing moon cycle
testDispersion(final_mod5) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod5, plot = F)
plot(simulationOutput) ## no assumptions violated

AIC(broad_abund_hab, final_mod2, final_mod3, final_mod4, final_mod5)
# habitat model 
summary(broad_abund_hab) # no significant structure variables

## save all AIC values for supplementary table
AIC_values <- data.frame(model=c("habitat", "arable 500m", "broadleaf 1500m", "conifer 3000m", "grassland 3000m"),
                         AIC=AIC(broad_abund_hab, final_mod2, final_mod3, final_mod4, final_mod5)$AIC, stringsAsFactors=FALSE) 
write.csv(AIC_values, file="Results/Broadleaf_only_abund_AIC.csv", row.names=FALSE)

## save lowest AIC model summary
broad_abund_hab_sum <- as.data.frame(coef(summary(broad_abund_hab))) #selecting full model coefficient averages
broad_abund_hab_sum$parameters <- row.names(broad_abund_hab_sum)
row.names(broad_abund_hab_sum) <- 1:nrow(broad_abund_hab_sum)
broad_abund_hab_sum$model <- "habitat"
broad_abund_hab_sum$AIC <- AIC(broad_abund_hab)
write.csv(broad_abund_hab_sum, file="Results/Broadleaf_only_abund_habitat.csv", row.names=FALSE)

###################################################################################################################################

# Conifer abundance

conifer_abund_arable1 <- glmer(conifer_abund ~ prop_arable.500 + (1|subcompartment) + (1|year),
                               data=moth_hab_data, family = poisson())
summary(conifer_abund_arable1)     
r.squaredGLMM(conifer_abund_arable1) # 0.1714292             
conifer_abund_arable2 <- glmer(conifer_abund ~ prop_arable.1500 + (1|subcompartment) + (1|year),
                               data=moth_hab_data, family = poisson())
summary(conifer_abund_arable2)                 
r.squaredGLMM(conifer_abund_arable2) # 0.01309323           
conifer_abund_arable3 <- glmer(conifer_abund ~ prop_arable.3000 + (1|subcompartment) + (1|year),
                               data=moth_hab_data, family = poisson())
summary(conifer_abund_arable3)                           
r.squaredGLMM(conifer_abund_arable3) # 0.000006308933          
#### ARABLE 500M #### (1500m previously)

conifer_abund_broad1 <- glmer(conifer_abund ~ prop_broadleaf.500 + (1|subcompartment) + (1|year),
                              data=moth_hab_data, family = poisson())
summary(conifer_abund_broad1)        
r.squaredGLMM(conifer_abund_broad1) # 0.2951056            
conifer_abund_broad2 <- glmer(conifer_abund ~ prop_broadleaf.1500 + (1|subcompartment) + (1|year),
                              data=moth_hab_data, family = poisson())
summary(conifer_abund_broad2)          
r.squaredGLMM(conifer_abund_broad2) # 0.1385571            
conifer_abund_broad3 <- glmer(conifer_abund ~ prop_broadleaf.3000 + (1|subcompartment) + (1|year),
                              data=moth_hab_data, family = poisson())
summary(conifer_abund_broad3)               
r.squaredGLMM(conifer_abund_broad3) # 0.05779810           
#### BROADLEAF 500M #### 

conifer_abund_conifer1 <- glmer(conifer_abund ~ prop_conifer.500 + (1|subcompartment) + (1|year),
                                data=moth_hab_data, family = poisson())
summary(conifer_abund_conifer1)                                       
r.squaredGLMM(conifer_abund_conifer1) # 0.04188325             
conifer_abund_conifer2 <- glmer(conifer_abund ~ prop_conifer.1500 + (1|subcompartment) + (1|year),
                                data=moth_hab_data, family = poisson())
summary(conifer_abund_conifer2)                   
r.squaredGLMM(conifer_abund_conifer2) # 0.0009349624             
conifer_abund_conifer3 <- glmer(conifer_abund ~ prop_conifer.3000 + (1|subcompartment) + (1|year),
                                data=moth_hab_data, family = poisson())
summary(conifer_abund_conifer3)                                                    
r.squaredGLMM(conifer_abund_conifer3) # 0.0004801142            
#### CONIFER 500M #### 

conifer_abund_grass1 <- glmer(conifer_abund ~ prop_grassland.500 + (1|subcompartment) + (1|year),
                              data=moth_hab_data, family = poisson())
summary(conifer_abund_grass1)                                 
r.squaredGLMM(conifer_abund_grass1) # 0.1029660              
conifer_abund_grass2 <- glmer(conifer_abund ~ prop_grassland.1500 + (1|subcompartment) + (1|year),
                              data=moth_hab_data, family = poisson())
summary(conifer_abund_grass2)                                            
r.squaredGLMM(conifer_abund_grass2) # 0.0005737802              
conifer_abund_grass3 <- glmer(conifer_abund ~ prop_grassland.3000 + (1|subcompartment) + (1|year),
                              data=moth_hab_data, family = poisson())
summary(conifer_abund_grass3)        
r.squaredGLMM(conifer_abund_grass3) # 0.007687833             
#### GRASSLAND 500M #### 

# save results
r2 <- c(r.squaredGLMM(conifer_abund_arable1)[1,1], r.squaredGLMM(conifer_abund_arable2)[1,1], r.squaredGLMM(conifer_abund_arable3)[1,1],
        r.squaredGLMM(conifer_abund_broad1)[1,1], r.squaredGLMM(conifer_abund_broad2)[1,1], r.squaredGLMM(conifer_abund_broad3)[1,1],
        r.squaredGLMM(conifer_abund_conifer1)[1,1], r.squaredGLMM(conifer_abund_conifer2)[1,1], r.squaredGLMM(conifer_abund_conifer3)[1,1],
        r.squaredGLMM(conifer_abund_grass1)[1,1], r.squaredGLMM(conifer_abund_grass2)[1,1], r.squaredGLMM(conifer_abund_grass3)[1,1])
land_cover <- c("Arable.500", "Arable.1500", "Arable.3000", "Broadleaf.500", "Broadleaf.1500", "Broadleaf.3000",
                "Conifer.500", "Conifer.1500", "Conifer.3000", "Grassland.500", "Grassland.1500", "Grassland.3000")
conifer_abund_land <- data.frame(Response=rep("Conifer moth abundance"), Land_cover=land_cover, R2=r2)
write.csv(conifer_abund_land, file="Results/Moths/Conifer_abund_land_cover.csv", row.names=FALSE)

## Now run models with habitat variables only 
conifer_abund_hab <- glmer.nb(conifer_abund ~ scale(dbh_average) + scale(basal_area) +
                                scale(per_broadleaf_canopy) + scale(complexity_score) + 
                                scale(canopy_openess) + scale(minimum_temperature) + 
                                scale(cloud_cover_.) + scale(dis_to_edge) + 
                                scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ (1|year) + 
                                (1|subcompartment), data = moth_hab_data, 
                              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
check_collinearity(conifer_abund_hab) # low correlations after removing moon cycle
testDispersion(conifer_abund_hab) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = conifer_abund_hab, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod2 <- glmer.nb(conifer_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_arable.500) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data, 
                       control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
check_collinearity(final_mod2) # low correlations after removing moon cycle
testDispersion(final_mod2) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod2, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod3 <- glmer.nb(conifer_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_broadleaf.500) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data,
                       control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
check_collinearity(final_mod3) # low correlations after removing moon cycle
testDispersion(final_mod3) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod3, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod4 <- glmer.nb(conifer_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_conifer.500) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data,
                       control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
check_collinearity(final_mod4) # low correlations after removing moon cycle
testDispersion(final_mod4) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod4, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod5 <- glmer.nb(conifer_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_grassland.500) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data,
                       control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
check_collinearity(final_mod5) # low correlations after removing moon cycle
testDispersion(final_mod5) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod5, plot = F)
plot(simulationOutput) ## no assumptions violated

AIC(conifer_abund_hab, final_mod2, final_mod3, final_mod4, final_mod5)
# arable model 
summary(final_mod2) # % broadleaf (positive) 

## save all AIC values for supplementary table
AIC_values <- data.frame(model=c("habitat", "arable 500m", "broadleaf 500m", "conifer 500m", "grassland 500m"),
                         AIC=AIC(conifer_abund_hab, final_mod2, final_mod3, final_mod4, final_mod5)$AIC, stringsAsFactors=FALSE) 
write.csv(AIC_values, file="Results/Conifer_abund_AIC.csv", row.names=FALSE)

## save lowest AIC model summary
final_mod2_sum <- as.data.frame(coef(summary(final_mod2))) 
final_mod2_sum$parameters <- row.names(final_mod2_sum)
row.names(final_mod2_sum) <- 1:nrow(final_mod2_sum)
final_mod2_sum$model <- "arable 500m"
final_mod2_sum$AIC <- AIC(final_mod2)
write.csv(final_mod2_sum, file="Results/Conifer_abund_habitat.csv", row.names=FALSE)


###################################################################################################################################

# Conifer ONLY abundance

conifer_abund_arable1 <- glmer(conifer_only_abund ~ prop_arable.500 + (1|subcompartment) + (1|year),
                               data=moth_hab_data, family = poisson())
summary(conifer_abund_arable1)     
r.squaredGLMM(conifer_abund_arable1) # 0.1488633              
conifer_abund_arable2 <- glmer(conifer_only_abund ~ prop_arable.1500 + (1|subcompartment) + (1|year),
                               data=moth_hab_data, family = poisson())
summary(conifer_abund_arable2)                 
r.squaredGLMM(conifer_abund_arable2) # 0.004906986            
conifer_abund_arable3 <- glmer(conifer_only_abund ~ prop_arable.3000 + (1|subcompartment) + (1|year),
                               data=moth_hab_data, family = poisson())
summary(conifer_abund_arable3)                           
r.squaredGLMM(conifer_abund_arable3) # 0.01331744           
#### ARABLE 500M #### 

conifer_abund_broad1 <- glmer(conifer_only_abund ~ prop_broadleaf.500 + (1|subcompartment) + (1|year),
                              data=moth_hab_data, family = poisson())
summary(conifer_abund_broad1)        
r.squaredGLMM(conifer_abund_broad1) # 0.1309578             
conifer_abund_broad2 <- glmer(conifer_only_abund ~ prop_broadleaf.1500 + (1|subcompartment) + (1|year),
                              data=moth_hab_data, family = poisson())
summary(conifer_abund_broad2)          
r.squaredGLMM(conifer_abund_broad2) # 0.04672491             
conifer_abund_broad3 <- glmer(conifer_only_abund ~ prop_broadleaf.3000 + (1|subcompartment) + (1|year),
                              data=moth_hab_data, family = poisson())
summary(conifer_abund_broad3)               
r.squaredGLMM(conifer_abund_broad3) # 0.06799532            
#### BROADLEAF 500M #### 

conifer_abund_conifer1 <- glmer(conifer_only_abund ~ prop_conifer.500 + (1|subcompartment) + (1|year),
                                data=moth_hab_data, family = poisson())
summary(conifer_abund_conifer1)                                       
r.squaredGLMM(conifer_abund_conifer1) # 0.008288149              
conifer_abund_conifer2 <- glmer(conifer_only_abund ~ prop_conifer.1500 + (1|subcompartment) + (1|year),
                                data=moth_hab_data, family = poisson())
summary(conifer_abund_conifer2)                   
r.squaredGLMM(conifer_abund_conifer2) # 0.01422989              
conifer_abund_conifer3 <- glmer(conifer_only_abund ~ prop_conifer.3000 + (1|subcompartment) + (1|year),
                                data=moth_hab_data, family = poisson())
summary(conifer_abund_conifer3)                                                    
r.squaredGLMM(conifer_abund_conifer3) # 0.01614609             
#### CONIFER 3000M #### 

conifer_abund_grass1 <- glmer(conifer_only_abund ~ prop_grassland.500 + (1|subcompartment) + (1|year),
                              data=moth_hab_data, family = poisson())
summary(conifer_abund_grass1)                                 
r.squaredGLMM(conifer_abund_grass1) # 0.0007032253               
conifer_abund_grass2 <- glmer(conifer_only_abund ~ prop_grassland.1500 + (1|subcompartment) + (1|year),
                              data=moth_hab_data, family = poisson())
summary(conifer_abund_grass2)                                            
r.squaredGLMM(conifer_abund_grass2) # 0.002594926               
conifer_abund_grass3 <- glmer(conifer_only_abund ~ prop_grassland.3000 + (1|subcompartment) + (1|year),
                              data=moth_hab_data, family = poisson())
summary(conifer_abund_grass3)        
r.squaredGLMM(conifer_abund_grass3) # 0.03585078              
#### GRASSLAND 3000M #### 

# save results
r2 <- c(r.squaredGLMM(conifer_abund_arable1)[1,1], r.squaredGLMM(conifer_abund_arable2)[1,1], r.squaredGLMM(conifer_abund_arable3)[1,1],
        r.squaredGLMM(conifer_abund_broad1)[1,1], r.squaredGLMM(conifer_abund_broad2)[1,1], r.squaredGLMM(conifer_abund_broad3)[1,1],
        r.squaredGLMM(conifer_abund_conifer1)[1,1], r.squaredGLMM(conifer_abund_conifer2)[1,1], r.squaredGLMM(conifer_abund_conifer3)[1,1],
        r.squaredGLMM(conifer_abund_grass1)[1,1], r.squaredGLMM(conifer_abund_grass2)[1,1], r.squaredGLMM(conifer_abund_grass3)[1,1])
land_cover <- c("Arable.500", "Arable.1500", "Arable.3000", "Broadleaf.500", "Broadleaf.1500", "Broadleaf.3000",
                "Conifer.500", "Conifer.1500", "Conifer.3000", "Grassland.500", "Grassland.1500", "Grassland.3000")
conifer_abund_land <- data.frame(Response=rep("Conifer only moth abundance"), Land_cover=land_cover, R2=r2)
write.csv(conifer_abund_land, file="Results/Moths/Conifer_only_abund_land_cover.csv", row.names=FALSE)

## Now run models with habitat variables only 
conifer_abund_hab <- glmer.nb(conifer_only_abund ~ scale(dbh_average) + scale(basal_area) +
                                scale(per_broadleaf_canopy) + scale(complexity_score) + 
                                scale(canopy_openess) + scale(minimum_temperature) + 
                                scale(cloud_cover_.) + scale(dis_to_edge) + 
                                scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ (1|year) + 
                                (1|subcompartment), data = moth_hab_data, 
                              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
check_collinearity(conifer_abund_hab) # low correlations after removing moon cycle
testDispersion(conifer_abund_hab) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = conifer_abund_hab, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod2 <- glmer.nb(conifer_only_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_arable.500) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data, 
                       control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
check_collinearity(final_mod2) # low correlations after removing moon cycle
testDispersion(final_mod2) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod2, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod3 <- glmer.nb(conifer_only_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_broadleaf.500) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data,
                       control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
check_collinearity(final_mod3) # low correlations after removing moon cycle
testDispersion(final_mod3) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod3, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod4 <- glmer.nb(conifer_only_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_conifer.3000) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data,
                       control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
check_collinearity(final_mod4) # low correlations after removing moon cycle
testDispersion(final_mod4) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod4, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod5 <- glmer.nb(conifer_only_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_grassland.3000) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data,
                       control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
check_collinearity(final_mod5) # low correlations after removing moon cycle
testDispersion(final_mod5) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod5, plot = F)
plot(simulationOutput) ## no assumptions violated

AIC(conifer_abund_hab, final_mod2, final_mod3, final_mod4, final_mod5)
# arable model 
summary(final_mod2) # dbh (negative) and openness (negative)

## save all AIC values for supplementary table
AIC_values <- data.frame(model=c("habitat", "arable 500m", "broadleaf 500m", "conifer 3000m", "grassland 3000m"),
                         AIC=AIC(conifer_abund_hab, final_mod2, final_mod3, final_mod4, final_mod5)$AIC, stringsAsFactors=FALSE) 
write.csv(AIC_values, file="Results/Conifer_only_abund_AIC.csv", row.names=FALSE)

## save lowest AIC model summary
final_mod2_sum <- as.data.frame(coef(summary(final_mod2))) 
final_mod2_sum$parameters <- row.names(final_mod2_sum)
row.names(final_mod2_sum) <- 1:nrow(final_mod2_sum)
final_mod2_sum$model <- "arable 500m"
final_mod2_sum$AIC <- AIC(final_mod2)
write.csv(final_mod2_sum, file="Results/Conifer_only_abund_habitat.csv", row.names=FALSE)


###################################################################################################################################

# Woodland abundance

wood_abund_arable1 <- glmer(woodland_abund ~ prop_arable.500 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(wood_abund_arable1)     
r.squaredGLMM(wood_abund_arable1) # 0.2083316              
wood_abund_arable2 <- glmer(woodland_abund ~ prop_arable.1500 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(wood_abund_arable2)                
r.squaredGLMM(wood_abund_arable2) # 0.04138059            
wood_abund_arable3 <- glmer(woodland_abund ~ prop_arable.3000 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(wood_abund_arable3)                                 
r.squaredGLMM(wood_abund_arable3) # 0.001923456           
#### ARABLE 500M #### (1500m previously)

wood_abund_broad1 <- glmer(woodland_abund ~ prop_broadleaf.500 + (1|subcompartment) + (1|year),
                           data=moth_hab_data, family = poisson())
summary(wood_abund_broad1)        
r.squaredGLMM(wood_abund_broad1) # 0.08863381             
wood_abund_broad2 <- glmer(woodland_abund ~ prop_broadleaf.1500 + (1|subcompartment) + (1|year),
                           data=moth_hab_data, family = poisson())
summary(wood_abund_broad2)           
r.squaredGLMM(wood_abund_broad2) # 0.07959963             
wood_abund_broad3 <- glmer(woodland_abund ~ prop_broadleaf.3000 + (1|subcompartment) + (1|year),
                           data=moth_hab_data, family = poisson())
summary(wood_abund_broad3)                 
r.squaredGLMM(wood_abund_broad3) # 0.1011072            
#### BROADLEAF 3000M #### (1500m previously)

wood_abund_conifer1 <- glmer(woodland_abund ~ prop_conifer.500 + (1|subcompartment) + (1|year),
                             data=moth_hab_data, family = poisson())
summary(wood_abund_conifer1)                                              
r.squaredGLMM(wood_abund_conifer1) # 0.006723381              
wood_abund_conifer2 <- glmer(woodland_abund ~ prop_conifer.1500 + (1|subcompartment) + (1|year),
                             data=moth_hab_data, family = poisson())
summary(wood_abund_conifer2)               
r.squaredGLMM(wood_abund_conifer2) # 0.0009393564              
wood_abund_conifer3 <- glmer(woodland_abund ~ prop_conifer.3000 + (1|subcompartment) + (1|year),
                             data=moth_hab_data, family = poisson())
summary(wood_abund_conifer3)                                                  
r.squaredGLMM(wood_abund_conifer3) # 0.06900126             
#### CONIFER 3000M #### 

wood_abund_grass1 <- glmer(woodland_abund ~ prop_grassland.500 + (1|subcompartment) + (1|year),
                           data=moth_hab_data, family = poisson())
summary(wood_abund_grass1)                             
r.squaredGLMM(wood_abund_grass1) # 0.01633504               
wood_abund_grass2 <- glmer(woodland_abund ~ prop_grassland.1500 + (1|subcompartment) + (1|year),
                           data=moth_hab_data, family = poisson())
summary(wood_abund_grass2)                                                
r.squaredGLMM(wood_abund_grass2) # 0.07442880               
wood_abund_grass3 <- glmer(woodland_abund ~ prop_grassland.3000 + (1|subcompartment) + (1|year),
                           data=moth_hab_data, family = poisson())
summary(wood_abund_grass3)     
r.squaredGLMM(wood_abund_grass3) # 0.00005575411              
#### GRASSLAND 1500M #### 

# save results
r2 <- c(r.squaredGLMM(wood_abund_arable1)[1,1], r.squaredGLMM(wood_abund_arable2)[1,1], r.squaredGLMM(wood_abund_arable3)[1,1],
        r.squaredGLMM(wood_abund_broad1)[1,1], r.squaredGLMM(wood_abund_broad2)[1,1], r.squaredGLMM(wood_abund_broad3)[1,1],
        r.squaredGLMM(wood_abund_conifer1)[1,1], r.squaredGLMM(wood_abund_conifer2)[1,1], r.squaredGLMM(wood_abund_conifer3)[1,1],
        r.squaredGLMM(wood_abund_grass1)[1,1], r.squaredGLMM(wood_abund_grass2)[1,1], r.squaredGLMM(wood_abund_grass3)[1,1])
land_cover <- c("Arable.500", "Arable.1500", "Arable.3000", "Broadleaf.500", "Broadleaf.1500", "Broadleaf.3000",
                "Conifer.500", "Conifer.1500", "Conifer.3000", "Grassland.500", "Grassland.1500", "Grassland.3000")
wood_abund_land <- data.frame(Response=rep("Woodland moth abundance"), Land_cover=land_cover, R2=r2)
write.csv(wood_abund_land, file="Results/Moths/Woodland_abund_land_cover.csv", row.names=FALSE)


## Now run models with habitat variables only 
wood_abund_hab <- glmer.nb(woodland_abund ~ scale(dbh_average) + scale(basal_area) +
                             scale(per_broadleaf_canopy) + scale(complexity_score) + 
                             scale(canopy_openess) + scale(minimum_temperature) + 
                             scale(cloud_cover_.) + scale(dis_to_edge) + 
                             scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ (1|year) + 
                             (1|subcompartment), data = moth_hab_data)
check_collinearity(wood_abund_hab) # low correlations after removing moon cycle
testDispersion(wood_abund_hab) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = wood_abund_hab, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod2 <- glmer.nb(woodland_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_arable.500) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data)
check_collinearity(final_mod2) # low correlations after removing moon cycle
testDispersion(final_mod2) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod2, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod3 <- glmer.nb(woodland_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_broadleaf.3000) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data)
check_collinearity(final_mod3) # low correlations after removing moon cycle
testDispersion(final_mod3) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod3, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod4 <- glmer.nb(woodland_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_conifer.3000) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data)
check_collinearity(final_mod4) # low correlations after removing moon cycle
testDispersion(final_mod4) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod4, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod5 <- glmer.nb(woodland_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_grassland.1500) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data)
check_collinearity(final_mod5) # low correlations after removing moon cycle
testDispersion(final_mod5) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod5, plot = F)
plot(simulationOutput) ## no assumptions violated

AIC(wood_abund_hab, final_mod2, final_mod3, final_mod4, final_mod5)
# habitat model 
summary(wood_abund_hab) # % broadleaf canopy (positive) and openness (negative)

## save all AIC values for supplementary table
AIC_values <- data.frame(model=c("habitat", "arable 500m", "broadleaf 3000m", "conifer 3000m", "grassland 1500m"),
                         AIC=AIC(wood_abund_hab, final_mod2, final_mod3, final_mod4, final_mod5)$AIC, stringsAsFactors=FALSE) 
write.csv(AIC_values, file="Results/Woodland_abund_AIC.csv", row.names=FALSE)

## save lowest AIC model summary
wood_abund_hab_sum <- as.data.frame(coef(summary(wood_abund_hab))) #selecting full model coefficient averages
wood_abund_hab_sum$parameters <- row.names(wood_abund_hab_sum)
row.names(wood_abund_hab_sum) <- 1:nrow(wood_abund_hab_sum)
wood_abund_hab_sum$model <- "habitat"
wood_abund_hab_sum$AIC <- AIC(wood_abund_hab)
write.csv(wood_abund_hab_sum, file="Results/Woodland_abund_habitat.csv", row.names=FALSE)


###################################################################################################################################

# Woodland ONLY abundance

wood_abund_arable1 <- glmer(woodland_only_abund ~ prop_arable.500 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(wood_abund_arable1)     
r.squaredGLMM(wood_abund_arable1) # 0.3797046               
wood_abund_arable2 <- glmer(woodland_only_abund ~ prop_arable.1500 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(wood_abund_arable2)                
r.squaredGLMM(wood_abund_arable2) # 0.06319592             
wood_abund_arable3 <- glmer(woodland_only_abund ~ prop_arable.3000 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(wood_abund_arable3)                                 
r.squaredGLMM(wood_abund_arable3) # 0.01911430            
#### ARABLE 500M #### 

wood_abund_broad1 <- glmer(woodland_only_abund ~ prop_broadleaf.500 + (1|subcompartment) + (1|year),
                           data=moth_hab_data, family = poisson())
summary(wood_abund_broad1)        
r.squaredGLMM(wood_abund_broad1) # 0.06248567              
wood_abund_broad2 <- glmer(woodland_only_abund ~ prop_broadleaf.1500 + (1|subcompartment) + (1|year),
                           data=moth_hab_data, family = poisson())
summary(wood_abund_broad2)           
r.squaredGLMM(wood_abund_broad2) # 0.0001638217              
wood_abund_broad3 <- glmer(woodland_only_abund ~ prop_broadleaf.3000 + (1|subcompartment) + (1|year),
                           data=moth_hab_data, family = poisson())
summary(wood_abund_broad3)                 
r.squaredGLMM(wood_abund_broad3) # 0.002042297             
#### BROADLEAF 500M #### 

wood_abund_conifer1 <- glmer(woodland_only_abund ~ prop_conifer.500 + (1|subcompartment) + (1|year),
                             data=moth_hab_data, family = poisson())
summary(wood_abund_conifer1)                                              
r.squaredGLMM(wood_abund_conifer1) # 0.004754402               
wood_abund_conifer2 <- glmer(woodland_only_abund ~ prop_conifer.1500 + (1|subcompartment) + (1|year),
                             data=moth_hab_data, family = poisson())
summary(wood_abund_conifer2)               
r.squaredGLMM(wood_abund_conifer2) # 0.01485690               
wood_abund_conifer3 <- glmer(woodland_only_abund ~ prop_conifer.3000 + (1|subcompartment) + (1|year),
                             data=moth_hab_data, family = poisson())
summary(wood_abund_conifer3)                                                  
r.squaredGLMM(wood_abund_conifer3) # 0.002703671              
#### CONIFER 1500M #### 

wood_abund_grass1 <- glmer(woodland_only_abund ~ prop_grassland.500 + (1|subcompartment) + (1|year),
                           data=moth_hab_data, family = poisson())
summary(wood_abund_grass1)                             
r.squaredGLMM(wood_abund_grass1) # 0.04651139                
wood_abund_grass2 <- glmer(woodland_only_abund ~ prop_grassland.1500 + (1|subcompartment) + (1|year),
                           data=moth_hab_data, family = poisson())
summary(wood_abund_grass2)                                                
r.squaredGLMM(wood_abund_grass2) # 0.1202149                
wood_abund_grass3 <- glmer(woodland_only_abund ~ prop_grassland.3000 + (1|subcompartment) + (1|year),
                           data=moth_hab_data, family = poisson())
summary(wood_abund_grass3)     
r.squaredGLMM(wood_abund_grass3) # 0.009689999               
#### GRASSLAND 1500M #### 

# save results
r2 <- c(r.squaredGLMM(wood_abund_arable1)[1,1], r.squaredGLMM(wood_abund_arable2)[1,1], r.squaredGLMM(wood_abund_arable3)[1,1],
        r.squaredGLMM(wood_abund_broad1)[1,1], r.squaredGLMM(wood_abund_broad2)[1,1], r.squaredGLMM(wood_abund_broad3)[1,1],
        r.squaredGLMM(wood_abund_conifer1)[1,1], r.squaredGLMM(wood_abund_conifer2)[1,1], r.squaredGLMM(wood_abund_conifer3)[1,1],
        r.squaredGLMM(wood_abund_grass1)[1,1], r.squaredGLMM(wood_abund_grass2)[1,1], r.squaredGLMM(wood_abund_grass3)[1,1])
land_cover <- c("Arable.500", "Arable.1500", "Arable.3000", "Broadleaf.500", "Broadleaf.1500", "Broadleaf.3000",
                "Conifer.500", "Conifer.1500", "Conifer.3000", "Grassland.500", "Grassland.1500", "Grassland.3000")
wood_abund_land <- data.frame(Response=rep("Woodland only moth abundance"), Land_cover=land_cover, R2=r2)
write.csv(wood_abund_land, file="Results/Moths/Woodland_only_abund_land_cover.csv", row.names=FALSE)


## Now run models with habitat variables only 
wood_abund_hab <- glmer.nb(woodland_only_abund ~ scale(dbh_average) + scale(basal_area) +
                             scale(per_broadleaf_canopy) + scale(complexity_score) + 
                             scale(canopy_openess) + scale(minimum_temperature) + 
                             scale(cloud_cover_.) + scale(dis_to_edge) + 
                             scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ (1|year) + 
                             (1|subcompartment), data = moth_hab_data)
check_collinearity(wood_abund_hab) # low correlations after removing moon cycle
testDispersion(wood_abund_hab) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = wood_abund_hab, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod2 <- glmer.nb(woodland_only_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_arable.500) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data)
check_collinearity(final_mod2) # low correlations after removing moon cycle
testDispersion(final_mod2) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod2, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod3 <- glmer.nb(woodland_only_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_broadleaf.500) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data)
check_collinearity(final_mod3) # low correlations after removing moon cycle
testDispersion(final_mod3) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod3, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod4 <- glmer.nb(woodland_only_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_conifer.1500) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data)
check_collinearity(final_mod4) # low correlations after removing moon cycle
testDispersion(final_mod4) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod4, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod5 <- glmer.nb(woodland_only_abund ~ scale(dbh_average) + scale(basal_area) +
                         scale(per_broadleaf_canopy) + scale(complexity_score) + 
                         scale(canopy_openess) + scale(minimum_temperature) + 
                         scale(cloud_cover_.) + scale(dis_to_edge) + 
                         scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_grassland.1500) + 
                         (1|year) + (1|subcompartment), data = moth_hab_data)
check_collinearity(final_mod5) # low correlations after removing moon cycle
testDispersion(final_mod5) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod5, plot = F)
plot(simulationOutput) ## no assumptions violated

AIC(wood_abund_hab, final_mod2, final_mod3, final_mod4, final_mod5)
# habitat model 
summary(wood_abund_hab) # no significant structural variables

## save all AIC values for supplementary table
AIC_values <- data.frame(model=c("habitat", "arable 500m", "broadleaf 500m", "conifer 1500m", "grassland 1500m"),
                         AIC=AIC(wood_abund_hab, final_mod2, final_mod3, final_mod4, final_mod5)$AIC, stringsAsFactors=FALSE) 
write.csv(AIC_values, file="Results/Woodland_only_abund_AIC.csv", row.names=FALSE)

## save lowest AIC model summary
wood_abund_hab_sum <- as.data.frame(coef(summary(wood_abund_hab))) #selecting full model coefficient averages
wood_abund_hab_sum$parameters <- row.names(wood_abund_hab_sum)
row.names(wood_abund_hab_sum) <- 1:nrow(wood_abund_hab_sum)
wood_abund_hab_sum$model <- "habitat"
wood_abund_hab_sum$AIC <- AIC(wood_abund_hab)
write.csv(wood_abund_hab_sum, file="Results/Woodland_only_abund_habitat.csv", row.names=FALSE)


###################################################################################################################################

# Grassland abundance

grass_abund_arable1 <- glmer(grass_abund ~ prop_arable.500 + (1|subcompartment) + (1|year),
                             data=moth_hab_data, family = poisson())
summary(grass_abund_arable1)   
r.squaredGLMM(grass_abund_arable1) # 0.09220429               
grass_abund_arable2 <- glmer(grass_abund ~ prop_arable.1500 + (1|subcompartment) + (1|year),
                             data=moth_hab_data, family = poisson())
summary(grass_abund_arable2)            
r.squaredGLMM(grass_abund_arable2) # 0.01702387             
grass_abund_arable3 <- glmer(grass_abund ~ prop_arable.3000 + (1|subcompartment) + (1|year),
                             data=moth_hab_data, family = poisson())
summary(grass_abund_arable3)                                      
r.squaredGLMM(grass_abund_arable3) # 0.005149757            
#### ARABLE 500M #### (1500m previously)

grass_abund_broad1 <- glmer(grass_abund ~ prop_broadleaf.500 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(grass_abund_broad1)     
r.squaredGLMM(grass_abund_broad1) # 0.08039054              
grass_abund_broad2 <- glmer(grass_abund ~ prop_broadleaf.1500 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(grass_abund_broad2)         
r.squaredGLMM(grass_abund_broad2) # 0.02534602              
grass_abund_broad3 <- glmer(grass_abund ~ prop_broadleaf.3000 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(grass_abund_broad3)                 
r.squaredGLMM(grass_abund_broad3) # 0.02755133             
#### BROADLEAF 500M #### (3000m previously)

grass_abund_conifer1 <- glmer(grass_abund ~ prop_conifer.500 + (1|subcompartment) + (1|year),
                              data=moth_hab_data, family = poisson())
summary(grass_abund_conifer1)                                                
r.squaredGLMM(grass_abund_conifer1) # 0.004300661               
grass_abund_conifer2 <- glmer(grass_abund ~ prop_conifer.1500 + (1|subcompartment) + (1|year),
                              data=moth_hab_data, family = poisson())
summary(grass_abund_conifer2)                
r.squaredGLMM(grass_abund_conifer2) # 0.001100853               
grass_abund_conifer3 <- glmer(grass_abund ~ prop_conifer.3000 + (1|subcompartment) + (1|year),
                              data=moth_hab_data, family = poisson())
summary(grass_abund_conifer3)                                                   
r.squaredGLMM(grass_abund_conifer3) # 0.006055737              
#### CONIFER 3000M #### (1500m previously)

grass_abund_grass1 <- glmer(grass_abund ~ prop_grassland.500 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(grass_abund_grass1)                          
r.squaredGLMM(grass_abund_grass1) # 0.01761178                
grass_abund_grass2 <- glmer(grass_abund ~ prop_grassland.1500 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(grass_abund_grass2)                                                      
r.squaredGLMM(grass_abund_grass2) # 0.05059977                
grass_abund_grass3 <- glmer(grass_abund ~ prop_grassland.3000 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(grass_abund_grass3)        
r.squaredGLMM(grass_abund_grass3) # 0.00003665149               
#### GRASSLAND 1500M #### 

# save results
r2 <- c(r.squaredGLMM(grass_abund_arable1)[1,1], r.squaredGLMM(grass_abund_arable2)[1,1], r.squaredGLMM(grass_abund_arable3)[1,1],
        r.squaredGLMM(grass_abund_broad1)[1,1], r.squaredGLMM(grass_abund_broad2)[1,1], r.squaredGLMM(grass_abund_broad3)[1,1],
        r.squaredGLMM(grass_abund_conifer1)[1,1], r.squaredGLMM(grass_abund_conifer2)[1,1], r.squaredGLMM(grass_abund_conifer3)[1,1],
        r.squaredGLMM(grass_abund_grass1)[1,1], r.squaredGLMM(grass_abund_grass2)[1,1], r.squaredGLMM(grass_abund_grass3)[1,1])
land_cover <- c("Arable.500", "Arable.1500", "Arable.3000", "Broadleaf.500", "Broadleaf.1500", "Broadleaf.3000",
                "Conifer.500", "Conifer.1500", "Conifer.3000", "Grassland.500", "Grassland.1500", "Grassland.3000")
grass_abund_land <- data.frame(Response=rep("Grassland moth abundance"), Land_cover=land_cover, R2=r2)
write.csv(grass_abund_land, file="Results/Moths/Grassland_abund_land_cover.csv", row.names=FALSE)

## Now run models with habitat variables only 
grass_abund_hab <- glmer(grass_abund ~ scale(dbh_average) + scale(basal_area) +
                           scale(per_broadleaf_canopy) + scale(complexity_score) + 
                           scale(canopy_openess) + scale(minimum_temperature) + 
                           scale(cloud_cover_.) + scale(dis_to_edge) + 
                           scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ (1|year) + 
                           (1|subcompartment), data = moth_hab_data, family = poisson(),
                         glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                         na.action = "na.fail")
check_collinearity(grass_abund_hab) # low correlations after removing moon cycle
testDispersion(grass_abund_hab) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = grass_abund_hab, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod2 <- glmer(grass_abund ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + 
                      scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_arable.500) + 
                      (1|year) + (1|subcompartment), data = moth_hab_data, family = poisson())
check_collinearity(final_mod2) # low correlations after removing moon cycle
testDispersion(final_mod2) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod2, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod3 <- glmer(grass_abund ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + 
                      scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_broadleaf.500) + 
                      (1|year) + (1|subcompartment), data = moth_hab_data, family = poisson())
check_collinearity(final_mod3) # low correlations after removing moon cycle
testDispersion(final_mod3) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod3, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod4 <- glmer(grass_abund ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + 
                      scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_conifer.3000) + 
                      (1|year) + (1|subcompartment), data = moth_hab_data, family = poisson())
check_collinearity(final_mod4) # low correlations after removing moon cycle
testDispersion(final_mod4) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod4, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod5 <- glmer(grass_abund ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + 
                      scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_grassland.1500) + 
                      (1|year) + (1|subcompartment), data = moth_hab_data, family = poisson())
check_collinearity(final_mod5) # low correlations after removing moon cycle
testDispersion(final_mod5) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod5, plot = F)
plot(simulationOutput) ## no assumptions violated

AIC(grass_abund_hab, final_mod2, final_mod3, final_mod4, final_mod5)
# habitat model 
summary(grass_abund_hab) # basal area (positive) and % broadleaf (positive)

## save all AIC values for supplementary table
AIC_values <- data.frame(model=c("habitat", "arable 500m", "broadleaf 500m", "conifer 3000m", "grassland 1500m"),
                         AIC=AIC(grass_abund_hab, final_mod2, final_mod3, final_mod4, final_mod5)$AIC, stringsAsFactors=FALSE) 
write.csv(AIC_values, file="Results/Grassland_abund_AIC.csv", row.names=FALSE)

## save lowest AIC model summary
grass_abund_hab_sum <- as.data.frame(coef(summary(grass_abund_hab))) 
grass_abund_hab_sum$parameters <- row.names(grass_abund_hab_sum)
row.names(grass_abund_hab_sum) <- 1:nrow(grass_abund_hab_sum)
grass_abund_hab_sum$model <- "habitat"
grass_abund_hab_sum$AIC <- AIC(grass_abund_hab)
write.csv(grass_abund_hab_sum, file="Results/Grassland_abund_habitat.csv", row.names=FALSE)



###################################################################################################################################

# Grassland ONLY abundance

grass_abund_arable1 <- glmer(grass_only_abund ~ prop_arable.500 + (1|subcompartment) + (1|year),
                             data=moth_hab_data, family = poisson())
summary(grass_abund_arable1)   
r.squaredGLMM(grass_abund_arable1) # 0.003647241                
grass_abund_arable2 <- glmer(grass_only_abund ~ prop_arable.1500 + (1|subcompartment) + (1|year),
                             data=moth_hab_data, family = poisson())
summary(grass_abund_arable2)            
r.squaredGLMM(grass_abund_arable2) # 0.04825056              
grass_abund_arable3 <- glmer(grass_only_abund ~ prop_arable.3000 + (1|subcompartment) + (1|year),
                             data=moth_hab_data, family = poisson())
summary(grass_abund_arable3)                                      
r.squaredGLMM(grass_abund_arable3) # 0.01035713             
#### ARABLE 1500M #### 

grass_abund_broad1 <- glmer(grass_only_abund ~ prop_broadleaf.500 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(grass_abund_broad1)     
r.squaredGLMM(grass_abund_broad1) # 0.03053912               
grass_abund_broad2 <- glmer(grass_only_abund ~ prop_broadleaf.1500 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(grass_abund_broad2)         
r.squaredGLMM(grass_abund_broad2) # 0.01571018               
grass_abund_broad3 <- glmer(grass_only_abund ~ prop_broadleaf.3000 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(grass_abund_broad3)                 
r.squaredGLMM(grass_abund_broad3) # 0.003379272              
#### BROADLEAF 500M #### 

grass_abund_conifer1 <- glmer(grass_only_abund ~ prop_conifer.500 + (1|subcompartment) + (1|year),
                              data=moth_hab_data, family = poisson())
summary(grass_abund_conifer1)                                                
r.squaredGLMM(grass_abund_conifer1) # 0.003518660                
grass_abund_conifer2 <- glmer(grass_only_abund ~ prop_conifer.1500 + (1|subcompartment) + (1|year),
                              data=moth_hab_data, family = poisson())
summary(grass_abund_conifer2)                
r.squaredGLMM(grass_abund_conifer2) # 0.0003669832                
grass_abund_conifer3 <- glmer(grass_only_abund ~ prop_conifer.3000 + (1|subcompartment) + (1|year),
                              data=moth_hab_data, family = poisson())
summary(grass_abund_conifer3)                                                   
r.squaredGLMM(grass_abund_conifer3) # 0.01836061               
#### CONIFER 3000M #### 

grass_abund_grass1 <- glmer(grass_only_abund ~ prop_grassland.500 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(grass_abund_grass1)                          
r.squaredGLMM(grass_abund_grass1) # 0.001016016                 
grass_abund_grass2 <- glmer(grass_only_abund ~ prop_grassland.1500 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(grass_abund_grass2)                                                      
r.squaredGLMM(grass_abund_grass2) # 0.01201995                 
grass_abund_grass3 <- glmer(grass_only_abund ~ prop_grassland.3000 + (1|subcompartment) + (1|year),
                            data=moth_hab_data, family = poisson())
summary(grass_abund_grass3)        
r.squaredGLMM(grass_abund_grass3) # 0.009056124                
#### GRASSLAND 1500M #### 

# save results
r2 <- c(r.squaredGLMM(grass_abund_arable1)[1,1], r.squaredGLMM(grass_abund_arable2)[1,1], r.squaredGLMM(grass_abund_arable3)[1,1],
        r.squaredGLMM(grass_abund_broad1)[1,1], r.squaredGLMM(grass_abund_broad2)[1,1], r.squaredGLMM(grass_abund_broad3)[1,1],
        r.squaredGLMM(grass_abund_conifer1)[1,1], r.squaredGLMM(grass_abund_conifer2)[1,1], r.squaredGLMM(grass_abund_conifer3)[1,1],
        r.squaredGLMM(grass_abund_grass1)[1,1], r.squaredGLMM(grass_abund_grass2)[1,1], r.squaredGLMM(grass_abund_grass3)[1,1])
land_cover <- c("Arable.500", "Arable.1500", "Arable.3000", "Broadleaf.500", "Broadleaf.1500", "Broadleaf.3000",
                "Conifer.500", "Conifer.1500", "Conifer.3000", "Grassland.500", "Grassland.1500", "Grassland.3000")
grass_abund_land <- data.frame(Response=rep("Grassland only moth abundance"), Land_cover=land_cover, R2=r2)
write.csv(grass_abund_land, file="Results/Moths/Grassland_only_abund_land_cover.csv", row.names=FALSE)

## Now run models with habitat variables only 
grass_abund_hab <- glmer(grass_only_abund ~ scale(dbh_average) + scale(basal_area) +
                           scale(per_broadleaf_canopy) + scale(complexity_score) + 
                           scale(canopy_openess) + scale(minimum_temperature) + 
                           scale(cloud_cover_.) + scale(dis_to_edge) + 
                           scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ (1|year) + 
                           (1|subcompartment), data = moth_hab_data, family = poisson(),
                         glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                         na.action = "na.fail")
check_collinearity(grass_abund_hab) # low correlations after removing moon cycle
testDispersion(grass_abund_hab) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = grass_abund_hab, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod2 <- glmer(grass_only_abund ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + 
                      scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_arable.1500) + 
                      (1|year) + (1|subcompartment), data = moth_hab_data, family = poisson(),
                      glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
check_collinearity(final_mod2) # low correlations after removing moon cycle
testDispersion(final_mod2) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod2, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod3 <- glmer(grass_only_abund ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + 
                      scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_broadleaf.500) + 
                      (1|year) + (1|subcompartment), data = moth_hab_data, family = poisson())
check_collinearity(final_mod3) # low correlations after removing moon cycle
testDispersion(final_mod3) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod3, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod4 <- glmer(grass_only_abund ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + 
                      scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_conifer.3000) + 
                      (1|year) + (1|subcompartment), data = moth_hab_data, family = poisson())
check_collinearity(final_mod4) # low correlations after removing moon cycle
testDispersion(final_mod4) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod4, plot = F)
plot(simulationOutput) ## no assumptions violated

final_mod5 <- glmer(grass_only_abund ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + scale(minimum_temperature) + 
                      scale(cloud_cover_.) + scale(dis_to_edge) + 
                      scale(wind_speed) + scale(dis_broadleaf) + scale(altitude)+ scale(prop_grassland.1500) + 
                      (1|year) + (1|subcompartment), data = moth_hab_data, family = poisson())
check_collinearity(final_mod5) # low correlations after removing moon cycle
testDispersion(final_mod5) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod5, plot = F)
plot(simulationOutput) ## no assumptions violated

AIC(grass_abund_hab, final_mod2, final_mod3, final_mod4, final_mod5)
# arable model 
summary(final_mod2) # no significant structural variables

## save all AIC values for supplementary table
AIC_values <- data.frame(model=c("habitat", "arable 1500m", "broadleaf 500m", "conifer 3000m", "grassland 1500m"),
                         AIC=AIC(grass_abund_hab, final_mod2, final_mod3, final_mod4, final_mod5)$AIC, stringsAsFactors=FALSE) 
write.csv(AIC_values, file="Results/Grassland_only_abund_AIC.csv", row.names=FALSE)

## save lowest AIC model summary
final_mod2_sum <- as.data.frame(coef(summary(final_mod2))) #selecting full model coefficient averages
final_mod2_sum$parameters <- row.names(final_mod2_sum)
row.names(final_mod2_sum) <- 1:nrow(final_mod2_sum)
final_mod2_sum$model <- "arable"
final_mod2_sum$AIC <- AIC(final_mod2)
write.csv(final_mod2_sum, file="Results/Grassland_only_abund_habitat.csv", row.names=FALSE)
