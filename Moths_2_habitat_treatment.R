##########################
#### user: Lisbeth Hordley
#### date: September 2021
#### info: Stourhead moths: differences in habitat between treatments

rm(list = ls())
options(scipen=999)

library(MuMIn)
library(DHARMa)
library(lme4)
library(glmmTMB)
library(dplyr)
library(multcomp)
library(broom)
library(ggplot2)
library(plyr)

moth_hab_data <- read.csv("Data/Stourhead_moths_final.csv", header=TRUE)

## change habitat codes
moth_hab_data$treatment <- recode(moth_hab_data$treatment, "Clear Fell" = "Stage 1", "Early Transitioning Irregular" = "Stage 2", "Irregular" = "Stage 3")

############# GLMS ##############
## Responses: Basal area, average dbh, canopy openness, % broadleaf canopy and complexity score
## Random effects: year + subcomparment
## Covariates: distance to edge, distance to broadleaved, min temp, 
## cloud cover, moon phase, wind speed, and aspect
## Predictors: treatment
## check correlation of fixed effects - probably don't need all of them

moth_hab_data <- moth_hab_data[,c("plot", "year", "treatment", "subcompartment", "basal_area", "dbh_average",
                          "per_broadleaf_canopy", "canopy_openess", "complexity_score")]

moth_hab_data$basal_area <- as.numeric(moth_hab_data$basal_area)
moth_hab_data$canopy_openess <- as.numeric(moth_hab_data$canopy_openess)
summary(moth_hab_data) # 1 NA in basal area, 11 NAs in canopy openness


########## 1. BASAL AREA ##########
## poisson model (log linear model had significant homogeneity of variance)
basal_mod <- glmer(basal_area ~ treatment + (1|year) + (1|subcompartment),
                   family="poisson", moth_hab_data)
summary(basal_mod)

sim1 <- DHARMa::simulateResiduals(fittedModel = basal_mod)
DHARMa::testDispersion(sim1)
plot(sim1) ## no assumptions violated

## Then, analyze post hoc comparisons as:
Treat.comp<-glht(basal_mod, mcp(treatment="Tukey"))
summary(Treat.comp)
# save results
res1 <- data.frame(broom::tidy(Treat.comp)[,c(1:2,4:5,7)])
res1$term <- "basal area"

letters <- data.frame(treatment = names(cld(Treat.comp)$mcletters$Letters),
                      L = cld(Treat.comp)$mcletters$Letters)
hdata1 <- merge(moth_hab_data, letters)
maxy_basal <- max(ddply(hdata1, "treatment", summarise, fivenum(basal_area)[5])[,2])*1.1

basal_area_p <- ggplot(hdata1, aes(x=treatment, y = basal_area)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA) + 
  geom_text(aes(x=treatment, y=maxy_basal, label = L), size=7) + 
  geom_point( position = position_jitter(w = 0.1, h = 0)) +
  ylab("Basal area")+
  xlab("Treatment") +
  theme_classic() +
  theme(text=element_text(size=20))
basal_area_p

########## 2. AVERAGE DBH ##########
## poisson model (log linear model had significant homogeneity of variance)
dbh_mod <- lmer(log(dbh_average) ~ treatment + (1|year) + (1|subcompartment), moth_hab_data)
summary(dbh_mod)

sim1 <- DHARMa::simulateResiduals(fittedModel = dbh_mod)
DHARMa::testDispersion(sim1)
plot(sim1) ## significant test for homogeneity of variance

## Then, analyze post hoc comparisons as:
Treat.comp<-glht(dbh_mod,mcp(treatment='Tukey'))
summary(Treat.comp)
# save results
res2 <- data.frame(broom::tidy(Treat.comp)[,c(1:2,4:5,7)])
res2$term <- "average dbh"
res_final <- rbind(res1, res2)

letters <- data.frame(treatment = names(cld(Treat.comp)$mcletters$Letters),
                      L = cld(Treat.comp)$mcletters$Letters)
hdata1 <- merge(moth_hab_data, letters)
maxy_dbh <- max(ddply(hdata1, "treatment", summarise, fivenum(dbh_average)[5])[,2])*1.1

average_dbh_p <- ggplot(hdata1, aes(x=treatment, y = dbh_average)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA) + 
  geom_text(aes(treatment, y=maxy_dbh, label = L), size=7) + 
  geom_point( position = position_jitter(w = 0.1, h = 0)) +
  ylab("Average DBH")+
  xlab("Treatment") +
  theme_classic() +
  theme(text=element_text(size=20))
average_dbh_p

########## 3. CANOPY OPENNESS ##########
## poisson model (log linear model had significant homogeneity of variance)
canopy_mod <- lmer(log(canopy_openess) ~ treatment + (1|year) + (1|subcompartment), moth_hab_data)
summary(canopy_mod)

sim1 <- DHARMa::simulateResiduals(fittedModel = canopy_mod)
DHARMa::testDispersion(sim1)
plot(sim1) ## all good

## Then, analyze post hoc comparisons as:
Treat.comp<-glht(canopy_mod,mcp(treatment='Tukey'))
summary(Treat.comp)
# save results
res3 <- data.frame(broom::tidy(Treat.comp)[,c(1:2,4:5,7)])
res3$term <- "canopy openness"
res_final <- rbind(res_final, res3)

letters <- data.frame(treatment = names(cld(Treat.comp)$mcletters$Letters),
                      L = cld(Treat.comp)$mcletters$Letters)
hdata1 <- merge(moth_hab_data, letters)
maxy_canopy <- max(ddply(hdata1, "treatment", summarise, fivenum(canopy_openess)[5])[,2])*1.1

canopy_p <- ggplot(hdata1, aes(x=treatment, y = canopy_openess)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA) + 
  geom_text(aes(treatment, y=maxy_canopy, label = L), size=7) + 
  geom_point(position = position_jitter(w = 0.1, h = 0)) +
  ylab("Canopy openness")+
  xlab("Treatment") +
  theme_classic() +
  theme(text=element_text(size=20))
canopy_p

########## 4. % BROADLEAF CANOPY ##########
# zero-inflated model
moth_hab_data$treatment = as.factor(moth_hab_data$treatment)
broadleaf_mod <- glmmTMB(per_broadleaf_canopy ~ treatment + (1|year) + (1|subcompartment),
                         data=moth_hab_data, ziformula=~1, family=poisson)
summary(broadleaf_mod)

sim1 <- DHARMa::simulateResiduals(fittedModel = broadleaf_mod)
DHARMa::testDispersion(sim1)
plot(sim1) ## looks good

## Then, analyze post hoc comparisons as:
Treat.comp<-glht(broadleaf_mod,mcp(treatment='Tukey'))
summary(Treat.comp)
# save results
res4<- data.frame(broom::tidy(Treat.comp)[,c(1:2,4:5,7)])
res4$term <- "percentage broadleaf"
res_final <- rbind(res_final, res4)

letters <- data.frame(treatment = names(cld(Treat.comp)$mcletters$Letters),
                      L = cld(Treat.comp)$mcletters$Letters)
hdata1 <- merge(moth_hab_data, letters)
maxy_broad <- max(ddply(hdata1, "treatment", summarise, fivenum(per_broadleaf_canopy)[5])[,2])*1.1

broadleaf_p <- ggplot(hdata1, aes(x=treatment, y = per_broadleaf_canopy)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA) + 
  geom_text(aes(treatment, y=maxy_broad, label = L), size=7) + 
  geom_point(position = position_jitter(w = 0.1, h = 0)) +
  ylab("Percentage \nbroadleaf canopy")+
  xlab("Treatment") +
  theme_classic()+
  theme(text=element_text(size=20))
broadleaf_p

########## 5. COMPLEXITY SCORE ##########
complex_mod <- lmer(complexity_score ~ treatment + (1|year) + (1|subcompartment), data=moth_hab_data)
summary(complex_mod) 

sim1 <- DHARMa::simulateResiduals(fittedModel = complex_mod)
DHARMa::testDispersion(sim1)
plot(sim1) ## looks good

## Then, analyze post hoc comparisons as:
Treat.comp<-glht(complex_mod,mcp(treatment='Tukey'))
summary(Treat.comp)
# save results
res5 <- data.frame(broom::tidy(Treat.comp)[,c(1:2,4:5,7)])
res5$term <- "complexity score"
res_final <- rbind(res_final, res5)
write.csv(res_final, file="Results/Habitat_treatment_tukey.csv", row.names=FALSE)

letters <- data.frame(treatment = names(cld(Treat.comp)$mcletters$Letters),
                      L = cld(Treat.comp)$mcletters$Letters)
hdata1 <- merge(moth_hab_data, letters)
maxy_complex <- max(ddply(hdata1, "treatment", summarise, fivenum(complexity_score)[5])[,2])*1.1

complex_p <- ggplot(hdata1, aes(x=treatment, y = complexity_score)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA) + 
  geom_text(aes(treatment, y=maxy_complex, label = L), size=7) + 
  geom_point(position = position_jitter(w = 0.1, h = 0)) +
  ylab("Complexity score")+
  xlab("Treatment") +
  theme_classic() +
  theme(text=element_text(size=20))
complex_p

###### Put all plots together
# 3 x 3 plots
library(ggpubr)
library(gridExtra)
library(grid)

layout_matrix <- matrix(c(1,1,2,2,3,3,6,4,4,5,5,6), nrow = 2, byrow = TRUE)
habitat_plots <- gridExtra::grid.arrange(basal_area_p, average_dbh_p, canopy_p, broadleaf_p,
                                         complex_p, layout_matrix = layout_matrix)
habitat_plots

ggsave(habitat_plots, file="Graphs/Habitat_treatment_moths.png", height=8, width=14)


