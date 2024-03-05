##########################
#### user: Lisbeth Hordley
#### date: September 2021
#### info: Stourhead moths: calculate richness and abundance indices to use in models

rm(list = ls())
options(scipen=999)

library(dplyr)
library(ggplot2)
library(fossil)
library(vegan)

## read in data
moth_data <- read.csv("Data/stourhead_moths.csv", header=TRUE)
hab_data <- read.csv("Data/Stourhead_hab_structure_draft.csv", header=TRUE)
moth_traits <- read.csv("Data/ecological_traits_stourhead_final.csv", header=TRUE)

## merge in trait info: species need to be classed into 4 guilds: woodland & grassland (habitat based) and conifer & broadleaf (food plant based)
# also strict categories for each guild too
moth_traits <- moth_traits[,c("scientific_name","broadleaf_trees", "broadleaf_only",
                              "coniferous_trees", "conifer_only", "Woodland", 
                              "Woodland_only", "grassland", "grassland_only")]
moth_traits <- moth_traits %>% mutate_if(is.integer,as.numeric)
## convert columns into binary 0/1 
## 1 and 2 = consistently and occasional use (change both these to 1)
## 3 and 4 = little evidence or only evidence from abroad (change these and blanks to 0)
moth_traits[moth_traits == 3] <- NA
moth_traits[moth_traits == 4] <- NA
moth_traits[is.na(moth_traits)] <- 0
moth_traits[moth_traits == 2] <- 1
str(moth_traits)

## merge in with moth data
moth_data <- merge(moth_data, moth_traits, by="scientific_name", all.x=TRUE)
moth_data$abundance <- as.numeric(moth_data$abundance)

########################################
### CALCULATE ABUNDANCE AND RICHNESS ###
########################################

moth_data_sum <- moth_data %>% 
  group_by(plot, year, minimum_temperature, wind_speed, moon_cycle, cloud_cover_., ) %>%
  summarise(richness = n_distinct(scientific_name),
            broadleaf_rich = n_distinct(scientific_name[broadleaf_trees==1], na.rm=TRUE),
            broadleaf_only_rich = n_distinct(scientific_name[broadleaf_only==1], na.rm=TRUE),
            conifer_rich = n_distinct(scientific_name[coniferous_trees==1], na.rm=TRUE),
            conifer_only_rich = n_distinct(scientific_name[conifer_only==1], na.rm=TRUE),
            woodland_rich = n_distinct(scientific_name[Woodland==1], na.rm=TRUE),
            woodland_only_rich = n_distinct(scientific_name[Woodland_only==1], na.rm=TRUE),
            grass_rich = n_distinct(scientific_name[grassland==1], na.rm=TRUE),
            grass_only_rich = n_distinct(scientific_name[grassland_only==1], na.rm=TRUE),
            tot_abund=sum(abundance),
            broadleaf_abund = sum(abundance[broadleaf_trees==1], na.rm=TRUE),
            broadleaf_only_abund = sum(abundance[broadleaf_only==1], na.rm=TRUE),
            conifer_abund = sum(abundance[coniferous_trees==1], na.rm=TRUE),
            conifer_only_abund = sum(abundance[conifer_only==1], na.rm=TRUE),
            woodland_abund = sum(abundance[Woodland==1], na.rm=TRUE),
            woodland_only_abund = sum(abundance[Woodland_only==1], na.rm=TRUE),
            grass_abund = sum(abundance[grassland==1], na.rm=TRUE),
            grass_only_abund = sum(abundance[grassland_only==1], na.rm=TRUE))
## remove NAs from host plant/habitat classifications - this is because many aggregate species
## don't have trait data, but they can be included in the total richness and total abundance estimates

## merge habitat variable data with moth data
moth_hab_data <- merge(moth_data_sum, hab_data, by.x=c("plot","year"), by.y=c("plot", "survey_year")) 
## remove unecessary columns 
moth_hab_data[,c("location","dbh_tree_1","dbh_tree_2","dbh_tree_3","dbh_tree_4","dbh_tree_5",
                 "tree_species_names", "canopy_.cover", "sub_canopy_.cover", "tall_understorey_cover",
                 "low_shrub_cover", "ground_layer_cover", "canopy_layer_score", "sub_canopy_layer_score",
                 "tall_understorey_score", "low_shrub_score", "ground_layer_score", "photo_number", "comments")] <- list(NULL)
moth_hab_data<-moth_hab_data[!(moth_hab_data$treatment=="Mid Transitioning Irregular"),] # remove mid-transitioning treatment
unique(moth_hab_data$treatment) ## 3 treatments: early transitioning irregular, irregular and clear fell 
length(unique(moth_hab_data$plot)) # 94 plots in total

## simple plots of abund and richness between treatments
ggplot(moth_hab_data, aes(treatment, tot_abund))+
  geom_boxplot()+ 
  theme_classic() + 
  ylab("Total abundance")

ggplot(moth_hab_data, aes(treatment, richness))+
  geom_boxplot()+ 
  theme_classic() + 
  ylab("Species richness")

## save file
write.csv(moth_hab_data, file="Data/Stourhead_moths_final.csv", row.names=FALSE)
