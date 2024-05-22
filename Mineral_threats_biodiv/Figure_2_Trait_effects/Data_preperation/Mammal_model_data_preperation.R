# Mammal data preperation RDS creation 

library(tidyverse)
# library(stringr); library(dplyr);library(readr);library(ggplot2);library(tidyr)
library(lqmm)
library(brms)

getwd()
# load mining threatened species
M_sp <- read.csv("Data/Mam_trait_model/Chordata_Mine_threatened_assessments.csv")

M_sp <- M_sp %>%
  pull(scientificName)
# full list of names including Name matches 
Names_final <- read_csv("Data/Mam_trait_model/Mammal_final_name_matches.csv")

# imputed amphibian traits 
imp_trait <- read_csv("Data/Mam_trait_model/Mam_imputed_traits_combine_wRange.csv")
glimpse(imp_trait)

unique(imp_trait$terrestrial)
sapply(imp_trait, function(y) sum(length(which(is.na(y))))) %>% 
  data.frame()

imp_trait <- imp_trait %>% 
  mutate(species = str_replace(phylacine_binomial, " ", "_"),
         activity_cycle = as.factor(activity_cycle),
         trophic_level = as.factor(trophic_level)) %>% 
  select(species, From_IUCN,
         adult_mass_g, 
         habitat_breadth_n, 
         dphy_plant, 
         dphy_vertebrate, 
         dphy_invertebrate,
         trophic_level,
         foraging_stratum,
         activity_cycle,
         litter_size_n,
         range_calculated,
         marine,
         freshwater,
         terrestrial_non_volant,
         terrestrial_volant)


# standardise all continuous variables 
idx <- base::sapply(imp_trait, class) == "numeric"
imp_trait[, idx] <- sapply(imp_trait[, idx], function(y) scale(y, center = TRUE, scale = TRUE))
glimpse(imp_trait)

imp_trait %>% 
  sapply(function(y) sum(length(which(is.na(y))))) %>% 
  data.frame()



# Read in Phylogenetic distance matrix
mam_phylodist <- read.table("Data/Mam_trait_model/Mammal_phylo_Distmatrix.txt")
# remove columns and rows of names not in trait database
mam_phylodist <- as.matrix(mam_phylodist)
mam_phylodist[1:10,1:10]

isSymmetric.matrix(mam_phylodist)
identical(colnames(mam_phylodist),rownames(mam_phylodist))
setdiff(colnames(mam_phylodist),rownames(mam_phylodist))


# Read in Distance correlation matrix 
mam_dist <- read.table("Data/Mam_trait_model/Mammal_centroid_Correlationmatrix.txt")
mam_dist[1:10, 1:10]
range(mam_dist)
mam_dist <- as.matrix(mam_dist)
identical(colnames(mam_dist),rownames(mam_dist))
isSymmetric(mam_dist)
isSymmetric.matrix(mam_dist)
is.positive.definite(mam_dist)


setdiff(colnames(mam_phylodist),colnames(mam_dist)) # should be empty

# data set of mined threatened species 
mam_mined <- Names_final %>%
  mutate(mine_thrnd = if_else(From_IUCN %in% M_sp, 1, 0)) %>%
  select(From_IUCN, mine_thrnd) %>%
  na.omit()

data <- left_join(mam_mined, imp_trait, by = "From_IUCN") %>%
  filter(species %in% colnames(mam_phylodist)) %>%
  distinct() %>% 
  select(-From_IUCN) %>%
  arrange(species) %>% 
  mutate(species_space  = species) 

# check correlation values of data 
idx <- sapply(data, FUN = class) == "numeric"
cor(data[,idx],data[,idx])

# match names of data to matrix 
match.names <- colnames(mam_phylodist) %in% (data$species)
mam_phylodist <- mam_phylodist[which(match.names == TRUE), which(match.names == TRUE)]
mam_dist <- mam_dist[which(match.names == TRUE), which(match.names == TRUE)]


write_rds(mam_phylodist, "Data/Mam_trait_model/Mammal_phylo_Distmatrix.rds")
write_rds(mam_dist, "Data/Mam_trait_model/Mammal_centroid_Correlationmatrix.rds")
write_csv(data, "Data/Mam_trait_model/Mammal_final_data.csv")


# species removed from dataset =====
# Remove Extinct or Extinct in the wild species 
Extinct_sp <- read_csv("../IUCN_data/Species_Pages/Outputs/Extinct_and_EW_sp2.csv")
extinct_list <- Extinct_sp$binomial

# make a list of assessed species at time of analysis
sp_Assessed <- Names_final %>% 
  filter(!From_IUCN %in% extinct_list ) %>% 
  pull(From_IUCN) %>% unique() 


# dataset of species in model 
sp_in_model <- Names_final %>% 
  mutate(mine_thrnd = if_else(From_IUCN %in% M_sp, 1, 0),
         species = str_replace(From_phylo, " ", "_")) %>% 
  filter(species %in% data$species) 


# species not used 
sp_not_used <- Names_final %>% 
  mutate(mine_thrnd = if_else(From_IUCN %in% M_sp, 1, 0),
         species = str_replace(From_phylo, " ", "_")) %>% 
  filter(!species %in% data$species)


sp_not_used_list <- sp_not_used %>%  pull(From_IUCN) %>% unique()

write_rds(sp_not_used_list, "Data/Mam_trait_model/Mam_species_not_used_IUCNnames.rds")
