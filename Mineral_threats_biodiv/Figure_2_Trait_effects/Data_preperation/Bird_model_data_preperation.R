# Bird model data preperation and rds creation


library(stringr); library(dplyr);library(readr);library(ggplot2);library(tidyr)
library(brms)
library(lqmm)

rm(list=ls())
gc()

getwd()
setwd("X:/edwards_lab1/User/bop21ipl/Chapter_One2")
start_time <- Sys.time()

# Data ---- 
M_sp <- read.csv("Data/Bird_trait_model/Chordata_Mine_threatened_assessments.csv")

M_sp <- M_sp %>%
  pull(scientificName)
# full list of names including Name matches 
Names_final <- read_csv("Data/Bird_trait_model/Bird_final_name_matches.csv")


# non imputed traits
NONimp_trait <- read_csv("Data/Bird_trait_model/Bird_traits_for_reduced_NO_imp_needed.csv")
NONimp_trait <- NONimp_trait %>% 
  # remove traits that are correlated
  select(-c(Marine.intertidal.or.coastal.supratidal, Artificial , Introduced.vegetation))

glimpse(NONimp_trait)

# imputed  traits 
imp_trait <- read_csv("Data/Bird_trait_model/Bird_imputed_traits2.csv")
# extract only remaining traits 
imp_trait <- imp_trait %>% 
  select(species, Generation_length_log_d, Carnivore, Herbivore)

# bird range data 
bird_range <- read_csv("Data/Bird_trait_model/Bird_ranges_ALL.csv") %>% 
  select(binomial, range_logstd)

# Load Phylo matrix
bird_phylodist <- read.table("Data/Bird_trait_model/Bird_phylo_Distmatrix.txt")
# remove columns and rows of names not in trait database
bird_phylodist <- as.matrix(bird_phylodist)
bird_phylodist[1:10,1:10]
isSymmetric.matrix(bird_phylodist)
identical(colnames(bird_phylodist),rownames(bird_phylodist))

# Load distance matrix 
bird_dist <- read.table("Data/Bird_trait_model/Bird_centroid_Correlationmatrix_FULL.txt")
bird_dist <- as.matrix(bird_dist)
bird_dist[1:10,1:10]
isSymmetric.matrix(bird_dist)
identical(colnames(bird_dist),rownames(bird_dist))

identical(colnames(bird_dist),colnames(bird_phylodist))
length(colnames(bird_dist))
length(colnames(bird_phylodist))

# match the phylo and spatial matricies
matrix_match <- colnames(bird_phylodist) %in% colnames(bird_dist)
bird_phylodist2 <- bird_phylodist[which(matrix_match==TRUE),which(matrix_match == TRUE)]
length(colnames(bird_phylodist2))


# data set of mined threatened species + range data
bird_mined <- Names_final %>%
  mutate(mine_thrnd = if_else(From_IUCN %in% M_sp, 1, 0),
         species = str_replace(From_phylo, " ", "_")) %>%
  left_join(bird_range, by = c("From_IUCN" = "binomial")) %>% 
  select(species, mine_thrnd, range_logstd) %>%
  distinct() %>% 
  na.omit()

# add imputed data
data <- left_join(NONimp_trait, bird_mined, by = "species") %>%
  left_join(imp_trait, by = "species") %>%
  filter(species %in% colnames(bird_phylodist2))%>%
  distinct()%>%
  arrange(species) %>% 
  mutate(species_space  = species,
         # rounding imputed binomial traits to 1 or 0
         across(.cols = Carnivore:Herbivore, round ),
         across(!starts_with("species"), .fns = function(x) round(x,digits = 6))) %>% 
  na.omit()


glimpse(data)

# check name matching 
identical((data$species), colnames(bird_phylodist2))
setdiff((data$species), colnames(bird_phylodist2))
setdiff(colnames(bird_phylodist2), (data$species))

# match the phylo and spatial matricies
sp_match <- colnames(bird_phylodist2) %in% data$species
bird_phylodist2 <- bird_phylodist2[which(sp_match==TRUE), which(sp_match == TRUE)]
length(colnames(bird_phylodist2))

bird_dist <- bird_dist[which(sp_match==TRUE), which(sp_match == TRUE)]
length(colnames(bird_dist))

identical((data$species), colnames(bird_phylodist2))
identical((data$species), colnames(bird_dist))

length(data$species)
nrow(bird_phylodist2)
nrow(bird_dist) #9022 

#clear some space
rm(bird_phylodist, bird_mined, bird_range, bird_trait); gc()

# data save ====
write_rds(bird_phylodist2, "Data/Bird_trait_model/Clean_data/Bird_phylo_Distmatrix.rds")
write_rds(bird_dist, "Data/Bird_trait_model/Clean_data/Bird_centroid_Correlationmatrix_FULL.rds")
write_csv(data, "Data/Bird_trait_model/Clean_data/Bird_final_data.csv")

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

write_rds(sp_not_used_list, "Data/Bird_trait_model/Bird_species_not_used_IUCNnames.rds")


