# Amphibian data preperation for trait modelling and RDS creation

library(stringr); library(dplyr);library(readr);library(ggplot2);library(tidyr)
library(brms)


# Data ---- 
# load species with minineral extraction threats Downloaded from IUCN Red List species pages 
M_sp <- read_csv("Data/Amph_trait_model/Chordata_Mine_threatened_assessments.csv")

M_sp <- M_sp %>%
  pull(scientificName)
# full list of names including Name matches 
Names_final <- read_csv("Data/Amph_trait_model/Amphibian_final_name_matches.csv")
Names_final$From_IUCN %>% length()

# imputed amphibian traits 
imp_trait <- read_csv("Data/Amph_trait_model/Amph_imputed_checked_final.csv")

# Read in Phylogenetic distance matrix
amph_phylodist <- read.table("Data/Amph_trait_model/Distance_matrix/Amphibian_phylo_Distmatrix.txt")
# remove columns and rows of names not in trait database
amph_phylodist <- amph_phylodist[-which(names(amph_phylodist) %in% c("Pristimantis_w.nigrum", "Scinax_v.signatus", "Scinax_x.signatus")),-which(names(amph_phylodist) %in% c("Pristimantis_w.nigrum", "Scinax_v.signatus", "Scinax_x.signatus"))]
amph_phylodist <- as.matrix(amph_phylodist)
amph_phylodist[1:10,1:10]
isSymmetric.matrix(amph_phylodist)
identical(colnames(amph_phylodist),rownames(amph_phylodist))

amph_dist <- read.table("Data/Amph_trait_model/Distance_matrix/Amphibian_Corr_matrix_SymPosDef.txt")
amph_dist[1:10, 1:10]
range(amph_dist)
amph_dist <- as.matrix(amph_dist)
identical(colnames(amph_dist),rownames(amph_dist))
isSymmetric(amph_dist)
isSymmetric.matrix(amph_dist)

setdiff(colnames(amph_phylodist),colnames(amph_dist)) # should be empty

# data set of mined threatened species 
amph_mined <- Names_final %>% 
  mutate(mine_thrnd = if_else(From_IUCN %in% M_sp, 1, 0),
         species = str_replace(From_phylo, " ", "_")) %>%
  select(species, mine_thrnd) %>%
  na.omit()

data <- left_join(amph_mined, imp_trait, by = "species") %>%
  filter(species %in% colnames(amph_phylodist))%>%
  distinct()%>%
  arrange(species) %>% 
  mutate(species_space  = species)


write_rds(amph_phylodist, "Data/Amph_trait_model/Clean_data/Distance_matrix/Amphibian_centroid_Correlationmatrix.rds")
write_rds(amph_dist, "Data/Amph_trait_model/Clean_data/Distance_matrix/Amphibian_Corr_matrix_SymPosDef.rds")
write_csv(data, "Data/Amph_trait_model/Clean_data/amphibian_final_data.csv")


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

sp_not_used_list <- sp_not_used %>%  pull(From_IUCN) %>%  unique()

write_rds(sp_not_used_list, "Data/Amph_trait_model/Amph_species_not_used_IUCNnames.rds")

