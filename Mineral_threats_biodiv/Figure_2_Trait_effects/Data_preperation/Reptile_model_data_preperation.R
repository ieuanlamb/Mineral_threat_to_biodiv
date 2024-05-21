# Reptile trait model script m2 bernoulli model

# Reptile trait model 1 to run on the HPC

library(stringr); library(dplyr);library(readr);library(ggplot2);library(tidyr)
library(lqmm)


# Data ---- 
M_sp <- read.csv("Data/Rep_trait_model/Chordata_Mine_threatened_assessments.csv")

M_sp <- M_sp %>%
  pull(scientificName)
# full list of names including Name matches 
Names_final <- read_csv("Data/Rep_trait_model/Reptile_Final_namematch2.csv")

# imputed  traits 
imp_trait <- read_csv("Data/Rep_trait_model/rep_loocv_imputed_sp.csv")

# Read in Phylogenetic distance matrix
rep_phylodist <- read.table("Data/Rep_trait_model/Rep_phylo_Distmatrix.txt")
# remove columns and rows of names not in trait database
rep_phylodist <- as.matrix(rep_phylodist)
rep_phylodist[1:10,1:10]
isSymmetric.matrix(rep_phylodist)
identical(colnames(rep_phylodist),rownames(rep_phylodist))
match <- colnames(rep_phylodist) != rownames(rep_phylodist)
colnames(rep_phylodist)[which(match == TRUE)] # "Phelsuma_v.nigra"
rownames(rep_phylodist)[which(match == TRUE)] # "Phelsuma_v-nigra"
rep_phylodist <- rep_phylodist[-which(rownames(rep_phylodist) == "Phelsuma_v-nigra"), -which(colnames(rep_phylodist) == "Phelsuma_v.nigra")]

colnames(rep_phylodist)[which(match == "Phelsuma_v-nigra")] 
setdiff(rownames(rep_phylodist),colnames(rep_phylodist))
identical(colnames(rep_phylodist),rownames(rep_phylodist))
ncol(rep_phylodist)

# read in distance correlation matrix 
rep_dist <- read.table("Data/Rep_trait_model/Reptile_centroid_Distmatrix.txt")
which(rep_dist != t(rep_dist), arr.ind = TRUE) # distance matrix is Symetrical 
identical(rownames(rep_dist), colnames(rep_dist))
match <- colnames(rep_dist) != rownames(rep_dist)

# # convert to a correlation matrix 
rep_dist <- (max(rep_dist) - rep_dist)/max(rep_dist) # scale of 0 : maxdist - 1
rep_dist <- as.matrix(rep_dist)
range(rep_dist)
nrow(rep_dist)
colnames(rep_dist)[which(match == TRUE)] 
rownames(rep_dist)[which(match == TRUE)] 
rep_dist <- rep_dist[-which(rownames(rep_dist) == "Phelsuma_v-nigra"), -which(colnames(rep_dist) == "Phelsuma_v.nigra")]
identical(rownames(rep_dist), colnames(rep_dist))

# check symmetry by selectinf which cells do not equal the transposed cell 
which(rep_dist != t(rep_dist), arr.ind = TRUE)
rep_dist <- as.matrix(rep_dist)

isSymmetric(rep_dist)
is.positive.definite(rep_dist)
rep_Dist <- make.positive.definite(rep_dist)
rep_Dist[1:10,1:10]
identical(colnames(rep_Dist), colnames(rep_dist))
identical(rownames(rep_Dist), rownames(rep_dist))
identical(rownames(rep_Dist), colnames(rep_Dist))

# match the phylo and spatial matricies
matrix_match <- colnames(rep_phylodist) %in% colnames(rep_Dist)
rep_phylodist2 <- rep_phylodist[which(matrix_match==TRUE),which(matrix_match == TRUE)]
identical(rownames(rep_phylodist2), colnames(rep_Dist))
identical(colnames(rep_phylodist2), rownames(rep_Dist))
setdiff(colnames(rep_phylodist2),colnames(rep_dist)) # should be empty

# data set of mined threatened species 
rep_mined <- Names_final %>%
  mutate(mine_thrnd = if_else(From_IUCN %in% M_sp, 1, 0),
         species = str_replace(From_Phylo, " ", "_")) %>%
  select(species, mine_thrnd) %>%
  na.omit()

data <- left_join(rep_mined, imp_trait, by = "species") %>%
  filter(species %in% colnames(rep_phylodist2))%>%
  distinct()%>%
  arrange(species) %>% 
  mutate(species_space  = species)
glimpse(data)
# check name matching 
identical((data$species), colnames(rep_phylodist2))
setdiff((data$species), colnames(rep_phylodist2))
setdiff(colnames(rep_phylodist2), (data$species))

# remove missing data from matricies 
rep_phylodist2 <- rep_phylodist2[which(rownames(rep_phylodist2) %in% data$species), which(colnames(rep_phylodist2) %in% data$species)]
rep_Dist <- rep_Dist[which(rownames(rep_Dist) %in% data$species), which(colnames(rep_Dist) %in% data$species)]


length(data$species)
nrow(rep_phylodist2)
nrow(rep_Dist)
identical((data$species), colnames(rep_Dist))
setdiff((data$species), colnames(rep_Dist))
identical(rownames(rep_Dist), colnames(rep_phylodist2))
identical(rownames(rep_phylodist2),colnames(rep_Dist))

# save data
write_rds(rep_phylodist2, "Data/Rep_trait_model/final_phylodist2.rds")
write_rds(rep_Dist, "Data/Rep_trait_model/final_Dist2.rds")
write_csv(data, "Data/Rep_trait_model/Reptile_final_data.csv")



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
         species = str_replace(From_Phylo, " ", "_")) %>% 
  filter(species %in% data$species) 


# species not used 
sp_not_used <- Names_final %>% 
  mutate(mine_thrnd = if_else(From_IUCN %in% M_sp, 1, 0),
         species = str_replace(From_Phylo, " ", "_")) %>% 
  filter(!species %in% data$species)


sp_not_used_list <- sp_not_used %>%  pull(From_IUCN) %>% unique() # 1338

write_rds(sp_not_used_list, "Data/Rep_trait_model/Rep_species_not_used_IUCNnames.rds")
