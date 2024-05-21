# Data preperation and matrix creation for fish
library(stringr); library(dplyr);library(readr);library(ggplot2);library(tidyr)
library(brms)
library(cmdstanr)
library(lqmm)
library(sf)
library(ape)
library(units)
getwd()

## Distance matrix ------
# read in centroids
fish_centroids <- st_read("Data/Fish_trait_model/Matrix/Fish_centroids_ALL.shp")
fish_centroids <- fish_centroids %>%
  st_transform(crs = "+proj=cea +datum=WGS84")

# read names 
Fish_names <- read_csv("Data/Fish_trait_model/Matrix/Fish_Name_match.csv")

# change names to match the phylo tree
fish_centsn <- left_join(fish_centroids, Fish_names, by = c("BINOMIAL"="From_IUCN")) %>%  
  select(species = From_phylo, geometry) %>% 
  mutate(species = str_replace(species, " ", "_")) %>% 
  na.omit()%>%
  distinct() %>%
  arrange(species) # 6615 species in phylotree of actinoperygii

#checks
multi_names <- as.data.frame(table(fish_centsn$species)) %>% # checking for name duplicates 
  filter(Freq > 1)
multi_names <- multi_names$Var1
# removes species names that are associated with multiple plylo names 
fish_centsn <- fish_centsn %>%
  filter(!species %in% multi_names)


# create a distance matrix
fish_Dist <- st_distance(fish_centsn, fish_centsn)

rownames(fish_Dist) <- fish_centsn$species
colnames(fish_Dist) <- fish_centsn$species

#check matrix
fish_Dist[1:10,1:10]


# check symmetry by selectinf which cells do not equal the transposed cell 
which(fish_Dist != t(fish_Dist), arr.ind = TRUE) # CHECK is distance matrix Symetrical?

# # convert to a correlation matrix 
fish_dist_cor <- (max(fish_Dist) - fish_Dist)/max(fish_Dist)

fish_dist_cor <- as.matrix(fish_dist_cor)
range(fish_dist_cor)

#check
which(fish_dist_cor != t(fish_dist_cor), arr.ind = TRUE) # CHECK is distance matrix Symetrical?
# fish_dist_cor <- fish_dist_cor[-c(4107,1653), -c(1653,4107)] # remove non matching species names 

rownames(fish_dist_cor)[rownames(fish_dist_cor)!= colnames(fish_dist_cor)]
colnames(fish_dist_cor)[rownames(fish_dist_cor)!= colnames(fish_dist_cor)]
# rename "Quietula_y-cauda"
# colnames(fish_dist_cor) <- str_replace(rownames(fish_dist_cor), "-",".")

fish_dist_cor <- set_units(fish_dist_cor, NULL )

isSymmetric(fish_dist_cor)
identical(rownames(fish_dist_cor), colnames(fish_dist_cor))

is.positive.definite(fish_dist_cor) 
fish_dist_cor <- make.positive.definite(fish_dist_cor)

# write.table(fish_dist_cor,"Data/FISH/Fishcentroid_corr_matrix2.tbl")

# phylo dist matrix ######################
# read in the MCC tree 
fish_tree <- read.tree("Data/Fish_trait_model/Matrix/actinopt_12k_raxml.tre.xz")

# read in names of species with range centroids 
names_list <- Fish_names %>% 
  filter(From_IUCN %in% fish_centroids$BINOMIAL) %>% 
  select(From_phylo)%>% 
  na.omit() %>% 
  pull(From_phylo) %>% 
  str_replace(" ", "_")

# 
# # tip labels 
labels <- fish_tree$tip.label

# difference between the two list (labels in the phylo that are not in IUCN names)
labelsdiff <- setdiff(labels,names_list)

#remove tips that DO NOT correspond to a IUCN species that have range centroids
fish_tree <- drop.tip(fish_tree, labelsdiff)

# for creating a distance matrix
phylo_dist_matrix <- vcv(fish_tree, corr = T)
glimpse(phylo_dist_matrix)
phylo_dist_matrix[1:10,1:10]

# reorder the columns
order <- sort(colnames(phylo_dist_matrix))
phylo_dist_matrixn <- phylo_dist_matrix[order,order]

#checks
phylo_dist_matrixn[1:10,1:10]
# check the diagonals are all equal
unique(phylo_dist_matrixn[col(phylo_dist_matrixn)==row(phylo_dist_matrixn)])
# number of speciesnames
length(colnames(phylo_dist_matrixn)) # 6621
unique(colnames(phylo_dist_matrixn) == rownames(phylo_dist_matrixn)) # TRUE
isSymmetric.matrix(phylo_dist_matrixn) # TRUE
is.positive.definite(phylo_dist_matrixn)

start_time <- Sys.time()
# Data ---- 
M_sp <- read_csv("Data/Fish_trait_model/Chordata_Mine_threatened_assessments.csv")

M_sp <- M_sp %>%
  pull(scientificName)
# full list of names including Name matches 
Names_final <- read_csv("Data/Fish_trait_model/Fish_Acti_final_names.csv")

# imputed  traits 
imp_trait <- read_csv("Data/Fish_trait_model/fish_imputed_traits.csv")
# range data 
range_data <- read_csv("Data/Fish_trait_model/Fish_range_calc_all.csv")
range_data <- range_data %>% 
  mutate(range_size_logst = scale(log(range_calc), center = TRUE, scale = TRUE)[,1]) %>% 
  left_join(Names_final, by = c("BINOMIAL" = "From_IUCN")) %>% 
  select(species = From_phylo, range_size_logst) %>% 
  mutate(species = str_replace(species, " ", "_")) %>% 
  na.omit()


# join two datasets
trait_data <- left_join(range_data, imp_trait, by = "species") 

# Load distance matrix
fish_phylodist <- phylo_dist_matrixn

# Load distance matrix 
fish_dist <- fish_dist_cor

identical(colnames(fish_dist),colnames(fish_phylodist))
length(colnames(fish_dist))
length(colnames(fish_phylodist))


# match the phylo and spatial matricies
matrix_match <- colnames(fish_phylodist) %in% colnames(fish_dist)
matrix_match2 <- colnames(fish_dist) %in% colnames(fish_phylodist)

fish_phylodist2 <- fish_phylodist[which(matrix_match==TRUE),which(matrix_match == TRUE)]
fish_dist2 <- fish_dist[which(matrix_match2==TRUE),which(matrix_match2 == TRUE)]
length(colnames(fish_phylodist2))
length(colnames(fish_dist2))

identical(colnames(fish_dist2),colnames(fish_phylodist2))
identical(rownames(fish_dist2),rownames(fish_phylodist2))

# data set of mined threatened species 
fish_mined <- Names_final %>%
  mutate(mine_thrnd = if_else(From_IUCN %in% M_sp, 1, 0),
         species = str_replace(From_phylo, " ", "_")) %>%
  select(species, mine_thrnd) %>%
  na.omit()

data <- left_join(fish_mined, trait_data, by = "species") %>%
  filter(species %in% colnames(fish_phylodist2))%>%
  distinct()%>%
  arrange(species) %>% 
  mutate(species_space  = species,
         # rounding imputed binomial traits to 1 or 0
         across(!starts_with("species"), .fns = function(x) round(x,digits = 6))) %>% 
  na.omit()

# check name matching 
identical((data$species), colnames(fish_phylodist2))
setdiff((data$species), colnames(fish_phylodist2))
setdiff(colnames(fish_phylodist2), (data$species))

# match the phylo and spatial matricies removing species with matched 
sp_match <- colnames(fish_phylodist2) %in% data$species
fish_phylodist2 <- fish_phylodist2[which(sp_match==TRUE), which(sp_match == TRUE)]
length(colnames(fish_phylodist2))

fish_dist2 <- fish_dist2[which(sp_match==TRUE), which(sp_match == TRUE)]
length(colnames(fish_dist2))


identical((data$species), colnames(fish_phylodist2))
identical((data$species), colnames(fish_dist2))
length(data$species)
nrow(fish_phylodist2)
nrow(fish_dist2)


# Save ====
write_csv(data, "Data/Fish_trait_model/final_data.csv")
write_rds(fish_dist2, "Data/Fish_trait_model/final_Dist2.rds")
write_rds(fish_phylodist2, "Data/Fish_trait_model/final_phylodist2.rds")



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

write_rds(sp_not_used_list, "Data/Fish_trait_model/Fish_species_not_used_IUCNnames.rds")
