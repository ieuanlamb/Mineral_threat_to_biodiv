# Fish Distance matrix 

# library(tidyverse)

library(tidyr)
library(dplyr)
library(readr)
library(stringr)
library(sf)
library(lqmm)
library(Rphylopars)
library(ape)
library(tidytree)

getwd()

# read in centroids
fish_centroids <- st_read("data/FISH/Centroids/Fish_centroids_ALL.shp")
fish_centroids <- fish_centroids %>%
  st_transform(crs = "+proj=cea +datum=WGS84")

# read names 
Fish_names <- read_csv("Ecological trait data/Synonyms/Fish_Name_match.csv")

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

write.table(fish_Dist, "data/FISH/Fish_centroid_Distmatrix.txt")
fish_Dist <- read.table("data/FISH/Fish_centroid_Distmatrix.txt")

# check symmetry by selectinf which cells do not equal the transposed cell 
which(fish_Dist != t(fish_Dist), arr.ind = TRUE) # CHECK is distance matrix Symetrical?

# # convert to a correlation matrix 
fish_dist_cor <- (max(fish_Dist) - fish_Dist)/max(fish_Dist)

fish_dist_cor <- as.matrix(fish_dist_cor)
range(fish_dist_cor)

#check
which(fish_dist_cor != t(fish_dist_cor), arr.ind = TRUE) # CHECK is distance matrix Symetrical?
fish_dist_cor <- fish_dist_cor[-c(4107,1653), -c(1653,4107)] # remove non matching species names 

rownames(fish_dist_cor)[rownames(fish_dist_cor)!= colnames(fish_dist_cor)]
colnames(fish_dist_cor)[rownames(fish_dist_cor)!= colnames(fish_dist_cor)]
# rename "Quietula_y-cauda"
rownames(fish_dist_cor) <- str_replace(rownames(fish_dist_cor), "-",".")

isSymmetric(fish_dist_cor)
identical(rownames(fish_dist_cor), colnames(fish_dist_cor))
is.positive.definite(fish_dist_cor) 
fish_dist_cor <- make.positive.definite(fish_dist_cor)

# save matrix 
write.table(fish_dist_cor, "data/FISH/Fishcentroid_Correlationmatrix_FULL.txt")

fish_dist_cor <- read.table("data/FISH/Fishcentroid_Correlationmatrix_FULL.txt")

# phylo dist matrix ######################

# read in the MCC tree 
fish_tree <- read.tree("Ecological trait data/Phylogenies/Fish/actinopt_12k_raxml.tre.xz")

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
unique(phylo_dist_matrix[col(phylo_dist_matrix)==row(phylo_dist_matrix)])
# number of speciesnames
length(colnames(phylo_dist_matrixn)) # 6621
unique(colnames(phylo_dist_matrixn) == rownames(phylo_dist_matrixn)) # TRUE
isSymmetric.matrix(phylo_dist_matrixn) # TRUE
is.positive.definite(phylo_dist_matrixn)

write.table(phylo_dist_matrixn, "data/FISH/Fish_phylo_Distmatrix.txt")



