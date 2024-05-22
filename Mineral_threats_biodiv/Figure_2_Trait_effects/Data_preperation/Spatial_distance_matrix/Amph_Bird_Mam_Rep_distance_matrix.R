# Spatial Matricies for use in ecological trait models 
# aim: To create a matricie of species proximity as ecological traits and mining threat are both linked to the
# spatial location of species

library(sf)
library(tidyverse)
library(tmap)
library(s2)
library(lqmm)


### Amphibian spatial distance matrix #########################################################################################

# full list of names including Name matches 
Names_final <- read_csv("Synonyms/Amphibian_final_name_matches.csv")

# load shapefiles  
world <- st_read("C:/Users/bop21ipl/My Drive/PhD/IUCN data/Species_Ranges/Data/Land/ne_50m_land_no_artic.shp")
world <- world %>%
  st_transform(crs = "+proj=cea +datum=WGS84")

amphs <- st_read("C:/Users/bop21ipl/My Drive/PhD/IUCN data/Species_Ranges/Data/AMPHS/AMPHIBIANS_ALL_09.12.21/AMPHIBIANS.shp")
amphs <- amphs %>%
  st_transform(crs = "+proj=cea +datum=WGS84")%>%
  filter(seasonal != 4)



# calculate a distance matrix from centroids of shape file
amphs_sp <- amphs %>%
  st_simplify(preserveTopology = TRUE, dTolerance = 1000) # simplify by a 1000m tolerance, could justify by the standard 10km radius of mining effect used in lit

# map plot checks 
tm_shape(world)+
  tm_polygons(border.col = NULL)+
  tm_shape(amphs_sp[["binomial" == "Agalychnis medinae"]]) +
    tm_polygons( col = "#8888ff", border.col = NULL)

spR1 <- amphs[1:20,]

tm_shape(world)+
  tm_polygons(border.col = NULL)+
  tm_shape(spR1) +
  tm_polygons("binomial")+
  tm_shape(amphs_sp[3202,]) +
  tm_polygons(col = "#ea7777", alpha = 0.5)

# make ranges valid through the 
amphs_sp <- amphs_sp %>%
  st_buffer(0)

# join the geometries by species to find the whole range
amphs_sp <- amphs_sp %>%
  group_by(binomial)%>%
  summarise(geometry = st_union(geometry))
# get centroids of species ranges
amph_cents <- amphs_sp %>%
  st_centroid()

st_write(amph_cents, "../IUCN data/Species_Ranges/Outputs/Amphibians/Amphibian_range_centroids.shp", layer = "Range_centroids")
amph_cents <- st_read("../IUCN data/Species_Ranges/Outputs/Amphibians/Amphibian_range_centroids.shp") %>% 
  st_transform(crs = "+proj=cea +datum=WGS84")

# change names to match the phylo tree
amph_centsn <- left_join(amph_cents, select(Names_final, binomial = From_IUCN, From_phylo), by = "binomial")%>%
  mutate(From_phylo = str_replace(From_phylo, " ", "_")) %>%
  select(species = From_phylo) %>%
  na.omit()%>%
  distinct() %>%
  arrange(species)

#checks
multi_names <- as.data.frame(table(amph_centsn$species)) %>% # checking for name duplicates 
  filter(Freq > 1)
multi_names <- multi_names$Var1
# removes species names that are associated with multiple plylo names 
amph_centsn <- amph_centsn %>%
  filter(!species %in% multi_names)

#list of names using phylo names that have range centroids
names_centroids <- (amph_centsn$species)
length(names_centroids) # 6645
 
write.list(names_centroids, "Distance_matrix/species_Wcentroids.csv", row.names = FALSE)

# plot a subset of centroids 
tm_shape(world) +
  tm_polygons(border.col = NULL)+
  tm_shape(amph_cents[1:100,]) +
  tm_bubbles(col = "#ea7777", alpha = 0.5, size = 0.1)

# create a distance matrix
amph_Dist <- st_distance(amph_centsn, amph_centsn) %>% 
  round(digits = 2)
# check symmety 
isSymmetric.matrix(amph_Dist)
#change row and column names to that of the associated species
rownames(amph_Dist) <- amph_centsn$species
colnames(amph_Dist) <- amph_centsn$species
range(amph_Dist)
# convert measure to correlation matrix
amph_dist_cor <- (max(amph_Dist) - amph_Dist)/max(amph_Dist) # scale of 0 : maxdist - 1
amph_dist_cor[1:10, 1:10]
range(amph_dist_cor)
isSymmetric.matrix(amph_dist_cor)

# save distance matrix 
write.table(amph_Dist, "Distance_matrix/Amphibian_centroid_Distmatrix.txt")

amph_Dist <- read.table("Distance_matrix/Amphibian_centroid_Distmatrix.txt")

write.table(amph_dist_cor, "Distance_matrix/Amphibian_centroid_Correlationmatrix.txt")

# distribution of mean correlation values 
plot <- amph_dist_cor %>% 
  tibble()
plot$Rowmean <- rowMeans(plot[,])

plot <- plot %>% select(Rowmean)

plot %>% 
  ggplot(aes(x = Rowmean)) +
  geom_density()


#
#### Mammal distance Matrix ##################################################################################################################
### Mammal phylo distance Matrix ##

# full list of names including Name matches 
mam_traits <- read_csv("Data_collated/Mam_imputed_traits_combine_wRange.csv")
mam_names <- mam_traits %>% 
  select(From_IUCN, phylacine_binomial)

# load shapefiles  
world <- st_read("C:/Users/bop21ipl/My Drive/PhD/IUCN data/Species_Ranges/Data/Land/ne_50m_land_no_artic.shp")
world <- world %>%
  st_transform(crs = "+proj=cea +datum=WGS84")

mams <- st_read("C:/Users/bop21ipl/My Drive/PhD/IUCN data/Species_Ranges/Data/MAMS/MAMMALS_ALL_09.12.21/MAMMALS.shp")
mams <- mams %>%
  st_transform(crs = "+proj=cea +datum=WGS84")%>%
  filter(seasonal != 4)


# map plot checks 
tm_shape(world)+
  tm_polygons(border.col = NULL)+
  tm_shape(mams[3,]) +
  tm_polygons( col = "#8888ff", border.col = NULL)

spR1 <- mams[1:20,]

tm_shape(spR1) +
  tm_polygons("binomial")+
  tm_shape( ) +
  tm_polygons(col = "#ea7777", alpha = 0.5)


# calculate a distance matrix from centroids of shape file
mams_sp <- mams %>%
  st_simplify(preserveTopology = TRUE, dTolerance = 1000) # simplify by a 1000m tolerance, could justify by the standard 10km radius of mining effect used in lit

# make ranges valid through the 
mams_sp <- mams_sp %>%
  st_buffer(0)

# join the geometries by species to find the whole range
mams_sp <- mams_sp %>%
  group_by(binomial)%>%
  summarise(geometry = st_union(geometry))
# get centroids of species ranges
mams_cents <- mams_sp %>%
  st_centroid()

# save mammal centroids layer
st_write(mams_cents, "../IUCN data/Species_Ranges/Outputs/Mammals/mammal_range_centroids.shp", layer = "Mammal_Range_centroids")

# change names to match the phylo tree
mams_centsn <- left_join(mams_cents,mam_names, by = c("binomial"="From_IUCN")) %>%  
  mutate(species = str_replace(phylacine_binomial, " ", "_")) %>%
  na.omit()%>%
  distinct() %>%
  select(species, geometry) %>% 
  arrange(species)

#checks
multi_names <- as.data.frame(table(mams_centsn$species)) %>% # checking for name duplicates 
  filter(Freq > 1)
multi_names <- multi_names$Var1
# removes species names that are associated with multiple plylo names 
mams_centsn <- mams_centsn %>%
  filter(!species %in% multi_names)


#list of names using phylo names that have range centroids
names_centroids <- mams_centsn %>% 
  tibble %>% 
  select(species)
length(names_centroids) # 5304

write_csv(names_centroids, "Distance_matrix/mam_species_Wcentroids.csv")

# plot a subset of centroids 
tm_shape(world) +
  tm_polygons(border.col = NULL)+
  tm_shape(mams_cents[1:100,]) +
  tm_bubbles(col = "#ea7777", alpha = 0.5, size = 0.1)

# create a distance matrix
mam_Dist <- st_distance(mams_centsn, mams_centsn)
#change row and column names to that of the associated species
rownames(mam_Dist) <- mams_centsn$species
colnames(mam_Dist) <- mams_centsn$species

#check matrix
mam_Dist[1:10,1:10]

# # save distance matrix 
write.table(mam_Dist, "Distance_matrix/Mammal_centroid_Distmatrix.txt")
# mam_Dist <- read.table("Distance_matrix/Mammal_centroid_Distmatrix.txt")
# 
# # convert to a correlation matrix 
mam_dist_cor <- (max(mam_Dist) - mam_Dist)/max(mam_Dist) # scale of 0 : maxdist - 1
mam_dist_cor[1:10, 1:10]
range(mam_dist_cor)
mam_dist_cor <- make.positive.definite(mam_dist_cor)

write.table(mam_dist_cor, "Distance_matrix/Mammal_centroid_Correlationmatrix.txt")
mam_dist_cor <- read.table("Distance_matrix/Mammal_centroid_Correlationmatrix.txt")

is.positive.definite(mam_dist_cor)

# distribution of mean correlation values 
plot <- mam_dist_cor %>% 
  tibble()
plot$Rowmean <- rowMeans(plot[,])

plot <- plot %>% select(Rowmean)

plot %>% 
  ggplot(aes(x = Rowmean)) +
  geom_density()


#### Reptile distance Matrix ##################################################################################################################
### Reptile phylo distance Matrix ##

# Rep names 
rep_names <- read_csv("Synonyms/Reptile_Final_namematch2.csv")

# Reptile species with traits - or imputed traits
rep_traits <- read_csv("Data_collated/rep_loocv_imputed_sp.csv") 
# full list of names including Name matches 

rep_sp <- rep_traits %>% pull(species)

# load shapefiles  
world <- st_read("C:/Users/bop21ipl/My Drive/PhD/IUCN data/Species_Ranges/Data/Land/ne_50m_land_no_artic.shp")
world <- world %>%
  st_transform(crs = "+proj=cea +datum=WGS84")

reps <- st_read("C:/Users/bop21ipl/My Drive/PhD/IUCN data/Species_Ranges/Data/REPS/REPTILES_ALL_09.12.21/REPTILES.shp")
reps <- reps %>%
  st_transform(crs = "+proj=cea +datum=WGS84")%>%
  filter(seasonal != 4)


# map plot checks 
tm_shape(world)+
  tm_polygons(border.col = NULL)+
  tm_shape(reps[3:10,]) +
  tm_polygons( col = "#8888ff", border.col = NULL)

spR1 <- reps[1:20,]

tm_shape(spR1) +
  tm_polygons("binomial")+
  tm_shape(world) +
  tm_polygons(col = "#ea7777", alpha = 0.5)


# turn off sf s2 
sf_use_s2(FALSE)

# calculate a distance matrix from centroids of shape file
reps_sp <- reps %>%
  st_simplify(preserveTopology = TRUE, dTolerance = 1000) # simplify by a 1000m tolerance, could justify by the standard 10km radius of mining effect used in lit

# make ranges valid through the 
reps_sp <- reps_sp %>%
  st_buffer(0)

# join the geometries by species to find the whole range
reps_sp <- reps_sp %>%
  group_by(binomial)%>%
  summarise(geometry = st_union(geometry))
# get centroids of species ranges
reps_cents <- reps_sp %>%
  st_centroid()

# save reptile centroids layer
st_write(reps_cents, "../IUCN data/Species_Ranges/Outputs/Reptiles/reptile_range_centroids.shp", layer = "Reptile_Range_centroids")
reps_cents <- st_read("../IUCN data/Species_Ranges/Outputs/Reptiles/reptile_range_centroids.shp", layer = "Reptile_Range_centroids")

# change names to match the phylo tree
reps_centsn <- left_join(reps_cents, rep_names, by = c("binomial"="From_IUCN")) %>%  
  select(species = From_Phylo, geometry) %>% 
  na.omit()%>%
  distinct() %>%
  arrange(species) # 7906

#checks
multi_names <- as.data.frame(table(reps_centsn$species)) %>% # checking for name duplicates 
  filter(Freq > 1)
multi_names <- multi_names$Var1
# removes species names that are associated with multiple plylo names 
reps_centsn <- reps_centsn %>%
  filter(!species %in% multi_names)


#list of names using phylo names that have range centroids
names_centroids <- reps_centsn %>% 
  tibble %>% 
  select(species)
length(names_centroids) # 7906

write_csv(names_centroids, "Distance_matrix/rep_species_Wcentroids.csv")

# plot a subset of centroids 
tm_shape(world) +
  tm_polygons(border.col = NULL)+
  tm_shape(reps_cents[1:100,]) +
  tm_bubbles(col = "#ea7777", alpha = 0.5, size = 0.1)

# create a distance matrix
rep_Dist <- st_distance(reps_centsn, reps_centsn)
#change row and column names to that of the associated species
rownames(rep_Dist) <- reps_centsn$species
colnames(rep_Dist) <- reps_centsn$species

#check matrix
rep_Dist[1:10,1:10]

# # save distance matrix 
write.table(rep_Dist, "Distance_matrix/Reptile_centroid_Distmatrix.txt")
rep_Dist <- read.table("Distance_matrix/Reptile_centroid_Distmatrix.txt")
which(rep_Dist != t(rep_Dist), arr.ind = TRUE) # distance matrix is Symetrical 

# # convert to a correlation matrix 
rep_dist_cor <- (max(rep_Dist) - rep_Dist)/max(rep_Dist) # scale of 0 : maxdist - 1

rep_dist_cor <- as.matrix(rep_dist_cor)
range(rep_dist_cor)

# check symmetry by selectinf which cells do not equal the transposed cell 
which(rep_dist_cor != t(rep_dist_cor), arr.ind = TRUE)
# which(rep_dist_cor_raw != t(rep_dist_cor_raw), arr.ind = TRUE)
# rep_dist_cor[c( 4485   , 1), c( 4485 ,   1)]
# rep_dist_cor_raw[c(4485  ,  1), c(4485  ,  1)]

# find rows where matching is FALSE 
match_F <- colnames(rep_dist_cor) != rownames(rep_dist_cor)
colnames(rep_dist_cor)[which(match_F == TRUE)] # Phelsuma_v.nigra
row.names(rep_dist_cor)[which(match_F == TRUE)] # Phelsuma_v-nigra
# Fix unmatching species names where possible
colnames(rep_dist_cor)[which(match_F == TRUE)] <- "Phelsuma_v_nigra"
rownames(rep_dist_cor)[which(match_F == TRUE)] <- "Phelsuma_v_nigra"

isSymmetric(rep_dist_cor)
is.positive.definite(rep_dist_cor) 
rep_dist_cor <- make.positive.definite(rep_dist_cor)

write.table(rep_dist_cor, "Distance_matrix/Reptile_centroid_Correlationmatrix.txt")
rep_dist_cor_raw <- read.table("Distance_matrix/Reptile_centroid_Correlationmatrix.txt") # Careful when reloading matrix in: can round rows and cause asymetry  
rep_dist_cor_raw <- round(rep_dist_cor_raw, digits = 6)
is.positive.definite(rep_dist_cor_raw) 

# test_Sym <- rep_dist_cor[rng,rng] %>% 
#   as.matrix()
# isSymmetric(test_Sym)
# test_Sym[1:10,1:10]
is.matrix(rep_dist_cor)
isSymmetric.matrix(rep_dist_cor)
isSymmetric(rep_dist_cor)



# distribution of mean correlation values 
plot <- rep_dist_cor %>% 
  tibble()
plot$Rowmean <- rowMeans(plot[,])

plot <- plot %>% select(Rowmean)

plot %>% 
  ggplot(aes(x = Rowmean)) +
  geom_density()


#### Bird distance Matrix ##################################################################################################################
### Bird distance Matrix ##

# bird names 
bird_names <- read_csv("Synonyms/Bird_final_name_matches.csv")

# Bird species with traits - or imputed traits
bird_traits <- read_csv("Data_collated/Bird_imputed_tratis.csv") 
# full list of names including Name matches 

bird_sp <- bird_traits %>% pull(species)

# load shapefiles  
world <- st_read("C:/Users/bop21ipl/My Drive/PhD/IUCN data/Species_Ranges/Data/Land/ne_50m_land_no_artic.shp")
world <- world %>%
  st_transform(crs = "+proj=cea +datum=WGS84")

birds <- st_read("C:/Users/bop21ipl/My Drive/PhD/IUCN data/Species_Ranges/Data/BIRDS/BOTW_Dec21/BOTW.gdb")
birds <- birds %>%
  st_transform(crs = "+proj=cea +datum=WGS84")%>%
  filter(seasonal != 4) %>% 
  if("geom" %in% colnames(birds)) {
    birds <- birds %>% 
      rename(geometry = geom) 
    print("renamed geom")
  }

  
# map plot checks 
tm_shape(world)+
  tm_polygons(border.col = NULL)+
  tm_shape(birds[3:10,]) +
  tm_polygons( col = "#8888ff", border.col = NULL)

spR1 <- birds[1:20,]

tm_shape(spR1) +
  tm_polygons("binomial")+
  tm_shape(world) +
  tm_polygons(col = "#ea7777", alpha = 0.5)


# turn off sf s2 
sf_use_s2(FALSE)

# calculate a distance matrix from centroids of shape file
birds_sp <- birds %>%
  st_simplify(preserveTopology = TRUE, dTolerance = 1000) # simplify by a 1000m tolerance, could justify by the standard 10km radius of mining effect used in lit

# make ranges valid through the 
birds_sp <- birds_sp %>%
  st_buffer(0)

# join the geometries by species to find the whole range
birds_sp <- birds_sp %>%
  group_by(binomial)%>%
  summarise(geometry = st_union(geometry))
# get centroids of species ranges
birds_cents <- birds_sp %>%
  st_centroid()

# save Bird centroids layer
st_write(birds_cents, "../IUCN data/Species_Ranges/Outputs/Birds/Bird_range_centroids.shp", layer = "Bird_Range_centroids")
birds_cents <- st_read("../IUCN data/Species_Ranges/Outputs/Birds/Bird_range_centroids.shp", layer = "Bird_Range_centroids")

# change names to match the phylo tree
birds_centsn <- left_join(birds_cents, bird_names, by = c("binomial"="From_IUCN")) %>%  
  select(species = From_Phylo, geometry) %>% 
  na.omit()%>%
  distinct() %>%
  arrange(species) # 7906

#checks
multi_names <- as.data.frame(table(birds_centsn$species)) %>% # checking for name duplicates 
  filter(Freq > 1)
multi_names <- multi_names$Var1
# removes species names that are associated with multiple plylo names 
birds_centsn <- birds_centsn %>%
  filter(!species %in% multi_names)


#list of names using phylo names that have range centroids
names_centroids <- birds_centsn %>% 
  tibble %>% 
  select(species)
length(names_centroids) # 7906

write_csv(names_centroids, "Distance_matrix/bird_species_Wcentroids.csv")

# plot a subset of centroids 
tm_shape(world) +
  tm_polygons(border.col = NULL)+
  tm_shape(birds_cents[1:100,]) +
  tm_bubbles(col = "#ea7777", alpha = 0.5, size = 0.1)

# create a distance matrix
bird_Dist <- st_distance(birds_centsn, birds_centsn)
#change row and column names to that of the associated species
rownames(bird_Dist) <- birds_centsn$species
colnames(bird_Dist) <- birds_centsn$species

#check matrix
bird_Dist[1:10,1:10]

# # save distance matrix 
write.table(bird_Dist, "Distance_matrix/Bird_centroid_Distmatrix.txt")
bird_Dist <- read.table("Distance_matrix/Bird_centroid_Distmatrix.txt")
which(bird_Dist != t(bird_Dist), arr.ind = TRUE) # distance matrix is Symetrical 

# # convert to a correlation matrix 
bird_dist_cor <- (max(bird_Dist) - bird_Dist)/max(bird_Dist) # scale of 0 : maxdist - 1

bird_dist_cor <- as.matrix(bird_dist_cor)
range(bird_dist_cor)

# check symmetry by selectinf which cells do not equal the transposed cell 
which(bird_dist_cor != t(bird_dist_cor), arr.ind = TRUE)
# which(bird_dist_cor_raw != t(bird_dist_cor_raw), arr.ind = TRUE)
# bird_dist_cor[c( 4485   , 1), c( 4485 ,   1)]
# bird_dist_cor_raw[c(4485  ,  1), c(4485  ,  1)]

# find rows where matching is FALSE 
match_F <- colnames(bird_dist_cor) != rownames(bird_dist_cor)
colnames(bird_dist_cor)[which(match_F == TRUE)] # Phelsuma_v.nigra
row.names(bird_dist_cor)[which(match_F == TRUE)] # Phelsuma_v-nigra
# Fix unmatching species names where possible
colnames(bird_dist_cor)[which(match_F == TRUE)] <- "Phelsuma_v_nigra"
rownames(bird_dist_cor)[which(match_F == TRUE)] <- "Phelsuma_v_nigra"

isSymmetric(bird_dist_cor)
is.positive.definite(bird_dist_cor) 
bird_dist_cor <- make.positive.definite(bird_dist_cor)

write.table(bird_dist_cor, "Distance_matrix/Bird_centroid_Correlationmatrix.txt")
bird_dist_cor_raw <- read.table("Distance_matrix/Bird_centroid_Correlationmatrix.txt") # Careful when reloading matrix in: can round rows and cause asymetry  
bird_dist_cor_raw <- round(bird_dist_cor_raw, digits = 6)
is.positive.definite(bird_dist_cor_raw) 

# test_Sym <- bird_dist_cor[rng,rng] %>% 
#   as.matrix()
# isSymmetric(test_Sym)
# test_Sym[1:10,1:10]
is.matrix(bird_dist_cor)
isSymmetric.matrix(bird_dist_cor)
isSymmetric(bird_dist_cor)



# distribution of mean correlation values 
plot <- bird_dist_cor %>% 
  tibble()
plot$Rowmean <- rowMeans(plot[,])

plot <- plot %>% select(Rowmean)

plot %>% 
  ggplot(aes(x = Rowmean)) +
  geom_density()

