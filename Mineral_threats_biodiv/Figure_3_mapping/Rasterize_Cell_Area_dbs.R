#### Rasterise the grids from Cell Area databases

library(sp)
library(raster)
library(tibble); library(stringr); library(dplyr); library(readr);library(ggplot2);library(tidyr)
library(sf)
library(viridis)
library(rgeos)
library(rgbif)
library(ggtext)
library(tmap)
library(rnaturalearth)
library(rnaturalearthdata)
library(stars)
library(lqmm)
library("RColorBrewer")

getwd()
setwd("X:/edwards_lab1/User/bop21ipl")
sf_use_s2(FALSE)


target_crs <- st_crs("+proj=moll +x_0=0 +y_0=0 +lat_0=0 +lon_0=0")

#Make bounding box to create a grid from ----
bbox <- st_bbox(c(xmin= -17596910, ymin= -6731548, xmax= 17596910, ymax= 8748765),
                crs = target_crs)

#Load global spatial Land file
# world <- st_read("Species_Ranges/Data/Land/ne_50m_land_no_artic.shp",layer = "ne_50m_land_no_artic")%>%
#   st_transform(crs ="+proj=cea +datum=WGS84") #Project to equal area
# extent(world)

#Load global spatial Land file
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_make_valid() %>%
  filter(continent != "Antarctica") %>% 
  st_transform(crs = target_crs)


#create a grid of 111000 size 
grid <- st_make_grid(bbox, cellsize = c(111000, 111000)) %>%
  st_sf() %>% 
  mutate(cell = 1:nrow(.))  %>% 
  st_transform(crs = target_crs)  


# Get coordinates or the center of each grid cell
world_grid_centroids <- grid %>% 
  st_centroid() %>% 
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>% 
  as.data.frame() %>% 
  dplyr::select(-geometry)

# Select only terrestrial cells for now
grid_test <- grid %>%
  st_intersection(world) %>%
  rowid_to_column(var = "cell2")

# Get coordinates or the center of each TERRESTRIAL grid cell
terst_grid_centroids <- grid_test %>%
  st_centroid() %>%
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>%
  as.data.frame() %>%
  dplyr::select(-geometry)

# check mapping 
tm_shape(world)+ 
  tm_polygons()+
  tm_shape(grid_test)+
  tm_borders()

####### load cell areas -----
# cell areas from the HPC run of above code

# Amphibians
amph_cell_Areas <- read_csv("IUCN_data/Species_Ranges/Data/amphs_minetht_cell_areas2.csv")
length(unique(amph_cell_Areas$species)) # 747
nrow(amph_cell_Areas)
amph_cell_Areas <- amph_cell_Areas %>% 
  distinct()
nrow(amph_cell_Areas)

#Birds
bird_cell_Areas <- read_csv("IUCN_data/Species_Ranges/Data/Bird_minetht_cell_areas2.csv")
length(unique(bird_cell_Areas$species)) # 556
nrow(bird_cell_Areas)
bird_cell_Areas[(duplicated(bird_cell_Areas)),]
bird_cell_Areas <- bird_cell_Areas %>% 
  distinct()
nrow(bird_cell_Areas)

# Mammals
mam_cell_Areas <- read_csv("IUCN_data/Species_Ranges/Data/Mams_minetht_cell_areas2.csv")
length(unique(mam_cell_Areas$species)) # 516
nrow(mam_cell_Areas)
mam_cell_Areas[(duplicated(mam_cell_Areas)),]
mam_cell_Areas <- mam_cell_Areas %>% 
  distinct()
nrow(mam_cell_Areas)


# Reptiles
rep_cell_Areas <- read_csv("IUCN_data/Species_Ranges/Data/reps_minetht_cell_areas2.csv")
length(unique(rep_cell_Areas$species)) # 737
nrow(rep_cell_Areas)
rep_cell_Areas[(duplicated(rep_cell_Areas)),]
rep_cell_Areas <- rep_cell_Areas %>% 
  distinct()
nrow(rep_cell_Areas)


# Fish
fish_cell_Areas <- read_csv("IUCN_data/Species_Ranges/Data/fish_minetht_cell_areas2.csv")
length(unique(fish_cell_Areas$species)) # 2015
nrow(fish_cell_Areas) 
fish_cell_Areas[(duplicated(fish_cell_Areas)),]
fish_cell_Areas <- fish_cell_Areas %>% 
  distinct()
nrow(fish_cell_Areas)

#### Rasterize ====

#### Amphibians ####
# sum any multiple distinct ranges within a single cell
amph_sp_cell_Values <- amph_cell_Areas %>%
  group_by(species,cell)%>%
  summarise(Areakm2 = sum(Areakm2))%>%
  # Add total species range column 
  group_by(species) %>%
  mutate(total_Areakm2 = sum(Areakm2))%>%
  # calculate proportion of sp. range within each grid cell 
  ungroup%>%
  mutate(Rng_prop_cell = Areakm2/total_Areakm2)


# Create dataframe of sum of Rng_proportions_cell
amph_cell_threat_val <- amph_sp_cell_Values %>%
  group_by(cell) %>%
  summarise(threat_val = sum(Rng_prop_cell)) %>%
  full_join(world_grid_centroids, by  = "cell")%>%
  mutate(threat_val = if_else(is.na(threat_val), 0, threat_val))%>%
  arrange(-x)%>%
  dplyr::select(cell,x,y, threat_val) %>% 
  arrange(cell)


# create grid shapefile with threat values 
amph_grid_ras <- left_join(grid, amph_cell_threat_val, by = "cell") %>% 
  dplyr::select(threat_val)

library("RColorBrewer")
# visualise raster
tm_shape(world)+
  tm_polygons()+
  tm_shape(amph_grid_ras)+
  tm_fill(col = "threat_val", alpha = 0.7, palette = "YlGnBu")

st_write(amph_grid_ras, "IUCN_data/Species_Ranges/Outputs/Amphibians/Amph_tht_cert_raster3.gpkg", append = FALSE)

#rasterize grid
amph_thrt_raster_cert <- amph_grid_ras %>% 
  st_rasterize()


tm_shape(world)+
  tm_polygons()+
  tm_shape(amph_thrt_raster_cert)+
  tm_raster(col = "threat_val", alpha = 0.7, palette = "YlGnBu")

tm_shape(world)+
  tm_fill(col = "#93ccdb")+
  tm_shape(amph_grid_ras)+
  tm_fill(col = "threat_val", palette = viridis(30, alpha = 0.6, option = "E"),
            style = "kmeans"
          )+
  tm_layout (frame = FALSE, bg.color = "transparent", main.title = "Amphibian threat value" , legend.title.color = "transparent")

# write the raster of threats to all amphibian species with certainty values of the proportion each species ranges within raster pixels
write_stars(amph_thrt_raster_cert, dsn = "IUCN_data/Species_Ranges/Outputs/Amphibians/Amph_tht_cert_raster3.tif", overwrite = TRUE)


#### Birds ####
# sum any multiple distinct ranges within a single cell
bird_sp_cell_Values <- bird_cell_Areas %>%
  group_by(species,cell)%>%
  summarise(Areakm2 = sum(Areakm2))%>%
  # Add total species range column 
  group_by(species) %>%
  mutate(total_Areakm2 = sum(Areakm2))%>%
  # calculate proportion of sp. range within each grid cell 
  ungroup%>%
  mutate(Rng_prop_cell = Areakm2/total_Areakm2)


# Create dataframe of sum of Rng_proportions_cell
bird_cell_threat_val <- bird_sp_cell_Values %>%
  group_by(cell) %>%
  summarise(threat_val = sum(Rng_prop_cell)) %>%
  full_join(world_grid_centroids, by  = "cell")%>%
  mutate(threat_val = if_else(is.na(threat_val), 0, threat_val))%>%
  arrange(-x)%>%
  dplyr::select(cell,x,y, threat_val) %>% 
  arrange(cell)


# create grid shapefile with threat values 
bird_grid_ras <- left_join(grid, bird_cell_threat_val, by = "cell") %>% 
  dplyr::select(threat_val)

library("RColorBrewer")
# visualise raster
tm_shape(world)+
  tm_polygons()+
  tm_shape(bird_grid_ras)+
  tm_fill(col = "threat_val", alpha = 0.7, palette = "YlGnBu")

st_write(bird_grid_ras, "IUCN_data/Species_Ranges/Outputs/Birds/bird_tht_cert_raster3.gpkg")

#rasterize grid
bird_thrt_raster_cert <- bird_grid_ras %>% 
  st_rasterize()

# write the raster of threats to all amphibian species with certainty values of the proportion each species ranges within raster pixels
write_stars(bird_thrt_raster_cert, dsn = "IUCN_data/Species_Ranges/Outputs/Birds/bird_tht_cert_raster3.tif", overwrite = TRUE)

#### Mammals ####
# sum any multiple distinct ranges within a single cell
mam_sp_cell_Values <- mam_cell_Areas %>%
  group_by(species,cell)%>%
  summarise(Areakm2 = sum(Areakm2))%>%
  # Add total species range column 
  group_by(species) %>%
  mutate(total_Areakm2 = sum(Areakm2))%>%
  # calculate proportion of sp. range within each grid cell 
  ungroup%>%
  mutate(Rng_prop_cell = Areakm2/total_Areakm2)


# Create dataframe of sum of Rng_proportions_cell
mam_cell_threat_val <- mam_sp_cell_Values %>%
  group_by(cell) %>%
  summarise(threat_val = sum(Rng_prop_cell)) %>%
  full_join(world_grid_centroids, by  = "cell")%>%
  mutate(threat_val = if_else(is.na(threat_val), 0, threat_val))%>%
  arrange(-x)%>%
  dplyr::select(cell,x,y, threat_val) %>% 
  arrange(cell)


# create grid shapefile with threat values 
mam_grid_ras <- left_join(grid, mam_cell_threat_val, by = "cell") %>% 
  dplyr::select(threat_val)

library("RColorBrewer")
# visualise raster
tm_shape(world)+
  tm_polygons()+
  tm_shape(mam_grid_ras)+
  tm_fill(col = "threat_val", alpha = 0.7, palette = "YlGnBu")

st_write(mam_grid_ras, "IUCN_data/Species_Ranges/Outputs/Mammals/mam_tht_cert_raster3.gpkg", append = F)

#rasterize grid
mam_thrt_raster_cert <- mam_grid_ras %>% 
  st_rasterize()

tm_shape(world)+
  tm_fill(col = "#93ccdb")+
  tm_shape(mam_grid_ras)+
  tm_fill(col = "threat_val", palette = viridis(30, alpha = 0.6, option = "E"),
          style = "kmeans")+
  tm_layout (frame = FALSE, bg.color = "transparent", main.title = "Mammals threat value" , legend.title.color = "transparent")


# write the raster of threats to all amphibian species with certainty values of the proportion each species ranges within raster pixels
write_stars(mam_thrt_raster_cert, dsn = "IUCN_data/Species_Ranges/Outputs/Mammals/mam_tht_cert_raster3.tif", overwrite = TRUE)


#### Reptiles ####
# sum any multiple distinct ranges within a single cell
rep_sp_cell_Values <- rep_cell_Areas %>%
  group_by(species,cell)%>%
  summarise(Areakm2 = sum(Areakm2))%>%
  # Add total species range column 
  group_by(species) %>%
  mutate(total_Areakm2 = sum(Areakm2))%>%
  # calculate proportion of sp. range within each grid cell 
  ungroup%>%
  mutate(Rng_prop_cell = Areakm2/total_Areakm2)


# Create dataframe of sum of Rng_proportions_cell
rep_cell_threat_val <- rep_sp_cell_Values %>%
  group_by(cell) %>%
  summarise(threat_val = sum(Rng_prop_cell)) %>%
  full_join(world_grid_centroids, by  = "cell")%>%
  mutate(threat_val = if_else(is.na(threat_val), 0, threat_val))%>%
  arrange(-x)%>%
  dplyr::select(cell,x,y, threat_val) %>% 
  arrange(cell)


# create grid shapefile with threat values 
rep_grid_ras <- left_join(grid, rep_cell_threat_val, by = "cell") %>% 
  dplyr::select(threat_val)

library("RColorBrewer")
# visualise raster
tm_shape(world)+
  tm_polygons()+
  tm_shape(rep_grid_ras)+
  tm_fill(col = "threat_val", alpha = 0.7, palette = "YlGnBu")

st_write(rep_grid_ras, "IUCN_data/Species_Ranges/Outputs/Reptiles/rep_tht_cert_raster3.gpkg", append = F)

#rasterize grid
rep_thrt_raster_cert <- rep_grid_ras %>% 
  st_rasterize()

tm_shape(world)+
  tm_fill(col = "#93ccdb")+
  tm_shape(rep_grid_ras)+
  tm_fill(col = "threat_val", palette = viridis(30, alpha = 0.6, option = "E"),
          style = "kmeans")+
  tm_layout (frame = FALSE, bg.color = "transparent", main.title = "Reptiles threat value" , legend.title.color = "transparent")


# write the raster of threats to all amphibian species with certainty values of the proportion each species ranges within raster pixels
write_stars(rep_thrt_raster_cert, dsn = "IUCN_data/Species_Ranges/Outputs/Reptiles/rep_tht_cert_raster3.tif", overwrite = TRUE)

#### Fish ####
# sum any multiple distinct ranges within a single cell
fish_sp_cell_Values <- fish_cell_Areas %>%
  group_by(species,cell)%>%
  summarise(Areakm2 = sum(Areakm2))%>%
  # Add total species range column 
  group_by(species) %>%
  mutate(total_Areakm2 = sum(Areakm2))%>%
  # calculate proportion of sp. range within each grid cell 
  ungroup%>%
  mutate(Rng_prop_cell = Areakm2/total_Areakm2)


# Create dataframe of sum of Rng_proportions_cell
fish_cell_threat_val <- fish_sp_cell_Values %>%
  group_by(cell) %>%
  summarise(threat_val = sum(Rng_prop_cell)) %>%
  full_join(world_grid_centroids, by  = "cell")%>%
  mutate(threat_val = if_else(is.na(threat_val), 0, threat_val))%>%
  arrange(-x)%>%
  dplyr::select(cell,x,y, threat_val) %>% 
  arrange(cell)


# create grid shapefile with threat values 
fish_grid_ras <- left_join(grid, fish_cell_threat_val, by = "cell") %>% 
  dplyr::select(threat_val)

library("RColorBrewer")
# visualise raster
tm_shape(world)+
  tm_polygons()+
  tm_shape(fish_grid_ras)+
  tm_fill(col = "threat_val", alpha = 0.7, palette = "YlGnBu")

st_write(fish_grid_ras, "IUCN_data/Species_Ranges/Outputs/Fish/fish_tht_cert_raster3.gpkg", append = F)

#rasterize grid
fish_thrt_raster_cert <- fish_grid_ras %>% 
  st_rasterize()

tm_shape(world)+
  tm_fill(col = "#93ccdb")+
  tm_shape(fish_grid_ras)+
  tm_fill(col = "threat_val", palette = viridis(30, alpha = 0.6, option = "E"),
          style = "kmeans")+
  tm_layout (frame = FALSE, bg.color = "transparent", main.title = "Reptiles threat value" , legend.title.color = "transparent")


# write the raster of threats to all amphibian species with certainty values of the proportion each species ranges within raster pixels
write_stars(fish_thrt_raster_cert, dsn = "IUCN_data/Species_Ranges/Outputs/Fish/fish_tht_cert_raster3.tif", overwrite = TRUE)
