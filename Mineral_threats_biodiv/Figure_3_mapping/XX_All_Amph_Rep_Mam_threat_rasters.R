# Script for standardised Maps Rasters 

# This script aims to calculate the area of each species range within each global 110x110km grid cell for every Amphibian, Mammal and Reptile that does
# not have any mineral extraction threats, and then create weighted scores of potential threat and convert then into a raster. The purpose of this is to 
# create a underlying value of biodiversity that is range size weighted to standardise the threat scores.
# The processes and aims are the same as in 01_Mine_thrt_Areas.R and 02_Rasterize_Cell_Areas_dbs.R. However, here - due to the size of opperation and 
# memory required for all species ranges - I have used for loops to save the opperation after running chunks of species at a time.
# This means some of the code is commented out and other parts have been added to accomodate when the loop has crashed and has been resarted from different
# save points.
# NOTE: FOR CLEARER CODE OF THE SAME PROCESSES BUT FOR LESS SPECIES REFER TO 01_Mine_thrt_Areas.R and 02_Rasterize_Cell_Areas_dbs.R.

library(raster)
library(tibble); library(stringr); library(dplyr); library(readr);library(ggplot2);library(tidyr)
library(sf)
library(viridis)
library(ggtext)
library(tmap)
library(rnaturalearth)
library(rnaturalearthdata)
library(stars)
library(lqmm)
library("RColorBrewer")

# ------ Set Up ------
sf_use_s2(FALSE)

# set mollwiede crs 
target_crs <- st_crs("+proj=moll +x_0=0 +y_0=0 +lat_0=0 +lon_0=0")

#Make bounding box to create a grid from ----
bbox <- st_bbox(c(xmin= -17596910, ymin= -6731548, xmax= 17596910, ymax= 8748765),
                crs = target_crs)

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

# read in list of mine threatened species and taxonomy
M_sp <- read_csv("Species_Pages/Raw_Data/CHORDATA_pg_Oil_Mining_threat/Chordata_Mine_threatened_assessments.csv")

# read list of extinct species 
extinct_sp <- read_csv("Species_Pages/Outputs/Extinct_and_EW_sp2.csv")
extinct_sp <- extinct_sp %>% 
  pull(binomial)

############## Ampnibians full non mining threatened species #######################

# load spatial data of species 
Amph_ranges <- st_read("Species_Ranges/Data/AMPHS/AMPHIBIANS_ALL_09.12.21/AMPHIBIANS.shp")

Amph_ranges <- st_transform(Amph_ranges, crs = target_crs)
Amph_ranges <- Amph_ranges %>% 
  rowid_to_column(var = "rowid")

# remove mining threatened species and extinct from the dataset
Amph_ranges <- Amph_ranges %>% 
  filter(!binomial %in% M_sp,
         seasonal != 4,
         !binomial %in% extinct_sp)

# reduce dataset to only necessary variables
Amph_ranges_T <- Amph_ranges  %>% 
  select(rowid, id_no, binomial)

# check mapping projection 
tm_shape(world)+ 
  tm_polygons()+
  tm_shape(grid)+
  tm_borders() +
  tm_shape(Amph_ranges_T[1:100,])+
  tm_polygons(col = "red")

# time running time for 
start_time <- Sys.time()
# save points
save <- c(seq(1000, length(unique(Amph_ranges_T$binomial)), by = 1000),length(unique(Amph_ranges_T$binomial)))
save_point <- Amph_ranges_T$binomial[save]

# Set up dataframe
cell_Areas <- tibble(binomial = character(),
                     cell2 = integer(), 
                     Areakm2 = double())

# first loop ran to save point 6000
Amph_ranges_T[6000:6001,]

# use for loop for each row
for(i in 6001:nrow(Amph_ranges_T)) {
  # find the cells that the range intersects 
  amph_cell_intersect <- st_intersects(Amph_ranges_T[i,], grid, sparse = TRUE) %>%
    unlist()
  
  # extract only grid cells that the species intersects with
  sp_grid <- grid %>%
    filter(cell %in% amph_cell_intersect)
  
  # find area of the species' range within the grids
  sp_cell_intersection <- st_intersection(st_make_valid(Amph_ranges_T[i,]), sp_grid)
  
  #calculate shape area within each cell
  sp_cell_Area  <- sp_cell_intersection %>%
    mutate(Areakm2 = units::set_units(sf::st_area(geometry),"km^2"),
           Areakm2 = as.numeric(Areakm2))%>%
    tibble()%>%
    dplyr::select( binomial, cell, Areakm2)
  
  cell_Areas <- rbind(cell_Areas, sp_cell_Area)
  
  if((sp_cell_Area$binomial %in% save_point)[1]) {
    write_csv(cell_Areas, paste0("Species_Ranges/Outputs/Backups/Amph_cellAreas_db4_",i,".csv"))
    print(i)
  }
  print(i)
}

write_csv(cell_Areas, paste0("Species_Ranges/Outputs/Backups/Amph_cellAreas_db4_6000_8571.csv"))
cell_Areas1 <- read_csv("Species_Ranges/Outputs/Backups/Amph_cellAreas_db4_6000.csv")

# join two datasets of file running
cell_Areas_full <- bind_rows(cell_Areas,cell_Areas1)
# remove any duplicate rows 
cell_Areas_full <- cell_Areas_full %>% 
  distinct()
#write
write_csv(cell_Areas_full, paste0("Species_Ranges/Outputs/Backups/Amph_cellAreas_db4_Full8571.csv"))
cell_Areas_full <- read_csv( "Species_Ranges/Outputs/Backups/Amph_cellAreas_db4_Full8571.csv")

# check that all species are included within the cell area database
length(unique(cell_Areas_full$binomial)) #7121
length(unique(Amph_ranges_T$binomial)) #7121

# sum any multiple distinct ranges within a single cell
amph_sp_cell_Values <- cell_Areas_full %>%
  group_by(binomial,cell)%>%
  summarise(Areakm2 = sum(Areakm2))%>%
  # Add total species range column 
  group_by(binomial) %>%
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
  select(threat_val)

library("RColorBrewer")
# visualise raster
tm_shape(world)+
  tm_polygons()+
  tm_shape(amph_grid_ras)+
  tm_fill(col = "threat_val", alpha = 0.7, palette = "YlGnBu")

# Save as Shapefile
st_write(amph_grid_ras, "Species_Ranges/Outputs/Amphibians/Amph_NONtht_cert_raster3.gpkg")

#rasterize grid
amph_Nonthrt_raster_cert <- amph_grid_ras %>% 
  st_rasterize()
# write the raster of threats to all amphibian species with certainty values of the proportion each species ranges within raster pixels
write_stars(amph_Nonthrt_raster_cert, dsn = "Species_Ranges/Outputs/Amphibians/Amph_NONtht_cert_raster3.tif", overwrite = TRUE)

end_time <- Sys.time()
end_time - start_time
# Time difference of 48.34403 mins

# Find SES raster by dividing threatened raster score by total score -------
Amph_threat_ras <- raster("Species_Ranges/Outputs/Amphibians/Amph_tht_cert_raster3.tif")
Amph_NONth_ras <- raster("Species_Ranges/Outputs/Amphibians/Amph_NONtht_cert_raster3.tif")

# visialise maps 
tm_shape(world)+
  tm_fill()+
  tm_shape(Amph_threat_ras)+
  tm_raster()

tm_shape(world)+
  tm_fill()+
  tm_shape(Amph_NONth_ras)+
  tm_raster()

# join rasters in stack  
amph_stack <- stack(Amph_threat_ras, Amph_NONth_ras)

# find total value of species certainty 
amph_stack$Amph_Sum_cert <- sum(amph_stack)

# divide threatened species certainty values by total
amph_stack$SES <- (amph_stack$Amph_tht_cert_raster/ amph_stack$Amph_Sum_cert)

tm_shape(world)+
  tm_fill()+
  tm_shape(amph_stack)+
  tm_raster()

tm_shape(world)+
  tm_fill()+
  tm_shape(amph_stack$SES)+
  tm_raster() +
  tm_layout(main.title = "Amphibians")
  

############## Reptiles full non mining threatened species #######################

# load spatial data of species 
Rep_ranges <- st_read("Species_Ranges/Data/REPS/REPTILES_ALL_09.12.21/REPTILES.shp")

Rep_ranges <- st_transform(Rep_ranges, crs =target_crs)
Rep_ranges <- Rep_ranges %>% 
  rowid_to_column(var = "rowid")

# remove mining threatened species from the dataset
Rep_ranges <- Rep_ranges %>% 
  filter(!binomial %in% M_sp,
         seasonal != 4,
         !binomial %in% extinct_sp)

# reduce dataset to only necessary variables
Rep_ranges_T <- Rep_ranges  %>% 
  select(rowid, id_no, binomial)

# check mapping projection 
tm_shape(world)+ 
  tm_polygons()+
  tm_shape(grid)+
  tm_borders() +
  tm_shape(Rep_ranges_T[1:100,])+
  tm_polygons(col = "red")

# save points
save <- c(seq(1000, length(unique(Rep_ranges_T$binomial)), by = 1000),length(unique(Rep_ranges_T$binomial)))
save_point <- Rep_ranges_T$binomial[save]

# Set up dataframe
cell_Areas <- tibble(binomial = character(),
                     cell = integer(), 
                     Areakm2 = double())

# use for loop for each row
# time running 
start_time <- Sys.time()
for(i in 1:nrow(Rep_ranges_T)) {
  # find the cells that the range intersects 
  rep_cell_intersect <- st_intersects(Rep_ranges_T[i,], grid, sparse = TRUE) %>%
    unlist()
  
  # extract only grid cells that the species intersects with
  sp_grid <- grid %>%
    filter(cell %in% rep_cell_intersect)
  
  # find area of the species' range within the grids
  sp_cell_intersection <- st_intersection(st_make_valid(Rep_ranges_T[i,]), sp_grid)
  
  #calculate shape area within each cell
  sp_cell_Area  <- sp_cell_intersection %>%
    mutate(Areakm2 = units::set_units(sf::st_area(geometry),"km^2"),
           Areakm2 = as.numeric(Areakm2))%>%
    tibble()%>%
    dplyr::select( binomial, cell, Areakm2)
  
  cell_Areas <- rbind(cell_Areas, sp_cell_Area)
  
  if(sp_cell_Area$binomial[1] %in% save_point) {
    write_csv(cell_Areas, paste0("Species_Ranges/Outputs/Backups/Rep_cellAreas_db4_",Rep_ranges_T$rowid[i],".csv"))
    print("save")
  }
  print(c(i, cell_Areas$binomial[i]))
}

write_csv(cell_Areas, paste0("Species_Ranges/Outputs/Backups/Rep_cellAreas_db3_full8484.csv"))
cell_Areas <- read_csv("Species_Ranges/Outputs/Backups/Rep_cellAreas_db3_full8484.csv")

length(unique(cell_Areas$binomial)) #8484
length(unique(Rep_ranges_T$binomial))
end_time <- Sys.time()
end_time - start_time

# sum any multiple distinct ranges within a single cell
rep_sp_cell_Values <- cell_Areas %>%
  group_by(binomial,cell)%>%
  summarise(Areakm2 = sum(Areakm2))%>%
  # Add total species range column 
  group_by(binomial) %>%
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
  dplyr::select(x,y, threat_val, cell)

# create grid shapefile with threat values 
rep_grid_ras <- left_join(grid, rep_cell_threat_val, by = "cell") %>% 
  select(threat_val)

# save as shapefile
st_write(rep_grid_ras, "Species_Ranges/Outputs/Reptiles/Rep_NONtht_cert_raster3.gpkg")

#rasterize grid
rep_Nonthrt_raster_cert <- rep_grid_ras %>% 
  st_rasterize()
# write the raster of threats to all amphibian species with certainty values of the proportion each species ranges within raster pixels
write_stars(rep_Nonthrt_raster_cert, dsn = "Species_Ranges/Outputs/Reptiles/Rep_NONtht_cert_raster3.tif", overwrite = TRUE)

# Find SES raster by dividing threatened raster score by total score -------
Rep_threat_ras <- raster("Species_Ranges/Outputs/Reptiles/Reps_tht_cert_raster3.tif")
Rep_NONth_ras <- raster("Species_Ranges/Outputs/Reptiles/Rep_NONtht_cert_raster3.tif")

# visialise maps 
tm_shape(world)+
  tm_fill()+
  tm_shape(Rep_threat_ras)+
  tm_raster()

tm_shape(world)+
  tm_fill()+
  tm_shape(Rep_NONth_ras)+
  tm_raster()

# join rasters in stack  
rep_stack <- stack(Rep_threat_ras, Rep_NONth_ras)

# find total value of species certainty 
rep_stack$Rep_Sum_cert <- sum(rep_stack)

# divide threatened species certainty values by total
rep_stack$SES <- (rep_stack$Reps_tht_cert_raster/ rep_stack$Rep_Sum_cert)
crs(rep_stack) <- crs(world)
rep_stack$SES_ter <-  crop(rep_stack$SES_ter, world)

tm_shape(world)+
  tm_fill()+
  tm_shape(rep_stack)+
  tm_raster()

tm_shape(world)+
  tm_fill()+
  tm_shape(rep_stack$SES_ter)+
  tm_raster() +
  tm_layout(main.title = "Reptiles")



############## Mammals full non mining threatened species #######################

# load spatial data of species 
mam_ranges <- st_read("Species_Ranges/Data/MAMS/MAMMALS_ALL_09.12.21/MAMMALS.shp")
mam_ranges <- st_transform(mam_ranges, crs =target_crs)
mam_ranges <- mam_ranges %>% 
  rowid_to_column(var = "rowid")

# remove mining threatened species from the dataset
mam_ranges <- mam_ranges %>% 
  filter(!binomial %in% M_sp,
         seasonal != 4,
         !binomial %in% extinct_sp)

# reduce dataset to only necessary variables
mam_ranges_T <- mam_ranges  %>% 
  select(rowid, id_no, binomial)

# check mapping projection 
tm_shape(world)+ 
  tm_polygons()+
  tm_shape(grid)+
  tm_borders() +
  tm_shape(mam_ranges_T[1:100,])+
  tm_polygons(col = "red")

# save points
save <- c(seq(1000, length(unique(mam_ranges_T$binomial)), by = 1000),length(unique(mam_ranges_T$binomial)))
save_point <- mam_ranges_T$binomial[save]

# check for valid geometries
check <- st_is_valid(mam_ranges_T)
unique(check)
# fix geometries
mam_ranges_T <- mam_ranges_T %>% 
  st_buffer(0)

# Set up empty dataframe
cell_Areas <- tibble(binomial = character(),
                     cell = integer(), 
                     Areakm2 = double())
# time running time for loop
start_time <- Sys.time()
for(i in 33:nrow(mam_ranges_T)) {
  # find the cells that the range intersects 
  mam_cell_intersect <- st_intersects(mam_ranges_T[i,], grid, sparse = TRUE) %>%
    unlist()
  
  # extract only grid cells that the species intersects with
  sp_grid <- grid %>%
    filter(cell %in% mam_cell_intersect)
  
  # find area of the species' range within the grids
  sp_cell_intersection <- st_intersection(mam_ranges_T[i,], sp_grid)
  
  #calculate shape area within each cell
  sp_cell_Area  <- sp_cell_intersection %>%
    mutate(Areakm2 = units::set_units(sf::st_area(geometry),"km^2"),
           Areakm2 = as.numeric(Areakm2))%>%
    tibble()%>%
    dplyr::select( binomial, cell, Areakm2)
  
  cell_Areas <- rbind(cell_Areas, sp_cell_Area)
  
  if(sp_cell_Area$binomial[1] %in% save_point) {
    write_csv(cell_Areas, paste0("Species_Ranges/Outputs/Backups/Mam_cellAreas_db3_", mam_ranges_T$rowid[i],".csv"))
    print("save")
  }
  print(c(i, cell_Areas$binomial[i]))
}

cell_Areas <- read_csv("Species_Ranges/Outputs/Backups/Mam_cellAreas_db3_Lastbackup.csv")

length(unique(cell_Areas$binomial)) #5845
length(unique(mam_ranges_T$binomial)) #5845

end_time <- Sys.time()
end_time - start_time

# sum any multiple distinct ranges within a single cell
mam_sp_cell_Values <- cell_Areas %>%
  group_by(binomial,cell)%>%
  summarise(Areakm2 = sum(Areakm2))%>%
  # Add total species range column 
  group_by(binomial) %>%
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
  dplyr::select(x,y, cell, threat_val)

# create grid shapefile with threat values 
mam_grid_ras <- left_join(grid, mam_cell_threat_val, by = "cell") %>% 
  select(threat_val)

# write as shape file 
st_write(mam_grid_ras, "Species_Ranges/Outputs/Mammals/Mam_NONtht_cert_raster3.gpkg")

#rasterize grid
mam_Nonthrt_raster_cert <- mam_grid_ras %>% 
  st_rasterize()
# write the raster of threats to all amphibian species with certainty values of the proportion each species ranges within raster pixels
write_stars(mam_Nonthrt_raster_cert, dsn = "Species_Ranges/Outputs/Mammals/Mam_NONtht_cert_raster3.tif", overwrite = TRUE)

end_time <- Sys.time()
end_time - start_time

# Find SES raster by dividing threatened raster score by total score -------
Mam_threat_ras <- raster("Species_Ranges/Outputs/Mammals/Mams_tht_cert_raster.tif")
crs(Mam_threat_ras) <- crs(world)
Mam_NONth_ras <- raster("Species_Ranges/Outputs/Mammals/Mam_NONtht_cert_raster3.tif")
crs(Mam_NONth_ras) <- crs(world)

# visialise maps 
tm_shape(world)+
  tm_fill()+
  tm_shape(Mam_threat_ras)+
  tm_raster()

tm_shape(world)+
  tm_fill()+
  tm_shape(Mam_NONth_ras)+
  tm_raster()

# join rasters in stack  
mam_stack <- stack(Mam_threat_ras, Mam_NONth_ras)
# find total value of species certainty 
mam_stack$Mam_Sum_cert <- sum(mam_stack)
# divide threatened species certainty values by total
mam_stack$SES <- (mam_stack$Mams_tht_cert_raster/ mam_stack$Mam_Sum_cert)

# clip to terrestrial landmass
mam_stack_ter <- crop(mam_stack, world)
mam_stack_ter <- mask(mam_stack, world)

tm_shape(world)+
  tm_fill()+
  tm_shape(mam_stack_ter)+
  tm_raster()

tm_shape(world)+
  tm_fill()+
  tm_shape(mam_stack_ter$SES)+
  tm_raster() +
  tm_layout(main.title = "Mammals")







############## Fish full non mining threatened species #######################
# RUNNING THE FOLLOWING CODE ON THE HPC DUE TO NOT BEING ABLE TO ALLOCATE MEMORY 
{
# # load spatial data of species 
# # DD 
# Fish_ranges_DD <- st_read("Species_Ranges/Data/FISH/Fish_Species_ranges_IUCN/Actinopterygii_16.11.22/Actinopterygii_DD_4581sp/data_0.shp",
#                           layer = 'data_0') %>% 
#   st_transform(crs = "+proj=cea +datum=WGS84")
# Fish_ranges_DD  <- Fish_ranges_DD %>% 
#   st_buffer(0)
# 
# # Redlist 
# Fish_ranges_RL <- st_read("Species_Ranges/Data/FISH/Fish_Species_ranges_IUCN/Actinopterygii_16.11.22/Actinopterygii_NT_VU_EN_CR_3844sp",
#                           layer = 'data_0') %>% 
#   st_transform(crs = "+proj=cea +datum=WGS84")
# Fish_ranges_RL  <- Fish_ranges_RL %>% 
#   st_buffer(0)
# 
# # LC1
# Fish_ranges_LC1 <- st_read("Species_Ranges/Data/FISH/Fish_Species_ranges_IUCN/Actinopterygii_16.11.22/Actinopterygii_LC_Charci_Cyptido_Pleur_Anguil_2247sp",
#                            layer = 'data_0') %>% 
#   st_transform(crs = "+proj=cea +datum=WGS84")
# Fish_ranges_LC1  <- Fish_ranges_LC1 %>% 
#   st_buffer(0)
# 
# # LC2
# Fish_ranges_LC2 <- st_read("Species_Ranges/Data/FISH/Fish_Species_ranges_IUCN/Actinopterygii_16.11.22/Actinopterygii_LC_perci_Cyprini_6214sp",
#                            layer = c('data_0','data_1', 'data_2')) %>% 
#   st_transform(crs = "+proj=cea +datum=WGS84")
# Fish_ranges_LC2  <- Fish_ranges_LC2 %>% 
#   st_buffer(0)
# 
# # LC3
# Fish_ranges_LC3 <- st_read("Species_Ranges/Data/FISH/Fish_Species_ranges_IUCN/Actinopterygii_16.11.22/Actinopterygii_LC_Siluri_Gobii_2350sp",
#                            layer = 'data_0') %>% 
#   st_transform(crs = "+proj=cea +datum=WGS84")
# Fish_ranges_LC3  <- Fish_ranges_LC3 %>% 
#   st_buffer(0)
# 
# # LC4
# Fish_ranges_LC4 <- st_read("Species_Ranges/Data/FISH/Fish_Species_ranges_IUCN/Actinopterygii_16.11.22/Actinopterygii_LC_Orders_exp_8largest_3682sp",
#                            layer =  c('data_0','data_1', 'data_2')) %>% 
#   st_transform(crs = "+proj=cea +datum=WGS84")
# Fish_ranges_LC4  <- Fish_ranges_LC4 %>% 
#   st_buffer(0)
# 
# # join all multiple potentially overlapping ranges
# fish_ranges <- rbind(Fish_ranges_DD, Fish_ranges_RL, Fish_ranges_LC1, Fish_ranges_LC2, Fish_ranges_LC3, Fish_ranges_LC4) #874
# 
# fish_ranges <- fish_ranges %>% 
#   select(binomial = BINOMIAL, geometry = '_ogr_geometry_') %>% 
#   rowid_to_column(var = "rowid")
# 
# glimpse(fish_ranges)
# 
# # remove species threatened by mining 
# # read in list of mine threatened species and taxonomy
# M_sp <- read.csv(here::here("Species_Pages", "Raw_Data/CHORDATA_pg_Oil_Mining_threat/Chordata_Mine_threatened/assessments.csv"))
# 
# # remove mining threatened species from the dataset
# fish_ranges <- fish_ranges %>% 
#   filter(!binomial %in% M_sp,
#          seasonal != 4)
# 
# # reduce dataset to only necessary variables
# fish_ranges <- fish_ranges  %>% 
#   select(rowid, binomial)
# 
# check mapping projection

# tm_shape(world)+
#   tm_polygons()+
#   tm_shape(grid)+
#   tm_borders() +
#   tm_shape(fish_ranges[10,])+
#   tm_polygons(col = "red")
#   
# 
# # save points
# save <- c(seq(1000, length(unique(fish_ranges$binomial)), by = 1000),length(unique(fish_ranges$binomial)))
# save_point <- fish_ranges$binomial[save]
# 
# 
# # Set up dataframe
# cell_Areas <- tibble(binomial = character(),
#                      cell = integer(), 
#                      Areakm2 = double())
# # use for loop for each row
# # Attempt to fix geometry
# library(lqmm)
# check <- st_is_valid(fish_ranges)
# unique(check)
# fish_ranges <- fish_ranges %>% 
#   st_make_valid()
# 
# # time running time for 
# start_time <- Sys.time()
# 
# for(i in 1:nrow(fish_ranges)) {
#   # find the cells that the range intersects 
#   mam_cell_intersect <- st_intersects(fish_ranges[i,], grid, sparse = TRUE) %>%
#     unlist()
#   
#   # extract only grid cells that the species intersects with
#   sp_grid <- grid %>%
#     filter(cell %in% mam_cell_intersect)
#   
#   # find area of the species' range within the grids
#   sp_cell_intersection <- st_intersection(fish_ranges[i,], sp_grid)
#   
#   #calculate shape area within each cell
#   sp_cell_Area  <- sp_cell_intersection %>%
#     mutate(Areakm2 = units::set_units(sf::st_area(geometry),"km^2"),
#            Areakm2 = as.numeric(Areakm2))%>%
#     tibble()%>%
#     dplyr::select( binomial, cell, Areakm2)
#   
#   cell_Areas <- rbind(cell_Areas, sp_cell_Area)
#   
#   if(sp_cell_Area$binomial %in% save_point) {
#     write_csv(cell_Areas, paste0("Species_Ranges/Outputs/Backups/fish_cellAreas_db_", fish_ranges$rowid[i],".csv"))
#     print("save")
#   }
#   print(c(i, cell_Areas$binomial[i]))
# }
# 
# # read backup 
# 
# # sum any multiple distinct ranges within a single cell
# fish_sp_cell_Values <- cell_Areas %>%
#   group_by(binomial,cell)%>%
#   summarise(Areakm2 = sum(Areakm2))%>%
#   # Add total species range column 
#   group_by(binomial) %>%
#   mutate(total_Areakm2 = sum(Areakm2))%>%
#   # calculate proportion of sp. range within each grid cell 
#   ungroup%>%
#   mutate(Rng_prop_cell = Areakm2/total_Areakm2)
# 
# 
# # Create dataframe of sum of Rng_proportions_cell
# fish_cell_threat_val <- fish_sp_cell_Values %>%
#   group_by(cell) %>%
#   summarise(threat_val = sum(Rng_prop_cell)) %>%
#   full_join(world_grid_centroids, by  = "cell")%>%
#   mutate(threat_val = if_else(is.na(threat_val), 0, threat_val))%>%
#   arrange(-x)%>%
#   dplyr::select(x,y, threat_val)
# 
# fish_Nonthrt_raster_cert <- rasterFromXYZ(fish_cell_threat_val, res = 111000, crs = "+proj=cea +datum=WGS84")
# 
# # write the raster of threats to all amphibian species with certainty values of the proportion each species ranges within raster pixels
# writeRaster(fish_Nonthrt_raster_cert, filename = here::here("Species_Ranges", "Outputs/Fish/Fish_NONtht_cert_raster.tif"), overwrite = TRUE)
# 
# end_time <- Sys.time()
# end_time - start_time
# 
# # visualise raster
# tm_shape(world)+
#   tm_polygons()+
#   tm_shape(mam_Nonthrt_raster_cert)+
#   tm_raster(alpha = 0.7)
# 
# 
}
# Find SES raster by dividing threatened raster score by total score -------
fish_threat_ras <- raster("Species_Ranges/Outputs/Fish/fish_tht_cert_raster.tif")
crs(fish_threat_ras) <- crs(world)
fish_NONth_ras <- raster("Species_Ranges/Outputs/Fish/fish_NONtht_cert_raster2.tif")
crs(fish_NONth_ras) <- crs(world)

# visialise maps 
tm_shape(world)+
  tm_fill()+
  tm_shape(fish_threat_ras)+
  tm_raster()

tm_shape(world)+
  tm_fill()+
  tm_shape(fish_NONth_ras)+
  tm_raster()

# join rasters in stack  
fish_stack <- stack(fish_threat_ras, fish_NONth_ras)


# find total value of species certainty 
fish_stack$fish_Sum_cert <- sum(fish_stack)

# divide threatened species certainty values by total
fish_stack$SES <- (fish_stack$fish_tht_cert_raster/ fish_stack$fish_Sum_cert)

# clip to terrestrial landmass
fish_stack_ter <- crop(fish_stack, world)
fish_stack_ter <- mask(fish_stack, world)

tm_shape(world)+
  tm_fill()+
  tm_shape(fish_stack_ter)+
  tm_raster()

tm_shape(world)+
  tm_fill()+
  tm_shape(fish_stack_ter$SES)+
  tm_raster() +
  tm_layout(main.title = "Fish")
