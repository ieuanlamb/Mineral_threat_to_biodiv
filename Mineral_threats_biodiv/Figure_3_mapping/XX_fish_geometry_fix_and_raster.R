# Fixing geometries and combining species cell area grids for fish

# NOTE: Due to the size and number of fish species range sizes this operation was conducted on HPC stanage.

library(readr)
library(dplyr)
library(tibble)
library(sf)
library(lqmm)
library(sf)
library(terra)
library(stars)

# ------- Set up -------- 
target_crs <- st_crs("+proj=moll +x_0=0 +y_0=0 +lat_0=0 +lon_0=0")

#Make bounding box to create a grid from ----
bbox <- st_bbox(c(xmin= -17596910, ymin= -6731548, xmax= 17596910, ymax= 8748765),
                crs = target_crs)

#Load global spatial Land file
world <- st_read("Species_Ranges/Data/Land/ne_50m_land_no_artic.shp",layer = "ne_50m_land_no_artic")%>%
  st_transform(crs = target_crs) #Project to equal area

#create a grid of 111000 size 
grid <- st_make_grid(bbox, cellsize = c(111000, 111000)) %>%
  st_sf() %>% 
  mutate(cell = 1:nrow(.)) 

# Get coordinates or the center of each grid cell
world_grid_centroids <- grid %>% 
  st_centroid() %>% 
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>% 
  as.data.frame() %>% 
  dplyr::select(-geometry)


############## Fish full non mining threatened species #######################
# load spatial data of species 
# DD 
Fish_ranges_DD <- st_read("Species_Ranges/Data/FISH/Fish_Species_ranges_IUCN/Actinopterygii_16.11.22/Actinopterygii_DD_4581sp/data_0.shp",
                          layer = 'data_0') %>% 
  st_transform(crs = target_crs)
Fish_ranges_DD  <- Fish_ranges_DD %>% 
  st_buffer(0)

glimpse(Fish_ranges_DD)

# Redlist 
Fish_ranges_RL <- st_read("Species_Ranges/Data/FISH/Fish_Species_ranges_IUCN/Actinopterygii_16.11.22/Actinopterygii_NT_VU_EN_CR_3844sp",
                          layer = 'data_0') %>% 
  st_transform(crs = target_crs)
Fish_ranges_RL  <- Fish_ranges_RL %>% 
  st_buffer(0)

# LC1
Fish_ranges_LC1 <- st_read("Species_Ranges/Data/FISH/Fish_Species_ranges_IUCN/Actinopterygii_16.11.22/Actinopterygii_LC_Charci_Cyptido_Pleur_Anguil_2247sp",
                           layer = 'data_0') %>% 
  st_transform(crs = target_crs)
Fish_ranges_LC1  <- Fish_ranges_LC1 %>% 
  st_buffer(0)

# LC2
Fish_ranges_LC2 <- st_read("Species_Ranges/Data/FISH/Fish_Species_ranges_IUCN/Actinopterygii_16.11.22/Actinopterygii_LC_perci_Cyprini_6214sp",
                           layer = c('data_0','data_1', 'data_2')) %>% 
  st_transform(crs = target_crs)
Fish_ranges_LC2  <- Fish_ranges_LC2 %>% 
  st_buffer(0)

# LC3
Fish_ranges_LC3 <- st_read("Species_Ranges/Data/FISH/Fish_Species_ranges_IUCN/Actinopterygii_16.11.22/Actinopterygii_LC_Siluri_Gobii_2350sp",
                           layer = 'data_0') %>% 
  st_transform(crs = target_crs)
Fish_ranges_LC3  <- Fish_ranges_LC3 %>% 
  st_buffer(0)

# LC4
Fish_ranges_LC4 <- st_read("Species_Ranges/Data/FISH/Fish_Species_ranges_IUCN/Actinopterygii_16.11.22/Actinopterygii_LC_Orders_exp_8largest_3682sp",
                           layer =  c('data_0','data_1', 'data_2')) %>% 
  st_transform(crs = target_crs)
Fish_ranges_LC4  <- Fish_ranges_LC4 %>% 
  st_buffer(0)

# LC4
Fish_ranges_NONActi <- st_read("Species_Ranges/Data/FISH/Fish_Species_ranges_IUCN/Non_Actinopterygii_24.04.23",
                           layer =  c('data_0')) %>% 
  st_transform(crs = target_crs) %>% 
  rename(BINOMIAL = SCI_NAME)

Fish_ranges_NONActi  <- Fish_ranges_NONActi %>% 
  st_buffer(0)

colnames(Fish_ranges_DD)  == colnames(Fish_ranges_NONActi) 

# join all multiple potentially overlapping ranges
fish_ranges <- rbind(Fish_ranges_DD, Fish_ranges_RL, Fish_ranges_LC1, Fish_ranges_LC2, Fish_ranges_LC3, Fish_ranges_LC4, Fish_ranges_NONActi) 

check <- st_is_valid(fish_ranges) 

fix_invalid_geom <- function(x) {
  print(x)
    fish_ranges_fix <- fish_ranges[x,] %>% 
      mutate(geometry == case_when(check[x]==FALSE ~ st_make_valid(geometry),
                                   TRUE ~ geometry))
    return(fish_ranges_fix)
}

sample <- fix_invalid_geom(3)

# check that length are matching 
length(check) == nrow(fish_ranges)

fish_ranges_fix <- bind_rows(lapply(1:nrow(fish_ranges), function = fix_invalid_geom ))

# CHECK that geometries are all TRUE
st_is_valid(fish_ranges_fix) %>% 
  unique()

# Reduce and rename cols 
fish_ranges_fix <- fish_ranges_fix %>% 
  select(binomial = BINOMIAL, seasonal = SEASONAL) %>% 
  rowid_to_column(var = "rowid")

# write into new dataframe 
st_write(fish_ranges_fix, file = "Species_Ranges/Data/FISH/Fish_Species_ranges_IUCN/Fish_ranges_fixed_geometries_moll.gpkg" )
glimpse(fish_ranges_fix)

## ***** Calcultate the Species cell areas on HPC stanage *****##### 
# SCRIPT; XX_All_Fish_cell_areas.R

# # Load all of the various data from back ups of the running loops. #####
SCA1 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_1000.csv")
SCA2 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_2000.csv")
SCA3 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_3000.csv")
SCA4 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_4000.csv")
SCA5 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_5000.csv")
SCA6 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_6000.csv")
SCA7 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_7000.csv")
SCA8 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_8000.csv")
SCA9 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_9000.csv")
SCA10 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_10000.csv")
SCA11 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_11000.csv")
SCA12 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_12000.csv")
SCA13 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_13000.csv")
SCA14 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_14000.csv")
SCA15 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_15000.csv")
SCA16 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_16000.csv")
SCA17 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_17000.csv")
SCA18 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_18000.csv")
SCA19 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_19000.csv")
SCA20 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_20000.csv")
SCA21 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_21000.csv")
SCA22 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_22000.csv")
SCA23 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_23000.csv")
SCA24 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_24000.csv")
SCA25 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_25000.csv")
SCA26 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_26000.csv")
SCA27 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_27000.csv")
SCAlast <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Fish_cellAreas_db3_17480.csv")

glimpse(SCA1)
glimpse(SCA2)
glimpse(SCA3)
glimpse(SCAlast)

length(unique(SCAlast$binomial))
range(SCA12$rowid)
length(unique(SCA13$binomial))
range(SCA13$rowid)
unique(setdiff(SCA12$binomial, SCA13$binomial))

# Bind all rows
SCA_all <- bind_rows(SCA1,
                     SCA2,
                     SCA3,
                     SCA4,
                     SCA5,
                     SCA6,
                     SCA7,
                     SCA8,
                     SCA9,
                     SCA10,
                     SCA11,
                     SCA12,
                     SCA13,
                     SCA14,
                     SCA15,
                     SCA16,
                     SCA17,
                     SCA18,
                     SCA19,
                     SCA20,
                     SCA21,
                     SCA22,
                     SCA23,
                     SCA24,
                     SCA25,
                     SCA26,
                     SCA27,
                     SCAlast)

# remove duplicates - these may have occured from starting of loops at different points repeating calculations for some species   
SCA_all <- SCA_all %>%
  distinct()

# save dataset of all cell areas 
write_csv(SCA_all, file = "IUCN_data/Species_Ranges/Outputs/Fish/fish_cellAreas_NONthrt_db3.csv")

# # sum any multiple distinct ranges within a single cell
fish_sp_cell_Values <- SCA_all %>%
  group_by(binomial,cell)%>%
  summarise(Areakm2 = sum(Areakm2))%>%
  # Add total species range column
  group_by(binomial) %>%
  mutate(total_Areakm2 = sum(Areakm2))%>%
  # calculate proportion of sp. range within each grid cell
  ungroup%>%
  mutate(Rng_prop_cell = Areakm2/total_Areakm2)
                                    
# # Create dataframe of sum of Rng_proportions_cell
fish_cell_threat_val <- fish_sp_cell_Values %>%
  group_by(cell) %>%
  summarise(threat_val = sum(Rng_prop_cell)) %>%
  full_join(world_grid_centroids, by  = "cell")%>%
  mutate(threat_val = if_else(is.na(threat_val), 0, threat_val))%>%
  arrange(-x)%>%
  dplyr::select(x,y,cell, threat_val)

write_csv(fish_cell_threat_val, "IUCN_data/Species_Ranges/Outputs/Fish/fish_cell_NONthrt_db3.csv")

# create grid shapefile with threat values 
fish_grid_ras <- left_join(grid, fish_cell_threat_val, by = "cell") %>% 
  dplyr::select(threat_val)

# save as shapefile 
st_write(fish_grid_ras, "IUCN_data/Species_Ranges/Outputs/Fish/Fish_NONtht_cert_raster3.gpkg")

#rasterize grid
fish_Nonthrt_raster_cert <- fish_grid_ras %>% 
  st_rasterize()
# write the raster of threats to all amphibian species with certainty values of the proportion each species ranges within raster pixels
write_stars(fish_Nonthrt_raster_cert, dsn =  "IUCN_data/Species_Ranges/Outputs/Fish/Fish_NONtht_cert_raster3.tif", overwrite = TRUE)

end_time <- Sys.time()
end_time - start_time

  tm_shape(fish_Nonthrt_raster_cert)+
  tm_raster(alpha = 0.7)


