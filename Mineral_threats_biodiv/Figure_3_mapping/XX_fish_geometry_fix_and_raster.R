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
# Downloaded from IUCN species ranges database https://www.iucnredlist.org/search in separate taxanomic classes; 
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


