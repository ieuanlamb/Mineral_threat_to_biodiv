#Additional impact area analysis. 

library(tidyverse)
library(terra)
library(sf)

# library(rgeos)
# library(rgbif)
library(rnaturalearth)
library(rnaturalearthdata)

set.seed(123)
getwd()
setwd("../IUCN_data")
sf_use_s2(FALSE)

# mollweide projection 
target_crs <- st_crs("+proj=moll +x_0=0 +y_0=0 +lat_0=0 +lon_0=0")


#Make bounding box to create a grid from ----
bbox <- st_bbox(c(xmin= -17596910, ymin= -6731548, xmax= 17596910, ymax= 8748765),
                crs = target_crs)


#Load global spatial Land file
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_make_valid() %>%
  filter(continent != "Antarctica")

#create a grid of 111000 size 
grid <- st_make_grid(bbox, cellsize = c(111000, 111000)) %>%
  st_sf() %>% 
  mutate(cell = 1:nrow(.))

# plot(grid, axes = T,pch = 3) # may look like a block if resolution is high


# Get coordinates or the center of each grid cell
world_grid_centroids <- grid %>% 
  st_centroid() %>% 
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>% 
  as.data.frame() %>% 
  dplyr::select(-geometry)


################################## All mining impacted bird species  ######################################
# read vertebrate taxonomy
tax <- read_csv("Species_Pages/Outputs/Chordata_taxonomy_01.03.22.csv")
glimpse(tax)
tax <- tax %>% 
  dplyr::select(internalTaxonId,scientificName, className) %>% 
  mutate(class = case_when(className == 'AMPHIBIA' ~ 'Amphibians',
                           className == 'AVES' ~ 'Birds',
                           className == 'MAMMALIA' ~ 'Mammals',
                           className == 'REPTILIA' ~ 'Reptiles',
                           TRUE ~ 'Fish'))

unique(tax$className)
unique(tax$class)

# read list of extinct species 
extinct_sp <- read_csv("Species_Pages/Outputs/Extinct_and_EW_sp2.csv")
extinct_sp <- extinct_sp %>% 
  pull(binomial)
# Read in the mine threat data
M_thrt <- st_read("Species_Ranges/Data/Mine_threatened_final/data_0.shp")

# read species that are impacted by mineral extraction
mining_impact_list <- read_csv("../Chapter_One2/Response_to_reviewers_CuBiol/mining_impact_list.csv")

impacted_birds <- mining_impact_list %>% 
  filter(className == "AVES")

# filter to only mining impacted species 
M_thrt <- M_thrt %>% 
  filter(SCI_NAME %in% impacted_birds$scientificName)
  
M_thrt <- M_thrt %>%
  st_transform(crs = target_crs)%>%
  rename(species = SCI_NAME, id_no = ID_NO)

M_thrt$SEASONAL %>% unique()

# join taxanomic group information
M_thrt <- M_thrt %>% 
  # remove species that are extinct
  filter(!species %in% extinct_sp,
         # remove passage ranges 
         !SEASONAL == 4)

# check the validity of shapes 
st_is_valid(M_thrt) %>% 
  unique()

#if not valid
M_thrt <- M_thrt %>%
  st_buffer(0)
i <- 3
# lapply function 
species_cell_area <- function(i){
  # find the cells that the range intersects 
  cell_intersect <- st_intersects(M_thrt[i,], grid, sparse = TRUE) %>%
    unlist()
  
  # extract only grid cells that the species intersects with
  sp_grid <- grid %>%
    filter(cell %in% cell_intersect)
  
  # find area of the species' range within the grids
  sp_cell_intersection <- st_intersection(M_thrt[i,], sp_grid)
  
  #calculate shape area within each cell
  sp_cell_Area  <- sp_cell_intersection %>%
    mutate(Areakm2 = units::set_units(sf::st_area(geometry),"km^2"),
           Areakm2 = as.numeric(Areakm2))%>%
    tibble()%>%
    dplyr::select(id_no, species,cell, Areakm2,cell)
  
  return(sp_cell_Area) 
  
}

start_t <- Sys.time()
Mine_thr_cell_A <- bind_rows(lapply(1:nrow(M_thrt), species_cell_area))
end_t <- Sys.time() - start_t

write_csv(Mine_thr_cell_A, "../Chapter_One2/Response_to_reviewers_CuBiol/Mining_impacted_sp_cell_areas.csv")
