# creating Rasters of where the ranges of species with mineral extraction threat 
# using a loop to reduce computational effort.

library(raster)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
# library(tmap)

set.seed(123)
getwd()
setwd("IUCN_data")
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

# read in list of extinct species 
extinct_sp <- read_csv("Species_Pages/Outputs/Extinct_and_EW_sp2.csv")
# create list to remove
extinct_sp <- extinct_sp %>% 
  pull(binomial)

################################## All mining threatened species ######################################
# read vertebrate taxonomy
tax <- read_csv("Species_Pages/Outputs/Chordata_taxonomy_01.03.22.csv")
# change class names and group fish classes 
tax <- tax %>% 
  dplyr::select(internalTaxonId,scientificName, className) %>% 
  mutate(class = case_when(className == 'AMPHIBIA' ~ 'Amphibians',
                           className == 'AVES' ~ 'Birds',
                           className == 'MAMMALIA' ~ 'Mammals',
                           className == 'REPTILIA' ~ 'Reptiles',
                           TRUE ~ 'Fish'))


# --------- Read in the mine threat spatial range data ---------
# Downloaded from IUCN Red List https://www.iucnredlist.org/search 
# Using search options to select for: vertebrates with any of the following threats
# 3.1 Oil & gas drilling , 3.2 Mining and quarrying, 9.2.1 Oil spills, 9.2.2 Seepage from mining
M_thrt <- st_read("Species_Ranges/Data/Mine_threatened_final/data_0.shp")
# reproject to mollweide
M_thrt <- M_thrt %>%
  st_transform(crs = target_crs)%>%
  rename(species = SCI_NAME, id_no = ID_NO)
# join taxanomic group information
M_thrt <- M_thrt %>% 
  left_join(tax, by = c("species" = "scientificName")) %>% 
  # remove species that are extinct
  filter(!species %in% extinct_sp
         # remove passage ranges 
         !SEASONAL == 4)
  
# check all species are matched to tax names
unique(M_thrt$class)
# find unmatched species if present 
M_thrt %>% 
  filter(is.na(class)) %>% 
  pull(species) %>% 
  unique()

# split species into taxanomic groups
# check the validity of shapes 
st_is_valid(M_thrt) %>% 
  unique()
# if not: fix geometries
M_thrt <- M_thrt %>%
  st_buffer(0)

# --------- function to calculate the area of each species range within each grid square ------------  
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
# bind results into table of all species, what grid cells their ranges overlap, and the area of their range within each grid cell
Mine_thr_cell_A <- bind_rows(lapply(1:nrow(M_thrt), species_cell_area))
end_t <- Sys.time() - start_t


