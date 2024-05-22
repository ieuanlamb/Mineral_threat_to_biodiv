# Calculate the area of each species range with in each global grid cell for fish species with mineral extraction threats
# Essentially the same as 01_Mine_thrt_Areas.R however was separated to run the code on the HPC Stanage
# creating Rasters using a loop to reduce computational effort.
# Note: done in HPC Stanage as shape files for fish ranges were too large to run on desktop

library(sp)
library(sf)
library(tidyverse)
# library(tmap)

set.seed(123)
# turn s2 off
sf_use_s2(FALSE)

# mollweide projection 
target_crs <- st_crs("+proj=moll +x_0=0 +y_0=0 +lat_0=0 +lon_0=0")

#Make bounding box to create a grid from ----
bbox <- st_bbox(c(xmin= -17596910, ymin= -6731548, xmax= 17596910, ymax= 8748765),
                crs = target_crs)


#create a grid of 111000 size 
grid <- st_make_grid(bbox, cellsize = c(111000, 111000)) %>%
  st_sf() %>% 
  mutate(cell = 1:nrow(.))


####### read data ####
# read vertebrate taxonomy
# read list of extinct species 
extinct_sp <- read_csv("data/Species_Pages/Outputs/Extinct_and_EW_sp2.csv")
extinct_sp <- extinct_sp %>% 
  pull(binomial)

# Read in the mine threat species spatial range data
Sp_ranges <- st_read("data/Species_Ranges/Data/FISH/Mine_threatened_03.04.23/data_0.shp")

Sp_ranges <- Sp_ranges %>%  
  st_transform(crs = target_crs)%>%
  # remove passage seasonal ranges
  filter(SEASONAL != 4) %>% 
  # remove species that are extinct
  rename(species = SCI_NAME, id_no = ID_NO) %>% 
  filter(!species %in% extinct_sp) 


# check the validity of shapes 
st_is_valid(Sp_ranges) %>% 
  unique()
#if not valid
Sp_ranges <- Sp_ranges %>%
  st_buffer(0)

#  function to calculate all cells each species range overlaps with ====
species_cell_area <- function(i){
  print(i)
  print(Sys.time())
   # find the cells that the range intersects 
  cell_intersect <- st_intersects(Sp_ranges[i,], grid, sparse = TRUE) %>%
    unlist()
  
  # extract only grid cells that the species intersects with
  sp_grid <- grid %>%
    filter(cell %in% cell_intersect)
  
  # find area of the species' range within the grids
  sp_cell_intersection <- st_intersection(Sp_ranges[i,], sp_grid)
  
  #calculate shape area within each cell
  sp_cell_Area  <- sp_cell_intersection %>%
    mutate(Areakm2 = units::set_units(sf::st_area(geometry),"km^2"),
           Areakm2 = as.numeric(Areakm2))%>%
    tibble()%>%
    dplyr::select(id_no, species,cell, Areakm2,cell)
  
  return(sp_cell_Area) 
  
}

start_t <- Sys.time()
Mine_thr_cell_A <- bind_rows(lapply(1:nrow(Sp_ranges), species_cell_area))
end_t <- Sys.time() - start_t

write_csv(Mine_thr_cell_A, "data/Species_Ranges/Data/Fish_minetht_cell_areas2.csv")


# check all species have run 
Mine_thr_cell_A <- read_csv("data/Species_Ranges/Data/Fish_minetht_cell_areas2.csv")

# number of threatened species with ranges lengths should match
length(unique(Sp_ranges$species)) 
length(unique(Mine_thr_cell_A$species))  

# IF NOT RUN THE FOLLOWING 
# 
# # Find any species that are missing from the cell areas database
# list <- unique(Sp_ranges$species)
# list2 <- unique(Mine_thr_cell_A$species)
# 
# missing_sp <- setdiff(list,list2)
# 
# # pull the row id of missing species 
# rowid <- Sp_ranges %>% rowid_to_column("ID") %>% 
#   filter(species %in% missing_sp) %>% 
#   pull(ID)
# 
# # calculate cell areas of missing values 
# Mine_thr_cell_A_missing <- bind_rows(lapply(rowid, species_cell_area))
# 
# #bind to rest of cell areas data
# Mine_thr_cell_A <- bind_rows(Mine_thr_cell_A, Mine_thr_cell_A_missing)

#overite
write_csv(Mine_thr_cell_A, "data/Species_Ranges/Data/Fish_minetht_cell_areas2.csv")

