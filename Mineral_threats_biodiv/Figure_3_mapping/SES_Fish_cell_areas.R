# SES Fish raster STANAGE

# fish SES script for HPC
library(readr)
library(dplyr)
library(tibble)
library(sf)
library(lqmm)

getwd()
# setwd("X:/edwards_lab1/User/bop21ipl/IUCN_data")
set.seed(1121995)

target_crs <- st_crs("+proj=moll +x_0=0 +y_0=0 +lat_0=0 +lon_0=0")

#Make bounding box to create a grid from ----
bbox <- st_bbox(c(xmin= -17596910, ymin= -6731548, xmax= 17596910, ymax= 8748765),
                crs = target_crs)


#Load global spatial Land file
# world <- st_read("data/Species_Ranges/Data/Land/ne_50m_land_no_artic.shp",layer = "ne_50m_land_no_artic")%>%
#   st_transform(crs = target_crs) #Project to equal area

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


# write into new dataframe 
# st_write(fish_ranges_fix, file = "data/Species_Ranges/Data/FISH/Fish_data/Species_Ranges_IUCN/Fish_ranges_fixed_geometries_moll.gpkg" )
# glimpse(fish_ranges_fix)

fish_ranges <- st_read(dsn = "data/Species_Ranges/Data/FISH/Fish_Species_ranges_IUCN/Fish_ranges_fixed_geometries_moll.gpkg" )


# remove species threatened by mining 
# read in list of mine threatened species and taxonomy
M_sp <- read_csv("data/Species_Pages/Chordata_Mine_threatened_assessments.csv")

# read list of extinct species 
extinct_sp <- read_csv("data/Species_Pages/Outputs/Extinct_and_EW_sp2.csv")
extinct_sp <- extinct_sp %>% 
  pull(binomial)

# remove mining threatened species from the dataset
fish_ranges <- fish_ranges %>% 
  filter(!binomial %in% M_sp,
         seasonal != 4,
         !binomial %in% extinct_sp,)

# reduce dataset to only necessary variables
fish_ranges <- fish_ranges  %>% 
  select(rowid, binomial, geometry = geom)

# reduce data to species remove species that have already been calculated
# fish_ranges <- fish_ranges[9001:10000,]

# remove species that area creating an error Gorgasia galzini
# fish_ranges <- fish_ranges %>% 
#   filter(binomial != "Gorgasia galzini")
  

# remove simplification as it causes an error 
{
# # simplyfy shapes 
# fish_ranges <- fish_ranges  %>% 
#   st_simplify(10000, preserveTopology = TRUE)
# 
# glimpse(fish_ranges)
# 
# check <- st_is_valid(fish_ranges)
# unique(check)
# fish_ranges_valid <- fish_ranges[check,]
# fish_ranges_invalid <-fish_ranges[which(check == FALSE),] %>% 
#   st_buffer(0)
# 
# fish_ranges <- rbind(fish_ranges_valid, fish_ranges_invalid) %>% 
#   arrange(rowid)
# st_is_valid(fish_ranges) %>% 
#   unique()
}

# save points
# save <- c(seq(1000, length(unique(fish_ranges$binomial)), by = 1000),length(unique(fish_ranges$binomial)))
# save_point <- fish_ranges$binomial[save]

# previous Run using for loop saved up to Altigena lippa: row = 2000, rowid = 2002
# remove these cells from the fish data


setdiff(fish_ranges$binomial, all_sp$binomial)

# function for calculating the cell values fore each species range 
sp_cell_area <- function(i) {
  species_data <- fish_ranges_grp[i,]
  
  print(c(i, species_data$binomial))
  # find the cells that the range intersects 
  fish_cell_intersect <- st_intersects(species_data, grid, sparse = TRUE) %>%
    unlist()
  
  # extract only grid cells that the species intersects with
  sp_grid <- grid %>%
    filter(cell %in% fish_cell_intersect)
  
  # find area of the species' range within the grids
  sp_cell_intersection <- st_intersection(species_data, sp_grid)
  
  #calculate shape area within each cell
  sp_cell_Area  <- sp_cell_intersection %>%
    mutate(Areakm2 = units::set_units(sf::st_area(geometry),"km^2"),
           Areakm2 = as.numeric(Areakm2))%>%
    tibble()%>%
    dplyr::select( binomial, rowid, cell, Areakm2)
  
  
  
  return(sp_cell_Area)
}

# process ranges in groups of species 
group_list <- seq(1000, nrow(fish_ranges), by = 1000)
# group_list <- seq(12000, nrow(fish_ranges), by = 1000) # #fist run on HPC saved to row 13000 but with erroneos looping  j-999:j not (j-999):j  
# group_list <- seq(7000, nrow(fish_ranges), by = 1000) # second run on HPC saved to row 9000
# group_list <- seq(1000, nrow(fish_ranges), by = 1000)

for(j in group_list) {
  fish_ranges_grp <- fish_ranges[(j-999):j,]
  
  # start from previous save point from the for loop
  cell_Areas_full <- bind_rows(lapply(1:nrow(fish_ranges_grp), sp_cell_area))
  write_csv(cell_Areas_full, paste0("data/Species_Ranges/Data/FISH/Fish_cellAreas_db3_10000.csv"))


}



