# SES threat raster for birds certainty raster calculation 

library(sp)
library(readr)
library(dplyr)
library(tibble)
library(sf)
library(lqmm)

getwd()
setwd("X:/edwards_lab1/User/bop21ipl")
set.seed(1121995)
sf_use_s2(FALSE)

target_crs <- st_crs("+proj=moll +x_0=0 +y_0=0 +lat_0=0 +lon_0=0")

#Make bounding box to create a grid from ----
bbox <- st_bbox(c(xmin= -17596910, ymin= -6731548, xmax= 17596910, ymax= 8748765),
                crs = target_crs)


#Load global spatial Land file
world <- st_read("IUCN_data/Species_Ranges/Data/Land/ne_50m_land_no_artic.shp",layer = "ne_50m_land_no_artic")%>%
  st_transform(crs =target_crs) #Project to equal area

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


############## Birds full non mining threatened species #######################

# time running time for 200 amphibians only terrestrial grids 
bird_db <- read_csv("IUCN_data/Species_Ranges/Data/BIRDS/Bird_ranges_noGeom.csv")

# read in list of mine threatened species and taxonomy
M_sp <- read_csv("IUCN_data/Species_Pages/Raw_Data/CHORDATA_pg_Oil_Mining_threat/Chordata_Mine_threatened_assessments.csv")

# read list of extinct species 
extinct_sp <- read.csv("IUCN_data/Species_Pages/Outputs/Extinct_and_EW_sp2.csv")
extinct_sp <- extinct_sp %>% 
  pull(binomial)

# use for loop for each row
# Attempt to fix geometry
# ensure multiple polygons
# Converting MULTISURFACE geometries to MULTIPOLYGON
# and make geometries valid
library(gdalUtilities)

ensure_multipolygons <- function(X) {
  tmp1 <- tempfile(fileext = ".gpkg")
  tmp2 <- tempfile(fileext = ".gpkg")
  st_write(X, tmp1)
  ogr2ogr(tmp1, tmp2, f = "GPKG", nlt = "MULTIPOLYGON")
  Y <- st_read(tmp2)
  st_sf(st_drop_geometry(X), geometry = st_geometry(Y))
}


# # Set up dataframe
# cell_Areas <- tibble(binomial = character(),
#                      cell = integer(),
#                      Areakm2 = double())

sp_cell_area <- function(i) {
  # find the cells that the range intersects 
  save <- c(20, seq(500, nrow(bird_ranges), by = 500), nrow(bird_ranges))
  species_data <- bird_ranges[i,]
  check <- st_is_valid(species_data)
  if(check == FALSE) {
    species_data <- st_make_valid(species_data)
  }
  
  bird_cell_intersect <- st_intersects(species_data, grid, sparse = TRUE) %>%
    unlist()
  
  # extract only grid cells that the species intersects with
  sp_grid <- grid %>%
    filter(cell %in% bird_cell_intersect)
  
  # find area of the species' range within the grids
  sp_cell_intersection <- st_intersection(species_data, sp_grid)
  
  #calculate shape area within each cell
  sp_cell_Area  <- sp_cell_intersection %>%
    mutate(Areakm2 = units::set_units(sf::st_area(geometry),"km^2"),
           Areakm2 = as.numeric(Areakm2))%>%
    tibble()%>%
    dplyr::select( binomial, rowid, cell, Areakm2)
  
  # cell_Areas <- bind_rows(cell_Areas, sp_cell_Area)
  
  if(i %in% save){
    # write_csv(cell_Areas, paste0("IUCN_data/Species_Ranges/Outputs/Backups/Bird_cellAreas2_grp",j,"db_",i,".csv"))
    print(i)
    end_time <- Sys.time()
    print("processing time to save point")
    print(end_time - start_time2)
  }
  
  return(sp_cell_Area)
}

# groups of Bird range data by file. 

groups <- c("1_3999" , "4000_7999", "8000_11999", "12000_15999", "16000_17572") # first HPC run saves groups 1 and 2
# groups <- c("8000_11999", "12000_15999", "16000_17572")
# groups <- c("12000_15999", "16000_17572")
for(j in groups){ 
  
  start_time <- Sys.time()
  print(paste0(start_time,"start loop ",j))
  # second attempt using ssubsets of bird ranges in Shapefiles
  bird_ranges <- read_sf(paste0("data/BIRDS/Bird_Ranges_",j,".gpkg" )) %>%
        st_transform(crs = target_crs)
  glimpse(bird_ranges)
 
  # bird_ranges <- bird_ranges[1:100,] # Time for 100  = 9 mins
  # bird_ranges <- bird_ranges[1:10,] # 
  bird_ranges <- bird_ranges %>% 
    rowid_to_column(var = "rowid")
  
  # remove mining threatened species from the dataset
  bird_ranges <- bird_ranges %>% 
    filter(!binomial %in% M_sp,
           seasonal != 4,
           !binomial %in% extinct_sp)
  bird_ranges <- bird_ranges %>%
    select(rowid, SISID, binomial, geometry = geom)

  bird_ranges <- ensure_multipolygons(bird_ranges) 
  
  bird_ranges <- bird_ranges %>% 
    st_simplify(10000, preserveTopology = TRUE)
  
  # data wrangling time
  end_time <- Sys.time()
  print("data wrangling time")
  print(end_time - start_time)
  
  # time running time for processing
  start_time2 <- Sys.time()
  
  cell_Areas_full <- bind_rows(lapply(1:nrow(bird_ranges), sp_cell_area))
  write.csv(cell_Areas_full, paste0("IUCN_data/Species_Ranges/Outputs/Backups/Bird_cellAreas2_db_",j,".csv"))
  
  # processing time
  end_time <- Sys.time()
  print("processing time for j")
  print(end_time - start_time2)

  print(j)
  end_time <- Sys.time()
  print("whole loop")
  print(end_time - start_time)
}


## compiling cell threat values data sets to create raster #####
# not for the HPC
library(raster)
library(tmap)

# load all data sets
bird_sp_cell_1 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Bird_cellAreas2_db_1_3999.csv")
bird_sp_cell_2 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Bird_cellAreas2_db_4000_7999.csv")
bird_sp_cell_3 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Bird_cellAreas2_db_8000_11999.csv")
bird_sp_cell_4 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Bird_cellAreas2_db_12000_15999.csv")
bird_sp_cell_5 <- read_csv("IUCN_data/Species_Ranges/Outputs/Backups/Bird_cellAreas2_db_16000_17572.csv")

cell_Areas <- bind_rows(bird_sp_cell_1, bird_sp_cell_2, bird_sp_cell_3, bird_sp_cell_4, bird_sp_cell_5) %>%
  select(-rowid) %>%
  distinct() # removes 1,060 rows

length(unique(cell_Areas$binomial)) #10990

# # sum any multiple distinct ranges within a single cell
bird_sp_cell_Values <- cell_Areas %>%
  group_by(binomial,cell) %>%
  summarise(Areakm2 = sum(Areakm2))%>%
  # Add total species range column
  group_by(binomial) %>%
  mutate(total_Areakm2 = sum(Areakm2))%>%
  # calculate proportion of sp. range within each grid cell
  ungroup %>%
  mutate(Rng_prop_cell = Areakm2/total_Areakm2)

# # Create dataframe of sum of Rng_proportions_cell
bird_cell_threat_val <- bird_sp_cell_Values %>%
  group_by(cell) %>%
  summarise(threat_val = sum(Rng_prop_cell)) %>%
  full_join(world_grid_centroids, by  = "cell")%>%
  mutate(threat_val = if_else(is.na(threat_val), 0, threat_val)) %>%
  arrange(-x)%>%
  dplyr::select(cell,x,y, threat_val) %>% 
  arrange(cell)

range(bird_cell_threat_val$threat_val)

write_csv(bird_cell_threat_val, paste0("IUCN_data/Species_Ranges/Outputs/Birds/bird_cell_NONthrt_db2_full.csv"))


# # create a raster of the coordinates
# bird_NONthreat_ras <- rasterFromXYZ(bird_cell_threat_val, res = 111000, crs = "+proj=cea +datum=WGS84")
# 
# writeRaster(bird_NONthreat_ras, filename = "IUCN_data/Species_Ranges/Outputs/Birds/bird_NONtht_cert_raster.tif", overwrite = TRUE)


bird_grid_ras <- left_join(grid, bird_cell_threat_val, by = "cell") %>% 
  dplyr::select(threat_val)

library("RColorBrewer")
# visualise raster
tm_shape(world)+
  tm_polygons()+
  tm_shape(bird_grid_ras)+
  tm_fill(col = "threat_val", alpha = 0.7, palette = "YlGnBu")

st_write(bird_grid_ras, "IUCN_data/Species_Ranges/Outputs/Birds/Bird_NONtht_cert_raster3.gpkg", append = FALSE)

#rasterize grid
bird_thrt_raster_cert <- bird_grid_ras %>% 
  st_rasterize()

tm_shape(world)+
  tm_polygons()+
  tm_shape(bird_thrt_raster_cert)+
  tm_raster(col = "threat_val", alpha = 0.7, palette = "YlGnBu")

tm_shape(world)+
  tm_fill(col = "#93ccdb")+
  tm_shape(bird_grid_ras)+
  tm_fill(col = "threat_val", palette = viridis(30, alpha = 0.6, option = "E"),
          style = "kmeans"
  )+
  tm_layout (frame = FALSE, bg.color = "transparent", main.title = "Birds threat value" , legend.title.color = "transparent")

# write the raster of threats to all amphibian species with certainty values of the proportion each species ranges within raster pixels
write_stars(bird_thrt_raster_cert, dsn = "IUCN_data/Species_Ranges/Outputs/Birds/Bird_NONtht_cert_raster3.tif", overwrite = TRUE)



