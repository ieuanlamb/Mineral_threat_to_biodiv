#### Rasterise the grids from Cell Area databases

library(tidyverse)
library(terra)
library(tidyterra)
library(sf)
library(viridis)
# library(rgeos)
# library(rgbif)
library(ggtext)
library(tmap)
library(rnaturalearth)
library(rnaturalearthdata)
library(stars)
library(lqmm)
library("RColorBrewer")

getwd()
setwd("X:/edwards_lab1/User/bop21ipl/IUCN_data")
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

write_sf(grid, "../Chapter_One2/Response_to_reviewers_CuBiol/methods_fig/global_grid.gpkg")

# Get coordinates or the center of each grid cell
world_grid_centroids <- grid %>% 
  st_centroid() %>% 
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>% 
  as.data.frame() %>% 
  dplyr::select(-geometry)


####### load cell areas -----
# cell areas from the HPC run of above code

#Birds
bird_cell_Areas <- read_csv("../Chapter_One2/Response_to_reviewers_CuBiol/Mining_impacted_sp_cell_areas.csv")
length(unique(bird_cell_Areas$species)) # 362
nrow(bird_cell_Areas)
bird_cell_Areas[(duplicated(bird_cell_Areas)),]
bird_cell_Areas <- bird_cell_Areas %>% 
  distinct()
nrow(bird_cell_Areas)


#### Rasterize ====

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
  tm_fill(fill = "threat_val", alpha = 0.7, palette = "YlGnBu"
    )

st_write(bird_grid_ras, "../Chapter_One2/Response_to_reviewers_CuBiol/bird_impacted_ras.gpkg")

#rasterize grid
bird_thrt_raster_cert <- bird_grid_ras %>% 
  st_rasterize()

# write the raster of threats to all amphibian species with certainty values of the proportion each species ranges within raster pixels
write_stars(bird_thrt_raster_cert, dsn = "../Chapter_One2/Response_to_reviewers_CuBiol/bird_impacted_ras.tif", overwrite = TRUE)

bird_thrt_raster_cert

# raster creation =====
three_col <- c( alpha("#5566AA",0.5), "#FFFFFF","#FFFF00","#CC2222")
#"#edcda4",

colfunc <- colorRampPalette(c("#93ccdb", #6f98a3",  # grey blue  # "#8b9cdd", # light blue
                                       "#20092e", # dark purple
                                       "#e73027", # red 
                                       "#ebe030", # yellow
                                       "#ffff9f"  # light yellow
))
colfunc1 <- colorRampPalette(c("#4575b5ff", "#f5fad2ff", "#d73027ff"))

# load as raster =====
bird_imp_raster_cert <- terra::rast( "../Chapter_One2/Response_to_reviewers_CuBiol/bird_impacted_ras.tif")

ter_bird <- terra::crop(bird_imp_raster_cert, world)
ter_bird <- terra::mask(ter_bird, world)
plot(ter_bird_n0)

ter_bird_n0 <- ter_bird
values(ter_bird_n0) <- ifelse(values(ter_bird_n0) == 0, NA, values(ter_bird_n0))

ter_bird_n0$bird_impacted_ras

Bird_thrt <- tm_shape(world)+
  tm_polygons(fill = "grey80",
              col = "#00000000",)+
  tm_shape(ter_bird_n0)+
  tm_raster(col = "bird_impacted_ras", 
            # col.scale =  tm_scale_continuous(values = viridis(200,direction = 1,  alpha = 1, begin = 0.2, option = "D"), value.na = "#000000"),
            col.scale =  tm_scale_continuous(values = viridis(200,direction = -1,  alpha = 1, begin = 0, option = "A"), value.na = "#00000000"),
            col.legend = tm_legend(title = "", 
                                   title.size = 0.9,
                                   reverse = T, 
                                   frame.lwd = 0,
                                   position = c(-0.03, 1)))+
  tm_layout(frame = F)
          
Bird_thrt
tmap_save(Bird_thrt, file = "../Chapter_One2/Response_to_reviewers_CuBiol/Bird_impacted_3.jpeg", width = 10,height = 4, dpi = 1000)

### SES map ====
bird_raster <- rast("Species_Ranges/Outputs/Birds/Bird_NONtht_cert_raster3.tif")
# # bird threat certaity raster
bird_thrt_raster <- rast("Species_Ranges/Outputs/Birds/Bird_tht_cert_raster3.tif")

bird_imp_raster_cert
bird_raster
bird_thrt_raster

bird_avrg <- (bird_imp_raster_cert/(bird_raster+bird_thrt_raster))
# give value to all cells even if no threat exists
values(bird_avrg)[is.na(values(bird_avrg))] <- 0 #<- if_else(is.na(values(bird_avrg)), 0 , values(bird_avrg))
values(bird_avrg)[is.na(values(bird_avrg))] #<- if_else(is.na(values(bird_avrg)), 0 , values(bird_avrg))

# reduce to terrestrial
ter_bird_avrg <- crop(bird_avrg, world)
ter_bird_avrg <- mask(ter_bird_avrg, world)

# Bird map
ter_birdSES_n0 <- ter_bird_avrg
values(ter_birdSES_n0) <- ifelse(values(ter_birdSES_n0) == 0, NA, values(ter_birdSES_n0))


bird_SES <- tm_shape(world)+
  tm_polygons(fill = "grey80",
              col = "#00000000",
              )+
  tm_shape(ter_birdSES_n0)+
  tm_raster(col = "bird_impacted_ras",
            # col.scale =  tm_scale_continuous(values = viridis(200,direction = 1,  alpha = 1, begin = 0.2, option = "D"), value.na = "#000000"),
            col.scale =  tm_scale_continuous(values = viridis(200,direction = -1,  alpha = 1, begin = 0, option = "A"), value.na = "#00000000"),
            col.legend = tm_legend(title = "", 
                                   title.size = 0.9,
                                   reverse = T, 
                                   frame.lwd = 0,
                                   position = c(-0.03, 1)
                                   )
  )+
  tm_layout(frame = F)
bird_SES

tmap_save(bird_SES, file = "../Chapter_One2/Response_to_reviewers_CuBiol/Scaled_bird_impact_map_2.jpeg", width = 10,height = 4, dpi = 1000)


