# Plotting the hotspot maps for species threatened by mining
#
# using rasters for Amphibians, Birds, Fish, Mammals, Reptiles from script Loop_for_cert_rasters 
#  ieuan lamb
# 
library(tidyverse)
library(sp)
library(raster)
library(sf)
library(cowplot)
library(viridis)
library(readr)
# library(sdmpredictors)
# library(rgdal)
library(tmap)
library(rnaturalearth)
library(rnaturalearthdata)
library(terra)

getwd()
setwd("X:/edwards_lab1/User/bop21ipl/IUCN_data/Species_Ranges/Outputs")

sf_use_s2(F)

# loading and collating can be bypassed by reading the All_tht_cert_Stack SEE plotting raster (not working yet)
# plotting colours ----
# three_col <- c( alpha("#6666AA",0.8), "#FFFFFF","#FFFF00","#FF0000")
three_col <- c( alpha("#5566AA",0.5), "#FFFFFF","#FFFF00","#CC2222")
#"#edcda4",
three_col <- c("#5566AA80", "#FFFFFF","#FFFF00","#CC2222")


colfunc <- colorRampPalette(c("#93ccdb", #6f98a3",  # grey blue  # "#8b9cdd", # light blue
                              "#20092e", # dark purple
                              "#e73027", # red 
                              "#ebe030", # yellow
                              "#ffff9f"  # light yellow
                              ))
colfunc1 <- colorRampPalette(c("#4575b5ff", "#f5fad2ff", "#d73027ff"))


#Load global spatial Land file ----
# world <- st_read("Species_Ranges/Data/Land/ne_50m_land_no_artic.shp",layer = "ne_50m_land_no_artic")%>%
#   st_transform(crs ="+proj=cea +datum=WGS84") #Project to equal area


target_crs <- st_crs("+proj=moll +x_0=0 +y_0=0 +lat_0=0 +lon_0=0")

#Load global spatial Land file
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_make_valid() %>%
  filter(continent != "Antarctica") %>% 
  st_transform(crs = target_crs)



# Load raster of threat certainty rasters ====
# where raster cell values = sum of the proportion of each species' range within the cell 
# amphibian threat certaity raster
amph_raster <- rast("Amphibians/Amph_tht_cert_raster3.tif")
ter_amph <- terra::crop(amph_raster, world)
ter_amph <- terra::mask(ter_amph, world)

# # bird threat certaity raster
bird_raster <- rast("Birds/Bird_tht_cert_raster3.tif")
ter_bird <- terra::crop(bird_raster, world)
ter_bird <- terra::mask(ter_bird, world)

# # fish threat certainty raster
fish_raster <- rast("Fish/Fish_tht_cert_raster3.tif")
ter_fish <- terra::crop(fish_raster, world)
ter_fish <- terra::mask(ter_fish, world)

# Mammal threat certainty raster
mam_raster <- rast("Mammals/Mam_tht_cert_raster3.tif")
ter_mam <- terra::crop(mam_raster, world)
ter_mam <- terra::mask(ter_mam, world)

# # Reptiles threat certainty raster
rep_raster <- rast("Reptiles/Rep_tht_cert_raster3.tif")
ter_rep <- terra::crop(rep_raster, world)
ter_rep <- terra::mask(ter_rep, world)


# Collate into one raster ====
all_sp_raster <- rast(list(amph_raster, 
                       bird_raster,
                       fish_raster,
                       mam_raster, 
                       rep_raster))
# 
# # sum the scores of each raster
# all_sp_raster$sum_all_sp_tht_cert_ras <- sum(amph_raster, 
#                                              bird_raster, 
#                                              fish_raster, 
#                                              mam_raster, 
#                                              rep_raster) 


sum_all_sp_tht_cert_ras <- sum(all_sp_raster$Amph_tht_cert_raster3, 
                               all_sp_raster$Bird_tht_cert_raster3, 
                               all_sp_raster$Fish_tht_cert_raster3, 
                               all_sp_raster$Mam_tht_cert_raster3, 
                               all_sp_raster$Rep_tht_cert_raster3 )

names(sum_all_sp_tht_cert_ras) <-  c("sum_all_sp_tht_cert_ras")

# write raster 
writeRaster(sum_all_sp_tht_cert_ras, filename = "All_tht_sum3.tiff", overwrite = TRUE)

# standardise
stnd_all_sp_thrt <- sum_all_sp_tht_cert_ras 
values(stnd_all_sp_thrt) <- scale(log(values(stnd_all_sp_thrt)+1), center = T,scale = T) 

# Terrestrial land mass
ter_Allsp <- terra::crop(stnd_all_sp_thrt, world)
ter_Allsp <- terra::mask(ter_Allsp, world)

#raw values 
ter_Allsp2 <- terra::crop(sum_all_sp_tht_cert_ras, world)
ter_Allsp2 <- terra::mask(ter_Allsp2, world)


#raw values 
all_sp_raster <- c(all_sp_raster,sum_all_sp_tht_cert_ras)
ter_all_sp_raster <- terra::crop(all_sp_raster, world)
ter_all_sp_raster <- terra::mask(ter_all_sp_raster, world)


tibble(values = values(stnd_all_sp_thrt)) %>% 
  ggplot(aes(x = values)) +
  geom_density()

# Plotting raster maps ====
ter_Allsp_n0 <- ter_Allsp2
# remove 0 values 
values(ter_Allsp_n0) <- ifelse(values(ter_Allsp_n0) == 0, NA, values(ter_Allsp_n0))

# All species
Allsp_thrt <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_Allsp_n0)+
  tm_raster(palette = viridis(200,direction = -1,  alpha = 0.6, option = "A"),
            style = "cont")+
  tm_layout (frame = FALSE, 
             legend.title.color = "transparent", legend.frame.lwd = 0
             )
Allsp_thrt
tmap_save(Allsp_thrt, file = "X:/edwards_lab1/User/bop21ipl/IUCN_data/Figures/Threat_value_viridis/All_species_logstd_viridis_ter_V3.jpeg", width = 8, height = 5, dpi = 1000)

# Amphibian threat value
ter_amph_n0 <- ter_amph
values(ter_amph_n0) <- ifelse(values(ter_amph_n0) == 0, NA, values(ter_amph_n0))

range(values(ter_amph_n0), na.rm = T)

tibble(values = values(ter_amph_logstd)) %>% 
  ggplot(aes(x = values)) +
  geom_density()

Amph_thrt <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_amph_n0)+
  tm_raster(palette = viridis(200,direction = -1,  alpha = 0.6, option = "A"),
            style = "cont")+
  tm_layout(frame = FALSE, 
             legend.title.color = "transparent"
            )

Amph_thrt
tmap_save(Amph_thrt, file = "X:/edwards_lab1/User/bop21ipl/IUCN_data/Figures/Threat_value_viridis/Amphibian_viridis_ter_V3.jpeg", width = 8,height = 5, dpi = 1000)

# Bird threat value
ter_bird_n0 <- ter_bird
values(ter_bird_n0) <- ifelse(values(ter_bird_n0) == 0, NA, values(ter_bird_n0))

Bird_thrt <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_bird_n0)+
  tm_raster(palette = viridis(200,direction = -1,  alpha = 0.6, option = "A"),
            style = "cont")+
  tm_layout(frame = FALSE, 
            legend.title.color = "transparent"
  )
Bird_thrt
tmap_save(Bird_thrt, file = "X:/edwards_lab1/User/bop21ipl/IUCN_data/Figures/Threat_value_viridis/Bird_viridis_ter_V3.jpeg", width = 8,height = 5, dpi = 1000)

# Fish threat value
ter_fish_n0 <- ter_fish
values(ter_fish_n0) <- ifelse(values(ter_fish_n0) == 0, NA, values(ter_fish_n0))

Fish_thrt <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_fish_n0)+
  tm_raster(palette = viridis(200,direction = -1,  alpha = 0.6, option = "A"),
            style = "cont")+
  tm_layout(frame = FALSE, 
            legend.title.color = "transparent"
  )
Fish_thrt
tmap_save(Fish_thrt, file = "X:/edwards_lab1/User/bop21ipl/IUCN_data/Figures/Threat_value_viridis/Fish_viridis_ter_V3.jpeg", width = 8,height = 5, dpi = 1000)

# Mammal threat value
ter_mam_n0 <- ter_mam
values(ter_mam_n0) <- ifelse(values(ter_mam_n0) == 0, NA, values(ter_mam_n0))

Mam_thrt <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_mam_n0)+
  tm_raster(palette = viridis(200,direction = -1,  alpha = 0.6, option = "A"),
            style = "cont")+
  tm_layout(frame = FALSE, 
            legend.title.color = "transparent")
Mam_thrt

tmap_save(Mam_thrt, file = "X:/edwards_lab1/User/bop21ipl/IUCN_data/Figures/Threat_value_viridis/Mammal_viridis_ter_V3.jpeg", width = 8,height = 5, dpi = 1000)

# Reptile threat value
ter_rep_n0 <- ter_rep
values(ter_rep_n0) <- ifelse(values(ter_rep_n0) == 0, NA, values(ter_rep_n0))

Rep_thrt <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_rep_n0)+
  tm_raster(palette = viridis(200,direction = -1,  alpha = 0.6, option = "A"),
            style = "cont")+
  tm_layout(frame = FALSE, 
            legend.title.color = "transparent", legend.frame.lwd = 0)
Rep_thrt
tmap_save(Rep_thrt, file = "X:/edwards_lab1/User/bop21ipl/IUCN_data/Figures/Threat_value_viridis/Reptile_viridis_ter_V3.jpeg", width = 8,height = 5, dpi = 1000)

# combing plots 
library(cowplot)

map_grid <- tmap_arrange(Allsp_thrt, Amph_thrt, Bird_thrt,Fish_thrt, Mam_thrt, Rep_thrt, ncol = 1 )
map_grid
tmap_save(map_grid, file = "X:/edwards_lab1/User/bop21ipl/IUCN_data/Figures/Threat_value_viridis/All_map_grid_viridis_ter_V3.jpeg", 
          width = 8, height = 20, dpi = 1000)



# METHOD 2 (USE THIS METHOD) taking 1% and 5% of terrestrial cells but keeping layer of threatened cells ======
# ALL SPECIES
ter_Allsp
ter_Allsp_H <- ter_Allsp

# method 2 taking 1% and 5% of terrestrial cells but keeping layer of threatened cells 
values( ter_Allsp_H ) <- ifelse(values( ter_Allsp_H ) == 0, NA, ifelse(values(ter_Allsp_H) <=  quantile(values(ter_Allsp_H) , prob = 0.95, na.rm = T),0, 
                                                                       ifelse(values(ter_Allsp_H) >  quantile( values(ter_Allsp_H), prob = 0.99, na.rm = T),2,1)))
plot(ter_Allsp_H)

writeRaster(ter_Allsp_H, filename = "All_threat_cert_hotspots.tiff")

# For amphibians 
ter_amph
ter_amph_H1 <- ter_amph
ter_amph_H2 <- ter_amph
# Method 1 example - taking 1% and 5% of threatened cells:
{values( ter_amph_H1 ) <- ifelse(values( ter_amph_H1 ) == 0, NA ,values( ter_amph_H1 ))
values(ter_amph_H1) <- ifelse(values(ter_amph_H1) <=  quantile( ter_amph_H1 , prob = 0.95),0, 
                              ifelse(values(ter_amph_H1) >  quantile( ter_amph_H1, prob = 0.99),2,1))

plot( ter_amph_H1 )
}
# method 2 taking 1% and 5% of terrestrial cells but keeping layer of threatened cells 
values( ter_amph_H2 ) <- ifelse(values( ter_amph_H2 ) == 0, NA, ifelse(values(ter_amph_H2) <=  quantile(values(ter_amph_H2) , prob = 0.95, na.rm = T),0, 
                                                                       ifelse(values(ter_amph_H2) >  quantile( values(ter_amph_H2), prob = 0.99, na.rm = T),2,1)))

plot( ter_amph_H2 )

# Birds 
ter_bird_H <- ter_bird
values( ter_bird_H ) <- ifelse(values( ter_bird_H ) == 0, NA, ifelse(values(ter_bird_H) <=  quantile(values(ter_bird_H) , prob = 0.95, na.rm = T),0, 
                                                                       ifelse(values(ter_bird_H) >  quantile( values(ter_bird_H), prob = 0.99, na.rm = T),2,1)))

plot( ter_bird_H )

# Fish 
ter_fish_H <- ter_fish
values( ter_fish_H ) <- ifelse(values( ter_fish_H ) == 0, NA, ifelse(values(ter_fish_H) <=  quantile(values(ter_fish_H) , prob = 0.95, na.rm = T),0, 
                                                                     ifelse(values(ter_fish_H) >  quantile( values(ter_fish_H), prob = 0.99, na.rm = T),2,1)))

plot( ter_fish_H )

# Mammals 
ter_mam_H <- ter_mam
values( ter_mam_H ) <- ifelse(values( ter_mam_H ) == 0, NA, ifelse(values(ter_mam_H) <=  quantile(values(ter_mam_H) , prob = 0.95, na.rm = T),0, 
                                                                     ifelse(values(ter_mam_H) >  quantile( values(ter_mam_H), prob = 0.99, na.rm = T),2,1)))

plot( ter_mam_H )

# Reptiles 
ter_rep_H <- ter_rep
values( ter_rep_H ) <- ifelse(values( ter_rep_H ) == 0, NA, ifelse(values(ter_rep_H) <=  quantile(values(ter_rep_H) , prob = 0.95, na.rm = T),0, 
                                                                     ifelse(values(ter_rep_H) >  quantile( values(ter_rep_H), prob = 0.99, na.rm = T),2,1)))

plot( ter_rep_H )

# Plotting terestrial hotspots METHOD 2  ----

break_label <- c("Species with MET","", "Top 5%","Top 1%")

# All species summed
allsphotspot <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_Allsp_H )+
  tm_raster(palette = three_col,
            labels = break_label,
            legend.reverse = TRUE)+
  tm_layout(frame = FALSE,
            # bg.color = "transparent",
            # main.title = "Vertebrates threat value hotspots",
            legend.title.color = "transparent",
            legend.show = FALSE
            )
allsphotspot
tmap_save(allsphotspot, file = "X:/edwards_lab1/User/bop21ipl/IUCN_data/Figures/Threat_hotspots_method2/All_sp_hotspot_ter2_no_lab.jpeg", width = 8,height = 5, dpi = 1000)

# Amphibians
amphhotspot <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_amph_H2 )+
  tm_raster(palette = three_col,
            labels = break_label,
            legend.reverse = TRUE)+
  tm_layout(frame = FALSE,
            # bg.color = "transparent",
            # main.title = "Amphibian threat value hotspots",
            legend.title.color = "transparent",
            legend.show = FALSE
            )
amphhotspot
# tmap_save(amphhotspot, file = "X:/edwards_lab1/User/bop21ipl/IUCN_data/Figures/Threat_hotspots_method2/Amph_sphotspot_ter2_no_lab.jpeg", width = 8,height = 5, dpi = 1000)
dev.off()
getwd()

# Birds 
birdsphotspot <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_bird_H )+
  tm_raster(palette = three_col,
            labels = break_label,
            legend.reverse = TRUE)+
  tm_layout(frame = FALSE,
            # bg.color = "transparent",
            # main.title = "Vertebrates threat value hotspots",
            legend.title.color = "transparent",
            legend.show = FALSE)
birdsphotspot
tmap_save(birdsphotspot, file = "X:/edwards_lab1/User/bop21ipl/IUCN_data/Figures/Threat_hotspots_method2/Bird_sphotspot_ter2_no_lab.jpeg", width = 8,height = 5, dpi = 1000)

# Fish
fishsphotspot <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_fish_H )+
  tm_raster(palette = three_col,
            labels = break_label,
            legend.reverse = TRUE)+
  tm_layout(frame = FALSE,
            # bg.color = "transparent",
            # main.title = "Vertebrates threat value hotspots",
            legend.title.color = "transparent",
            legend.show = FALSE)
fishsphotspot
tmap_save(fishsphotspot, file = "X:/edwards_lab1/User/bop21ipl/IUCN_data/Figures/Threat_hotspots_method2/Fish_sphotspot_ter2_no_lab.jpeg", width = 8,height = 5, dpi = 1000)

# Reptiles
repsphotspot <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_rep_H )+
  tm_raster(palette = three_col,
            labels = break_label,
            legend.reverse = TRUE)+
  tm_layout(frame = FALSE,
            # bg.color = "transparent",
            # main.title = "Vertebrates threat value hotspots",
            legend.title.color = "transparent",
            legend.show = FALSE)
repsphotspot
tmap_save(repsphotspot, file = "X:/edwards_lab1/User/bop21ipl/IUCN_data/Figures/Threat_hotspots_method2/Reptile_sphotspot_ter2_no_lab.jpeg", width = 8,height = 5, dpi = 1000)

# Mammals
mamsphotspot <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_mam_H )+
  tm_raster(palette = three_col,
            labels = break_label,
            legend.reverse = TRUE)+
  tm_layout(frame = FALSE,
            # bg.color = "transparent",
            # main.title = "Vertebrates threat value hotspots",
            legend.title.color = "transparent",
            legend.show = FALSE
            )
mamsphotspot
tmap_save(mamsphotspot, file = "X:/edwards_lab1/User/bop21ipl/IUCN_data/Figures/Threat_hotspots_method2/Mammal_sphotspot_ter2_no_lab.jpeg", width = 8,height = 5, dpi = 1000)

# plot only the legend horizontally
legend_only <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_mam_H )+
  tm_raster(palette = three_col,
            labels = break_label,
            legend.reverse = TRUE,
            legend.is.portrait = FALSE
            )+
  tm_layout(frame = FALSE, legend.frame.lwd = 0,
            legend.title.size = 0,
            legend.title.color = "transparent",
            legend.stack = "horizontal",
            legend.only = TRUE,
  )

legend_only

tmap_save(legend_only, file = "X:/edwards_lab1/User/bop21ipl/IUCN_data/Figures/Threat_hotspots_method2/Legend_only_Horz2.jpeg", width = 3,height = 1, dpi = 1000)




# Regional heatmaps ----
#South America
S_Amerbbox <-  st_bbox(c(xmin = -11000000, ymin = -4000000,xmax = -3800000, ymax = 3000000))

tm_shape(world, bbox= st_bbox(S_Amerbbox))+
  tm_fill(col = "black")

S_Amer_allsp <- tm_shape(world, bbox= st_bbox(S_Amerbbox))+
  tm_fill(col = "#93ccdb")+
  # tm_borders(col = "#555555")+
  tm_shape(ter_all_sp_raster$sum_all_sp_tht_cert_ras, bbox= st_bbox(S_Amerbbox))+
  tm_raster(palette =  colfunc(25),
            style = "cont",
            title = "Threat Value")+  
  tm_layout(bg.color = "white")

tmap_save(S_Amer_allsp, file = "Figures/Final_maps_tmap/S_Amer_All_sp_value.jpeg", width = 7,height = 7, dpi = 1000)

# Africa 
Africabbox <-  st_bbox(c( xmin = -2000000, xmax = 6000000, ymin = -4000000, ymax = 3000000))

tm_shape(world , bbox= st_bbox(Africabbox)
)+
  tm_fill(col = "black")+  
  tm_layout(bg.color = "white")

Africa_allsp <- tm_shape(world, bbox= st_bbox(Africabbox))+
  tm_fill(col = "#93ccdb")+
  # tm_borders(col = "#555555")+
  tm_shape(ter_all_sp_raster$sum_all_sp_tht_cert_ras, bbox= st_bbox(Africabbox))+
  tm_raster(palette =  colfunc(10),
            style = "cont",
            title = "Threat Value")+  
  tm_layout(bg.color = "white")

tmap_save(Africa_allsp, file = "Figures/Final_maps_tmap/Africa_All_sp_value.jpeg", width = 7,height = 7, dpi = 1000)

# South East Asia 
SE_Asiabbox <-  st_bbox(c( xmin = 8500000, xmax = 20000000, ymin = -4000000, ymax = 3000000))

tm_shape(world , bbox= st_bbox(SE_Asiabbox))+
  tm_fill(col = "black")+  
  tm_layout(bg.color = "white")

SE_Asia_allsp <- tm_shape(world, bbox= st_bbox(SE_Asiabbox))+
  tm_fill(col = "#93ccdb")+
  # tm_borders(col = "#555555")+
  tm_shape(ter_all_sp_raster$sum_all_sp_tht_cert_ras, bbox= st_bbox(SE_Asiabbox))+
  tm_raster(palette =  colfunc(25),
            style = "cont",
            title = "Threat Value")+  
  tm_layout(bg.color = "white")

tmap_save(SE_Asia_allsp, file = "Figures/Final_maps_tmap/SE_Asia_All_sp_value.jpeg", width = 8,height = 7, dpi = 1000)

# plot(ter_all_sp_raster$Bird_tht_cert_raster)
# plot(world$geometry, col = "grey97", border = "grey70",axes = T, main = "bird sp terestrial richness")
# plot(ter_all_sp_raster$sum_all_sp_tht_cert_ras, col = alpha(rev(cividis(35)),0.95), add = T)
# 

# Creating hotspots of the certainty maps METHOD 1  taking 1% and 5% of threatened cells: -----
#  Hotspots: 5% quantile and 1% quantile
ter_all_sp_H <- ter_all_sp_raster
# Set 0's to NA to allow the quantile to calculate from cells with value >0 - max
values( ter_all_sp_H ) <- ifelse(values( ter_all_sp_H ) == 0, NA ,values( ter_all_sp_H ))

# for sum of all species 
# set top 1% of cells to value 2, top 5% to value 1
values( ter_all_sp_H$sum_all_sp_tht_cert_ras ) <- ifelse(values( ter_all_sp_H$sum_all_sp_tht_cert_ras ) <=  quantile( ter_all_sp_H$sum_all_sp_tht_cert_ras , prob = 0.95),0, 
                                                         ifelse(values( ter_all_sp_H$sum_all_sp_tht_cert_ras ) >  quantile( ter_all_sp_H$sum_all_sp_tht_cert_ras , prob = 0.99),2,1))

plot( ter_all_sp_H$sum_all_sp_tht_cert_ras )


# For amphibians
# set top 1% of cells to value 2, top 5% to value 1
values( ter_all_sp_H$Amph_tht_cert_raster ) <- ifelse(values( ter_all_sp_H$Amph_tht_cert_raster ) <=  raster::quantile( ter_all_sp_H$Amph_tht_cert_raster , prob = 0.95),0, 
                                                      ifelse(values( ter_all_sp_H$Amph_tht_cert_raster ) >  quantile( ter_all_sp_H$Amph_tht_cert_raster , prob = 0.99),2,1))

plot( ter_all_sp_H$Amph_tht_cert_raster )


# repeat for birds
# set top 1% of cells to value 2, top 5% to value 1
values( ter_all_sp_H$Bird_tht_cert_raster ) <- ifelse(values( ter_all_sp_H$Bird_tht_cert_raster ) <=  quantile( ter_all_sp_H$Bird_tht_cert_raster , prob = 0.95),0, 
                                                      ifelse(values( ter_all_sp_H$Bird_tht_cert_raster ) >  quantile( ter_all_sp_H$Bird_tht_cert_raster , prob = 0.99),2,1))

plot( ter_all_sp_H$Bird_tht_cert_raster )

# repeat for fish
# set top 1% of cells to value 2, top 5% to value 1
values( ter_all_sp_H$fish_tht_cert_raster ) <- ifelse(values( ter_all_sp_H$fish_tht_cert_raster ) <=  quantile( ter_all_sp_H$fish_tht_cert_raster , prob = 0.95),0, 
                                                      ifelse(values( ter_all_sp_H$fish_tht_cert_raster ) >  quantile( ter_all_sp_H$fish_tht_cert_raster , prob = 0.99),2,1))

plot( ter_all_sp_H$fish_tht_cert_raster )

# repeat for mammals
# set top 1% of cells to value 2, top 5% to value 1
values( ter_all_sp_H$Mams_tht_cert_raster ) <- ifelse(values( ter_all_sp_H$Mams_tht_cert_raster ) <=  quantile( ter_all_sp_H$Mams_tht_cert_raster , prob = 0.95),0, 
                                                      ifelse(values( ter_all_sp_H$Mams_tht_cert_raster ) >  quantile( ter_all_sp_H$Mams_tht_cert_raster , prob = 0.99),2,1))

plot( ter_all_sp_H$Mams_tht_cert_raster )

# repeat for reptiles
# set top 1% of cells to value 2, top 5% to value 1
values( ter_all_sp_H$Reps_tht_cert_raster ) <- ifelse(values( ter_all_sp_H$Reps_tht_cert_raster ) <=  quantile( ter_all_sp_H$Reps_tht_cert_raster , prob = 0.95),0, 
                                                      ifelse(values( ter_all_sp_H$Reps_tht_cert_raster ) >  quantile( ter_all_sp_H$Reps_tht_cert_raster , prob = 0.99),2,1))

plot( ter_all_sp_H$Reps_tht_cert_raster )


# Plotting terestrial hotspots METHOD 1  ----

break_label <- c("Threat Present","", "Top 5%","Top 1%")

# All species summed
allsphotspot <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_all_sp_H$sum_all_sp_tht_cert_ras )+
  tm_raster(title = "All species",
            # breaks = c(0, 1, 1.5, 2),
            palette = three_col,
            labels = break_label,
            legend.reverse = TRUE)

# tmap_save(allsphotspot, file = "Figures/Final_maps_tmap/Final_map_All_sphotspot.jpeg", width = 8,height = 5, dpi = 1000)

# Amphibians
amphhotspot <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_all_sp_H$Amph_tht_cert_raster3 )+
  tm_raster(palette = three_col,
            title = "Amphibians",
            labels = break_label,
            legend.reverse = TRUE)

# tmap_save(amphhotspot, file = "Figures/Final_maps_tmap/Final_map_Amph_sphotspot.jpeg", width = 8,height = 5, dpi = 1000)

# Birds 
birdsphotspot <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_all_sp_H$Bird_tht_cert_raster )+
  tm_raster(palette = three_col,
            title = "Birds",
            labels = break_label,
            legend.reverse = TRUE)

tmap_save(birdsphotspot, file = "Figures/Final_maps_tmap/Final_map_Bird_sphotspot.jpeg", width = 8,height = 5, dpi = 1000)

# Fish
fishsphotspot <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_all_sp_H$fish_tht_cert_raster )+
  tm_raster(palette = three_col,
            title = "Fish",
            labels = break_label,
            legend.reverse = TRUE)

tmap_save(fishsphotspot, file = "Figures/Final_maps_tmap/Final_map_Fish_sphotspot.jpeg", width = 8,height = 5, dpi = 1000)

# Reptiles
repsphotspot <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_all_sp_H$Reps_tht_cert_raster )+
  tm_raster(palette = three_col,
            title = "Reptiles",
            labels = break_label,
            legend.reverse = TRUE)

tmap_save(repsphotspot, file = "Figures/Final_maps_tmap/Final_map_Reptile_sphotspot.jpeg", width = 8,height = 5, dpi = 1000)

# Mammals
mamsphotspot <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_all_sp_H$Mams_tht_cert_raster )+
  tm_raster(palette = three_col,
            title = "Mammals",
            labels = break_label,
            legend.reverse = TRUE)

tmap_save(mamsphotspot, file = "Figures/Final_maps_tmap/Final_map_Mammal_sphotspot.jpeg", width = 8,height = 5, dpi = 1000)


