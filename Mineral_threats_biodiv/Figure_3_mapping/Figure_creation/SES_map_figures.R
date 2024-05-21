# Standardised hotspot maps. 

library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(tibble)
library(ggplot2)
library(scales)
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

getwd()
setwd( "X:/edwards_lab1/User/bop21ipl")
sf_use_s2(FALSE)
# loading and collating can be bypassed by reading the All_tht_cert_Stack SEE plotting raster (not working yet)
# plotting colours ----
# three_col <- c( alpha("#6666AA",0.2), "#FFFFFF","#FFFF00","#FF0000")
three_col <- c( alpha("#5566AA",0.5), "#FFFFFF","#FFFF00","#CC2222")
# three_col <- c( "#5566AA77", "#FFFFFF","#FFFF00","#CC2222")
# three_col <- c( "#86C96C", "#FFFFFF","#FFFF00","#CC2222")
# 
# three_col <- c( "#e8be7c", "#FFFFFF","#7596d8","#345462")
# three_col <- c( "#B8DFA9", "#FFFFFF","#8F6CC9","#4F1D4D")
# three_col <- c( "#440154", "#FFFFFF","#5ec962","#fde725")

#Load global spatial Land file ----
target_crs <- st_crs("+proj=moll +x_0=0 +y_0=0 +lat_0=0 +lon_0=0")


world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_make_valid() %>%
  filter(continent != "Antarctica") %>% 
  st_transform(crs = target_crs)


# Load raster of threat certainty rasters ====
# where raster cell values = sum of the proportion of each species' range within the cell 
# amphibian threat certaity raster
amph_raster <- raster("IUCN_data/Species_Ranges/Outputs/Amphibians/Amph_NONtht_cert_raster3.tif")

# # bird threat certaity raster
bird_raster <- raster("IUCN_data/Species_Ranges/Outputs/Birds/Bird_NONtht_cert_raster3.tif")

# # fish threat certainty raster
fish_raster <- raster("IUCN_data/Species_Ranges/Outputs/Fish/fish_NONtht_cert_raster3.tif")

# # Mammal threat certainty raster
mam_raster <- raster("IUCN_data/Species_Ranges/Outputs/Mammals/Mam_NONtht_cert_raster3.tif")
# 
# # Reptiles threat certainty raster
rep_raster <- raster("IUCN_data/Species_Ranges/Outputs/Reptiles/Rep_NONtht_cert_raster3.tif")


# Collate into one raster ====
all_sp_raster <- stack(amph_raster,
                       bird_raster,
                       fish_raster,
                       mam_raster,
                       rep_raster)
# 
# # sum the scores of each raster
all_sp_raster$sum_all_sp_tht_cert_ras <- sum(amph_raster, bird_raster, 
                                             fish_raster,
                                             mam_raster, rep_raster)

# Load raster of threat certainty rasters ====
# where raster cell values = sum of the proportion of each species' range within the cell 
# amphibian threat certaity raster
amph_thrt_raster <- raster("IUCN_data/Species_Ranges/Outputs/Amphibians/Amph_tht_cert_raster3.tif")

# # bird threat certaity raster
bird_thrt_raster <- raster("IUCN_data/Species_Ranges/Outputs/Birds/Bird_tht_cert_raster3.tif")

# # fish threat certainty raster
fish_thrt_raster <- raster("IUCN_data/Species_Ranges/Outputs/Fish/fish_tht_cert_raster3.tif")
# 
# # Mammal threat certainty raster
mam_thrt_raster <- raster("IUCN_data/Species_Ranges/Outputs/Mammals/mam_tht_cert_raster3.tif")
# 
# # Reptiles threat certainty raster
rep_thrt_raster <- raster("IUCN_data/Species_Ranges/Outputs/Reptiles/Rep_tht_cert_raster3.tif")

# tm_shape(amph_thrt_raster)+
#   tm_raster(palette = viridis(30, alpha = 0.6, option = "E"))+
#   tm_shape(world)+
#   tm_borders(col = "white", lwd = 1)+
#   tm_layout(bg.color = "grey90")
# 
# tm_shape(amph_raster)+
#   tm_raster(palette = viridis(30, alpha = 0.6, option = "E"))+
#   tm_shape(world)+
#   tm_borders(col = "white", lwd = 1)+
#   tm_layout(bg.color = "grey90")


# averaged rasters ====
amph_avrg <- amph_thrt_raster/(amph_raster+amph_thrt_raster)
# give value to all cells even if no threat exists
values(amph_avrg) <- if_else(is.na(values(amph_avrg)), 0 , values(amph_avrg))
# reduce to terrestrial
ter_amph_avrg <- crop(amph_avrg, world)
ter_amph_avrg <- mask(ter_amph_avrg, world)

bird_avrg <- bird_thrt_raster/(bird_raster+bird_thrt_raster)
# give value to all cells even if no threat exists
values(bird_avrg) <- if_else(is.na(values(bird_avrg)), 0 , values(bird_avrg))
# reduce to terrestrial
ter_bird_avrg <- crop(bird_avrg, world)
ter_bird_avrg <- mask(ter_bird_avrg, world)

fish_avrg <- fish_thrt_raster/(fish_raster+fish_thrt_raster)
# give value to all cells even if no threat exists
values(fish_avrg) <- if_else(is.na(values(fish_avrg)), 0 , values(fish_avrg))
# reduce to terrestrial
ter_fish_avrg <- crop(fish_avrg, world)
ter_fish_avrg <- mask(ter_fish_avrg, world)

mam_avrg  <- mam_thrt_raster/(mam_raster+mam_thrt_raster)
# give value to all cells even if no threat exists
values(mam_avrg) <- if_else(is.na(values(mam_avrg)), 0 , values(mam_avrg))
# reduce to terrestrial
ter_mam_avrg <- crop(mam_avrg, world)
ter_mam_avrg<- mask(ter_mam_avrg, world)

rep_avrg  <- rep_thrt_raster/(rep_raster+rep_thrt_raster)
values(rep_avrg) <- if_else(is.na(values(rep_avrg)), 0 , values(rep_avrg))
# reduce to terrestrial
ter_rep_avrg<- crop(rep_avrg, world)
ter_rep_avrg <- mask(ter_rep_avrg, world)

# All species SES
ALL_sp_SES <- sum(amph_thrt_raster, bird_thrt_raster, fish_thrt_raster, mam_thrt_raster, rep_thrt_raster)/ 
  sum(amph_thrt_raster, bird_thrt_raster, fish_thrt_raster, mam_thrt_raster, rep_thrt_raster,
      amph_raster, bird_raster,fish_raster, mam_raster, rep_raster)

values(ALL_sp_SES) <- if_else(is.na(values(ALL_sp_SES)), 0 , values(ALL_sp_SES))
# reduce to terrestrial
ter_ALL_sp_SES<- crop(ALL_sp_SES, world)
ter_ALL_sp_SES <- mask(ter_ALL_sp_SES, world)

all_SES_stk <-  stack(amph_avrg, bird_avrg, fish_avrg, mam_avrg, rep_avrg, ALL_sp_SES)
names <- c("amph_avrg", "bird_avrg", "fish_avrg", "mam_avrg", "rep_avrg", "ALL_sp_SES")
names(all_SES_stk)  <- names

#
writeRaster(all_SES_stk, filename = "IUCN_data/Species_Ranges/Outputs/All_SES_Stack3.stk", overwrite = TRUE)

all_SES_stk <- stack("IUCN_data/Species_Ranges/Outputs/All_SES_Stack3.tif")
names <- c("amph_avrg", "bird_avrg", "fish_avrg", "mam_avrg", "rep_avrg", "ALL_sp_SES")
names(all_SES_stk)  <- names

# plot heat maps of SES =====
# All species
# amphibian map 

ter_AllspSES_n0 <- ter_ALL_sp_SES
values(ter_AllspSES_n0) <- ifelse(values(ter_AllspSES_n0) == 0, NA, values(ter_AllspSES_n0))

Allsp_SES <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_AllspSES_n0)+
  tm_raster(palette = viridis(200,direction = -1,  alpha = 0.6, option = "A"),
            style = "cont")+
  tm_layout (frame = FALSE, 
             legend.title.color = "transparent"
  )
Allsp_SES

tmap_save(Allsp_SES, file = "IUCN_data/Figures/SES_viridis/All_species_SES_viridis_cont_ter5.jpeg", width = 8,height = 5, dpi = 1000)

# amphibian map 
ter_amphSES_n0 <- ter_amph_avrg
values(ter_amphSES_n0) <- ifelse(values(ter_amphSES_n0) == 0, NA, values(ter_amphSES_n0))

amph_SES <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_amphSES_n0)+
  tm_raster(palette = viridis(200,direction = -1,  alpha = 0.6, option = "A"),
            style = "cont")+
  tm_layout (frame = FALSE, 
             legend.title.color = "transparent"
  )

amph_SES
tmap_save(amph_SES, file = "IUCN_data/Figures/SES_viridis/Amphibian_SES_viridis_cont_ter5.jpeg", width = 8,height = 5, dpi = 1000)

# Bird map
ter_birdSES_n0 <- ter_bird_avrg
values(ter_birdSES_n0) <- ifelse(values(ter_birdSES_n0) == 0, NA, values(ter_birdSES_n0))


bird_SES <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_birdSES_n0)+
  tm_raster(palette = viridis(200,direction = -1,  alpha = 0.6, option = "A"),
            style = "cont")+
  tm_layout (frame = FALSE, 
             legend.title.color = "transparent"
  )
bird_SES
tmap_save(bird_SES, file = "IUCN_data/Figures/SES_viridis/Bird_SES_viridis_cont_ter5.jpeg", width = 8,height = 5, dpi = 1000)

# fish map
ter_fishSES_n0 <- ter_fish_avrg
values(ter_fishSES_n0) <- ifelse(values(ter_fishSES_n0) == 0, NA, values(ter_fishSES_n0))

fish_SES <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_fishSES_n0)+
  tm_raster(palette = viridis(200,direction = -1,  alpha = 0.6, option = "A"),
            style = "cont")+
  tm_layout (frame = FALSE, 
             legend.title.color = "transparent"
  )
fish_SES
tmap_save(fish_SES, file = "IUCN_data/Figures/SES_viridis/Fish_SES_viridis_cont_ter5.jpeg", width = 8,height = 5, dpi = 1000)

# Mammal map
ter_mamSES_n0 <- ter_mam_avrg
values(ter_mamSES_n0) <- ifelse(values(ter_mamSES_n0) == 0, NA, values(ter_mamSES_n0))

mam_SES <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_mamSES_n0)+
  tm_raster(palette = viridis(200,direction = -1,  alpha = 0.6, option = "A"),
            style = "cont")+
  tm_layout (frame = FALSE, 
             legend.title.color = "transparent")
mam_SES
tmap_save(mam_SES, file = "IUCN_data/Figures/SES_viridis/Mammal_SES_viridis_cont_ter5.jpeg", width = 8,height = 5, dpi = 1000)

# Reptile map 
ter_repSES_n0 <- ter_rep_avrg
values(ter_repSES_n0) <- ifelse(values(ter_repSES_n0) == 0, NA, values(ter_repSES_n0))

rep_SES <-  tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_repSES_n0)+
  tm_raster(palette = viridis(200,direction = -1,  alpha = 0.6, option = "A"),
            style = "cont")+
  tm_layout (frame = FALSE, 
             legend.title.color = "transparent")
rep_SES
tmap_save(rep_SES, file = "IUCN_data/Figures/SES_viridis/Reptile_SES_viridis_cont_ter5.jpeg", width = 8,height = 5, dpi = 1000)


# combing plots 
library(cowplot)

ses_map_grid <- tmap_arrange(Allsp_SES, amph_SES, bird_SES, fish_SES, mam_SES, rep_SES, ncol = 1 )
ses_map_grid
tmap_save(ses_map_grid, file = "X:/edwards_lab1/User/bop21ipl/IUCN_data/Figures/SES_viridis/AllSES_map_grid_viridis_ter5.jpeg", 
          width = 8, height = 20, dpi = 1000)


# Creating hotspots of the certainty maps -----
#  Hotspots: 5% quantile and 1% quantile
ter_all_sp_H <- ter_ALL_sp_SES
# for sum of all species 
# set top 1% of cells to value 2, top 5% to value 1
values( ter_all_sp_H) <- ifelse(values( ter_all_sp_H ) ==0 , NA, 
                                       ifelse(values( ter_all_sp_H) <  quantile( ter_all_sp_H , prob = 0.95),0, 
                                                         ifelse(values( ter_all_sp_H) >=  quantile( ter_all_sp_H, prob = 0.99),2,1)))

plot( ter_all_sp_H )
writeRaster(ter_all_sp_H, filename = "All_SES_hotspots.tiff")


# For amphibians
# set top 1% of cells to value 2, top 5% to value 1
ter_amph_avrg_H <- ter_amph_avrg
values( ter_amph_avrg_H ) <- ifelse(values( ter_amph_avrg_H ) ==0 , NA, 
                                    ifelse(values( ter_amph_avrg_H ) <  quantile( ter_amph_avrg_H , prob = 0.95),0, 
                                                      ifelse(values( ter_amph_avrg_H ) >=  quantile( ter_amph_avrg_H , prob = 0.99),2,1)))

plot( ter_amph_avrg_H )

# repeat for birds
# set top 1% of cells to value 2, top 5% to value 1
ter_bird_avrg_H <- ter_bird_avrg
values( ter_bird_avrg_H ) <- ifelse(values( ter_bird_avrg_H ) ==0 , NA, 
                                    ifelse(values( ter_bird_avrg_H ) <  quantile( ter_bird_avrg_H , prob = 0.95),0, 
                                                      ifelse(values( ter_bird_avrg_H ) >=  quantile( ter_bird_avrg_H , prob = 0.99),2,1)))

plot( ter_bird_avrg_H )

# repeat for fish
# set top 1% of cells to value 2, top 5% to value 1
ter_fish_avrg_H <- ter_fish_avrg
values(ter_fish_avrg_H) <- ifelse(values( ter_fish_avrg_H ) ==0 , NA, 
                                  ifelse(values( ter_fish_avrg_H ) <  quantile( ter_fish_avrg_H , prob = 0.95),0,
                                         # 0.99 quantile == 1 so need to change code
                                    ifelse(values( ter_fish_avrg_H ) >= quantile( ter_fish_avrg_H , prob = 0.99), 2,1)))

quantile( ter_fish_avrg , prob = c(0.95,0.99))

hist(values(ter_fish_avrg),)
plot(ter_fish_avrg_H)
plot( ter_fish_avrg )

# repeat for mammals
# set top 1% of cells to value 2, top 5% to value 1
ter_mam_avrg_H <- ter_mam_avrg
values( ter_mam_avrg_H ) <- ifelse(values( ter_mam_avrg_H ) ==0 , NA, 
                                   ifelse(values( ter_mam_avrg_H ) <  quantile( ter_mam_avrg_H , prob = 0.95),0, 
                                    ifelse(values( ter_mam_avrg_H ) >=  quantile( ter_mam_avrg_H , prob = 0.99),2,1)))


plot(ter_mam_avrg_H)
plot( ter_mam_avrg )

# repeat for reptiles
# set top 1% of cells to value 2, top 5% to value 1
ter_rep_avrg_H <- ter_rep_avrg
values( ter_rep_avrg_H ) <- ifelse(values( ter_rep_avrg_H ) ==0 , NA, 
                                   ifelse(values( ter_rep_avrg_H ) <  quantile( ter_rep_avrg_H , prob = 0.95),0, 
                                    ifelse(values( ter_rep_avrg_H ) >=  quantile( ter_rep_avrg_H , prob = 0.99),2,1)))

plot(ter_rep_avrg_H)

# Plotting terestrial hotspots ----

break_label <- c("Threat Present","", "Top 5%","Top 1%")

# All species summed
allsphotspot <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_all_sp_H)+
  tm_raster(palette = three_col,
            labels = break_label,
            legend.reverse = TRUE)+
  tm_layout(frame = FALSE,
            legend.show = F)
allsphotspot
tmap_save(allsphotspot, file = "IUCN_data/Figures/SES_hotspots/Final_SES_All_sp_hotspot_V4.jpeg", width = 8,height = 5, dpi = 1000)

# Amphibians
amphhotspot <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_amph_avrg_H)+
  tm_raster(palette = three_col,
            labels = break_label,
            legend.reverse = TRUE)+
  tm_layout(frame = FALSE,
            legend.show = F)

amphhotspot
tmap_save(amphhotspot, file = "IUCN_data/Figures/SES_hotspots/Final_SES_Amph_sphotspot_V4.jpeg", width = 8,height = 5, dpi = 1000)

# Birds 
birdsphotspot <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_bird_avrg_H)+
  tm_raster(palette = three_col,
            labels = break_label,
            legend.reverse = TRUE)+
  tm_layout(frame = FALSE,
            legend.show = F)
birdsphotspot
tmap_save(birdsphotspot, file = "IUCN_data/Figures/SES_hotspots/Final_SES_Bird_sphotspot_V4.jpeg", width = 8, height = 5, dpi = 1000)

# Fish
fishsphotspot <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_fish_avrg_H)+
  tm_raster(palette = three_col,
            labels = break_label,
            legend.reverse = TRUE)+
  tm_layout(frame = FALSE,
            legend.show = F)
fishsphotspot
tmap_save(fishsphotspot, file = "IUCN_data/Figures/SES_hotspots/Final_SES_Fish_sphotspot_V4.jpeg", width = 8,height = 5, dpi = 1000)

#mammals
mamsphotspot <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_mam_avrg_H)+
  tm_raster(palette = three_col,
            title = "Mammals",
            labels = break_label,
            legend.reverse = TRUE)+
  tm_layout(frame = FALSE,
            legend.show = F)
mamsphotspot
tmap_save(mamsphotspot, file = "IUCN_data/Figures/SES_hotspots/Final_SES_Mammal_sphotspot_V4.jpeg", width = 8,height = 5, dpi = 1000)

# Reptiles
repsphotspot <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_rep_avrg_H)+
  tm_raster(palette = three_col,
            labels = break_label,
            legend.reverse = TRUE)+
  tm_layout(frame = FALSE,
            legend.show = F)
repsphotspot
tmap_save(repsphotspot, file = "IUCN_data/Figures/SES_hotspots/Final_SES_Reptile_sphotspot_V4.jpeg", width = 8,height = 5, dpi = 1000)

# plot only the legend horizontally
legend_only <- tm_shape(world)+
  tm_fill(col = "grey80")+
  tm_shape(ter_rep_avrg_H )+
  tm_raster(palette = three_col,
            labels = break_label,
            legend.reverse = TRUE,
            legend.is.portrait = FALSE)+
  tm_layout(frame = FALSE,
            legend.title.color = "transparent",
            legend.stack = "horizontal",
            legend.only = TRUE,
  )

legend_only

tmap_save(legend_only, file = "IUCN_data/Figures/SES_hotspots/Legend_only_SES_Horz_V2.jpeg", width = 3, height = 1, dpi = 1000)

