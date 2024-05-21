# Additional mining maps comparison 

library(tidyverse)
library(terra)
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
library(data.table)
library(patchwork)
library(rmarkdown)

getwd()
setwd("X:/edwards_lab1/User/bop21ipl")
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





# load cell areas and create raster for mining and quarrying threat and seepage threat =====

# list of species only threatened by mining and quarrying and mine seepage 
threatened_sp <- read_csv("IUCN_data/Species_Pages/Raw_Data/CHORDATA_pg_Oil_Mining_threat/threats.csv")

threatened_sp <- threatened_sp %>% 
  filter(code %in% c("3.2","9.2.2")) %>% 
  pull(internalTaxonId)

threatened_sp %>% unique() %>% length() # 4116 species threatened by Mining and quarrying or seepage from mining


# Amphibians
amph_cell_Areas <- read_csv("IUCN_data/Species_Ranges/Data/amphs_minetht_cell_areas2.csv")
length(unique(amph_cell_Areas$species)) # 747
nrow(amph_cell_Areas)
amph_cell_Areas <- amph_cell_Areas %>% 
  distinct()
nrow(amph_cell_Areas)

#Birds
bird_cell_Areas <- read_csv("IUCN_data/Species_Ranges/Data/Bird_minetht_cell_areas2.csv")
length(unique(bird_cell_Areas$species)) # 556
nrow(bird_cell_Areas)
bird_cell_Areas[(duplicated(bird_cell_Areas)),]
bird_cell_Areas <- bird_cell_Areas %>% 
  distinct()
nrow(bird_cell_Areas)

# Mammals
mam_cell_Areas <- read_csv("IUCN_data/Species_Ranges/Data/Mams_minetht_cell_areas2.csv")
length(unique(mam_cell_Areas$species)) # 516
nrow(mam_cell_Areas)
mam_cell_Areas[(duplicated(mam_cell_Areas)),]
mam_cell_Areas <- mam_cell_Areas %>% 
  distinct()
nrow(mam_cell_Areas)


# Reptiles
rep_cell_Areas <- read_csv("IUCN_data/Species_Ranges/Data/reps_minetht_cell_areas2.csv")
length(unique(rep_cell_Areas$species)) # 737
nrow(rep_cell_Areas)
rep_cell_Areas[(duplicated(rep_cell_Areas)),]
rep_cell_Areas <- rep_cell_Areas %>% 
  distinct()
nrow(rep_cell_Areas)


# Fish
fish_cell_Areas <- read_csv("IUCN_data/Species_Ranges/Data/fish_minetht_cell_areas2.csv")
length(unique(fish_cell_Areas$species)) # 2015
nrow(fish_cell_Areas) 
fish_cell_Areas[(duplicated(fish_cell_Areas)),]
fish_cell_Areas <- fish_cell_Areas %>% 
  distinct()
nrow(fish_cell_Areas)


# bind all cell areas 
all_cell_areas <- bind_rows(amph_cell_Areas,
                            bird_cell_Areas,
                            mam_cell_Areas,
                            rep_cell_Areas,
                            fish_cell_Areas)


# reduce to only species with mining threat
all_cell_areas <- all_cell_areas %>% 
  filter(id_no %in% threatened_sp)

# sum any multiple distinct ranges within a single cell
sp_cell_Values <- all_cell_areas %>%
  group_by(species,cell)%>%
  summarise(Areakm2 = sum(Areakm2))%>%
  # Add total species range column 
  group_by(species) %>%
  mutate(total_Areakm2 = sum(Areakm2))%>%
  # calculate proportion of sp. range within each grid cell 
  ungroup%>%
  mutate(Rng_prop_cell = Areakm2/total_Areakm2)


# Create dataframe of sum of Rng_proportions_cell
cell_threat_val <- sp_cell_Values %>%
  group_by(cell) %>%
  summarise(threat_val = sum(Rng_prop_cell)) %>%
  full_join(world_grid_centroids, by  = "cell")%>%
  mutate(threat_val = if_else(is.na(threat_val), 0, threat_val))%>%
  arrange(-x)%>%
  dplyr::select(cell,x,y, threat_val) %>% 
  arrange(cell)


# create grid shapefile with threat values 
grid_ras <- left_join(grid, cell_threat_val, by = "cell") %>% 
  dplyr::select(threat_val)

library("RColorBrewer")
# visualise raster
tm_shape(world)+
  tm_polygons()+
  tm_shape(grid_ras)+
  tm_fill(col = "threat_val", alpha = 0.7, palette = "YlGnBu")

st_write(grid_ras, "Chapter_One2/Response_to_reviewers_CuBiol/All_tax_MQ_SEEP_raster.gpkg", append = FALSE)

#rasterize grid
thrt_raster_cert <- grid_ras %>% 
  st_rasterize()

tm_shape(world)+
  tm_polygons()+
  tm_shape(thrt_raster_cert)+
  tm_raster(col = "threat_val", alpha = 0.7, palette = "YlGnBu")

tm_shape(world)+
  tm_fill(col = "#93ccdb")+
  tm_shape(bird_grid_ras)+
  tm_fill(col = "threat_val", palette = viridis(30, alpha = 0.6, option = "E"),
          style = "kmeans"
  )+
  tm_layout (frame = FALSE, bg.color = "transparent", main.title = "Birds threat value" , legend.title.color = "transparent")


tm_shape(world)+
  tm_polygons(col = "grey80")+
  tm_shape(thrt_raster_cert)+
  tm_raster(#col = "bird_impacted_ras",
            # col.scale =  tm_scale_continuous(values = viridis(200,direction = 1,  alpha = 1, begin = 0.2, option = "D"), value.na = "#000000"),
            col.scale =  tm_scale_continuous(values = viridis(200,direction = -1,  alpha = 1, begin = 0, option = "A"), value.na = "#00000000"),
            col.legend = tm_legend(title = "mining threat",
                                   reverse = T,
                                   
            ))


# write the raster of threats to all  species with certainty values of the proportion each species ranges within raster pixels
write_stars(thrt_raster_cert, dsn = "Chapter_One2/Response_to_reviewers_CuBiol/All_tax_MQ_SEEP_raster.tif", overwrite = TRUE)

# All species rasters of global diversity weighted by range size ==== 
# Load raster of threat certainty rasters ===
# where raster cell values = sum of the proportion of each species' range within the cell 
# amphibian  certaity raster
amph_raster <- raster("IUCN_data/Species_Ranges/Outputs/Amphibians/Amph_NONtht_cert_raster3.tif")

# # bird  certaity raster
bird_raster <- raster("IUCN_data/Species_Ranges/Outputs/Birds/Bird_NONtht_cert_raster3.tif")

# # fish  certainty raster
fish_raster <- raster("IUCN_data/Species_Ranges/Outputs/Fish/fish_NONtht_cert_raster3.tif")

# # Mammal  certainty raster
mam_raster <- raster("IUCN_data/Species_Ranges/Outputs/Mammals/Mam_NONtht_cert_raster3.tif")
# 
# # Reptiles certainty raster
rep_raster <- raster("IUCN_data/Species_Ranges/Outputs/Reptiles/Rep_NONtht_cert_raster3.tif")



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

# All species SES
ALL_sp_sum <- sum(amph_thrt_raster, bird_thrt_raster, fish_thrt_raster, mam_thrt_raster, rep_thrt_raster,
      amph_raster, bird_raster,fish_raster, mam_raster, rep_raster)

# change NA values to 0 
terra::values(ALL_sp_sum) <- if_else(is.na(values(ALL_sp_sum)), 0 , terra::values(ALL_sp_sum))

terra::writeRaster(ALL_sp_sum, overwrite = T, "IUCN_data/Species_Ranges/Outputs/All_species_diversity_range_weighted.tiff")

ALL_sp_sum %>% plot
## combine mining raster and threat raster =====
# load mining raster from Maus 2022
mining_rast2 <- terra::rast("data/Maus_mine_data/global_miningarea_v2_30arcminute.tif")
# mining_rast2 <- terra::rast("data/Maus_mine_data/global_miningarea_v2_30arcsecond.tif")

mining_rast

bird_imp_raster_cert <- terra::rast( "Chapter_One2/Response_to_reviewers_CuBiol/bird_impacted_ras.tif")
mining_rast2 <- project(mining_rast2, bird_imp_raster_cert)

mining_rast2 %>% plot

# load threat data
thrt_raster_cert <- terra::rast( "Chapter_One2/Response_to_reviewers_CuBiol/All_tax_MQ_SEEP_raster.tif")
# load total diversity raster 
all_sp_div <- terra::rast("IUCN_data/Species_Ranges/Outputs/All_species_diversity_range_weighted.tiff")

par(mfrow = c(1,1))
plot1 <- log10(all_sp_div$All_species_diversity_range_weighted+1) %>% plot

#transform raster to the same projection
thrt_raster_cert <- project(thrt_raster_cert, mining_rast2)
all_sp_div <- project(all_sp_div, mining_rast2)

#crop threat map and diversity map to terrestrial areas
ter_thrt_raster_cert <- terra::crop(thrt_raster_cert, mining_rast2)
ter_thrt_raster_cert <- terra::mask(ter_thrt_raster_cert, mining_rast2)

# reduce to terrestrial
ter_all_sp_div<- terra::crop(all_sp_div, mining_rast2)
ter_all_sp_div<- terra::mask(ter_all_sp_div, mining_rast2)

ter_all_sp_div %>% plot
mining_rast2 %>% plot
ter_thrt_raster_cert %>%  plot

# calculate standardised threat map
ter_SES_threat <- ter_thrt_raster_cert / ter_all_sp_div

par(mfrow = c(3,1))
plot(ter_thrt_raster_cert)
plot(ter_SES_threat)
plot(mining_rast2)


# check extents match 
ext(ter_thrt_raster_cert) == ext(mining_rast2)
ext(ter_all_sp_div) == ext(mining_rast2)

ext(ter_SES_threat) == ext(mining_rast2)

# write rasters with matching extents to use 
terra::writeRaster(ter_thrt_raster_cert, overwrite = T, "Chapter_One2/Response_to_reviewers_CuBiol/transformed_TER_All_tax_MQ_SEEP_raster.tif")
terra::writeRaster(ter_SES_threat, overwrite = T,  "Chapter_One2/Response_to_reviewers_CuBiol/transformed_TER_All_tax_SES_MQ_SEEP_raster.tif")
terra::writeRaster(mining_rast2, overwrite = T,  "Chapter_One2/Response_to_reviewers_CuBiol/transformed_Maus_global_miining_areaV2_30arcminute.tif")
# terra::writeRaster(mining_rast2, overwrite = T,  "Chapter_One2/Response_to_reviewers_CuBiol/transformed_Maus_global_miining_areaV2_5arcminute.tif")


##### load in using raster package to use bivariate plot funciton ====
ter_thrt_raster <- raster::raster( "Chapter_One2/Response_to_reviewers_CuBiol/transformed_TER_All_tax_MQ_SEEP_raster.tif")
ter_SES_threat <-  raster::raster( "Chapter_One2/Response_to_reviewers_CuBiol/transformed_TER_All_tax_SES_MQ_SEEP_raster.tif")
# convert to RasterLayer
ter_SES_threat <- ter_SES_threat$All_tax_MQ_SEEP_raster

mining_rast <- raster::raster("Chapter_One2/Response_to_reviewers_CuBiol/transformed_Maus_global_miining_areaV2_30arcminute.tif")
# mining_rast <- raster::raster("Chapter_One2/Response_to_reviewers_CuBiol/transformed_Maus_global_miining_areaV2_5arcminute.tif")

par(mfrow = c(3,1))
ter_thrt_raster %>% plot
ter_SES_threat %>% plot
mining_rast %>%  plot

crs(ter_thrt_raster)
crs(mining_rast)
ext(ter_thrt_raster) == ext(mining_rast)
ext(ter_SES_threat) == ext(mining_rast)

ter_thrt_raster_log <- log10(ter_thrt_raster +1)
mining_rast_log <- log10(mining_rast+1)

# basic bivariate plot =====
nBreaks = 5

# Distribution plots ====
mining_val <- tibble(val = values(mining_rast_log) %>% na.omit())
threat_val <- tibble(val = values(ter_thrt_raster_log) %>% na.omit())
SES_val <- tibble(val = values(ter_SES_threat) %>% na.omit())

ggplot(mining_val,aes(x = val)) +
  geom_density()+
ggplot(threat_val,aes(x = val)) +
  geom_density()+
ggplot(SES_val,aes(x = val)) +
  geom_density()

library("spData")

data(jenks71, package="spData")
pal1 <- c("wheat1", "red3")
opar <- par(mfrow=c(2,3))

range(mining_val$val)
n <- 5

par(mfrow= c(3,1))
pal1 <- c("#dfdfd6", "#0096FB")
pal2 <- c("#dfdfd6", "#d76700")
mine_f <-plot(classIntervals(mining_val$val, n=n, style="fisher"), pal=pal1, main="mine area")
mining_val$val %>% range
mine_f <-plot(classIntervals(mining_val$val, n=n, style="fixed",
                             fixedBreaks = c(0, 0.2, 0.8, 1.5, 2.532184)), pal=pal1, main="mine area")

threat_f <-plot(classIntervals(threat_val$val, n=n, style="fisher"), pal=pal2, main="threat score")
SES_f <-plot(classIntervals(SES_val$val, n=n, style="fisher"), pal=pal2, main="SES threat")

threat_f <-plot(classIntervals(threat_val$val, n=n, style="quantile"), pal=pal2, main="threat score")
SES_f <-plot(classIntervals(SES_val$val, n=n, style="quantile"), pal=pal2, main="SES threat")


# create the bivariate  for  threat values raster FISHER BREAKS =====
breaks <- c(0, 0.2, 0.8, 1.5, 2.532184)
col.matrixF <- colmat(nbreaks = nBreaks, breakstyle = "fisher", 
                      xlab = "Mine area", ylab = "Threat score", 
                      # upperleft = "#f0f022", upperright = "#111111", bottomleft = "#E8E8E8", bottomright = "#7787D4",
                      # bottomright = "#0096EB", upperright = "#820050", bottomleft= "#BEBEBE", upperleft = "#FFE60F",
                      bottomright = "#0096FB", upperright = "#420050", bottomleft= "#dfdfd0", upperleft = "#d76700",
                      saveLeg = FALSE, plotLeg = TRUE )


biv_plot <- bivariate.map(rasterx = mining_rast_log, 
                          rastery = ter_thrt_raster_log,
                          export.colour.matrix = FALSE,
                          colourmatrix = col.matrixF
)

# Convert to dataframe for plotting with ggplot
bivMapDFF <- setDT(as.data.frame(biv_plot, xy = TRUE))
colnames(bivMapDFF)[3] <- "BivValue"
bivMapDFF <- melt(bivMapDFF, id.vars = c("x", "y"),
                  measure.vars = "BivValue",
                  value.name = "bivVal",
                  variable.name = "Variable")

ext(biv_plot)

# Make the map using ggplot
map_f <- ggplot(bivMapDFF, aes(x = x, y = y)) +
  geom_raster(aes(fill = bivVal)) +
  # scale_y_continuous(breaks = seq(-20, 60, by = 10),
  # labels = paste0(seq(-20, 60, 10), "°")) +
  scale_fill_gradientn(colours = col.matrixF, na.value = "transparent") +
  theme_bw() +
  theme(text = element_text(size = 10, colour = "black")) +
  borders(colour = "black", size = 0.5) +
  coord_sf(expand = FALSE, xlim =c(xmin(biv_plot),xmax(biv_plot)), ylim = c(ymin(biv_plot),ymax(biv_plot))) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        strip.text = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black")) +
  # labs(x = "Longitude", y = "Latitude") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        # panel.background = element_rect(fill = "#F0F0F0"),
  )
map_f

ggsave(plot = map_f , filename =  "Chapter_One2/Response_to_reviewers_CuBiol/mining_threat_val_bivariateFISHER.jpg", width = 8, height = 5, dpi = 300)

# change legend text size to 9 
fig2 <- {map_f + ggtitle("A \nSpecies Threat")}  + 
  inset_element(BivLegend + theme(plot.background = element_rect(fill = "white",
                                                                 colour = NA,
                                                                 ),
                                  axis.title = element_text(size = 9)), 
                left = 0, bottom = -0.1, right = 0.2, top = 0.2,
                align_to = "full") 

fig2
ggsave(plot = fig2 , filename =  "Chapter_One2/Response_to_reviewers_CuBiol/mining_threat_val_bivariate_FISHER_wlegend.jpg", 
       units = "in", 
       width = 8, height = 5, dpi = 300)


# create the bivariate  for SES threat (Community sensitivity score) values raster FISHER BREAKS =====

col.matrixF2 <- colmat(nbreaks = nBreaks, breakstyle = "fisher",
                      xlab = "Mine area", ylab = "CS score", 
                      # upperleft = "#f0f022", upperright = "#111111", bottomleft = "#E8E8E8", bottomright = "#7787D4",
                      # bottomright = "#0096EB", upperright = "#820050", bottomleft= "#BEBEBE", upperleft = "#FFE60F",
                      bottomright = "#0096FB", upperright = "#420050", bottomleft= "#dfdfd0", upperleft = "#d76700",
                      saveLeg = FALSE, plotLeg = TRUE)


biv_plot2 <- bivariate.map(rasterx = mining_rast_log, 
                          rastery = ter_SES_threat,
                          export.colour.matrix = FALSE,
                          colourmatrix = col.matrixF2
)

# Convert to dataframe for plotting with ggplot
bivMapDFF2 <- setDT(as.data.frame(biv_plot2, xy = TRUE))
colnames(bivMapDFF2)[3] <- "BivValue"
bivMapDFF2 <- melt(bivMapDFF2, id.vars = c("x", "y"),
                  measure.vars = "BivValue",
                  value.name = "bivVal",
                  variable.name = "Variable")

ext(biv_plot2)

# Make the map using ggplot
map_SESf <- ggplot(bivMapDFF2, aes(x = x, y = y)) +
  geom_raster(aes(fill = bivVal)) +
  # scale_y_continuous(breaks = seq(-20, 60, by = 10),
  # labels = paste0(seq(-20, 60, 10), "°")) +
  scale_fill_gradientn(colours = col.matrixF, na.value = "transparent") +
  theme_bw() +
  theme(text = element_text(size = 10, colour = "black")) +
  borders(colour = "black", size = 0.5) +
  coord_sf(expand = FALSE, xlim =c(xmin(biv_plot),xmax(biv_plot)), ylim = c(ymin(biv_plot),ymax(biv_plot))) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        strip.text = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black")) +
  # labs(x = "Longitude", y = "Latitude") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        # panel.background = element_rect(fill = "#F0F0F0"),
  )
map_SESf

# ggsave(plot = map_SESf , filename =  "Chapter_One2/Response_to_reviewers_CuBiol/mining_SESthreat_val_bivariateFISHER.jpg", width = 8, height = 5, dpi = 300)

# change legend text size to 9 
fig3 <- {map_SESf + ggtitle("B \nCommunity Sensitivity")}  + 
  inset_element(BivLegend + theme(plot.background = element_rect(fill = "white",
                                                                 colour = NA,
  ),
  axis.title = element_text(size = 9)), 
  left = 0, bottom = -0.1, right = 0.2, top = 0.2,
  align_to = "full") 

fig3

# ggsave(plot = fig3 , filename =  "Chapter_One2/Response_to_reviewers_CuBiol/mining_SESthreat_val_bivariate_FISHER_wlegend.jpg", 
#        units = "in", 
#        width = 8, height = 5, dpi = 300)



# join both plots together ====
fig4 <- {fig2 / fig3 }  
fig4
ggsave(plot = fig4 , filename =  "Chapter_One2/Response_to_reviewers_CuBiol/Figure_S1_both_threat_mining_bivplots_FISHER.jpg",
       units = "in",
       width = 8, height = 10, dpi = 300)

#
# create the bivariate threat raster QUANTILE BREAK ====

nBreaks <- 10
col.matrixQ <- colmat(nbreaks = nBreaks, breakstyle = "quantile",
                      xlab = "Mine density", ylab = "Threat score", 
                      # upperleft = "#f0f022", upperright = "#111111", bottomleft = "#E8E8E8", bottomright = "#7787D4",
                      # bottomright = "#0096EB", upperright = "#820050", bottomleft= "#BEBEBE", upperleft = "#FFE60F",
                      bottomright = "#0096FB", upperright = "#420050", bottomleft= "#dfdfd6", upperleft = "#d76700",
                      saveLeg = FALSE, plotLeg = TRUE)

biv_plot <- bivariate.map(rasterx = mining_rast_log, 
                          rastery = ter_thrt_raster_log,
                          export.colour.matrix = FALSE,
                          colourmatrix = col.matrixQ
)

# Convert to dataframe for plotting with ggplot
bivMapDFQ <- setDT(as.data.frame(biv_plot, xy = TRUE))
colnames(bivMapDFQ)[3] <- "BivValue"
bivMapDFQ <- melt(bivMapDFQ, id.vars = c("x", "y"),
                  measure.vars = "BivValue",
                  value.name = "bivVal",
                  variable.name = "Variable")

ext(biv_plot)

# Make the map using ggplot
map_q <- ggplot(bivMapDFQ, aes(x = x, y = y)) +
  geom_raster(aes(fill = bivVal)) +
  # scale_y_continuous(breaks = seq(-20, 60, by = 10), 
                     # labels = paste0(seq(-20, 60, 10), "°")) +
  # scale_x_continuous(breaks = seq(50,175,25), 
                     # labels = paste0(seq(50,175,25), "°")) +
  scale_fill_gradientn(colours = col.matrixQ, na.value = "transparent") + 
  theme_bw() +
  theme(text = element_text(size = 10, colour = "black")) +
  borders(colour = "black", size = 0.5) +
  coord_sf(expand = FALSE, xlim =c(xmin(biv_plot),xmax(biv_plot)), ylim = c(ymin(biv_plot),ymax(biv_plot))) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        strip.text = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black")) +
  # labs(x = "Longitude", y = "Latitude") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        # panel.background = element_rect(fill = "#F0F0F0"),
        )
map_q

# figrue
# change legend text size to 9 
fig5 <- {map_q + ggtitle("A \nSpecies Threat")}  + 
  inset_element(BivLegend + theme(plot.background = element_rect(fill = "white",
                                                                 colour = NA,
  ),
  axis.title = element_text(size = 9)), 
  left = 0, bottom = -0.1, right = 0.2, top = 0.2,
  align_to = "full") 

fig5
ggsave(plot = fig5 , filename =  "Chapter_One2/Response_to_reviewers_CuBiol/mining_threat_val_bivariate_QUANTILE_wlegend.jpg", 
       units = "in", 
       width = 8, height = 5, dpi = 300)



# # create the bivariate SES raster QUANTILE BREAK ====

col.matrixSESQ <- colmat(nbreaks = nBreaks, breakstyle = "quantile",
                      xlab = "Mine density", ylab = "CS score", 
                      # upperleft = "#f0f022", upperright = "#111111", bottomleft = "#E8E8E8", bottomright = "#7787D4",
                      # bottomright = "#0096EB", upperright = "#820050", bottomleft= "#BEBEBE", upperleft = "#FFE60F",
                      bottomright = "#0096FB", upperright = "#420050", bottomleft= "#dfdfd6", upperleft = "#d76700",
                      saveLeg = FALSE, plotLeg = TRUE)

biv_plot2 <- bivariate.map(rasterx = mining_rast_log, 
                          rastery = ter_SES_threat,
                          export.colour.matrix = FALSE,
                          colourmatrix = col.matrixSESQ
)

# Convert to dataframe for plotting with ggplot
bivMapDFSESQ <- setDT(as.data.frame(biv_plot2, xy = TRUE))
colnames(bivMapDFSESQ)[3] <- "BivValue"
bivMapDFSESQ <- melt(bivMapDFSESQ, id.vars = c("x", "y"),
                  measure.vars = "BivValue",
                  value.name = "bivVal",
                  variable.name = "Variable")

ext(biv_plot2)

# Make the map using ggplot
map_SESq <- ggplot(bivMapDFSESQ, aes(x = x, y = y)) +
  geom_raster(aes(fill = bivVal)) +
  # scale_y_continuous(breaks = seq(-20, 60, by = 10), 
  # labels = paste0(seq(-20, 60, 10), "°")) +
  # scale_x_continuous(breaks = seq(50,175,25), 
  # labels = paste0(seq(50,175,25), "°")) +
  scale_fill_gradientn(colours = col.matrixQ, na.value = "transparent") + 
  theme_bw() +
  theme(text = element_text(size = 10, colour = "black")) +
  borders(colour = "black", size = 0.5) +
  coord_sf(expand = FALSE, xlim =c(xmin(biv_plot),xmax(biv_plot)), ylim = c(ymin(biv_plot),ymax(biv_plot))) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        strip.text = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black")) +
  # labs(x = "Longitude", y = "Latitude") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        # panel.background = element_rect(fill = "#F0F0F0"),
  )
map_SESq

fig6 <- {map_SESq + ggtitle("B \nCommunity Sensitivity")}  + 
  inset_element(BivLegend + theme(plot.background = element_rect(fill = "white",
                                                                 colour = NA,
  ),
  axis.title = element_text(size = 9)), 
  left = 0, bottom = -0.1, right = 0.2, top = 0.2,
  align_to = "full") 

fig6
ggsave(plot = fig6 , filename =  "Chapter_One2/Response_to_reviewers_CuBiol/mining_SES_val_bivariate_QUANTILE_wlegend.jpg", 
       units = "in", 
       width = 8, height = 5, dpi = 300)

# join both figures together =====
fig7 <- {fig5 / fig6 }  
fig7
ggsave(plot = fig7 , filename =  "Chapter_One2/Response_to_reviewers_CuBiol/Figure_S1_both_threat_mining_bivplots_QUANTILE.jpg",
       units = "in",
       width = 8, height = 10, dpi = 300)

#
# create the bivariate raster COMPLETE CLUSTER BREAKS =====
classInt::classIntervals()
col.matrixP <- colmat(nbreaks = nBreaks, breakstyle = "hclust",
                      xlab = "Mine density log10", ylab = "Threat score log 10",
                      # upperleft = "#f0f022", upperright = "#111111", bottomleft = "#E8E8E8", bottomright = "#7787D4",
                      # bottomright = "#0096EB", upperright = "#820050", bottomleft= "#BEBEBE", upperleft = "#FFE60F",
                      bottomright = "#0096FB", upperright = "#420050", bottomleft= "#dfdfd0", upperleft = "#d76700",
                      saveLeg = FALSE, plotLeg = TRUE )


biv_plot <- bivariate.map(rasterx = mining_rast_log, 
                          rastery = ter_thrt_raster_log,
                          export.colour.matrix = FALSE,
                          colourmatrix = col.matrixP
)


# Convert to dataframe for plotting with ggplot
bivMapDFP <- setDT(as.data.frame(biv_plot, xy = TRUE))
colnames(bivMapDFP)[3] <- "BivValue"
bivMapDFP <- melt(bivMapDFP, id.vars = c("x", "y"),
                  measure.vars = "BivValue",
                  value.name = "bivVal",
                  variable.name = "Variable")

ext(biv_plot)

# Make the map using ggplot
map_p <- ggplot(bivMapDFP, aes(x = x, y = y)) +
  geom_raster(aes(fill = bivVal)) +
  # scale_y_continuous(breaks = seq(-20, 60, by = 10), 
  # labels = paste0(seq(-20, 60, 10), "°")) +
  # scale_x_continuous(breaks = seq(50,175,25), 
  # labels = paste0(seq(50,175,25), "°")) +
  scale_fill_gradientn(colours = col.matrixQ, na.value = "transparent") + 
  theme_bw() +
  theme(text = element_text(size = 10, colour = "black")) +
  borders(colour = "black", size = 0.5) +
  coord_sf(expand = FALSE, xlim =c(xmin(biv_plot),xmax(biv_plot)), ylim = c(ymin(biv_plot),ymax(biv_plot))) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        strip.text = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black")) +
  # labs(x = "Longitude", y = "Latitude") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        # panel.background = element_rect(fill = "#F0F0F0"),
  )
map_p

# Compare breaks ====
library(patchwork)

fig <- {{map_q + ggtitle("Quantile breaks")} / {map_f + ggtitle("Fisher breaks")} / {map_p + ggtitle("kmeans breaks")} } + 
  inset_element(BivLegend + theme(plot.background = element_rect(fill = "white",
                                                                 colour = NA)), 
                left = 0.5, bottom = 0.5, right = 1.2, top = 0.90,
                align_to = "full") +
  plot_annotation(caption = "Both plots are made on the same data, but have breaks defined differently")
fig

fig <- { {map_f + ggtitle("Fisher breaks")} / {map_p + ggtitle("kmeans breaks")} } + 
  inset_element(BivLegend + theme(plot.background = element_rect(fill = "white",
                                                                 colour = NA)), 
                left = 0.5, bottom = 0.5, right = 1.2, top = 0.90,
                align_to = "full") +
  plot_annotation(caption = "Both plots are made on the same data, but have breaks defined differently")
fig


fig2 <- {map_q}  + 
  inset_element(BivLegend + theme(plot.background = element_rect(fill = "white",
                                                                 colour = NA)), 
                left = -0.8, bottom = 0.7, right = 1.1, top = 0.98,
                align_to = "left",) 
fig2 <- {map_f + ggtitle("Fisher breaks")}  + 
  inset_element(BivLegend + theme(plot.background = element_rect(fill = "white",
                                                                 colour = NA)), 
                left = 0, bottom = -0.1, right = 0.2, top = 0.2,
                align_to = "full") 
fig2
