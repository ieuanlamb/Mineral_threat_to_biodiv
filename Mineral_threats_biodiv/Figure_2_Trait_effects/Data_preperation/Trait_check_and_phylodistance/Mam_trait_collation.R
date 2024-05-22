# Mammal trait collation 
##########################################
# check coverage for etard traits
# join Etard traits and range data in single dataset
# check correlation 
# match names and add phylo tree
# check phylogenetic signal
# format for imputation 
##########################################

library(tidyverse)
library(ape)
library(treeplyr)
library(Rphylopars)
library(phangorn)
library(picante)

getwd()
set.seed(1121995)

setwd("Ecological trait data")
# load Etard traits dataset 
Mams_Etard <- read_csv("Etard_etal_2020/Vert_ABMR_gaps/Mammals.csv")
glimpse(Mams_Etard)

# load extinct species list 
extinct <- read_csv("../IUCN_data/Species_Pages/Outputs/Extinct_and_EW_sp.csv") %>% 
  pull(scientificName)

# Add synonyms to select species that have phylo and IUCN data also 
Mam_names <- read_csv("Synonyms/Mammal_final_name_matches.csv") 

data <- Mam_names %>% 
  left_join(Mams_Etard, by = c("From_trait" = "Best_guess_binomial"))%>%
  arrange(From_IUCN) %>% 
  # calculate value for habitat specialists. NOTE: probably redundant due to colinearity
  rowwise() %>%
  mutate(count = sum(c(Forest, Savanna, Shrubland, Grassland, Wetland,Rocky.areas,Caves.and.subterranean,Desert,
                       Marine,Marine.intertidal.or.coastal.supratidal,Artificial,Introduced.vegetation,Other.Unknown)
                     , na.rm = TRUE)) %>%
  ungroup()%>%
  mutate(specialist = case_when(count > 1 ~ 0, # non specialist
                                count == 1 ~ 1, # specialist
                                TRUE ~ NaN),
         specialist = na_if(specialist, NaN),
         Forest.specialist = case_when(specialist == 1 & Forest == 1 ~ 1, # forest specialist
                                       specialist == NA ~ NaN, # specialist
                                       TRUE ~ 0),
         Forest.specialist = na_if(specialist, NaN)
  )



  
data[1:10,c(33,34)]
colnames(data)
# Check Coverage 
glimpse(data)
# coverage of data = 1 - number of NA/ number of rows
coverage <- 1 - apply(data[,c(7:26,33,34)], 2, function(x) sum(is.na(x)))/nrow(data)
# Keep traits with coverage above 60%
coverage.keep <- names(coverage[which(coverage > 0.6)])

data <- data %>% 
  select(From_IUCN:From_phylo, paste(coverage.keep))

# Check correlation 
cor.check <- data %>%
  select(where(is.numeric))%>%
  na.omit()

cor(cor.check, cor.check) > 0.5



#### COMBINE DATA ###############################################################################################
# load synonyms,  imputed trait
combine_names <- read_csv("COMBINE/COMBINE/COMBINE_archives/taxonomy_crosswalk.csv")
# trait data
combine_data1 <- read_csv("COMBINE/COMBINE/COMBINE_archives/imputation_phylo_1.csv")
# imputed trait
combine_reported <- read_csv("COMBINE/COMBINE/COMBINE_archives/trait_data_reported.csv")

synonyms <- read_csv("Synonyms/Mam_IU_etard_syn.csv") %>%
  select(-internalTaxonId)

# Phylogenetic VCV
{
  # calculate maximum credibility tree from 1000 phylacine trees
  # phylacine <- read.nexus("Phylogenies/Mammals/PHYLACINE/DATA/DataS1/Complete_phylogeny.nex")

  # Mam_MCC_phylacine <- maxCladeCred(phylacine)
  # write.nexus(Mam_MCC_phylacine, file = "Phylogenies/Mammals/PHYLACINE/mam_MCC_phylacine.nex")
}
Mam_MCC_phylacine <- read.nexus("Phylogenies/Mammals/PHYLACINE/mam_MCC_phylacine.nex")

data2  <- combine_data1 %>%
  left_join(synonyms,by = c("iucn2020_binomial" = "Synonyms"))

length(unique(data2$From_IUCN)) + sum(is.na(data2$From_IUCN))

# join synonyms with Phylacine names of rows that did not match a current IUCN name
match2_data2 <- data2 %>% 
  filter(is.na(From_IUCN)) %>% 
  select(-From_IUCN) %>% 
  left_join(synonyms,by = c("phylacine_binomial" = "Synonyms")) %>% 
  filter(!is.na(From_IUCN))
  
# add synonyms matched through matching with Phylacine name and remove unmatched data
data2 <- data2 %>% 
  filter(!is.na(From_IUCN)) %>% 
  bind_rows(.,match2_data2)
  
  
# select traits that had >60% coverage before imputation or > 90% without imputation
# adult_mass_g             
# adult_body_length_mm     
# dphy_plant
# dphy_vertebrate          
# dphy_invertebrate        
# trophic_level   
# foraging_stratum         
# activity_cycle 
# dispersal_km
# habitat_breadth_n 
# litter_size_n
glimpse(data2)

data2 <- data2 %>% 
  select(From_IUCN, iucn2020_binomial, phylacine_binomial, order, adult_mass_g, adult_body_length_mm, habitat_breadth_n, dphy_plant, dphy_vertebrate, dphy_invertebrate, trophic_level, foraging_stratum, activity_cycle, dispersal_km, litter_size_n)

# Data that wasn't imputed but has coverage
# freshwater
# marine
# terrestrial_non-volant
# terrestrial_volant
# biogeographic_realm
glimpse(combine_reported)

r <- combine_reported %>% 
  select(iucn2020_binomial, freshwater, marine, terrestrial_non_volant = "terrestrial_non-volant", terrestrial_volant, biogeographical_realm) %>% 
  mutate(terrestrial =  case_when(terrestrial_non_volant == 1 ~ 1,
                                  terrestrial_volant == 1 ~ 1,
                                  terrestrial_volant == 0 &  terrestrial_non_volant == 0 ~ 0,
                                  TRUE ~ NA ))
unique(r$terrestrial)

# join impted and non imputed datasets
data2 <- left_join(data2, r, by = "iucn2020_binomial")

glimpse(data2)
# Save Final Dataset for Mammal traits

write_csv(data2, "Data_collated/Mam_imputed_traits_combine.csv")

# Adding range sizes to dataset 
data2 <- read_csv("Data_collated/Mam_imputed_traits_combine.csv")
mam_ranges <- read_csv("X:/edwards_lab1/User/bop21ipl/data/MAMS/Mammal_rangesizes_calc_ALL.csv")
mam_ranges <- mam_ranges %>% 
  select(From_IUCN = binomial, range_calculated)
data3 <- left_join(data2, mam_ranges, by = "From_IUCN") %>% 
  # remove extinct species 
  filter(!From_IUCN %in% extinct)

nrow(data3)  
nrow(mam_ranges)  

# count missing values 
data3 %>% 
  sapply(function(y) sum(length(which(is.na(y))))) %>% 
  data.frame()

data3 <- data3 %>% 
  filter(!is.na(range_calculated),
         !is.na(freshwater))

write_csv(data3, "Data_collated/Mam_imputed_traits_combine_wRange.csv")


# drop tips of phylotree that don't have range centroids
names <- data3 %>%
  mutate(phylacine_binomial = str_replace(phylacine_binomial, " ", "_"))  %>% 
  pull(phylacine_binomial)
length(names)

# read in list of names that have range centroids
mam_names_centroids <- read.csv("Distance_matrix/mam_species_Wcentroids.csv")
mam_names_centroids <- mam_names_centroids$species

# list of phylogenetic tree tip labels 
labels <- Mam_MCC_phylacine$tip.label
length(labels) # 5831
length(mam_names_centroids) #5304
# find list of species to drop 
labelsdiff <- setdiff(labels, mam_names_centroids)

# create tree for the species that have trait data 
tree <- drop.tip(Mam_MCC_phylacine, labelsdiff)

# create a variance covariance matrix for species with trait data
mam_vcv <- vcv(tree, corr = T)
mam_vcv[1:10,1:10]

# reorder the columns 
order <- sort(colnames(mam_vcv))
phylo_dist_matrixn <- mam_vcv[order,order]

#checks   
phylo_dist_matrixn[1:10,1:10]
# check the diagonals are all equal 
unique(phylo_dist_matrixn[col(phylo_dist_matrixn)==row(phylo_dist_matrixn)]) 
# number of speciesnames 
length(colnames(phylo_dist_matrixn)) # 5304

# make symetrical 
range(phylo_dist_matrixn)
phylo_dist_matrixn <- round(phylo_dist_matrixn, digits = 6)
isSymmetric(phylo_dist_matrixn)
library(lqmm)
is.positive.definite(phylo_dist_matrixn)
write.table(phylo_dist_matrixn, "Distance_matrix/Mammal_phylo_Distmatrix.txt")




