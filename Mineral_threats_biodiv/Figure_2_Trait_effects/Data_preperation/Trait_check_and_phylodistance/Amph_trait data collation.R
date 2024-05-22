# Trait_data 
# Collation 

library(tidyverse)
library(readr)
library(ape)
library(taxize)
library(forcats)
library(Rphylopars)

getwd()
# setwd("/Users/ieuan/Google Drive/My Drive/PhD")
setwd("X:/edwards_lab1/User/bop21ipl")

# Mine threatened species 
M_sp <- read.csv(file = "C:/Users/bop21ipl/My Drive/PhD/IUCN data/Species_Pages/Raw_Data/CHORDATA_pg_Oil_Mining_threat/Chordata_Mine_threatened/assessments.csv")
M_sp <- pull(M_sp, scientificName)

# taxonomy list form IUCN Species pages
# tax <- read.csv(file = ("IUCN data/Species_Pages/Outputs/Chordata_taxonomy_01.03.22.csv"), header = TRUE)
tax <- read.csv("C:/Users/bop21ipl/My Drive/PhD/IUCN data/Species_Pages/Outputs/Chordata_taxonomy_01.03.22.csv")
# remove columns 
tax <- tax[,2:10]
# Create new data for Taxa groups
tax <- tax %>%
  mutate(Class = as.factor(className))
# rewrite the levels of class so that the Class are fish = fish
# resultant levels are "FISH"     "AMPHIBIA" "BIRDS"     "MAMMALIA" "REPTILIA"
levels(tax$Class)[c(1,4:5,7,9)] <- "Fish"
levels(tax$Class)[c(2,3,4,5)] <- c("Amphibians", "Birds","Mammals", "Reptiles")
#reorder the levels
tax <- tax %>%
  mutate(Class = fct_relevel(Class,"Amphibians", "Birds","Fish","Mammals", "Reptiles"))


# Read in the Ecological trait datasets collated by Etard et al. 2020 +++++++++++++++++++++++++++++
Amphs_Etard <- read_csv("Ecological trait data/Etard_etal_2020/Vert_ABMR_gaps/Amphibians.csv")

glimpse(Amphs_Etard)

# read 100 amphibian trees from vert life 
{# amph_tree <- read.nexus("Phylogenies/Amphibians/output.nex")

# Get the MCC tree 
# MCC_tree <- maxCladeCred(amph_tree)

# write the MCC tree 
# write.nexus(MCC_tree, file = "Phylogenies/Amphibians/100amphMCC.nex")
}
# read in the MCC tree 
amph_MCC_tree <- read.nexus("Ecological trait data/Phylogenies/Amphibians/100amphMCC.nex")

# species names for amphibians in phylo 
tree_tip_amph <- amph_MCC_tree$tip.label

# Read in name matches
# Names_final <- read_csv("Ecological trait data/Synonyms/Amphibian_final_name_matches.csv") %>%
#   select( -(...1))

# final names list. species with ranges and with a phylogenetic position
Names_final <- read_csv("Ecological trait data/Synonyms/Amphibian_name_sp_w_range_and_phy.csv") %>%
  select( -(...1))

#calculate area of species ranges 
{
  # NOW IN Sp_ranges_calculated.R SCRIPT!!!!!
  # library(sf)
  # # 1) find species with multiple seasonal ranges
  # Amph_ranges <- Amph_ranges %>% mutate( seasonal = as.factor(seasonal)) 
  # levels(Amph_ranges$seasonal) # all amphibians only have resident ranges therefore we can assume there is no overlap of shape files
  # 
  # # sum the Shape_Area
  # Amph_ranges_calc <- Amph_ranges %>%
  #   tibble() %>%
  #   group_by(binomial) %>%
  #   # units used by the IUCN are km^2
  #   summarise(range_calculated = sum(SHAPE_Area))
  # 
  # write.csv(Amph_ranges_calc, "../IUCN data/Species_Ranges/Outputs/Amphibians/Amph_RangeSizes.csv")
  }

Amph_ranges_calc <- read_csv("IUCN_data/Species_Ranges/Outputs/Amphibians/Amphs_rangesizes_calc.csv")

# Read species synonyms list 
IU_amphiweb_etard_syn <- read_csv("Ecological trait data/Synonyms/IU_amphiweb_etard_syn.csv")
# Extract species range sizes from range_calculated col 
{
  Amph_Rsize <- tibble(From_IUCN = Amph_ranges_calc$binomial,
                           Range_km2 = Amph_ranges_calc$range_calculated) 
  # correct any names that are not matching 
  Amph_Rsize_wNames <- left_join(Amph_Rsize, IU_amphiweb_etard_syn, by = c("From_IUCN" = "Synonyms")) %>%
    filter(!(From_IUCN != From_IUCN.y & From_IUCN.y %in% .$From_IUCN))
  
  # Checks
  multi_names <- Amph_Rsize_wNames %>% # checking for name duplicates 
    group_by(From_IUCN)%>%
    count()%>%
    filter(n > 1) # 0 
  Amph_Rsize_wNames %>%   # check for non matching names between datasets
    filter(From_IUCN %in% multi_names$From_IUCN,
           !From_IUCN.y %in% Amph_Rsize_wNames$From_IUCN)  # 0 
  Amph_Rsize_wNames%>% # check for species missing in iucn taxonomy 
    filter(is.na(From_IUCN.y)) # 1 = Indirana tenuilingua no population found and not assessed (DD)
  
  # CLEAn the range dataset
  Amph_Rsize_wNames <- Amph_Rsize_wNames %>%
    select(From_IUCN = From_IUCN.y, Range_km2)
}


# prepare species data for phylogenetic imputation 
Amphs_trait <- Amphs_Etard %>%
  # select traits that have >40% coverage 
  rename(From_trait = Best_guess_binomial) %>%
  # calculate value for habitat speciealists. NOTE: probably redundant due to colinearity
  rowwise() %>%
  mutate(count = sum(c(Forest, Savanna, Shrubland, Grassland, Wetland,Rocky.areas,Caves.and.subterranean,Desert,
                       Marine,Marine.intertidal.or.coastal.supratidal,Artificial,Introduced.vegetation,Other.Unknown)
                     , na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(specialist = case_when(count > 1 ~ 0, # non specialist
                                count == 1 ~ 1, # specialist
                                TRUE ~ NaN),
         specialist = na_if(specialist, "NaN") ) %>%
  # join data to the list of names (ie all names remain in dataset)
  right_join(Names_final, by = "From_trait") %>%
  # add species ranges to the data set
  left_join(Amph_Rsize_wNames, by = "From_IUCN") %>% 
  # remove duplicated rows
  distinct() %>%
  filter(!is.na(From_phylo)) %>%
  mutate(From_phylo = str_replace(From_phylo," ","_"),
         # Standardise Body length and Range size around the mean
         BLlog_std = scale(log(Body_length_mm), center = TRUE, scale = TRUE)[,1],
         # log range size then standardise
         Rangelogkm2_std = scale(log(Range_km2), center = TRUE, scale = TRUE)[,1],
         Hab_bread_IU_std = scale(Habitat_breadth_IUCN, center = TRUE, scale = TRUE)[,1],)
  # rename(species = From_phylo)
  # select(species = From_phylo, BLlog_std, Body_length_mm, Rangelogkm2_std, Range_km2, Hab_bread_IU_std, Habitat_breadth_IUCN, specialist)
  
unique(Amphs_trait$Trophic_level)
unique(Amphs_trait$Diel_activity)

Amphs_trait <- Amphs_trait %>% 
  mutate(Carnivore = case_when(Trophic_level == "Carnivore" ~ 1,
                             is.na(Trophic_level) ~ NaN,
                             TRUE ~ 0),
       Herbivore = case_when(Trophic_level == "Herbivore" ~ 1,
                             is.na(Trophic_level) ~ NaN,
                             TRUE ~ 0),
       Omnivore = case_when(Trophic_level == "Omnivore" ~ 1,
                            is.na(Trophic_level) ~ NaN,
                            TRUE ~ 0),
       Nocturnal = case_when(Diel_activity  == "Nocturnal" ~ 1,
                             is.na(Diel_activity) ~ NaN,
                             TRUE ~ 0),) %>% 
  select(-c(Trophic_level, Diel_activity, Artificial_habitat_use))

# reorder to suit rphylopars structure
Amphs_trait <- Amphs_trait %>% 
  mutate(species = From_phylo) %>% 
  relocate(species) %>% 
  select(-c(Order, Family,Genus, From_trait,From_IUCN, From_phylo, Note, Other.Unknown))
  

# species that have multiple rows due to splitting synonyms
multiname <- Amphs_trait %>%
  group_by(species) %>%
  count()%>%
  filter(n>1) %>% 
  pull(species) # 29 species names 
# remove speceis with multiple rows for one phylo species name, 
# NOTE if multiple rows are included phylopars will treat rows as multiple observations per species and repredict all values based on phylogenetic signal 
Amphs_trait <- Amphs_trait %>%
  filter(!species %in% multiname)
Amphs_trait #  6,639 

glimpse(Amphs_trait)
# check coverage 
trait_na_count <- sapply(Amphs_trait, FUN = function(y) sum(length(which(is.na(y))))) %>% 
  data.frame()
colnames(trait_na_count) <- "NA_count"
trait_na_count <- trait_na_count %>% 
  mutate(coverage = (nrow(Amphs_trait)- NA_count)/nrow(Amphs_trait),
         above60 = coverage > 0.6) %>% 
  rownames_to_column(var = "Trait") 

write_csv(trait_na_count, "Ecological trait data/Data_collated/Amph_trait_na_count.csv")

# list of traits cith coverage over 60% 
good_coverage <- trait_na_count %>% 
  filter(above60 == TRUE) %>% 
  pull("Trait")

traits <- Amphs_trait %>% 
  select(good_coverage)

glimpse(traits)

# check overdispersion =====
# overdispersed data can lead to poor imputation 
over_dispersed <- function(i) {
  trait_val <- traits[,i] %>% pull(.)
  dispersion <- var(trait_val, na.rm = TRUE) - mean(trait_val, na.rm = TRUE)
  tbl <- tibble(trait = colnames(traits[,i]),
                dispersion = dispersion)
  return(tbl)
}
over_dispersed(ncol(traits))

trait_dispersion <- bind_rows(lapply(2:ncol(traits), FUN = over_dispersed)) %>% 
  mutate(overdispersed = if_else(dispersion > 0 ,T, F))

write_csv(trait_dispersion, "Ecological trait data/Data_collated/Amph_trait_data_overdispersion.csv")

# Checking phylogenetic signal ====
# clip the phylo tree to the list of species that have at least one trait

# read tree
amph_MCC_tree <- read.nexus("Ecological trait data/Phylogenies/Amphibians/100amphMCC.nex")
amph_MCC_tree <- drop.tip(amph_MCC_tree, "Homo_sapiens")

# Assign tree
species_list <- traits$species
tree <- keep.tip(amph_MCC_tree, tip = species_list)

# check for multiple names
traits %>%  
  group_by(species) %>% 
  count() %>% 
  filter(n > 1)

lambda_check <- function (i) {
  trait_temp <- traits[,c(1,i)]
  p_temp <- phylopars(trait_temp, tree, model = "lambda")
  lambda <- tibble(trait = colnames(traits[,i]),
                   lambda = p_temp$model$lambda)
  
  return(lambda)
}

trait_lambda <- bind_rows(lapply(2:ncol(traits), FUN = lambda_check))

# lambda > 0.6 considered reasonable for imputation 
trait_lambda <- trait_lambda %>% 
  mutate(Above.6 = if_else(lambda > 0.6, TRUE, FALSE))

write_csv(trait_lambda , "Ecological trait data/Data_collated/Amph_trait_lambda.csv")

good_lambda <- trait_lambda %>%
  filter(Above.6 == TRUE) %>%
  pull(trait)
 {#[1] "Body_length_mm"         "Habitat_breadth_IUCN"   "Forest"                 "Savanna"                "Wetland"                "Rocky.areas"            "Caves.and.subterranean"
#  [8] "count"                  "BLlog_std"              "Rangelogkm2_std"        "Hab_bread_IU_std"    
}

# reduce to traits with suitable lambda  and coverage for imputation 
glimpse(traits)
traits <- traits %>% 
  select(species, all_of(good_lambda), -Body_length_mm, -Habitat_breadth_IUCN )

# correlation matrix 
cor(na.omit(traits[,-1]))
cor(na.omit(traits[,-1])) > 0.5

# remove count for correlation with multiple other variables
traits <- traits %>% 
  select(-count)

write_csv(traits, "Ecological trait data/Data_collated/Amph_traits_final_ver3.csv")

