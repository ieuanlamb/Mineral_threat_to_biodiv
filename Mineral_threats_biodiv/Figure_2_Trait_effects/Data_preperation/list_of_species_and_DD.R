#list of species used in models 

library(tidyverse)

getwd()
setwd("X:/edwards_lab1/User/bop21ipl/Chapter_One2/Data")

# list of all vertebrates assessed by the IUCN 
# mine threatened species 
M_sp <- read.csv("X:/edwards_lab1/User/bop21ipl/IUCN_data/Species_Pages/Raw_Data/CHORDATA_pg_Oil_Mining_threat/Chordata_Mine_threatened_assessments.csv")

glimpse(M_sp)

# list species 
M_sp <- M_sp %>%
  pull(scientificName)

# Dataset for all Chordata and add taxonomy ====
Chordata <- read_csv("X:/edwards_lab1/User/bop21ipl/IUCN_data/Species_Pages/Outputs/Chordata_pg_ALL_03.05.23.csv")
glimpse(Chordata)

tax <- read_csv("X:/edwards_lab1/User/bop21ipl/IUCN_data/Species_Pages/Outputs/Chordata_taxonomy_03.05.23.csv")

# remove columns 
tax <- tax[,1:9]
glimpse(tax)

# Change internalTaxonId to a factor vector to join the Chrodata and tax dataframes by
tax <- tax %>%
  mutate(ID = as.factor(internalTaxonId))
Chordata <- Chordata %>%
  mutate(ID = as.factor(internalTaxonId))

# full joing of taxonomy and assessment data sets
Chordata <- full_join(Chordata, tax, by = "ID" ) 

# Create new data for Taxa groups
Chordata <- Chordata %>%
  mutate(Class = as.factor(className))
# rewrite the levels of class so that the Class are fish = fish
# resultant levels are "FISH"     "AMPHIBIA" "BIRDS"     "MAMMALIA" "REPTILIA"
levels(Chordata$Class)[c(1,4:5,7,9)] <- "Fish"
levels(Chordata$Class)[c(2,3,4,5)] <- c("Amphibians", "Birds","Mammals", "Reptiles")
#reorder the levels
Chordata <- Chordata %>%
  mutate(Class = fct_relevel(Class,"Amphibians", "Birds","Fish","Mammals", "Reptiles"))

# Check that the data sets match shoulf have 0 rows
Chordata %>%
  filter(scientificName.x != scientificName.y)

glimpse(Chordata)

# remove duplicate cols, adjust col names
Chordata <- Chordata %>%
  dplyr::select(-c(internalTaxonId.y , scientificName.y)) %>%
  rename(internalTaxonId =  internalTaxonId.x ,
         scientificName =  scientificName.x)


# Add Mine threatened variable -----  
# Use M_sp to create a BINOMIAL 1= threatened by mining, 0 = not threatened
Chordata <- Chordata %>%
  mutate(Minethrt = as.factor(if_else(scientificName %in% M_sp,1,0)),
         Category_abr = as.factor( case_when(redlistCategory == "Least Concern" | redlistCategory == "Lower Risk/least concern"  ~ "LC",
                                             redlistCategory ==  "Near Threatened"| redlistCategory == "Lower Risk/near threatened" ~ "NT",
                                             redlistCategory ==  "Vulnerable" ~ "VU",  
                                             redlistCategory == "Endangered"~ "EN", 
                                             redlistCategory == "Critically Endangered" ~ "CR", 
                                             redlistCategory == "Extinct" ~ "EX",  
                                             redlistCategory == "Data Deficient" ~ "DD" ,
                                             TRUE ~ "Other")))

glimpse(Chordata)


vert_cats <- Chordata %>% 
  select(From_IUCN = scientificName, Category_abr)

# amphibians in model ====
# full list of names including Name matches 
Names_final <- read_csv("Amph_trait_model/Amphibian_final_name_matches.csv")

# imputed amphibian traits 
imp_trait <- read_csv("Amph_trait_model/Amph_imputed_checked_final.csv")

# new species list
new_syn <- read_csv("X:/edwards_lab1/User/bop21ipl/IUCN_data/Species_Pages/Raw_Data/AMPHS_MAM_REP_pg_ALL_03.05.23/synonyms.csv")
new_syn <- new_syn %>% 
  mutate(species = str_c(genusName, speciesName, sep = " ")) %>% 
  select(scientificName, species) %>% 
  distinct()


# Read in Phylogenetic distance matrix
amph_phylodist <- read.table("Amph_trait_model/Distance_matrix/Amphibian_phylo_Distmatrix.txt")
# remove columns and rows of names not in trait database
amph_phylodist <- amph_phylodist[-which(names(amph_phylodist) %in% c("Pristimantis_w.nigrum", "Scinax_v.signatus", "Scinax_x.signatus")),-which(names(amph_phylodist) %in% c("Pristimantis_w.nigrum", "Scinax_v.signatus", "Scinax_x.signatus"))]
amph_phylodist <- as.matrix(amph_phylodist)
amph_phylodist[1:10,1:10]
isSymmetric.matrix(amph_phylodist)
identical(colnames(amph_phylodist),rownames(amph_phylodist))

amph_dist <- read.table("Amph_trait_model/Distance_matrix/Amphibian_Corr_matrix_SymPosDef.txt")
amph_dist[1:10, 1:10]
range(amph_dist)
amph_dist <- as.matrix(amph_dist)
identical(colnames(amph_dist),rownames(amph_dist))
isSymmetric(amph_dist)
isSymmetric.matrix(amph_dist)

setdiff(colnames(amph_phylodist),colnames(amph_dist)) # should be empty

# data set of mined threatened species 
amph_mined <- Names_final %>%
  mutate(mine_thrnd = if_else(From_IUCN %in% M_sp, 1, 0),
         species = str_replace(From_phylo, " ", "_")) %>%
  select(From_IUCN, species, mine_thrnd) %>%
  na.omit()

data <- left_join(amph_mined, imp_trait, by = "species") %>%
  filter(species %in% colnames(amph_phylodist))%>%
  distinct()%>%
  arrange(species) %>% 
  mutate(species_space  = species)

# add species threat category and summarise
data_cat <- left_join(data, vert_cats) %>% 
  select(From_IUCN, species, Category_abr, mine_thrnd)

# write_csv(data_cat,"Amph_trait_model/amph_splist_and_category.csv")

# after being checked manually 
data_cat <- read_csv("Amph_trait_model/amph_splist_and_category.csv") 

sum <- data_cat %>% 
  group_by(Category_abr) %>% 
  summarise(n = n()) %>% 
  # filter(Category_abr %in% c("CR", "DD", "EN", "LC", "NT", "VU")) %>% 
  ungroup() %>% 
  mutate(proportion = n/sum(n),
         total = sum(n))
 
# DD n = 932 (14%)


# Birds ====
# run in tandem with Bird_model_data_preperation
# data set of mined threatened species + range data
bird_mined <- Names_final %>%
  mutate(mine_thrnd = if_else(From_IUCN %in% M_sp, 1, 0),
         species = str_replace(From_phylo, " ", "_")) %>%
  left_join(bird_range, by = c("From_IUCN" = "binomial")) %>% 
  select(From_IUCN, species, mine_thrnd, range_logstd) %>%
  distinct() %>% 
  na.omit()

# add imputed data
data <- left_join(NONimp_trait, bird_mined, by = "species") %>%
  left_join(imp_trait, by = "species") %>%
  filter(species %in% colnames(bird_phylodist2))%>%
  distinct()%>%
  arrange(species) %>% 
  mutate(species_space  = species,
         # rounding imputed binomial traits to 1 or 0
         across(.cols = Carnivore:Herbivore, round ),
         across(!starts_with("species") & !From_IUCN, .fns = function(x) round(x,digits = 6))
         ) %>%
  na.omit()


glimpse(data)

# add species threat category and summarise
data_cat <- left_join(data, vert_cats) %>% 
  select(From_IUCN, species, Category_abr, mine_thrnd)


write_csv(data_cat, "Data/Bird_trait_model/bird_splist_and_category.csv")

