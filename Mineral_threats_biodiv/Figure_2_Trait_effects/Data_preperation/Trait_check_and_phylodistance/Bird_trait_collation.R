# Bird trait collation
# this scripts using the synonyms data from Bird_Name_fixing to match the species from the IUCN database to the trait datasets and the phylogenetic tree
# it also selects the traits for bird species that are suitable for imputation 
# and selects the species that have centroids and so can matched to the spatial distance matrix used at the analysis stage
# check correlation of traits 
library(tidyverse)
library(tibble); library(stringr); library(dplyr); library(readr);library(ggplot2);library(tidyr)
library(Rphylopars)
library(ape)
getwd()
setwd("X:/edwards_lab1/User/bop21ipl/Ecological trait data")
# load traits
Birds_Etard <- read_csv("Etard_etal_2020/Vert_ABMR_gaps/Birds.csv")
glimpse(Birds_Etard)

# load range data 
Bird_ranges <- read_csv("../IUCN_data/Species_Ranges/Outputs/Birds/Bird_ranges_ALL.csv")
Bird_ranges <- Bird_ranges %>% 
  select(binomial,range_logstd)

# load extinct species list 
extinct <- read_csv("../IUCN_data/Species_Pages/Outputs/Extinct_and_EW_sp.csv") %>% 
  pull(scientificName)

# load synonyms 
Name_final <- read_csv("Synonyms/Bird_final_name_matches.csv")
Name_final <- Name_final %>% 
  select(From_IUCN, From_trait, From_phylo)


# bind names to trait database 
Bird_traits <- Name_final %>% 
  left_join(Birds_Etard, by = c("From_trait" = "Best_guess_binomial")) %>% 
  # joining only species that have range data
  right_join(Bird_ranges, by = c("From_IUCN" = "binomial" )) %>% 
  filter(!is.na(From_phylo),
         !From_IUCN %in% extinct) %>% 
  select(-c(Order,Family,Genus,Note,Other.Unknown))
glimpse(Bird_traits) # 9636

# check the coverage =====
trait_na_count <- sapply(Bird_traits, FUN = function(y) sum(length(which(is.na(y))))) %>% 
  data.frame()
colnames(trait_na_count) <- "NA_count"
trait_na_count <- trait_na_count %>% 
  mutate(coverage = (nrow(Bird_traits)- NA_count)/nrow(Bird_traits),
         above60 = coverage > 0.6) %>% 
  rownames_to_column(var = "Trait")  

write_csv(trait_na_count, "Data_collated/bird_trait_na_count.csv")

good_coverage <- trait_na_count %>% 
  filter(above60 == TRUE) %>% 
  pull("Trait")

# filter traits to ones with suitable coverage for imputation 
Bird_traits_db <- Bird_traits %>% 
  select(good_coverage)

# join range data
Bird_traits_db <- Bird_traits_db %>% 
  rowwise() %>%
  mutate(count = sum(c(Forest, Savanna, Shrubland, Grassland, Wetland,Rocky.areas,Caves.and.subterranean,Desert,
                       Marine, Marine.intertidal.or.coastal.supratidal, Artificial,Introduced.vegetation)
                     , na.rm = TRUE)) %>%
  ungroup()%>%
  mutate(specialist = case_when(count > 1 ~ 0, # non specialist
                                count == 1 ~ 1, # specialist
                                TRUE ~ NaN),
         specialist = na_if(specialist, "NaN"),
         Forest.specialist = case_when(specialist == 1 & Forest == 1 ~ 1, # forest specialist
                                       specialist == NA ~ NaN, # specialist
                                       TRUE ~ 0),
         Forest.specialist = na_if(specialist, "NaN"))  %>% 
  select(-Artificial_habitat_use ) %>% 
  distinct()

glimpse(Bird_traits_db)


# check correlation ====
idx <- sapply(Bird_traits_db, FUN =  class) == "numeric"
trait_cor <- cor(na.omit(Bird_traits_db[,idx]),na.omit(Bird_traits_db[,idx])) > 0.5
trait_cor_num <- cor(na.omit(Bird_traits_db[,idx]),na.omit(Bird_traits_db[,idx])) 
# correlation of:
# marine and marine.intertidal 
# count and artificial
# forest and specialist
trait_cor_num[trait_cor]
Bird_traits_db <-  Bird_traits_db %>% 
  select(-c(specialist)) %>% 
  distinct()

# species with multiple phylo names for one IUCN name
multi <- Bird_traits_db %>%  
  group_by(From_phylo) %>% 
  count() %>% 
  filter(n > 1) %>% 
  pull(From_phylo)

# remove multi species 
Bird_traits_db <- Bird_traits_db %>% 
  filter(!From_phylo %in% multi)

write_csv(Bird_traits_db, "Data_collated/Bird_traits_final2.csv")


# Cleaning traits data base =====
Bird_traits_db <- read_csv("Data_collated/Bird_traits_final2.csv")
Bird_MCC <- read.nexus("Phylogenies/Birds/100birdMCC.nex")

# clean name columns
traits <- Bird_traits_db %>% 
  rename(species = From_phylo) %>%
  mutate(species = str_replace(species, " ", "_")) %>% 
  filter(!is.na(species)) %>% 
  select(-c(From_IUCN,From_trait)) %>% 
  distinct()

# clean Trophic level and diel activity into binary traits]
unique(traits$Trophic_level)
unique(traits$Diel_activity)

traits <- traits %>%  
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
  select(-c(Trophic_level, Diel_activity))

glimpse(traits)

# checking overdispersion ====
# examining the count variable for overdispersion: important for accurate imputation 
# if positive then data is overdisperesed.
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

write_csv(trait_dispersion, "Data_collated/Bird_trait_data_overdispersion.csv")


# Checking phylogenetic signal ====
# clip the phylo tree to the list of species that have at least one trait
names <- traits$species 
tree <- keep.tip(Bird_MCC, tip = names)

# Trait data for imputation Checking lambda values

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

trait_lambda <- bind_rows(lapply(2:23, FUN = lambda_check))
trait_lambda <- trait_lambda %>% 
  mutate(Above.6 = if_else(lambda > 0.6, TRUE, FALSE))

write_csv(trait_lambda , "Data_collated/Bird_trait_lambda.csv")

good_lambda <- trait_lambda %>%
  filter(Above.6 == TRUE) %>%
  pull(trait)
{# [1] "Body_mass_g"                             "Generation_length_d"                     "Litter_clutch_size"                     
# [4] "Forest"                                  "Grassland"                               "Wetland"                                
# [7] "Rocky.areas"                             "Caves.and.subterranean"                  "Marine"                                 
# [10] "Marine.intertidal.or.coastal.supratidal" "Carnivore"                               "Herbivore"     
}

# exclude Marine.intertidal.or.coastal.supratidal for correlation with Marine. see above.

# Clean data to add the only trait that need imputing ====
# remove traits with lambda values < 0.6
for_imputation <- traits %>% 
  select(species,
         Body_mass_g,
         Generation_length_d,
         Litter_clutch_size,
         Forest,
         Grassland,
         Wetland,
         Rocky.areas,
         Caves.and.subterranean,
         Marine,
         Carnivore,
         Herbivore,
         count) %>% 
  # center traits
  mutate(Body_mass_log_g = scale(log(Body_mass_g), center = TRUE, scale = TRUE)[,1],
         Generation_length_log_d = scale(log(Generation_length_d), center = TRUE, scale = TRUE)[,1],
         count_logstd = scale(log(1 + count), center = TRUE, scale = TRUE)[,1],
         Litter_clutch_size_log = scale(log(Litter_clutch_size), center = TRUE, scale = TRUE)[,1]) %>% 
  select(-c(Body_mass_g,
            Generation_length_d,
            Litter_clutch_size,
            count)) %>% 
  distinct()

length(for_imputation$species) # 9617
length(unique(for_imputation$species)) #9617
# multiple rows per species 
multi <- for_imputation %>% 
  group_by(species) %>% 
  count() %>% 
  filter(n > 1) %>% 
  pull(species)

for_imputation %>% 
  filter(species %in% multi)

glimpse(for_imputation)

for_imputation %>%
  ggplot(aes(x = Litter_clutch_size_log  )) +
  geom_density()


glimpse(for_imputation)

write_csv(for_imputation, "Data_collated/Bird_traits_for_impuation3.csv")

# Phylo_distance matrix #####
tree
# phylogenetic correlation matrix 
bird_vcv <- vcv(tree, corr = TRUE)
bird_vcv[1:10,1:10]

# reorder the columns 
order <- sort(colnames(bird_vcv))
phylo_dist_matrixn <- bird_vcv[order,order]

#checks   
phylo_dist_matrixn[1:10,1:10]
# check the diagonals are all equal 
unique(phylo_dist_matrixn[col(phylo_dist_matrixn)==row(phylo_dist_matrixn)]) 
# number of speciesnames 
length(colnames(phylo_dist_matrixn)) # 9758

# make symetrical 
range(phylo_dist_matrixn)
phylo_dist_matrixn <- round(phylo_dist_matrixn, digits = 6)
isSymmetric(phylo_dist_matrixn)
library(lqmm)
is.positive.definite(phylo_dist_matrixn)
write.table(phylo_dist_matrixn, "Distance_matrix/Bird_phylo_Distmatrix_for_imputation_model.txt")




# Clean data for model type 2 with more traits and less imputaion ======
glimpse(traits)
for_reduced_imp <- traits %>% 
  filter(!is.na(Habitat_breadth_IUCN),
         !is.na(Nocturnal)) # reduced by 125 species

trait_na_count <- sapply(for_reduced_imp, FUN = function(y) sum(length(which(is.na(y))))) %>% 
  data.frame()

# traits that still need some imputation even with reduced model
for_reduced_imputation <- for_reduced_imp %>% 
  mutate(Generation_length_log_d = scale(log(Generation_length_d), center = TRUE, scale = TRUE)[,1],
         Litter_clutch_size_log = scale(log(Litter_clutch_size), center = TRUE, scale = TRUE)[,1]) %>% 
  select(species, Generation_length_log_d,Litter_clutch_size_log, Carnivore, Herbivore)

write_csv(for_reduced_imputation, "Data_collated/Bird_traits_for_reduced_imputation.csv")

# traits for model that have full coverage
no_impitation <- for_reduced_imp %>% 
  mutate(Body_mass_log_g = scale(log(Body_mass_g), center = TRUE, scale = TRUE)[,1],
         count_logstd = scale(log(1 + count), center = TRUE, scale = TRUE)[,1]) %>% 
  select(-c(Generation_length_d,Litter_clutch_size, Carnivore, Herbivore, Omnivore)) # 9,492 

no_impitation %>% na.omit() #9,492 
  
write_csv(no_impitation, "Data_collated/Bird_traits_for_reduced_NO_imp_needed.csv")
