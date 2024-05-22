# Traits for Fish 

library(tidyverse)
library(readr)
library(ape)
library(taxize)
library(rfishbase)
rfishbase::

getwd()
# setwd("/Users/ieuan/Google Drive/My Drive/PhD/Ecological trait data")
setwd("/Volumes/shared/edwards_lab1/User/bop21ipl/Ecological trait data/")

# Name matching database Actino
Names_fish <- read_csv("Synonyms/Fish_Name_match.csv") %>% 
  select(-(...1))

# Read range values
Fish_ranges <- read_csv("../IUCN_data/Species_Ranges/Outputs/Fish/Fish_range_calc_all.csv")

fish_ranges <- Names_fish %>% 
  right_join(Fish_ranges, by = c("From_IUCN" = "BINOMIAL"))

fb_tables()

# "countecosystem", "diet","ecosystem", "ecology", "fecundity", "food", "foodtroph", "maturity", "morphdat", "predatortroph","reproduc"
ecology <- fb_tbl("ecology", server = "fishbase") %>% 
  select(SpecCode, DietTroph, DietSeTroph, FoodTroph, FoodSeTroph )
View(ecology)

species <- fb_tbl("species", server = "fishbase") %>% 
  mutate(sci_name = str_c(Genus, Species, sep = " ")) %>% 
  left_join(ecology, by = "SpecCode") %>% 
  right_join(fish_ranges, by = c("sci_name" = "fishbase_name")) %>% 
  filter(!is.na(From_phylo)) %>% 
  select(From_IUCN, From_phylo, sci_name, Length, Fresh, Brack, Saltwater, DietTroph, DietSeTroph, FoodTroph, FoodSeTroph ) %>% 
  distinct()

glimpse(species)

# check the coverage ====
trait_na_count <- sapply(species, FUN = function(y) sum(length(which(is.na(y))))) %>% 
  data.frame()
colnames(trait_na_count) <- "NA_count"
trait_na_count <- trait_na_count %>% 
  mutate(coverage = (nrow(species)- NA_count)/nrow(species),
         above60 = coverage > 0.6) %>% 
  rownames_to_column(var = "trait")

write_csv(trait_na_count,"Data_collated/Fish_trait_na_count.csv")

good_coverage <- trait_na_count %>% 
  filter(above60 == TRUE) %>% 
  rownames_to_column() %>% 
  pull("rowname")
#  "Length"    "Fresh"     "Brack"     "Saltwater"

# write species that fields 
write_csv(species, "Data_collated/Fish_fishbase_traits2.csv")

# join IUCN Names and Phylo names 
Fish_traits <- species %>% 
  select(species = From_phylo, Length, Fresh, Brack, Saltwater) %>% 
  mutate(species = str_replace(species, " ", "_")) %>% 
  distinct() # 6441

# check multinames
multi <- Fish_traits %>% 
  group_by(species) %>% 
  count() %>% 
  filter(n > 1) 

# remove species with multiple phylogenetic names relating to one IUCN name
Fish_traits <- Fish_traits %>% 
  filter(!species %in% multi$species) #6616

# checking overdispersion ====
# examining the count variable for overdispersion: important for accurate imputation 
# if positive then data is overdisperesed.
over_dispersed <- function(i) {
  trait_val <- Fish_traits[,i] %>% pull(.)
  dispersion <- var(trait_val, na.rm = TRUE) - mean(trait_val, na.rm = TRUE)
  tbl <- tibble(trait = colnames(Fish_traits[,i]),
                dispersion = dispersion)
  return(tbl)
}
over_dispersed(ncol(Fish_traits))

trait_dispersion <- bind_rows(lapply(2:ncol(Fish_traits), FUN = over_dispersed)) %>% 
  mutate(overdispersed = if_else(dispersion > 0 ,T, F))

write_csv(trait_dispersion, "Data_collated/Fish_trait_data_overdispersion.csv")

# log to help overdispersion 
Fish_traits <- Fish_traits %>% 
  mutate(Length_logst = scale(log(Length), scale= T, center = T)[,1]) %>% 
  select(-Length)

# check correlation ====
cor(na.omit(Fish_traits[,-1]), na.omit(Fish_traits[,-1]))
cor(na.omit(Fish_traits[,-1]), na.omit(Fish_traits[,-1])) > 0.5 | cor(na.omit(Fish_traits[,-1]), na.omit(Fish_traits[,-1])) < -0.5

# freshwater and salt water negatively correlated

# Checking phylogenetic signal ====
# clip the phylo tree to the list of species that have at least one trait
names <- Fish_traits$species 
length(names)
length(unique(names))

tree <- read.tree("Phylogenies/Fish/actinopt_12k_raxml.tre.xz")
tree <- keep.tip(tree, tip = names)

# Trait data for imputation Checking lambda values
# check for multiple names
Fish_traits %>%  
  group_by(species) %>% 
  count() %>% 
  filter(n > 1)

lambda_check <- function (i) {
  trait_temp <- Fish_traits[,c(1,i)]
  p_temp <- phylopars(trait_temp, tree, model = "lambda")
  lambda <- tibble(trait = colnames(Fish_traits[,i]),
                   lambda = p_temp$model$lambda)
  
  return(lambda)
}
lambda_check(2)

trait_lambda <- bind_rows(lapply(2:ncol(Fish_traits), FUN = lambda_check))
trait_lambda <- trait_lambda %>% 
  mutate(Above.6 = if_else(lambda > 0.6, TRUE, FALSE))

write_csv(trait_lambda , "Data_collated/Fish_trait_lambda.csv")

good_lambda <- trait_lambda %>%
  filter(Above.6 == TRUE) %>%
  pull(trait)

# Clean for imputation ====
glimpse(Fish_traits)

write_csv(Fish_traits, "Data_collated/Fish_traits_for_imputation2.csv")

