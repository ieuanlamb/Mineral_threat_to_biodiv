# Reptile trait collation 
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

# Load reptile trait data 
Reps_Etard <- read_csv("Etard_etal_2020/Vert_ABMR_gaps/Reptiles.csv")
glimpse(Reps_Etard)

# Load phylotree 
Rep_MCCtree <- read.tree("Phylogenies/Reptiles/Reptile_100_MCCtree.nexus")

# Load name matching database 
Rep_names <- read_csv("Synonyms/Reptile_Final_namematch2.csv")

# species that need have trait information for two species but only one tip.label on the phylogentetic tree 
# *should average these trait data or remove them from the tree and NOT impute their trait data. 
Rep_names %>% 
  group_by(From_Phylo) %>% 
  count() %>% 
  filter(n >1)

rep_range_data <- read_csv("../IUCN_data/Species_Ranges/Outputs/Reptiles/Rep_rangesizes_calc_sp_groups_ALL.csv") %>% 
  select(-SHAPE_area)

# combine range data to etard traits and calculate a specialist trait

data <- Rep_names %>% 
  left_join(Reps_Etard, by = c("From_Trait" = "Best_guess_binomial"))%>%
  arrange(From_IUCN) %>% 
  # join and reduce to species that have range data
  right_join(rep_range_data, by = c("From_IUCN" = "binomial")) %>% 
  # calculate value for habitat specialists. NOTE: probably redundant due to colinearity
  rowwise() %>%
  mutate(count = sum(c(Forest, Savanna, Shrubland, Grassland, Wetland,Rocky.areas,Caves.and.subterranean,Desert,
                       Marine,Marine.intertidal.or.coastal.supratidal,Artificial,Introduced.vegetation,Other.Unknown)
                     , na.rm = FALSE)) %>%
  ungroup()%>%
  mutate(specialist = case_when(count > 1 ~ 0, # non specialist
                                count == 1 ~ 1, # specialist
                                TRUE ~ NaN),
         specialist = na_if(specialist, "NaN"),
         Forest.specialist = case_when(specialist == 1 & Forest == 1 ~ 1, # forest specialist
                                       specialist == NA ~ NaN, # specialist
                                       TRUE ~ 0),
         Forest.specialist = na_if(specialist, "NaN"),
         Carnivore = case_when(Trophic_level == "Carnivore" ~ 1,
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
                                TRUE ~ 0)) %>% 
           select(-c(Trophic_level, Diel_activity, Artificial_habitat_use, Order, Family, Genus, Other.Unknown, Note,Introduced.vegetation)) %>% 
    filter(!is.na(From_Phylo)) %>% 
    distinct()
  
glimpse(data)


# Check Coverage ======
glimpse(data)

# check the coverage
trait_na_count <- sapply(data, FUN = function(y) sum(length(which(is.na(y))))) %>% 
  data.frame()
colnames(trait_na_count) <- "NA_count"
trait_na_count <- trait_na_count %>% 
  mutate(coverage = (nrow(data)- NA_count)/nrow(data),
         above60 = coverage > 0.6) %>% 
  rownames_to_column(var = "Trait")  

write_csv(trait_na_count, "Data_collated/Rep_trait_na_count.csv")

good_coverage <- trait_na_count %>% 
  filter(above60 == TRUE) %>% 
  pull("Trait")


# reduce dataset to good coverage 
data <- data %>% 
  select(species = From_Phylo, Body_mass_g, range_calculated) 

# check correlation ====
cor(na.omit(data[,2:3]),na.omit(data[,2:3])) 
cor(na.omit(data[,2:3]),na.omit(data[,2:3])) > 0.5
#                   Body_mass_g range_calculated
# Body_mass_g        1.0000000        0.1474162
# range_calculated   0.1474162        1.0000000

# checking overdispersion ====
# examining the count variable for overdispersion: important for accurate imputation 
# if positive then data is overdisperesed.
over_dispersed <- function(i) {
  trait_val <- data[,i] %>% pull(.)
  dispersion <- var(trait_val, na.rm = TRUE) - mean(trait_val, na.rm = TRUE)
  tbl <- tibble(trait = colnames(data[,i]),
                dispersion = dispersion)
  return(tbl)
}
over_dispersed(ncol(data))

trait_dispersion <- bind_rows(lapply(2:ncol(data), FUN = over_dispersed)) %>% 
  mutate(overdispersed = if_else(dispersion > 0 ,T, F))

write_csv(trait_dispersion, "Data_collated/Rep_trait_data_overdispersion.csv")

# log and standardise to help with over dispersion
data <- data %>% 
  mutate(Body_mass_log_st = scale(log(data$Body_mass_g), center = TRUE, scale = TRUE)[,1],
         Range_log_st = scale(log(data$range_calculated), center = TRUE, scale = TRUE)[,1]) %>% 
  distinct()

# Checking phylogenetic signal ====
# clip the phylo tree to the list of species that have at least one trait
names <- data$species # 7907
tree <- keep.tip(Rep_MCCtree, tip = names)

# Trait data for imputation Checking lambda values
# check for multiple names
data %>%  
  group_by(species) %>% 
  count() %>% 
  filter(n > 1)
# Anolis_luciae somehow duplicated but not removed with distinct remove with top_n
data <- data %>%  
  group_by(species) %>% 
  slice_head(n = 1)

lambda_check <- function (i) {
  trait_temp <- data[,c(1,i)]
  p_temp <- phylopars(trait_temp, tree, model = "lambda")
  lambda <- tibble(trait = colnames(data[,i]),
                   lambda = p_temp$model$lambda)
  
  return(lambda)
}

trait_lambda <- bind_rows(lapply(2:ncol(data), FUN = lambda_check))
trait_lambda <- trait_lambda %>% 
  mutate(Above.6 = if_else(lambda > 0.6, TRUE, FALSE))

write_csv(trait_lambda , "Data_collated/Rep_trait_lambda.csv")

good_lambda <- trait_lambda %>%
  filter(Above.6 == TRUE) %>%
  pull(trait)

# keep traits thhat have lambda > 0.6 
data <- data %>% 
  select(species, all_of(good_lambda)) %>% 
  ungroup()

glimpse(data) # 7906
write_csv(data, "Data_collated/Reptile_trait_final2.csv")
data <- read_csv("Data_collated/Reptile_trait_final2.csv")

# traits for imputation
# impute bodysize 
for_imputation <- data %>% 
  select(species, Body_mass_log_st)

write_csv(for_imputation, "Data_collated/Reptile_trait_for_imputaion2.csv")

# calculate phylogenetic correlation matrix ########
rep_vcv <- vcv(tree, corr = T)
rep_vcv[1:10,1:10]

# reorder the columns 
order <- sort(colnames(rep_vcv))
phylo_dist_matrixn <- rep_vcv[order,order]

#checks   
phylo_dist_matrixn[1:10,1:10]
# check the diagonals are all equal 
unique(phylo_dist_matrixn[col(phylo_dist_matrixn)==row(phylo_dist_matrixn)]) 
# number of speciesnames 
length(colnames(phylo_dist_matrixn)) # 7906

# make symetrical 
range(phylo_dist_matrixn)
phylo_dist_matrixn <- round(phylo_dist_matrixn, digits = 6)
isSymmetric(phylo_dist_matrixn)
library(lqmm)
is.positive.definite(phylo_dist_matrixn)
write.table(phylo_dist_matrixn, "Distance_matrix/Rep_phylo_Distmatrix2.txt")



