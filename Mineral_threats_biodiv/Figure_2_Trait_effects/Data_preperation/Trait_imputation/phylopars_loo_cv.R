# phylopars LOO Cross Validation 
# based on a script provided by Oscar Morton 

library(Rphylopars)
library(tidyverse)
library(stringr); library(dplyr);library(readr);library(ggplot2);library(tidyr);library(tibble)
library(ape)
library(phangorn)
library(treeplyr)
library(ggtree)

getwd()
# Reptile Trait data ######################
trait_data <- read_csv("Data_collated/Reptile_trait_for_imputaion2.csv")
trait_data <- trait_data %>% na.omit() %>% 
  select(species, contains("_st")) %>% 
  distinct()# 7722

names <-trait_data$species

# read in the RepMCC tree 
tree <- read.tree("Phylogenies/Reptiles/Reptile_100_MCCtree.nexus")
tree <- keep.tip(tree, names)

setdiff(names,tree$tip.label)
length(unique(names))
length(names)

## Get cov matrix
p <- phylopars(trait_data = trait_data,tree = tree)
## make empty matrices
imputed_var <- imputed_ind <- imputed_sp <- trait_data
imputed_var[,-1] <- imputed_sp[,-1] <- NA
# save points
save <- c(seq(1000,7000, by = 1000),nrow(trait_data))
## run LOO CV
for(i in 1:nrow(trait_data))
{
  print(i)
    temp_trait_data <- trait_data
    temp_trait_data[i,2] <- NA
    p_temp <- phylopars(trait_data = temp_trait_data,tree = tree, phylocov_fixed = 
                          p$pars$phylocov)
    # create dataset of just imputed values when left out
    imputed_sp[i,2] <- left_join(select(temp_trait_data, species), 
                                   rownames_to_column(as.data.frame(p_temp$anc_recon)), 
                                   by = c("species" = "rowname"))[i,2]
    imputed_var[i,2] <- left_join(select(temp_trait_data, species), 
                                    rownames_to_column(as.data.frame(p_temp$anc_var)), 
                                    by = c("species" = "rowname"))[i,2]

  if(i %in% save) {
    # variants of imputed datapoints
    write_csv(imputed_var, paste0("Data_collated/Rep_loocv_imputed/rep_loocv_",i,"imputed_var.csv"))
    # imputed data points using loo 
    write_csv(imputed_sp, paste0("Data_collated/Rep_loocv_imputed/rep_loocv_",i,"imputed_sp.csv"))
    }
}

# # Data frame of full imputation of missing data with variance terms 
# trait_p_var <- as.data.frame(p$anc_var) 
# colnames(trait_p_var) <- paste(colnames(trait_p_var), "var", sep = "_")
# 
# # imputed data from phylopars in a dataframe  
# trait_data_test <- p$anc_recon %>% 
#   data.frame() %>% 
#   rownames_to_column(var = "species") %>% 
#   tibble() %>% 
#   filter(species %in% imputed_sp$species) %>% 
#   arrange(species)

# variants of imputed datapoints
write_csv(imputed_var, "Data_collated/rep_loocv2_imputed_var.csv")
# imputed data points using loo 
write_csv(imputed_sp, "Data_collated/rep_loocv2_imputed_sp.csv")
# trait data before imputation (ie. raw traits)
write_csv(imputed_ind, "Data_collated/rep_loocv2_imputed_ind.csv")

write_csv(imputed_var, "Google Drive/My Drive/PhD/Back_ups/Data_backups/rep_loocv2_imputed_var.csv")
# imputed data points using loo 
write_csv(imputed_sp, "Google Drive/My Drive/PhD/Back_ups/Data_backups/rep_loocv2_imputed_sp.csv")
# trait data before imputation (ie. raw traits)
write_csv(imputed_ind, "Google Drive/My Drive/PhD/Back_ups/Data_backups/rep_loocv2_imputed_ind.csv")

imputed_var <- read_csv( "Data_collated/rep_loocv_imputed_var.csv")
imputed_sp <- read_csv( "Data_collated/rep_loocv_imputed_sp.csv")
imputed_ind <- read_csv("Data_collated/rep_loocv_imputed_ind.csv")

glimpse(imputed_var)
glimpse(imputed_sp)
glimpse(imputed_ind)

# check NA's
trait_na <- imputed_ind %>% sapply(FUN = function(y) sum(length(which(is.na(y)))))


## Summarise error
# MAE absolute error 
sp_mean_abs_error <- colMeans(imputed_sp[,-1] - trait_data_test[,-1],na.rm=TRUE)
# Body_mass_log_st     Range_log_st 
# -0.0003005101    -0.0138353329 

## p squared prediction coefficient Guenard et al 2013.
## P2 = 1 when all predictions perfectly match the observations, whereas values below 1 indicate imperfect predictions. Values P2 close to 0 
## (negative or positive) indicate that predictions have poor accuracy, being no better thanwhat would be expected from chance alone.
##http://adn.biol.umontreal.ca/~numericalecology/Reprints/Guenard_et_al_MEE_2013.pdf
# sp_mean_sq_error <- colMeans((Mean_imp_sp[,-1] - Actual[,-1])^2,na.rm=TRUE)
# Sam_var <- Actual[,-1] %>% summarise_all(~ var(.x, na.rm = TRUE))
# 1 - sp_mean_sq_error/Sam_var

# mean squared prediction error 
# P^2 = 1 - (sum(yi - predicted_yi)^2) / (true variance) 
# # P^2 = 1 - MSE mean squared error / (true variance) 
# 
MPSE <- colMeans((imputed_sp[,-1] - imputed_ind[,-1])^2,na.rm=TRUE) 
Sam_var <- imputed_ind[,-1] %>% summarise_all(~ var(.x, na.rm=TRUE))
p_squared <- 1 - MPSE/Sam_var
# Body_mass_log_st 0.7421824
# Range_log_st     0.1337892

# long form
p_squared <- 1 - colMeans((imputed_sp[,-1] - imputed_ind[,-1])^2, na.rm=TRUE)/var(imputed_ind[,-1], na.rm = TRUE) %>% 
  diag() %>% data.frame()

# Evaluate the imputation 
# compare real data to imputed loo data
trait_real_comp1 <- trait_data %>% 
  filter(!is.na(Body_mass_log_st)) %>% 
  mutate(across( Body_mass_g:Range_log_st, function(x) round(x,digits = 6))) %>% 
  distinct()

multi <- trait_real_comp1 %>%  
  group_by(species) %>% 
  
  count() %>% 
  filter(n > 1) %>% 
  pull(species)

trait_real_comp1 %>% 
  filter(species == multi)

# loo imputed traits of real data 
impute_loo_comp1 <- imputed_sp %>% 
  filter(species %in% trait_real_comp1$species) %>% 
  tibble() %>% 
  distinct()

length(unique(trait_real_comp1$species))
setdiff(trait_real_comp1$species, impute_loo_comp1$species)

# transform back to real scale 
# real mean and sd 
# ?scale to find the values that were used to standardise the data before imputation 
# to find whether the imputation is within reasonable boundaries 
trait_sumry <- trait_data %>% 
  summarise(BM_mean = mean(Body_mass_g, na.rm = TRUE),
            BM_sd = sd(Body_mass_g, na.rm = TRUE),
            RS_mean = mean(range_calculated, na.rm = TRUE),
            RS_sd = sd(range_calculated, na.rm = TRUE)) 

body_mass <- tibble(BM_real = trait_real_comp1$Body_mass_log_st,
                    BM_loo_imp = impute_loo_comp1$Body_mass_log_st) %>% 
  mutate(BM_real_t = ((BM_real + trait_sumry$BM_mean)*trait_sumry$BM_sd),
         BM_imp_t = ((BM_loo_imp+ trait_sumry$BM_mean)*trait_sumry$BM_sd) ,
          difference = (BM_real_t - BM_imp_t)^2) 


body_mass %>% 
  ggplot(aes(x = BM_real_t, y = BM_imp_t, colour =  difference )) + 
  geom_point(size = 3, alpha = 0.5)+
  geom_abline(slope = 1, intercept = 0) 
  # geom_abline(slope = 1, intercept = trait_sumry$BM_mean) + 
  # geom_abline(slope = 1, intercept = -trait_sumry$BM_mean)  
  # 


#compare real range values to loo imputed values 
trait_real_comp2 <- trait_data %>% 
  filter(!is.na(Range_log_st)) %>% 
  mutate(across( Body_mass_g:Range_log_st, function(x) round(x,digits = 6))) %>% 
  distinct()

multi <- trait_real_comp2 %>%  
  group_by(species) %>% 
  count() %>% 
  filter(n > 1) %>% 
  pull(species)

trait_real_comp2 <- trait_real_comp2 %>% 
  filter(species == multi)

impute_loo_comp2 <- imputed_sp %>% 
  filter(species %in% trait_real_comp2$species) %>% 
  tibble() %>% 
  distinct()

length(unique(trait_real_comp2$species))
setdiff(trait_real_comp2$species, impute_loo_comp2$species)

## transform back to real scale 
# real mean and sd 
range_size <- tibble(RS_real = trait_real_comp2$Range_log_st,
                    RS_loo_imp = impute_loo_comp2$Range_log_st) %>% 
  mutate(RS_real_t = ((RS_real + trait_sumry$RS_mean)*trait_sumry$RS_sd),
         RS_imp_t = ((RS_loo_imp+ trait_sumry$RS_mean)*trait_sumry$RS_sd) ,
         difference = (RS_real_t - RS_imp_t)^2) 


range_size %>% 
  ggplot(aes(x = RS_real_t, y = RS_imp_t, colour =  difference )) + 
  geom_point(size = 3, alpha = 0.5)+
  geom_abline(slope = 1, intercept = 0) 


# Apmh Trait data  ##############################################
trait_data <- read_csv("Data_collated/Amphibian_for_imputation_ver2.csv")

names <- trait_data$species

# read in the amphMCC tree 
tree <- read.nexus("Phylogenies/Amphibians/100amphMCC.nex")
tree <- keep.tip(tree, names)


# traits 1 ======
traits1 <- trait_data %>% 
  select(species,BLlog_std:count_std) %>% 
  distinct()

length(traits1$species)
length(unique(traits1$species))

# multiple rows per species 
multi <- traits1 %>% 
  group_by(species) %>% 
  count() %>% 
  filter(n > 1) %>% 
  pull(species)

# traits1 %>% 
#   filter(species %in% multi) %>% 
#   group_by(species) %>% 
#   summarise(across(.cols = c(range_logstd:count_log), .fns = function(y) mean(y, na.rm = TRUE )))

# remove species with multiple rows
traits1 <- traits1 %>% 
  filter(!species %in% multi)
# check
length(traits1$species)
length(unique(traits1$species))

# imputation
p1 <- phylopars(trait_data = traits1, tree = tree, phylo_correlated = TRUE)

## make empty matrices
imputed_var <- imputed_ind <- imputed_sp <- traits1
imputed_var[,-1] <- imputed_sp[,-1] <- NA
save <- c( seq(1000, 6000, by = 1000), nrow(traits1))
## run LOO CV
for(i in 1:nrow(traits1))
{
  print(i)
  for(j in 1:(ncol(traits1)-1))
  {
    print(j)
    temp_trait_data <- traits1
    temp_trait_data[i,j+1] <- NA
    p_temp <- phylopars(trait_data = temp_trait_data,tree = tree, phylocov_fixed = 
                          p1$pars$phylocov)
    
    imputed_sp[i,j+1] <- left_join(select(temp_trait_data, species), 
                                   rownames_to_column(as.data.frame(p_temp$anc_recon)), 
                                   by = c("species" = "rowname"))[i,j+1]
    imputed_var[i,j+1] <- left_join(select(temp_trait_data, species), 
                                    rownames_to_column(as.data.frame(p_temp$anc_var)), 
                                    by = c("species" = "rowname"))[i,j+1]
  }
  if(i %in% save){
    write_csv(imputed_var, paste0("Data_collated/AmphG_loocv_back_up/amphG1_loocv_sav",i,"_var.csv"))
    write_csv(imputed_sp, paste0("Data_collated/AmphG_loocv_back_up/amphG1_loocv_sav",i,"_sp.csv"))
  }
}


write_csv(imputed_var, "Data_collated/AmphG_loocv_back_up/amphG1_loocv_imputed_var.csv")
write_csv(imputed_sp, "Data_collated/AmphG_loocv_back_up/amphG1_loocv_imputed_sp.csv")
write_csv(imputed_ind, "Data_collated/AmphG_loocv_back_up/amphG1_loocv_imputed_ind.csv")

imputed_var <- read_csv("Data_collated/AmphG_loocv_back_up/amphG1_loocv_imputed_var.csv")
imputed_sp <- read_csv("Data_collated/AmphG_loocv_back_up/amphG1_loocv_imputed_sp.csv")
imputed_ind <- read_csv( "Data_collated/AmphG_loocv_back_up/amphG1_loocv_imputed_ind.csv")

## p squared prediction coefficient Guenard et al 2013.
## P2 = 1 when all predictions perfectly match the observations, whereas values below 1 indicate imperfect predictions. Values P2 close to 0 
## (negative or positive) indicate that predictions have poor accuracy, being no better thanwhat would be expected from chance alone.
##http://adn.biol.umontreal.ca/~numericalecology/Reprints/Guenard_et_al_MEE_2013.pdf
# sp_mean_sq_error <- colMeans((Mean_imp_sp[,-1] - Actual[,-1])^2,na.rm=TRUE)
# Sam_var <- Actual[,-1] %>% summarise_all(~ var(.x, na.rm = TRUE))
# 1 - sp_mean_sq_error/Sam_var

p_squared <- 1 - colMeans((imputed_sp[,-1] - imputed_ind[,-1])^2, na.rm=TRUE)/ imputed_ind[,-1] %>% summarise_all(~ var(.x, na.rm = TRUE)) 
# BLlog_std Hab_bread_IU_std count_std
# 1 0.7759939        0.5458554 0.4813179

# Data frame of full imputation of missing data with variance terms 
trait_p_var <- as.data.frame(p1$anc_var) 
colnames(trait_p_var) <- paste(colnames(trait_p_var), "var", sep = "_")

trait_data_test <- p1$anc_recon %>% 
  data.frame() %>% 
  rownames_to_column(var = "species") %>% 
  tibble() %>% 
  filter(species %in% imputed_sp$species) %>% 
  arrange(species)

glimpse(imputed_var)
glimpse(imputed_sp)
glimpse(imputed_ind)


# trait group 2 =====
# impute variables that are binomial occupancy of habitats 
traits2 <- trait_data %>% 
  select(-c(BLlog_std, Hab_bread_IU_std, count_std)) %>% 
  # na.omit() %>% 
  distinct()

# identical(for_imputation$Forest.specialist, for_imputation$specialist)

# multiple rows per species 
multi <- traits2 %>% 
  group_by(species) %>% 
  count() %>% 
  filter(n > 1) %>% 
  pull(species)

# remove species with multiple rows
traits2 <- traits2 %>% 
  filter(!species %in% multi)

# list of species in trait group
names2 <- traits2$species

tree2 <- keep.tip(tree,names2)

length(unique(names2))
length(names2)


# check
length(traits2$species)
length(unique(traits2$species))
tree2

# impute
p2 <-  phylopars(trait_data = traits2, tree = tree,  phylo_correlated = TRUE)
## make empty matrices
imputed_var <- imputed_ind <- imputed_sp <- traits2
imputed_var[,-1] <- imputed_sp[,-1] <- NA
## run LOO CV
save <- c( seq(1000, 9000, by = 1000), nrow(traits2))
for(i in 1:nrow(traits2))
{
  print(i)
  for(j in 1:(ncol(traits2)-1))
  {
    print(j)
    temp_trait_data <- traits2
    temp_trait_data[i,j+1] <- NA
    p_temp <- phylopars(trait_data = temp_trait_data,tree = tree, phylocov_fixed = 
                          p2$pars$phylocov)
    
    imputed_sp[i,j+1] <- left_join(select(temp_trait_data, species), 
                                   rownames_to_column(as.data.frame(p_temp$anc_recon)), 
                                   by = c("species" = "rowname"))[i,j+1]
    imputed_var[i,j+1] <- left_join(select(temp_trait_data, species), 
                                    rownames_to_column(as.data.frame(p_temp$anc_var)), 
                                    by = c("species" = "rowname"))[i,j+1]
  }
  if(i %in% save) {
    write_csv(imputed_var, paste0("Data_collated/AmphG_loocv_back_up/amphG2_loocv_sav",i,"_var.csv"))
    write_csv(imputed_sp, paste0("Data_collated/AmphG_loocv_back_up/amphG2_loocv_sav",i,"_sp.csv"))
  }
}
# load loo_cv imputed traits 
imputed_var <- read_csv( "Data_collated/AmphG_loocv_back_up/amphG2_loocv_sav6639_var.csv")
imputed_sp <- read_csv( "Data_collated/AmphG_loocv_back_up/amphG2_loocv_sav6639_sp.csv")

write_csv(imputed_var, "Data_collated/AmphG_loocv_back_up/amphG2_loocv_imputed_var.csv")
write_csv(imputed_sp, "Data_collated/AmphG_loocv_back_up/amphG2_loocv_imputed_sp.csv")
write_csv(imputed_ind, "Data_collated/AmphG_loocv_back_up/amphG2_loocv_imputed_ind.csv")

imputed_var <- read_csv("Data_collated/AmphG_loocv_back_up/amphG2_loocv_imputed_var.csv")
imputed_sp <- read_csv("Data_collated/AmphG_loocv_back_up/amphG2_loocv_imputed_sp.csv")
imputed_ind <- read_csv( "Data_collated/AmphG_loocv_back_up/amphG2_loocv_imputed_ind.csv")


glimpse(imputed_var)
glimpse(imputed_sp)
glimpse(imputed_ind)

# round binomial imputed traits
imputed_sp <- imputed_sp %>% 
  mutate(across(.cols = c(2:ncol(imputed_sp)), .fns = function(x) if_else(x >= 0.5, 1, 0)))
# reduce dataset to only rows where real values area available 
# find NA in trait dataset 
trait_na_count <- sapply(imputed_ind, FUN = function(y) sum(length(which(is.na(y))))) %>% 
  data.frame() 
# remove non known values to true 
imputed_ind <- imputed_ind %>% 
  na.omit() 

# reduce to only real values 
imputed_var <- imputed_var %>% filter(species %in% imputed_ind$species)
imputed_sp <- imputed_sp %>% filter(species %in% imputed_ind$species)

write_csv(imputed_var, "Data_collated/amphG2_loocv_imputed_var.csv")
write_csv(imputed_sp, "Data_collated/amphG2_loocv_imputed_sp.csv")
write_csv(imputed_ind, "Data_collated/amphG2_loocv_imputed_ind.csv")

imputed_var <- read_csv( "Data_collated/amphG2_loocv_imputed_var.csv")
imputed_sp <- read_csv( "Data_collated/amphG2_loocv_imputed_sp.csv")
imputed_ind <- read_csv("Data_collated/amphG2_loocv_imputed_ind.csv")

# Directly compare the correct imputation to real known data.
# Imputed data
glimpse(imputed_sp)
# check variables bounded by 0 : 1
range(imputed_sp$Forest)
range(imputed_sp$Wetland)

# real data
glimpse(imputed_ind)

check <- imputed_sp == imputed_ind 
check_count <- check %>% data.frame() %>% sapply(FUN = function(y) sum(length(which(y == TRUE)))/nrow(check)) %>% 
  data.frame()

# proportion of correct predictions for 
# species                1.0000000
# Forest                 0.8354922
# Savanna                0.9150629
# Wetland                0.8776832
# Rocky.areas            0.9591044
# Caves.and.subterranean 0.9801999






# Bird Trait data ############## 
# use Bird_trait_imputaion.R script for reference
getwd()
trait_data <- read_csv("Data_collated/Bird_traits_for_impuation3.csv")

glimpse(trait_data)
# load phylo tree 
tree <- read.nexus("Phylogenies/Birds/100birdMCC.nex")
# tree <- keep.tip(tree, names)

## impute variables that are standardised around 0 ## GROUP 1 ################################################-
traits1 <- trait_data %>% 
  select(species,range_logstd,Body_mass_log_g , Generation_length_log_d ,  Litter_clutch_size_log, count_logstd ) %>% 
  # na.omit() %>% 
  distinct()
  
# multiple rows per species 
multi <- traits1 %>% 
  group_by(species) %>% 
  count() %>% 
  filter(n > 1) %>% 
  pull(species)

# remove species with multiple rows
traits1 <- traits1 %>% 
  filter(!species %in% multi)

# check
length(traits1$species)
length(unique(traits1$species))

# list of names in trait database
names1 <- traits1$species

# subset tree to trait groups
tree1 <- keep.tip(tree, names1)

length(unique(traits1$species))

# imputation
p1 <- phylopars(trait_data = traits1, tree = tree, phylo_correlated = TRUE)

## make empty matrices
imputed_var <- imputed_ind <- imputed_sp <- traits1
imputed_var[,-1] <- imputed_sp[,-1] <- NA
save <- c( seq(1000, 9000, by = 1000), nrow(trait1))
## run LOO CV
for(i in 1:nrow(traits1))
{
  print(i)
  for(j in 1:(ncol(traits1)-1))
  {
    print(j)
    temp_trait_data <- traits1
    temp_trait_data[i,j+1] <- NA
    p_temp <- phylopars(trait_data = temp_trait_data,tree = tree, phylocov_fixed = 
                          p1$pars$phylocov)
    
    imputed_sp[i,j+1] <- left_join(select(temp_trait_data, species), 
                                   rownames_to_column(as.data.frame(p_temp$anc_recon)), 
                                   by = c("species" = "rowname"))[i,j+1]
    imputed_var[i,j+1] <- left_join(select(temp_trait_data, species), 
                                    rownames_to_column(as.data.frame(p_temp$anc_var)), 
                                    by = c("species" = "rowname"))[i,j+1]
  }
  if(i %in% save){
    write_csv(imputed_var, paste0("Data_collated/birdG1_loocv_sav",i,"_var.csv"))
    write_csv(imputed_sp, paste0("Data_collated/birdG1_loocv_sav",i,"_sp.csv"))
  }
}

write_csv(imputed_var, "Data_collated/birdG1_loocv_imputed_var.csv")
write_csv(imputed_sp, "Data_collated/birdG1_loocv_imputed_sp.csv")
write_csv(imputed_ind, "Data_collated/birdG1_loocv_imputed_ind.csv")


# second imputation check with larger set of traits due to better coverage from only using species with range data
imputed_var <- read_csv("Data_collated/birdG1_loocv2_imputed_var.csv")
imputed_sp <- read_csv("Data_collated/birdG1_loocv2_imputed_sp.csv")
imputed_ind <- read_csv("Data_collated/birdG1_loocv2_imputed_ind.csv")


## Summarise error
# MAE absolute error 
sp_mean_abs_error <- colMeans(imputed_sp[,-1] - imputed_ind[,-1],na.rm=TRUE) %>% 
  data.frame()
# Body_mass_log_g         0.0007045750
# Generation_length_log_d 0.0008963551
# count_logstd            0.0107057038
# Litter_clutch_size_log  0.0055348277


# MPSE MEAN SQUARE PREDICTED ERROR
# MPSE as sum(((true - imputed)^2)/n_spp))
MSPE <- colMeans((imputed_sp[,-1] - imputed_ind[,-1])^2, na.rm = TRUE) 

## p squared prediction coefficient Guenard et al 2013.
## P2 = 1 when all predictions perfectly match the observations, whereas values below 1 indicate imperfect predictions. Values P2 close to 0 
## (negative or positive) indicate that predictions have poor accuracy, being no better thanwhat would be expected from chance alone.
##http://adn.biol.umontreal.ca/~numericalecology/Reprints/Guenard_et_al_MEE_2013.pdf
# sp_mean_sq_error <- colMeans((Mean_imp_sp[,-1] - Actual[,-1])^2,na.rm=TRUE)
# Sam_var <- Actual[,-1] %>% summarise_all(~ var(.x, na.rm = TRUE))
# 1 - sp_mean_sq_error/Sam_var

# mean squared prediction error 
# P^2 = 1 - (mean(predicted_yi - yi)^2) / (true variance) 
# # P^2 = 1 - MPSE / (true variance) 
Sam_var <- imputed_ind[,-1] %>% summarise_all(~ var(.x, na.rm = TRUE))
p_squared <- 1 - MSPE/Sam_var
# Body_mass_log_g   Generation_length_log_d   count_logstd    Litter_clutch_size_log
#   0.938827               0.9454201           0.05348601                0.55835

p_squared <- 1 - colMeans((imputed_sp[,-1] - imputed_ind[,-1])^2, na.rm=TRUE)/Sam_var 


# Data frame of full imputation of missing data with variance terms 
trait_p_var <- as.data.frame(p1$anc_var) 
colnames(trait_p_var) <- paste(colnames(trait_p_var), "var", sep = "_")

trait_data_test <- p1$anc_recon %>% 
  data.frame() %>% 
  rownames_to_column(var = "species") %>% 
  tibble() %>% 
  filter(species %in% imputed_sp$species) %>% 
  arrange(species)

glimpse(imputed_var)
glimpse(imputed_sp)
glimpse(imputed_ind)




# trait group 2 ====
# impute variables that are binomial occupancy of habitats 
traits2 <- trait_data %>% 
  mutate(Forest.specialist = as.numeric(Forest.specialist)) %>%  # saving errorin script. Output documents are fine but some changes are missing 
  select(species, Forest:Artificial, Forest.specialist) %>% 
  na.omit() %>% 
  distinct()

identical(for_imputation$Forest.specialist, for_imputation$specialist)

# multiple rows per species 
multi <- traits2 %>% 
  group_by(species) %>% 
  count() %>% 
  filter(n > 1) %>% 
  pull(species)

# remove species with multiple rows
traits2 <- traits2 %>% 
  filter(!species %in% multi)

# list of species in trait group
names2 <- traits2$species

tree2 <- keep.tip(tree,names2)

length(unique(names2))
length(names2)


# check
length(traits2$species)
length(unique(traits2$species))
tree2
unique(traits2$Forest.specialist)
unique(traits2$specialist)

# impute
p2 <-  phylopars(trait_data = traits2, tree = tree,  phylo_correlated = TRUE)
## make empty matrices
imputed_var <- imputed_ind <- imputed_sp <- traits2
imputed_var[,-1] <- imputed_sp[,-1] <- NA
## run LOO CV
save <- c( seq(1000, 9000, by = 1000), nrow(traits2))
for(i in 1:nrow(traits2))
{
  print(i)
  for(j in 1:(ncol(traits2)-1))
  {
    print(j)
    temp_trait_data <- traits2
    temp_trait_data[i,j+1] <- NA
    p_temp <- phylopars(trait_data = temp_trait_data,tree = tree, phylocov_fixed = 
                          p2$pars$phylocov)
    
    imputed_sp[i,j+1] <- left_join(select(temp_trait_data, species), 
                                   rownames_to_column(as.data.frame(p_temp$anc_recon)), 
                                   by = c("species" = "rowname"))[i,j+1]
    imputed_var[i,j+1] <- left_join(select(temp_trait_data, species), 
                                    rownames_to_column(as.data.frame(p_temp$anc_var)), 
                                    by = c("species" = "rowname"))[i,j+1]
  }
  if(i %in% save) {
    write_csv(imputed_var, paste0("Data_collated/birdG2_loocv_sav",i,"_var.csv"))
    write_csv(imputed_sp, paste0("Data_collated/birdG2_loocv_sav",i,"_sp.csv"))
  }
}

write_csv(imputed_var, paste0("Data_collated/birdG2_loocv2_var.csv"))
write_csv(imputed_sp, paste0("Data_collated/birdG2_loocv2_sp.csv"))


# load loo_cv imputed traits 
imputed_var <- read_csv( "Data_collated/BirdG2_loocv_back_up/birdG2_loocv2_sav8371_var.csv")
imputed_sp <- read_csv( "Data_collated/BirdG2_loocv_back_up/birdG2_loocv2_sav8371_sp.csv")

glimpse(imputed_var)
glimpse(imputed_sp)
glimpse(imputed_ind)

# round binomial imputed traits
imputed_sp <- imputed_sp %>% 
  mutate(across(.cols = c(2:ncol(imputed_sp)), .fns = function(x) round(x,digits = 0)))
# reduce dataset to only rows where real values area available 
identical(imputed_sp$species, traits2$species)
# find NA in trait dataset 
trait_na_count <- sapply(traits2, FUN = function(y) sum(length(which(is.na(y))))) %>% 
  data.frame() # none in this dataset

write_csv(imputed_var, "Data_collated/birdG2_loocv_imputed_var.csv")
write_csv(imputed_sp, "Data_collated/birdG2_loocv_imputed_sp.csv")
write_csv(imputed_ind, "Data_collated/birdG2_loocv_imputed_ind.csv")

imputed_var <- read_csv( "Data_collated/birdG2_loocv_imputed_var.csv")
imputed_sp <- read_csv( "Data_collated/birdG2_loocv_imputed_sp.csv")
imputed_ind <- read_csv("Data_collated/birdG2_loocv_imputed_ind.csv")

glimpse(imputed_var)
glimpse(imputed_sp)
glimpse(imputed_ind)


# manual check of number of correct predictions for binomial traits
check <- imputed_sp == imputed_ind 
check_count <- check %>% data.frame() %>% lapply(FUN = function(y) sum(length(which(y == TRUE)))/nrow(check)) %>% 
 bind_rows() %>% 
  gather()

write_csv(check_count, "Data_collated/BirdG2_trait_loocv_accuray.csv")

# Proportion of correctly imputed binomial values 
#
# species                1.0000000
# Forest                 0.8678772
# Grassland              0.8383706
# Wetland                0.8380122
# Rocky.areas            0.9388365
# Caves.and.subterranean 0.9948632
# Marine                 0.9800502
# Carnivore              0.9499462
# Herbivore              0.9504241

# Fish Trait data -----
trait_data <- read_csv("Data_collated/fish_traits_for_imputation.csv")
trait_data <- read_csv("Data_collated/Fish_traits_for_imputation2.csv")

names <- trait_data$species
head(names)

# check for duplicate names 
trait_data %>% 
  group_by(species) %>% 
  count() %>% 
  filter(n > 1)

length(unique(trait_data$species))
length(trait_data$species)

# phylogentic tree 
tree <- read.tree("Phylogenies/Fish/actinopt_12k_raxml.tre.xz")
tree <- keep.tip(tree, names)

traits1 <- trait_data %>% 
  select(species, Length) %>% 
  mutate(Length = scale(log(Length), center = TRUE, scale = TRUE)[,1])

traits1 %>% 
  summarise(max = max(Length, na.rm = T),
            min = min(Length, na.rm = T))

## Get cov matrix
p <- phylopars(trait_data = traits1, tree = tree)
## make empty matrices
imputed_var <- imputed_ind <- imputed_sp <- traits1
imputed_var[,-1] <- imputed_sp[,-1] <- NA
# save points
save <- c(seq(1000, 8000, by = 2000),nrow(traits1))
## run LOO CV
for(i in 1:nrow(traits1))
{
  print(i)
  for(j in 1:(ncol(traits1)-1))
  {
    print(j)
    temp_trait_data <- traits1
    temp_trait_data[i,j+1] <- NA
    p_temp <- phylopars(trait_data = temp_trait_data,tree = tree, phylocov_fixed = 
                          p$pars$phylocov)
    # create dataset of just imputed values when left out
    imputed_sp[i,j+1] <- left_join(select(temp_trait_data, species), 
                                   rownames_to_column(as.data.frame(p_temp$anc_recon)), 
                                   by = c("species" = "rowname"))[i,j+1]
    imputed_var[i,j+1] <- left_join(select(temp_trait_data, species), 
                                    rownames_to_column(as.data.frame(p_temp$anc_var)), 
                                    by = c("species" = "rowname"))[i,j+1]
  }
  if(i %in% save) {
    # variants of imputed datapoints
    write_csv(imputed_var, paste0("Data_collated/fishG1_loocv_",i,"imputed_var.csv"))
    # imputed data points using loo 
    write_csv(imputed_sp, paste0("Data_collated/fishG1_loocv_",i,"imputed_sp.csv"))
  }
}
write_csv(imputed_ind, "Data_collated/fishG1_loocv_imputed_ind.csv")

imputed_sp <- read_csv("Data_collated/fishG1_loocv_8131imputed_sp.csv")
imputed_var <- read_csv("Data_collated/fishG1_loocv_8131imputed_var.csv")
imputed_ind <- read_csv("Data_collated/fishG1_loocv_imputed_ind.csv")

glimpse(imputed_sp)
glimpse(imputed_var)
glimpse(imputed_ind)

## p squared prediction coefficient Guenard et al 2013.
## P2 = 1 when all predictions perfectly match the observations, whereas values below 1 indicate imperfect predictions. Values P2 close to 0 
## (negative or positive) indicate that predictions have poor accuracy, being no better than what would be expected from chance alone.
#ical#http://adn.biol.umontreal.ca/~numerecology/Reprints/Guenard_et_al_MEE_2013.pdf
# sp_mean_sq_error <- colMeans((Mean_imp_sp[,-1] - Actual[,-1])^2,na.rm=TRUE)
# Sam_var <- Actual[,-1] %>% summarise_all(~ var(.x, na.rm = TRUE))
# 1 - sp_mean_sq_error/Sam_var

# mean squared prediction error 
# P^2 = 1 - (sum(yi - predicted_yi)^2) / (true variance) 
# # P^2 = 1 - MSE mean squared error / (true variance) 
# 
MPSE <- colMeans((imputed_sp[,-1] - imputed_ind[,-1])^2,na.rm=TRUE) 
Sam_var <- imputed_ind[,-1] %>% summarise_all(~ var(.x, na.rm=TRUE))
p_squared <- 1 - MPSE/Sam_var
# Length  0.7720206

# long form
p_squared <- 1 - colMeans((imputed_sp[,-1] - imputed_ind[,-1])^2, na.rm=TRUE)/var(imputed_ind[,-1], na.rm = TRUE) %>% 
  diag() %>% data.frame() 


# Binomial traits of habitat
traits2 <- trait_data %>% 
  select(-Length)

## Get cov matrix
p <- phylopars(trait_data = traits2, tree = tree)
## make empty matrices
imputed_var <- imputed_ind <- imputed_sp <- traits2
imputed_var[,-1] <- imputed_sp[,-1] <- NA
# save points
save <- c(seq(1000, 8000, by = 2000),nrow(traits2))
## run LOO CV
for(i in 1:nrow(traits2))
{
  print(i)
  for(j in 1:(ncol(traits2)-1))
  {
    print(j)
    temp_trait_data <- traits2
    temp_trait_data[i,j+1] <- NA
    p_temp <- phylopars(trait_data = temp_trait_data,tree = tree, phylocov_fixed = 
                          p$pars$phylocov)
    # create dataset of just imputed values when left out
    imputed_sp[i,j+1] <- left_join(select(temp_trait_data, species), 
                                   rownames_to_column(as.data.frame(p_temp$anc_recon)), 
                                   by = c("species" = "rowname"))[i,j+1]
    imputed_var[i,j+1] <- left_join(select(temp_trait_data, species), 
                                    rownames_to_column(as.data.frame(p_temp$anc_var)), 
                                    by = c("species" = "rowname"))[i,j+1]
  }
  if(i %in% save) {
    # variants of imputed datapoints
    write_csv(imputed_var, paste0("Data_collated/fishG2_loocv_",i,"imputed_var.csv"))
    # imputed data points using loo 
    write_csv(imputed_sp, paste0("Data_collated/fishG2_loocv_",i,"imputed_sp.csv"))
  }
}
write_csv(imputed_ind, "Data_collated/fishG2_loocv_imputed_ind.csv")

imputed_var <- read_csv("Data_collated/fishG2_loocv_8131imputed_var.csv")
imputed_sp <- read_csv("Data_collated/fishG2_loocv_8131imputed_sp.csv")
imputed_ind <- read_csv("Data_collated/fishG2_loocv_imputed_ind.csv")

glimpse(imputed_var)
glimpse(imputed_sp)
glimpse(imputed_ind)

# round to binomial
imputed_sp <- imputed_sp %>% 
  mutate(across(.cols = c(2:4), .fns = function(x) round(x, digits = 0)))

write_csv(imputed_sp, "Data_collated/fishG2_loocv_imputed_sp.csv")
write_csv(imputed_var, "Data_collated/fishG2_loocv_imputed_var.csv")

# Prediction accuracy 
# remove species where we don't know the real values
missing <- which(is.na(imputed_ind$Fresh))

imputed_sp_check <- imputed_sp[-c(missing),] 
imputed_ind_check <- imputed_ind[-c(missing),] 

# proportion of correct estimates
check <- imputed_sp_check == imputed_ind_check 
check_count <- check %>% data.frame() %>% sapply(FUN = function(y) sum(length(which(y == TRUE)))/nrow(check)) %>%  
  data.frame()
# species    1.0000000
# Fresh      0.9662043
# Brack      0.8682219
# Saltwater  0.9720050


