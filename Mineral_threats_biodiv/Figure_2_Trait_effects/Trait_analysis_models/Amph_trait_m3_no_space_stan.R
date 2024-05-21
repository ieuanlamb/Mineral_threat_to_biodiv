# Amphibian trait model 3 to run on the HPC
# more traits ie habitat traits with bernoulli 

library(stringr); library(dplyr);library(readr);library(ggplot2);library(tidyr)
library(brms)

getwd()
# setwd("X:/edwards_lab1/User/bop21ipl")

# Data Load ====
data <- read_csv("Data/Amph_trait_model/amphibian_final_data.csv")
amph_dist <- readRDS("Data/Amph_trait_model/Distance_matrix/Amphibian_Corr_matrix_SymPosDef.rds")
amph_phylodist <- readRDS("Data/Amph_trait_model/Distance_matrix/Amphibian_centroid_Correlationmatrix.rds")

# check name matching 
identical((data$species), colnames(amph_phylodist))
setdiff((data$species), colnames(amph_phylodist))
setdiff(colnames(amph_phylodist),(data$species))
identical((data$species), colnames(amph_dist))
setdiff((data$species), colnames(amph_dist))
identical(rownames(amph_dist), colnames(amph_dist))
identical(rownames(amph_phylodist), colnames(amph_phylodist))

# model =====

m3 <- brm(mine_thrnd  ~ BLlog_std  * Rangelogkm2_std  + Forest + Savanna + Wetland + Rocky.areas + Caves.and.subterranean + 
            (1|gr(species, cov = A)), 
          family = bernoulli,
          prior = c(prior(normal(0,1), class = "Intercept"),
                    prior(normal(0,0.5), class = "b"),
                    prior(normal(0,1), class = "sd")),
          data = data, 
          data2 = list(A = amph_phylodist),
          file = "models/Amph_trait_m3_no_space_stan.rds",
          output_dir = "logs/cmdstanr",
          output_basename = "Amph_trait_m3_nospace_stan",
          # control = list(adapt_delta = 0.9),
          backend = "cmdstanr", chains = 4, cores = 4, threads = threading(2), thin = 1,
          warmup = 1000, iter = 2000,
          seed = 9)

