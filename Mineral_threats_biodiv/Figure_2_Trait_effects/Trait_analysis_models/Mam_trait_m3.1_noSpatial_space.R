# Mammal Trait model noSpatial
# using COMBINE phylo tree and imputed traits

# library(tidyverse)
library(stringr); library(dplyr);library(readr);library(ggplot2);library(tidyr)
library(lqmm)
library(brms)

getwd()

# Data load ====

data <- read_csv("Data/Mam_trait_model/Mammal_final_data.csv")
mam_dist <- readRDS("Data/Mam_trait_model/Mammal_centroid_Correlationmatrix.rds")
mam_phylodist <- readRDS("Data/Mam_trait_model/Mammal_phylo_Distmatrix.rds")

# check name matching 
identical((data$species), colnames(mam_phylodist))
identical((data$species), colnames(mam_dist))
identical(rownames(mam_dist), colnames(mam_dist))
identical(rownames(mam_phylodist), colnames(mam_phylodist))
setdiff(colnames(mam_phylodist), (data$species))

# model ====

m3 <- brm(mine_thrnd  ~   adult_mass_g * range_calculated + litter_size_n +
            habitat_breadth_n +
            dphy_plant * dphy_vertebrate * dphy_invertebrate +
            activity_cycle +
            trophic_level +
            foraging_stratum +
            terrestrial_non_volant +
            terrestrial_volant +
            marine +
            freshwater +
            (1|gr(species, cov = A)),
          family = bernoulli,
          prior = c(prior(normal(0,1), class = "Intercept"),
                    prior(normal(0,0.5), class = "b"),
                    prior(normal(0,1), class = "sd")),
          data = data,
          data2 = list(A = mam_phylodist),
          file = "models/Mam_trait_m3_no_space_stan.rds",
          output_dir = "logs/cmdstanr",
          output_basename = "Mam_trait_m3_nospace_stan",
          backend = "cmdstanr",
          chains = 4, cores = 4, threads = threading(2), thin = 1,
          warmup = 1000, iter = 2000,
          # stanvars = stanvars,  # note our `stanvars`
          seed = 9)

# LOO Compare spatial element of amphibian full model -----
# mam1 <- read_rds("models/Mam_trait_m3.1.rds")
# 
# mam_upd <- update(mam1, 
#                   formula. = .  ~  adult_mass_g * range_calculated + litter_size_n +
#                     habitat_breadth_n +
#                     dphy_plant * dphy_vertebrate * dphy_invertebrate +
#                     activity_cycle +
#                     trophic_level +
#                     foraging_stratum +
#                     (1|gr(species, cov = A)),
#                   newdata = data,
#                   data2 = list(A = mam_phylodist),
#                   seed = 1
# )
# 
# saveRDS(mam_upd, "models/Mam_trait_m3.2_no_spatial.rds")
# 
# # mam1 <- add_criterion(mam1, "loo") 
# # 
# # mam_upd <- add_criterion(mam_upd, "loo") 
# # 
# # loo_compare(mam1, mam_upd, criterion = "loo")
# 
# # mcmc_plot(m2, type = "trace")
# # mcmc_plot(m1, type = "areas", prob = 0.95)
# # m1
