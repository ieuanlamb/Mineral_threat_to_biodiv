# Bird traits model 2 reduced and not imputed 

library(stringr); library(dplyr);library(readr);library(ggplot2);library(tidyr)
library(brms)
library(lqmm)

rm(list=ls())
gc()

getwd()
start_time <- Sys.time()

# Data load ====

data <- read_csv("Data/Bird_trait_model/Bird_final_data.csv")
bird_dist <- readRDS("Data/Bird_trait_model/Bird_centroid_Correlationmatrix_FULL.rds")
bird_phylodist2 <- readRDS("Data/Bird_trait_model/Bird_phylo_Distmatrix.rds")

#check
identical((data$species), colnames(bird_phylodist2))
identical((data$species), colnames(bird_dist))

length(data$species)
nrow(bird_phylodist2)
nrow(bird_dist) #9022 

# model ====

# reduced datasets
data_t <- data[1:200,]

phylo_t <- bird_phylodist2[1:200,1:200]

space_t <- bird_dist[1:200,1:200]
glimpse(data_t)


# m3_test <- brm(mine_thrnd  ~ range_logstd.x * Body_mass_log_g + Habitat_breadth_IUCN + Generation_length_log_d +
#                  Forest + Savanna + Grassland + Wetland + Rocky.areas + Caves.and.subterranean + Desert  + Marine +
#                  specialist +
#                  count_logstd +
#                  Carnivore + Herbivore,
#                (1|gr(species, cov = A)) +
#                (1|gr(species_space, cov = B)),
#                family = bernoulli,
#                prior = c(prior(normal(0,1), class = "Intercept"),
#                          prior(normal(0,0.5), class = "b"),
#                          prior(normal(0,1), class = "sd")
#                ),
#                data = data_t,
#                data2 = list(A = phylo_t, B = space_t),
#                # file = "models/bird_trait_t200_m3..rds",
#                # control = list(adapt_delta = 0.9),
#                backend = "cmdstanr",
#                chains = 4, cores = 4, threads = threading(2), thin = 1,
#                warmup = 1000, iter = 2000)

# summary(m3_test)


# model 2 with full dataset
m3 <- brm(mine_thrnd  ~ range_logstd.x * Body_mass_log_g + 
            Forest + Savanna + Grassland + Wetland + Rocky.areas + Caves.and.subterranean + Desert  + Marine +
            count_logstd +
            Carnivore + 
            (1|gr(species, cov = A)) +
            (1|gr(species_space, cov = B)),
          family = bernoulli,
          prior = c(prior(normal(0,1), class = "Intercept"),
                    prior(normal(0,0.5), class = "b"),
                    prior(normal(0,1), class = "sd")),
          data = data,
          data2 = list(A = bird_phylodist2, B = bird_dist),
          file = "models/bird_reduced_stan_chain2.rds",
          # control = list(adapt_delta = 0.9),
          # output_dir = "logs/cmdstanr",
          # output_basename = "Bird_m3_reduced",
          # backend = "cmdstanr",
          chains = 1, cores = 1,  #threads = threading(2), thin = 1,
          warmup = 1000, iter = 2000,
          seed = 2)

# mcmc_plot(m2,type = "trace")
# mcmc_plot(m2,type = "rhat")
# mcmc_plot(m2,type = "acf_bar")
# mcmc_plot(m2,type = "neff")
# mcmc_plot(m2,type = "areas")
