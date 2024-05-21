# Reptile trait model script m2 bernoulli model

# Reptile trait model 1 to run on the HPC

library(stringr); library(dplyr);library(readr);library(ggplot2);library(tidyr)
library(brms)
library(lqmm)


getwd()

# Data load ====
data <- read_csv("Data/Rep_trait_model/final_data.csv")
rep_phylodist2 <- readRDS("Data/Rep_trait_model/final_phylodist2.rds")
rep_Dist <- readRDS("Data/Rep_trait_model/final_Dist2.rds")

# checks 
length(data$species)
nrow(rep_phylodist2)
nrow(rep_Dist)
identical((data$species), colnames(rep_Dist))
setdiff((data$species), colnames(rep_Dist))
identical(rownames(rep_Dist), colnames(rep_phylodist2))
identical(rownames(rep_phylodist2),colnames(rep_Dist))

# model ====

m2 <- brm(mine_thrnd  ~ Body_mass_log_st * Range_log_st +
            (1|gr(species, cov = A)) +
            (1|gr(species_space, cov = B)),
          family = bernoulli,
          prior = c(prior(normal(0,1), class = "Intercept"),
                    prior(normal(0,0.5), class = "b"),
                    prior(normal(0,1), class = "sd")),
          data = data,
          data2 = list(A = rep_phylodist2, B = rep_Dist),
          file = "models/Rep_trait_full_m2_stan.rds",
          output_dir = "logs/cmdstanr",
          output_basename = "Rep_trait_full_m2_stan",
          backend = "cmdstanr", chains = 4, cores = 4, threads = threading(2), thin = 1,
          warmup = 1000, iter = 2000,
          seed = 9)


m2
# 
# mcmc_plot(m1,type = "trace")
# mcmc_plot(m1,type = "rhat")
# mcmc_plot(m1,type = "acf_bar")
# mcmc_plot(m1,type = "neff")
# mcmc_plot(m1,type = "areas")
