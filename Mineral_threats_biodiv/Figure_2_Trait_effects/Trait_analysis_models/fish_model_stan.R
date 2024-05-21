# fish model for stanage

library(tidyverse)
library(brms)



# load data inputs ====
data <- read_csv("Data/Fish_trait_model/final_data.csv")
fish_phylodist2 <- readRDS("Data/Fish_trait_model/final_phylodist2.rds")
fish_dist2 <- readRDS("Data/Fish_trait_model/final_Dist2.rds")


# model run with full dataset ====
m1 <- brm(mine_thrnd  ~ range_size_logst * Length_logstd + Fresh + Brack + Saltwater,
          family = bernoulli,
          prior = c(prior(normal(0,1), class = "Intercept"),
                    prior(normal(0,0.5), class = "b")),
          data = data,
          backend = "cmdstanr",
          chains = 4, cores = 4,  threads = threading(2), thin = 1,
          warmup = 1000, iter = 2000,
          seed = 9)


# # model 1 with full dataset
# m1 <- brm(mine_thrnd  ~ range_size_logst * Length_logstd + Fresh * Brack * Saltwater +
#             (1|gr(species, cov = A)),
#           family = bernoulli,
#           prior = c(prior(normal(0,1), class = "Intercept"),
#                     prior(normal(0,0.5), class = "b"),
#                     prior(normal(0,1), class = "sd")),
#           data = data,
#           data2 = list(A = fish_phylodist2),
#           file = "models/fish_trait_m2_no_space.rds",
#           # control = list(adapt_delta = 0.9),
#           # output_dir = "logs/cmdstanr",
#           # output_basename = "fish_trait_basic",
#           backend = "cmdstanr",
#           chains = 4, cores = 4,  threads = threading(2), thin = 1,
#           warmup = 1000, iter = 2000,
#           seed = 9)

# 
# mcmc_plot(m1,type = "trace")
# mcmc_plot(m1,type = "rhat")
# mcmc_plot(m1,type = "acf_bar")
# mcmc_plot(m1,type = "neff")
# mcmc_plot(m1,type = "areas")


# model 2 with full dataset
m2 <- brm(mine_thrnd  ~  range_size_logst * Length_logstd + Fresh + Brack + Saltwater +
            (1|gr(species, cov = A)) +
            (1|gr(species_space, cov = B)),
          family = bernoulli,
          prior = c(prior(normal(0,1), class = "Intercept"),
                    prior(normal(0,0.5), class = "b"),
                    prior(normal(0,1), class = "sd")),
          data = data,
          data2 = list(A = fish_phylodist2, B = fish_dist2),
          file = "models/fish_trait_m2.rds",
          control = list(adapt_delta = 0.9),
          output_dir = "logs/cmdstanr",
          output_basename = "fish_trait_m2",
          backend = "cmdstanr",
          chains = 4, cores = 4,  threads = threading(2), thin = 1,
          warmup = 1000, iter = 2000,
          seed = 9)
