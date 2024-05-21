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


# m2 <- brm(mine_thrnd  ~ Body_mass_log_st * Range_log_st +
#             (1|gr(species, cov = A)) +
#             (1|gr(species_space, cov = B)),
#           family = bernoulli,
#           prior = c(prior(normal(0,1), class = "Intercept"),
#                     prior(normal(0,0.5), class = "b"),
#                     prior(normal(0,1), class = "sd")),
#           data = data,
#           data2 = list(A = rep_phylodist2, B = rep_Dist),
#           file = "models/Rep_trait_full_m2.rds",
#           # control = list(adapt_delta = 0.9),
#           backend = "cmdstanr", chains = 4, cores = 4, threads = threading(2), thin = 1,
#           warmup = 1000, iter = 3000,
#           seed = 9)
# 
# 
# m2


# LOO Compare spatial element of amphibian full model -----
# rep1 <- read_rds("models/Rep_trait_m2.rds") 

# rep_upd <- update(rep1, 
#                    formula. = .  ~ Body_mass_log_st * Range_log_st + (1|gr(species, cov = A)),
#                    newdata = data,
#                    data2 = list(A = rep_phylodist2)
# )
# saveRDS(rep_upd, "models/Rep_trait_m2_noSpace.rds")

# full model may not have run yet so run no_spatial model in its intierety instead of updating the full model
rep_upd <- brm(mine_thrnd  ~ Body_mass_log_st * Range_log_st +
                (1|gr(species, cov = A)),
              family = bernoulli,
              prior = c(prior(normal(0,1), class = "Intercept"),
                        prior(normal(0,0.5), class = "b"),
                        prior(normal(0,1), class = "sd")),
              data = data,
              data2 = list(A = rep_phylodist2),
              file = "models/Rep_trait_m2_noSpace.rds",
              # control = list(adapt_delta = 0.9),
              backend = "cmdstanr", chains = 4, cores = 4, threads = threading(2), thin = 1,
              warmup = 1000, iter = 3000,
              seed = 9)


# if full model has not run then this may cause error >>>>
rep1 <- read_rds("models/Rep_trait_full_m2.rds")
rep_upd <- read_rds("models/Rep_trait_m2_noSpace.rds")

rep1 <- add_criterion(rep1, "loo") 

rep_upd <- add_criterion(rep_upd, "loo") 

loo_compare(rep1, rep_upd, criterion = "loo")
# elpd_diff se_diff
# rep1       0.0       0.0 
# rep_upd -183.3      19.2 

# mcmc_plot(m1,type = "trace")
# mcmc_plot(m1,type = "rhat")
# mcmc_plot(m1,type = "acf_bar")
# mcmc_plot(m1,type = "neff")
# mcmc_plot(m1,type = "areas")
