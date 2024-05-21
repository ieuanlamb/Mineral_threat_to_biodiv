# models with one chain check

library(brms)
library(tidyverse)

getwd()
set.seed(2)

# read models 
chain1 <- readRDS("models/Bird_trait_m3_stan_chain1.rds")
chain2 <- readRDS("models/Bird_trait_m3_stan_chain2.rds")
chain3 <- readRDS("models/Bird_trait_m3_stan_chain3.rds")
chain4 <- readRDS("models/Bird_trait_m3_stan_chain4.rds")
chain15 <- readRDS("models/Bird_trait_m3_stan_chain15.rds")
chain16 <- readRDS("models/Bird_trait_m3_stan_chain16.rds")

# summary outputs of each chain 
summary(chain1)
summary(chain2)
summary(chain3)
summary(chain4)
summary(chain15)
summary(chain16)


chains_comb <- combine_models(chain1,chain2,chain3,chain4, chain16, chain16, check_data = TRUE)

# diagnostic checks 
mcmc_plot(chains_comb, type = 'trace')
mcmc_plot(chains_comb, type = 'rhat')
mcmc_plot(chains_comb,type = "acf_bar")
mcmc_plot(chains_comb,type = "neff")
mcmc_plot(chains_comb,type = "areas")

summary(chains_comb)

# save combined model 
write_rds(chains_comb, file = "models/Bird_trait_m3_combinedchains.rds")


# bird model m3 with reduced number of species but more traits =====
m2_c1 <- readRDS("models/bird_reduced_stan_chain1.rds")
m2_c2 <- readRDS("models/bird_reduced_stan_chain2.rds")
m2_c3 <- readRDS("models/bird_reduced_stan_chain3.rds")
m2_c4 <- readRDS("models/bird_reduced_stan_chain4.rds")

# combine all individual chains into one model 
m2_comb <- combine_models(m2_c1,m2_c2,m2_c3,m2_c4, check_data = TRUE)

# diagnostic checks 
mcmc_plot(m2_comb, type = 'trace')
mcmc_plot(m2_comb, type = 'rhat')
mcmc_plot(m2_comb,type = "acf_bar")
mcmc_plot(m2_comb,type = "neff")
mcmc_plot(m2_comb,type = "areas")

# save combined model 
write_rds(m2_comb, file = "models/Bird_reduced_stan_combinedchains.rds")


# bird model m2 without correlation with reduced number of species but more traits =====
m2_c1 <- readRDS("models/bird_reduced_stan_alt_chain1.rds")
m2_c2 <- readRDS("models/bird_reduced_stan_alt_chain2.rds")
m2_c3 <- readRDS("models/bird_reduced_stan_alt_chain3.rds")
m2_c4 <- readRDS("models/bird_reduced_stan_alt_chain4.rds")

# combine all individual chains into one model 
m2_comb <- combine_models(m2_c1,m2_c2,m2_c3,m2_c4, check_data = TRUE)

# diagnostic checks 
mcmc_plot(m2_comb, type = 'trace')
mcmc_plot(m2_comb, type = 'rhat')
mcmc_plot(m2_comb,type = "acf_bar")
mcmc_plot(m2_comb,type = "neff")
mcmc_plot(m2_comb,type = "areas")

# save combined model 
write_rds(m2_comb, file = "models/Bird_reduced_stan_alt_combinedchains.rds")


