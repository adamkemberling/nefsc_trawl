#############################################################################
############################# Von Bertalanffy Sheepshead  ###################
#############################################################################
# https://github.com/grantdadams/Spatial-Growth-Models
# Analysis from Adams, G.D., Leaf, R.T., Ballenger, J.C., Arnott, S.A., Mcdonough, C.J., 2018. Spatial variability in the growth of Sheepshead (Archosargus probatocephalus) in the Southeast US : Implications for assessment and management. Fish. Res. 206, 35–43. doi:10.1016/j.fishres.2018.04.023

# clear the workspace
rm(list = ls(all = T))

# load packages
library(rstan)
library(coda)
library(tidyr)
library(dplyr)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores() - 1)

##### READ AND PREPARE DATA #####
haddock <- readRDS(file = here::here("data/single_species/haddock_bio.rds"))

# Set up design matrix
design_mat <- matrix(c(rep(1, nrow(haddock))), nrow = nrow(haddock))

# If you were also relating growth paramters to enviro data
# design_mat <- model.matrix(~Regime + SST + Climate, data = haddock) 

# Put into list to pass to STAN
dataList = list(
  n_i          = nrow(haddock), # Sample size
  n_r          = length(unique(haddock$yearclass)), # Number of year classes
  n_pred       = ncol(design_mat), # Number of predictors
  log_length_i = log(haddock$length_cm), # Vector of log fork length of fish i
  age_i        = haddock$age, # Vector of age of fish i
  design_mat   = design_mat, # n_i * n_pred design matrix of linear predictors
  r_i          = as.numeric(haddock$yearclass) # Integer vector for cohort of fish i
)

plot(y = haddock$length_cm, x = haddock$age)


##### MCMC DIMENSIONS #####
ni = 5000
burn = 100
nChains = 2

##### RUN THE BASE MODEL IN STAN #####
StanFitBase <- stan(file = here::here('R/size_at_age/stan/Stan models/VBGF_Model_1.stan'), 
                    data = dataList, 
                    iter = ni, 
                    chains = nChains, 
                    cores = 2, 
                    verbose = TRUE, 
                    warmup = burn, 
                    control = list(max_treedepth = 14, adapt_delta = 0.9))



##### RUN THE MODEL IN STAN WITH AND RANDOM COHORT EFFECTS #####
StanFitHierarchical <- stan(file = here::here('R/size_at_age/stan/Stan models/VBGF_Model_2.stan'), 
                            data = dataList, 
                            iter = ni, 
                            chains = nChains, 
                            cores = 2, 
                            verbose = TRUE, 
                            warmup = burn, control = list(max_treedepth = 14, adapt_delta = 0.9))

##### SAVE #####
mod_list <- list(StanFitBase, StanFitHierarchical)
file_name <- paste("Stan models/VBGF_Stan_models_",gsub("-","_",as.Date(trunc(Sys.time(),"day"))),".RData", sep = "")
save(mod_list, file = file_name)


