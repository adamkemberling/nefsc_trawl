####
# Cohort Von Bert Testing
# Following Code setup from "hierarchical_size_at_age.R" ~ Grant Adams UW
# Objective: JAGS VBGF + Hierarchical Effect for things like Year/Cohort


####  Packages  ####
library(targets)
library(gmRi)
#library(patchwork)
library(tidyverse)
library(here)
library(runjags)
library(coda)
library(bayesplot)
library(tidybayes)


####  Data  ####

# Load the list of data for each of the species
tar_load("vonbert_species_bio")

# Set one as guinea pig that has good data
haddock <- vonbert_species_bio$haddock %>% 
  mutate(
    decade = floor_decade(est_year),
    age = as.numeric(age),
    year = factor(est_year),
    cohort_year = yearclass) %>% 
  # filter(cohort_year >= 1970) %>% 
  mutate(yearclass = factor(yearclass),
         yearclass = fct_drop(yearclass)) %>% 
  drop_na(length_cm, age)

# Plot age/length data
ggplot(haddock, aes(age, length_cm, color = yearclass)) +
  geom_point(alpha = 0.3) + 
  facet_wrap(~decade, ncol = 1) +
  theme(legend.position = "right") +
  labs(x = "age", y = "length (cm)", subtitle = "Haddock", color = "Cohort")

# Plot Cohort Abundances & Number of Ages Present in Data
yearclass_counts <- haddock %>% count(yearclass) 
cohort_ages <- haddock %>% count(yearclass, age)
ggplot(cohort_ages, aes(x = age, y = yearclass)) + 
  geom_text(data = yearclass_counts, 
            aes(y = yearclass, x = -1.5, 
                label = str_pad(str_c("N = ", n), width = 10, side = "right", pad = " ")), 
            size = 2.25) +
  geom_point(aes(size = n)) + 
  theme_classic() + 
  theme(legend.position = "bottom") +
  labs(x = "Age", y = "Cohort Year")

# Check that lengths are good
# range(haddock$length_cm)

#Save the data so Grant can test
# saveRDS(haddock, file = here::here("data/single_species/haddock_bio.rds"))


# Assign Data to List 
dat <- list(
  age = haddock$age,                      # Age
  length = haddock$length_cm,             # Length
  group = haddock$yearclass,              # Group Hierarchy
  N = length(haddock$age),                # Number of Fish
  G = length(unique(haddock$yearclass))   # Number of year classes
)



####  JAGS Model  ####
jags.mod <-
  "model{

    # Likelihood Statement / Model
    for(i in 1:N){
      length[i] ~ dnorm(y.hat[i], tau.y)
      y.hat[i] = Linf[group[i]] * (1-exp(-k[group[i]] * (age[i] - t0[group[i]] )))
  }
 
  # SD
  tau.y = pow(sigma, -2)
  sigma ~ dunif(0,100) 
   
  # Level-2 parameters for each group (location, cohort, etc.)
  for(j in 1:G){
    # T(L, R) used to truncate left and right limits
    Linf[j]  ~ dnorm(mu.Linf, tau.Linf)
    k[j]     ~ dnorm(mu.k, tau.k) 
    t0[j]    ~ dnorm(mu.t0, tau.t0) 
  }
   
  # Get hyperparameters on untransformed scale (if transformed)
  mu.Linf = exp(log.mu.Linf) 
  mu.k    = exp(log.mu.k)
  # mu.t0

  # Priors for level-2 parameters
  log.mu.Linf ~ dnorm(0, 0.0001)
  log.mu.k    ~ dnorm(0, 0.0001)
  mu.t0       ~ dnorm(0, 0.0001)
   
  # Precision
  tau.Linf = pow(sig.Linf, -2)
  tau.k    = pow(sig.k, -2)
  tau.t0   = pow(sig.t0, -2)
   
  # SD of parameters
  sig.Linf ~ dunif(0,10)
  sig.k    ~ dunif(0,10)
  sig.t0   ~ dunif(0,10)
}"

# write model to a text file
writeLines(jags.mod, here::here("R/size_at_age/yearclass_jags_model.txt"))



##### PARAMETERS TO MONITOR #####
params = c("Linf", 
           "k", 
           "t0", 
           "log.mu.Linf", 
           "log.mu.k", 
           "mu.Linf", 
           "mu.k", 
           "mu.t0", 
           "sig.Linf",
           "sig.k",
           "sig.t0",
           "sigma" 
)

##### MCMC DIMENSIONS #####
ni = 1500
nb = 4000
na = 1000
nt = 10
nc = 3
n.iter = ni + nb

##### RUN THE MODEL IN JAGS #####
runJagsOut <- run.jags(model = here::here("R/size_at_age/yearclass_jags_model.txt") , 
                       monitor = params ,
                       data = dat ,
                       n.chains = nc ,
                       adapt = na ,
                       burnin = nb ,
                       sample = ni ,
                       thin = nt ,
                       summarise = FALSE ,
                       plots = FALSE )

# Converged?
# plot(runJagsOut)

# Summary
jags_summ <- summary(runJagsOut)

# Clean up using {tidybayes}
out_tidy <- tidy_draws(runJagsOut)

# Pulling a key for index:yearclass pairs
yearclass_id_key <- data.frame(yearclass = sort(unique(dat$group))) %>% 
  mutate(yearclass_num = row_number(),
         yearclass_num = str_c("[", yearclass_num, "]"))

# Pulling Linf
linf_post <- out_tidy %>% select(starts_with("Linf")) %>%
  rename_all(.funs = ~ str_replace(.x, "Linf", "")) %>% 
  pivot_longer(cols = everything(), names_to = "yearclass_num", values_to = "Linf") %>% 
  left_join(yearclass_id_key, by = "yearclass_num")

ggplot(linf_post, aes(x = Linf, y = fct_rev(yearclass))) +
  geom_boxplot() + 
  labs(x = "Posterior Distribution: Linf (cm)", y = "Cohort Year")

# Pulling K
K_post <- out_tidy %>% select(starts_with("k")) %>%
  rename_all(.funs = ~ str_replace(.x, "k", "")) %>% 
  pivot_longer(cols = everything(), names_to = "yearclass_num", values_to = "K") %>% 
  left_join(yearclass_id_key, by = "yearclass_num")

ggplot(K_post, aes(x = K, y = fct_rev(yearclass))) +
  geom_boxplot() + 
  labs(y = "Posterior Distribution: K")


# Pulling t0
t0_post <- out_tidy %>% select(starts_with("t0")) %>%
  rename_all(.funs = ~ str_replace(.x, "t0", "")) %>% 
  pivot_longer(cols = everything(), names_to = "yearclass_num", values_to = "t0") %>% 
  left_join(yearclass_id_key, by = "yearclass_num")

ggplot(t0_post, aes(x = t0, y = fct_rev(yearclass))) +
  geom_boxplot() + labs(y = "Posterior Distribution: t0 ")



####  Distribution Plotting  ####
bayesplot::mcmc_areas(runJagsOut$mcmc, pars = str_c("Linf[", c(1:51), "]"))
