####  Spectra BRT's  ####
# Goal is to try to be more creative with evaluating how much drivers matter




####  Packages  ####
library(targets)
library(here)
library(gmRi)
library(patchwork)
library(scales)
library(dismo)
library(gbm)
library(tidyverse)
library(tidymodels)
library(parsnip)
library(lightgbm)
library(bonsai)

# Package Conflicts
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")

# # Support functions
# source(here("R/support/sizeSpectra_support.R"))

# Resource Path
res_path <- gmRi::cs_path("res")



# Set a seed for reproducing stochastic elements
set.seed(123)

####  Online Resources  ####

# Good rundown of how BRT's and gradient descent work:
# https://bradleyboehmke.github.io/HOML/gbm.html

# Using the gbm package and tuning information
# https://uc-r.github.io/gbm_regression




####  Load Modeling Data  ####



#### Load from `spectra_lm_data_prep.R`
spectra_results <- read_csv(here::here("data/size_and_spectra_model_data.csv")) %>% mutate(
  survey_area = str_replace_all(survey_area, "-", "_"),
  survey_area = str_replace_all(survey_area, " ", "_"),
  survey_area = factor(survey_area))



####__________________####

####  BRT from Dismo  ####

# Can't even get them to fit without tuning separately
# Solved this by scaling the response so units aren't so small

# Details on diagnostics:
# https://rspatial.org/raster/sdm/9_sdm_brt.html


####_________####
#### Median Length Model  ####


# Subset what we're trying to model
length_model_dat <- spectra_results %>% 
  select(med_len, year, survey_area, landings,sst, zp_large)  %>% 
  drop_na(sst) %>% 
  as.data.frame()


# Pick what one we want to try with
model_data <- length_model_dat %>% 
  select(-year)

# Set splits
splits <- initial_split(model_data, prop = 7/10)

# Prepare testing and training
data_other <- training(splits)
data_test  <- testing(splits)


# Train the Model
gbm_mod <- gbm.step(
  data             = model_data,  
  gbm.x            = c(2:ncol(model_data)), 
  gbm.y            = 1, 
  n.folds          = 10,
  cv.folds         = 5,
  family           = "gaussian", 
  n.trees          = 200, 
  max.trees        = 15000,
  tree.complexity  = 5,
  step.size        = 100,
  learning.rate    = 0.001,
  bag.fraction     = 1,
  tolerance.method = "auto")



# Summary variable importance
gbm_mod$contributions %>% 
  mutate(var =  fct_reorder(var, rel.inf)) %>% 
  arrange(desc(var)) %>% 
  ggplot(aes(y = var, x = rel.inf)) +
  geom_col() +
  scale_x_continuous(labels = label_number(suffix = "%")) +
  labs(y = "Independent Variable", x = "Relative Information", 
       title = "Gulf of Maine Spectra - BRT Variable Importance")




# Whats the fitted etc
results_df <- data.frame(
  "year" = length_model_dat$year,
  "survey_area" = model_data$survey_area,
  "actual" = model_data$med_len,
  "fitted" = predict(gbm_mod, model_data))  %>%
  mutate(
    # set_type = ifelse(year %in% data_other$Year, "Training", "Testing"),
    residual = fitted - actual,
    resid_direction = ifelse(residual < 0, "darkred", "royalblue"))




# Plot the Fitted and Actual Values
ggplot(results_df) +
  geom_line(aes(year, actual, color = "Actual"), linewidth = 1) +
  geom_line(aes(year, fitted, color = "BRT Predicted"), linewidth = 1) +
  scale_color_gmri() + 
  facet_wrap(~survey_area)



# What are the Residuals
ggplot(results_df) +
  geom_hline(yintercept = 0, color = "gray70") +
  geom_segment(aes(x = year, xend = year, y = 0, yend = residual, color= I(resid_direction))) +
  geom_point(aes(year, residual, color= I(resid_direction))) +
  labs(x = "Year", y = "BRT Residuals", color = "Testing/Training",
       title = "Median Length BRT Residuals")+ 
  facet_wrap(~survey_area)


# Partial Dependence Plots
gbm_mod %>%
  partial(
    pred.var = "sst", 
    n.trees = gbm_mod$n.trees, 
    grid.resolution = 100) %>%
  autoplot(
    rug = TRUE, 
    train = data_other)



####_________####

#### Spectra Intercept Model  ####

spectra_model_dat <- spectra_results %>% 
  select(b, year, survey_area,  landings,sst, zp_large) %>% 
  drop_na(sst) %>% 
  as.data.frame()


# Pick what one we want to try with
model_data <- spectra_model_dat

# Set splits
splits <- initial_split(model_data, prop = 7/10)

# Prepare testing and training
data_other <- training(splits)
data_test  <- testing(splits)


# Train the Model
gbm_mod <- gbm.step(
  data             = model_data,  
  gbm.x            = c(2:ncol(model_data)), 
  gbm.y            = 1, 
  n.folds          = 10,
  cv.folds         = 5,
  family           = "gaussian", 
  n.trees          = 200, 
  max.trees        = 15000,
  tree.complexity  = 5,
  step.size        = 100,
  learning.rate    = 0.001,
  bag.fraction     = 1,
  tolerance.method = "auto")



# Summary variable importance
gbm_mod$contributions %>% 
  mutate(var =  fct_reorder(var, rel.inf)) %>% 
  arrange(desc(var)) %>% 
  ggplot(aes(y = var, x = rel.inf)) +
  geom_col() +
  scale_x_continuous(labels = label_number(suffix = "%")) +
  labs(y = "Independent Variable", x = "Relative Information", 
       title = "Gulf of Maine Spectra - BRT Variable Importance")




# Whats the fitted etc
results_df <- data.frame(
  "year" = model_data$year,
  "survey_area" = model_data$survey_area,
  "actual" = model_data$b,
  "fitted" = predict(gbm_mod, model_data))  %>%
  mutate(
    # set_type = ifelse(year %in% data_other$Year, "Training", "Testing"),
    residual = fitted - actual,
    resid_direction = ifelse(residual < 0, "darkred", "royalblue"))




# Plot the Fitted and Actual Values
ggplot(results_df) +
  geom_line(aes(year, actual, color = "Actual"), linewidth = 1) +
  geom_line(aes(year, fitted, color = "BRT Predicted"), linewidth = 1) +
  scale_color_gmri() + 
  facet_wrap(~survey_area)



# What are the Residuals
ggplot(results_df) +
  geom_hline(yintercept = 0, color = "gray70") +
  geom_segment(aes(x = year, xend = year, y = 0, yend = residual, color= I(resid_direction))) +
  geom_point(aes(year, residual, color= I(resid_direction))) +
  labs(x = "Year", y = "BRT Residuals", color = "Testing/Training",
       title = "Spectra Slope BRT Residuals")+ 
  facet_wrap(~survey_area)





####_________####

####  Regional Models  ####



####  GB BRT  ####


# Set region
region  <- "GB"


# Prepare Starting Data
model_dat <- region_indices %>%
  filter(survey_area == region) %>%
  dplyr::select(Year, b) %>%
  mutate(Year = as.numeric(Year)) %>%
  left_join(drivers_scaled, by = join_by(Year == year)) %>%
  drop_na()

# Scale the response 
region_mu <- mean(model_dat$b)
region_sd <- sd(model_dat$b)
model_dat$b <- scale(model_dat$b)




# Train the Model
gb_gbm <- gbm.step(
  data             = model_dat,  
  gbm.x            = c(3:ncol(model_dat)), 
  gbm.y            = 2, 
  n.folds          = 10,
  family           = "gaussian", 
  n.trees          = 200, 
  max.trees        = 20000,
  tree.complexity  = 5,
  step.size        = 100,
  learning.rate    = 0.00015,
  bag.fraction     = 0.8,
  tolerance.method = "auto"
)

# Assign so we can recycle code
region_gbm <- gb_gbm

# Summary variable importance
region_gbm$contributions %>% 
  mutate(var =  fct_reorder(var, rel.inf)) %>% 
  arrange(desc(var)) %>% 
  ggplot(aes(y = var, x = rel.inf)) +
  geom_col() +
  scale_x_continuous(labels = label_number(suffix = "%")) +
  labs(y = "Predictor", x = "Relative Information", 
       title = "Georges Bank Spectra - BRT Variable Importance")




# Whats the fitted etc
results_df <- data.frame(
  "year" = model_dat$Year,
  "actual" = (model_dat$b * region_sd) + region_mu,
  "fitted" = (predict(region_gbm, model_dat) * region_sd) + region_mu) %>%
  mutate(
    set_type = ifelse(year %in% data_other$Year, "Training", "Testing"),
    residual = fitted - actual,
    resid_direction = ifelse(residual < 0, "darkred", "royalblue"))



# Plot the Fitted and Actual Values
ggplot(results_df) +
  geom_line(aes(year, actual, color = "Actual"), linewidth = 1) +
  geom_line(aes(year, fitted, color = "BRT Predicted"), linewidth = 1) +
  scale_color_gmri()



# What are the Residuals
ggplot(results_df) +
  geom_hline(yintercept = 0, color = "gray70") +
  geom_segment(aes(x = year, xend = year, y = 0, yend = residual, color= I(resid_direction))) +
  geom_point(aes(year, residual, color= I(resid_direction))) +
  labs(x = "Year", y = "BRT Residuals", color = "Testing/Training", 
       title = "Georges Bank Spectra BRT Residuals")


# # What about the different variables curves
# gb_gbm$gbm.call$predictor.names
# gbm::plot.gbm(gb_gbm, i.var = 1:2)



####  SNE BRT  ####


# Set region
region  <- "SNE"


# Prepare Starting Data
model_dat <- region_indices %>%
  filter(survey_area == region) %>%
  dplyr::select(Year, b) %>%
  mutate(Year = as.numeric(Year)) %>%
  left_join(drivers_scaled, by = join_by(Year == year)) %>%
  drop_na() %>%
  dplyr::select(-ends_with("slope"))

# Scale the response 
region_mu <- mean(model_dat$b)
region_sd <- sd(model_dat$b)
model_dat$b <- scale(model_dat$b)

# Does SNE even have a trend - slightly...
lm(b ~ Year, data = model_dat)
ggplot(model_dat, aes(Year, b)) +
  geom_line() +
  geom_smooth(method = "lm")


# Train the Model
sne_gbm <- gbm.step(
  # data             = data_other,  
  # gbm.x            = c(3:ncol(data_other)), 
  data             = model_dat,  
  gbm.x            = c(3:ncol(model_dat)), 
  gbm.y            = 2, 
  n.folds          = 5,
  family           = "gaussian", 
  n.trees          = 10, 
  max.trees        = 1000,
  tree.complexity  = 5,
  step.size        = 1,
  learning.rate    = 0.0001,
  bag.fraction     = .8)


# Assign so we can recycle code
region_gbm <- sne_gbm


# Summary variable importance
region_gbm$contributions %>% 
  mutate(var =  fct_reorder(var, rel.inf)) %>% 
  arrange(desc(var)) %>% 
  ggplot(aes(y = var, x = rel.inf)) +
  geom_col() +
  scale_x_continuous(labels = label_percent(scale = .01)) +
  labs(y = "Predictor", x = "Relative Information", 
       title = "Southern New England Spectra - BRT Variable Importance")




# Whats the fitted etc
results_df <- data.frame(
  "year" = model_dat$Year,
  "actual" = (model_dat$b * region_sd) + region_mu,
  "fitted" = (predict(region_gbm, model_dat) * region_sd) + region_mu) %>%
  mutate(
    #set_type = ifelse(year %in% data_other$Year, "Training", "Testing"),
    residual = fitted - actual,
    resid_direction = ifelse(residual < 0, "darkred", "royalblue"))




# Plot the Fitted and Actual Values
ggplot(results_df) +
  geom_line(aes(year, actual, color = "Actual"), linewidth = 1) +
  geom_line(aes(year, fitted, color = "BRT Predicted"), linewidth = 1) +
  scale_color_gmri()



# What are the Residuals
ggplot(results_df) +
  geom_hline(yintercept = 0, color = "gray70") +
  geom_segment(aes(x = year, xend = year, y = 0, yend = residual, color= I(resid_direction))) +
  geom_point(aes(year, residual, color= I(resid_direction))) +
  labs(x = "Year", y = "BRT Residuals", color = "Testing/Training", 
       title = "Southern New England Spectra BRT Residuals")







####  MAB BRT  ####


# Set region
region  <- "MAB"


# Prepare Starting Data
model_dat <- region_indices %>%
  filter(survey_area == region) %>%
  dplyr::select(Year, b) %>%
  mutate(Year = as.numeric(Year)) %>%
  left_join(drivers_scaled, by = join_by(Year == year)) %>%
  drop_na()

# Scale the response 
region_mu <- mean(model_dat$b)
region_sd <- sd(model_dat$b)
model_dat$b <- scale(model_dat$b)



# Train the Model
mab_gbm <- gbm.step(
  data             = model_dat,  
  gbm.x            = c(3:ncol(model_dat)), 
  # data             = data_other,  
  # gbm.x            = c(3:ncol(data_other)), 
  gbm.y            = 2, 
  n.folds          = 10,
  family           = "gaussian", 
  n.trees          = 10, 
  max.trees        = 18000,
  tree.complexity  = 5,
  step.size        = 1,
  learning.rate    = 0.0001,
  bag.fraction     = 0.9,
  tolerance.method = "auto")

# What things 

# Assign so we can recycle code
region_gbm <- mab_gbm

# Plot the predictor curves
gbm.plot(gbm.object = region_gbm, variable.no = 2, write.title = FALSE, plot.layout = c(1,1))
gbm.plot(gbm.object = region_gbm, variable.no = 5, write.title = FALSE, plot.layout = c(1,1))
gbm.plot(gbm.object = region_gbm, variable.no = 7, write.title = FALSE, plot.layout = c(1,1))
gbm.plot(gbm.object = region_gbm, variable.no = 8, write.title = FALSE, plot.layout = c(1,1))
gbm.plot(gbm.object = region_gbm, variable.no = 9, write.title = FALSE, plot.layout = c(1,1))
gbm.plot(gbm.object = region_gbm, variable.no = 10, write.title = FALSE, plot.layout = c(1,1))
gbm.plot(gbm.object = region_gbm, variable.no = 11, write.title = FALSE, plot.layout = c(1,1))

# Summary variable importance
region_gbm$contributions %>% 
  mutate(var =  fct_reorder(var, rel.inf)) %>% 
  arrange(desc(var)) %>% 
  ggplot(aes(y = var, x = rel.inf)) +
  geom_col() +
  scale_x_continuous(labels = label_number(suffix = "%")) +
  labs(y = "Predictor", x = "Relative Information", 
       title = "Mid-Atlantic Bight Spectra - BRT Variable Importance")




# Whats the fitted etc
results_df <- data.frame(
  "year" = model_dat$Year,
  "actual" = (model_dat$b * region_sd) + region_mu,
  "fitted" = (predict(region_gbm, model_dat) * region_sd) + region_mu)  %>% 
  mutate(
    residual = fitted - actual,
    resid_direction = ifelse(residual > 0, "royalblue", "darkred")
  )




# Plot the Fitted and Actual Values
ggplot(results_df) +
  geom_line(aes(year, actual, color = "Actual"), linewidth = 1) +
  geom_line(aes(year, fitted, color = "BRT Predicted"), linewidth = 1) +
  scale_color_gmri()



# What are the Residuals
ggplot(results_df) +
  geom_hline(yintercept = 0, color = "gray70") +
  geom_segment(aes(x = year, xend = year, y = 0, yend = residual, color= I(resid_direction))) +
  geom_point(aes(year, residual, color= I(resid_direction))) +
  scale_size_manual(values = c(2, 1)) +
  scale_alpha_manual(values = c(1, 0.4)) +
  labs(x = "Year", y = "BRT Residuals", color = "Testing/Training", 
       title = "Mid-Atlantic Bight Spectra BRT Residuals")







####_______________####
####  Doesn't Work - Tidymodel BRT  ####
####  Boosted Regression Trees - Tidymodels  ####

# Can BRT's resolve the changes in spectra
# better than segmented regressions



# Set seed
set.seed(123)


# Pick the Region
region = "GoM"


# Prepare Starting Data
model_dat <- region_indices %>%
  filter(survey_area == region) %>%
  dplyr::select(Year, b) %>%
  mutate(Year = as.numeric(Year)) %>%
  left_join(drivers_scaled, by = join_by(Year == year)) %>%
  drop_na() 

corrplot::corrplot(cor(model_dat))

# Scale the response 
region_mu <- mean(model_dat$b)
region_sd <- sd(model_dat$b)
model_dat$b <- scale(model_dat$b)


# Prepare Splits
splits     <- initial_split(model_dat, prop = 4/5)
data_other <- training(splits)
data_test  <- testing(splits)


# Prepare the validation data
val_set <- validation_split(data_other, prop = 0.80)

# Or prepare an n-fold cross-validation set
folds <- vfold_cv(model_dat, v = 10)



####  Build the Model  ####
brt_mod <- boost_tree(trees = tune())  %>%
  set_mode(mode = "regression") %>% 
  set_engine("lightgbm", seed = 63233)


####  Create the Recipe  ####
# this is where you'd do feature engineering: scaling, dummy coding, dates
mod_recipe <- recipe(b ~ ., data = data_other) %>% 
  step_lag(is.numeric, lag = 2:3)


####  Create the Workflow  ####

mod_workflow <- workflow() %>% 
  add_model(brt_mod) %>% 
  add_recipe(mod_recipe) 


####  Create Tuning Grid  ####

# BRT is one main argument "trees"
tune_args(brt_mod)

# We can use a regular grid to pick values for testing
tree_grid <- grid_regular(
  trees(range = c(50,2000)), 
  levels = 30)



####  Fit/Train the model  ####




# Fit it with tuning
mod_res <- mod_workflow %>% 
  tune_grid(
    val_set,
    grid = tree_grid,
    control = control_grid(save_pred = T))



# Plot some results from tuning
# Doesn't look like its learning...
mod_res %>% 
  collect_metrics() %>% 
  filter(.metric == "rmse") %>% 
  ggplot(aes(x = trees, y = mean)) + 
  geom_point() + 
  geom_line() + 
  ylab("rmse") 





####_______________####

####  GAM  ####

# More coefficients than data
 library(mgcv)

# Pick the Region
region = "GoM"


# Prepare Starting Data
model_dat <- region_indices %>%
  filter(survey_area == region) %>%
  dplyr::select(Year, b) %>%
  mutate(Year = as.numeric(Year)) %>%
  left_join(drivers_scaled, by = join_by(Year == year)) %>%
  drop_na() 


ggplot(model_dat, aes(Gulf_of_Maine_landings,b)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(model_dat, aes(All_sst_anom, b)) +
  geom_point() +
  geom_smooth(method = "lm")


# Global model
gom_gam <- gam(
  formula = b ~ s(All_sst_anom, k = 4) + s(Gulf_of_Maine_zp_small, k = 4) +
    s(Gulf_of_Maine_zp_large, k = 4) + s(Gulf_of_Maine_landings, k = 4) + s(All_gulf_stream_index, k = 4),
  data = model_dat,
  method = "REML")



# Pick the Region
region = "GB"


# Prepare Starting Data
model_dat <- region_indices %>%
  filter(survey_area == region) %>%
  dplyr::select(Year, b) %>%
  mutate(Year = as.numeric(Year)) %>%
  left_join(drivers_scaled, by = join_by(Year == year)) %>%
  drop_na() 




ggplot(model_dat, aes(Georges_Bank_landings,b)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(model_dat, aes(All_sst_anom, b)) +
  geom_point() +
  geom_smooth(method = "lm")

drivers_scaled %>% 
  pivot_longer(ends_with("landings"), names_to = "region", values_to = "landings") %>% 
  ggplot(aes(year, landings)) +
  geom_line() +
  facet_wrap(~region) +
  geom_vline(xintercept = 1982, linetype = 2)



# Extract



# Global model
gb_gam <- gam(
  formula = b ~ s(All_sst_anom, k = 4) + s(Georges_Bank_zp_small, k = 4) +
    s(Georges_Bank_zp_large, k = 4) + s(Georges_Bank_landings, k = 4) + s(All_gulf_stream_index, k = 4),
  data = model_dat,
  method = "REML")





# Pick the Region
region = "MAB"


# Prepare Starting Data
model_dat <- region_indices %>%
  filter(survey_area == region) %>%
  dplyr::select(Year, b) %>%
  mutate(Year = as.numeric(Year)) %>%
  left_join(drivers_scaled, by = join_by(Year == year)) %>%
  drop_na()  %>% 
  rename_all(~str_replace(.x, "-", "_"))


# Global model
mab_gam <- gam(
  formula = b ~ 
    s(All_sst_anom, k = 4) + 
    s(`Mid_Atlantic_Bight_zp_small`, k = 4) +
    s(`Mid_Atlantic_Bight_zp_large`, k = 4) + 
    s(`Mid_Atlantic_Bight_landings`, k = 4) + 
    s(All_gulf_stream_index, k = 4),
  data = model_dat,
  method = "REML")

# Pick the Region
region = "SNE"


# Prepare Starting Data
model_dat <- region_indices %>%
  filter(survey_area == region) %>%
  dplyr::select(Year, b) %>%
  mutate(Year = as.numeric(Year)) %>%
  left_join(drivers_scaled, by = join_by(Year == year)) %>%
  drop_na()  %>% 
  rename_all(~str_replace(.x, "-", "_"))


# Global model
sne_gam <- gam(
  formula = b ~ 
    s(All_sst_anom, k = 4) + 
    s(Mid_Atlantic_Bight_zp_small, k = 4) +
    s(Mid_Atlantic_Bight_zp_large, k = 4) + 
    s(Southern_New_England_landings, k = 4) + 
    s(All_gulf_stream_index, k = 4),
  data = model_dat,
  method = "REML")






####  GAM Dredge  ####
