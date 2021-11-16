#### Size at Age Testing
#### 11/11/2021



####  Packages  ####
library(FSAdata)
library(FSA)
library(car)
library(tidyverse)
library(gmRi)





#####  Load Data  ####

# Access Biological Data with {gmRi}
nmfs_bio <- gmRi::gmri_survdat_prep(survdat_source = "bio")

# Pull Cod
cod <- filter(nmfs_bio, comname == "atlantic cod")




####  EDA  ####
cod %>% ggplot(aes(length_cm, age)) +
  geom_point(alpha = 0.1, color = "gray60") +
  facet_grid(season~catchsex)


# pull relevant columns for length-age and weight-age relationship
cod_prep <- cod %>% 
  select(tow_id = id, est_towdate, survey_area, comname, 
         indid, length_cm, indwt, age, sex, maturity)


#### VBert Example  ####

# Check the model structure/code
( vb <- vbFuns(param="Typical") )

# Get reasonable starting points:
cod_starts <- vbStarts(length_cm ~ age, data = cod_prep)


# Use nls to estimate VBGF parameters using starting points
cod_fit <- nls(length_cm ~ vb(age, Linf, K, t0), data = cod_prep, start = cod_starts)

# Access parameters with coef()
coef(cod_fit)


# Use boot() to get bootstrapped confidence intervals
cod_boot1 <- Boot(cod_fit)  # Be patient! Be aware of some non-convergence
# confint(cod_boot1)



####  Plotting Predicted Values  ####

# Get prediction and confidence intervals over range of ages
ages <- seq(-1,12,by=0.2)
preds1 <- data.frame(ages,
                     predict(cod_fit, data.frame(age = ages)),
                     confint(cod_boot1))
names(preds1) <- c("age","fit","LCI","UCI")



# Restrict predictions to only ages found in the data
agesum <- group_by(cod_prep, sex) %>%
  summarize(minage = min(age, na.rm = T),
            maxage = max(age, na.rm = T))

# filter
preds2 <- filter(preds1, 
                 age >= min(agesum$minage), 
                 age <= max(agesum$maxage))


# Plot the overall age-length relationship
ggplot() + 
  #geom_ribbon(data= filter(preds2, LCI > 2), aes(x=age, ymin=LCI, ymax=UCI), fill="gray40") +
  geom_point(data = cod_prep, aes(y = length_cm, x = age), size=2, alpha=0.1) +
  geom_line(data=preds1,aes(y = fit, x = age), size = 1, linetype = 3, color = "gray70") +
  geom_line(data=preds2, aes(x=age, y = fit), size = 1.5, linetype = 1, color = gmri_cols("gmri blue")) +
  labs(x = "age", y = "Total Length (cm)", subtitle = "Atlantic Cod")






