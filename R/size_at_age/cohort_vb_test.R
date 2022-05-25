# This was removed from size at age EDA when it became clear that data might not be available to do this:



# Cohort Based Von Bertalannffy
"
Fish born as part of the same cohort are likely to exhibit similar patterns in
growth due to their exposure to similar environments at the same points in their development.

For this reason and others it may be meaningful to capture that cohort effect 
when tracking changes in growth coefficients.

One way to do this is with mixed effects models.
"

# Grab one species to use as tester
test_species <- species_data$`winter flounder` %>% 
  mutate(
    age = as.numeric(age),
    year = factor(est_year),
    yearclass = factor(yearclass))


# Formula for nlme model
vbert_form <- ~ Linf * (1 - exp(-K * (t - t0)))

# Attempt to include random effect for cohort with "yearclass"
#vbert_form <-  ~ Linf * (1 - exp(-K * (t - t0))) + (1 | yearclass)


# use derive to construct function
nfun <- deriv(expr = vbert_form,
              namevec = c("Linf", "K", "t0"),
              function.arg = c("Linf", "K", "t0"))

# Set starting points vector
test_starts <-  vbStarts(length_cm ~ as.numeric(age), data = test_species)

# Play around until it runs
nlmer(length_cm ~ nfun(age, Linf, K, t0)  ~ (1 | yearclass),
      data = test_species, 
      start = test_starts)


