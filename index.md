# Github Pages Navigation

 1. [NEFSC size spectrum build EDA](https://adamkemberling.github.io/nefsc_trawl/R/qaqc_reports/03_nefsc_eda.html)
 
A look at the overall data trends following "normal" data prep for the nefsc bottom trawl data. Includes overall abundance and biomass as well as the catch composition and stratified area abundances.
 
 2. [Aggregation Level validation](https://adamkemberling.github.io/nefsc_trawl/R/qaqc_reports/02_survdat_stratification_validation.html)

Report confirming the validity of getting area stratified abundance and biomass at an individual or size level, and then aggregating by year or season later. Basically, can area stratified catch rates be applied for numbers at length and will they compare well to overall abundance and biomass.

 3. [Identifying Possible Issues with Fall Haddock Numbers](https://adamkemberling.github.io/nefsc_trawl/R/qaqc_reports/Haddock_check.html)

When digging into the data for Haddock it appeared that the biomass for haddock in the Henry Bigelow Era, as derived from length-weight relationships was distinctly higher than for in the spring as well as during the Albatross years. This report dug into whether the relationship was tied to a specific survdat pull.


 4. [Size spectrum slope sensitivity testing](https://adamkemberling.github.io/nefsc_trawl/R/qaqc_reports/09_ss_sensitivity.html)

When calculating size spectrum slopes the question arose of whether to truncate the data by minimum body size to account for size-selectivity in the survey gear. This report looked into the sensitivity of results to adjustments to the size cutoff threshold.

 5. [SURVDAT Pull Diagnostics](https://adamkemberling.github.io/nefsc_trawl/R/qaqc_reports/survdat_pull_check.html)
 
Side by Side comparison of all the survdat pulls we have in inventory. Aim was to identify whether or not the abundance and biomass conversions accounting for the survey vessel switch were being applied as intended and across different data sets.
 
 
 6. [Comparing Albatross to Bigelow Conversions using Haddock Data](https://adamkemberling.github.io/nefsc_trawl/R/qaqc_reports/albatross_bigelow_conversions.html)
 
This report seeks to check the accuracy of the spring and fall conversions of haddock biomass from the nefsc groundfish survey data. This really identifies the inconsistent application of conversion factors for fall haddock.
 
 7. [Abundance and Biomass Retrospective Check](https://adamkemberling.github.io/nefsc_trawl/R/qaqc_reports/survdat_ab_bio_check.html)
 
 This report seeks to identify for which species of which we have published literature there are notable differences in abundance and biomass from using the newest survdat pull which has corrected conversion issues.