# Github Pages Navigation

---

## QA/QC Reports

The following markdown documents were made to verify various aspects of the research workflow. These were made to validate integral steps like the validity of length-weight relationships, the application of area-stratification, and the investigation of patterns that stand out within these reports.

**These include:**   


 1. [NEFSC size spectrum build EDA](https://adamkemberling.github.io/nefsc_trawl/R/qaqc_reports/nefsc_eda.html)
 
A look at the overall data trends following "normal" data prep for the nefsc bottom trawl data. Includes overall abundance and biomass as well as the catch composition and stratified area abundances.

 2. [Aggregation Level Validation & LW Relationship Fits](https://adamkemberling.github.io/nefsc_trawl/R/qaqc_reports/stratification_validation.html)

Report confirming the validity of getting area stratified abundance and biomass at an individual or size level, and then aggregating by year or season later. 

Basically, can area stratified catch rates be applied for numbers at length, and will they compare well to the aggregate abundance and biomass of all lengths caught when area-stratified.

End of report details how aggregate biomass compares to length-weight relationships from wigley 2006 paper and fishbase coefficients. Concludes with timeline that omits species that don't meet cutoff.

 3. [Identifying Possible Issues with Fall Haddock Numbers](https://adamkemberling.github.io/nefsc_trawl/R/qaqc_reports/Haddock_check.html)

When digging into the data for Haddock it appeared that the biomass for haddock in the Henry Bigelow Era, as derived from length-weight relationships was distinctly higher than for in the spring as well as during the Albatross years. This report dug into whether the relationship was tied to a specific survdat pull or not.

The results from this markdown indicated that there was indeed an incorrect application of catch adjustment for some species, which was then corrected resulting in a new survdat pull.

 4. [SURVDAT Pull Diagnostics](https://adamkemberling.github.io/nefsc_trawl/R/qaqc_reports/survdat_pull_check.html)
 
Because on the investigation into problems seen in the haddock data it was decided to go back and do comparisons across different data pulls.

This side-by-side comparison of all the survdat pulls we have in inventory was done to identify whether or not the abundance and biomass conversions accounting for the survey vessel switch were being applied as intended and across different data sets.


 5. [Comparing Albatross to Bigelow Conversions using Haddock Data](https://adamkemberling.github.io/nefsc_trawl/R/qaqc_reports/albatross_bigelow_conversions.html)
 
This report seeks to check the accuracy of the spring and fall conversions of haddock biomass from the nefsc groundfish survey data. This really identifies the inconsistent application of conversion factors for fall haddock.

---

##  Summary of Retrospective Problems Across Different SURVDAT Resources

[2019 Abundance and Biomass Retrospective Check](https://adamkemberling.github.io/nefsc_trawl/R/qaqc_reports/abundance_biomass_check_2019.html)

&

[2020 Abundance and Biomass Retrospective Check](https://adamkemberling.github.io/nefsc_trawl/R/qaqc_reports/abundance_biomass_check_2020.html) 
 
 These two reports seeks to identify for which species of which we have published literature there are notable differences in abundance and biomass from using the newest survdat pull which has corrected conversion issues.



 ---
 
## Core Size Spectrum Analysis


[Size Spectrum Methods Comparison](https://adamkemberling.github.io/nefsc_trawl/R/qaqc_reports/sizespectra_methods_comparison.html)

This document was done to directly compare the results of calculating size-spectrum characteristics using different methods. 

 
[Size Spectrum Results and Biomass Allocation Across Body  Size](https://adamkemberling.github.io/nefsc_trawl/R/nmfs_size_spectra/bodymass_allocation.html)

This report digs into what the outcomes of the size spectrum estimations look like, and how their results break out into different stanzas. There's also timelines of each region showing relative prevalence of differently sized individuals.

---

## Preliminary Exploration - Outdated

The following reports are no longer relevant as either the input data has changed or their premise is no longer relevant.

[Size spectrum slope sensitivity testing](https://adamkemberling.github.io/nefsc_trawl/R/qaqc_reports/ss_sensitivity.html)

When calculating size spectrum slopes the question arose of whether to truncate the data by minimum body size to account for size-selectivity in the survey gear. This report looked into the sensitivity of results to adjustments to the size cutoff threshold.
