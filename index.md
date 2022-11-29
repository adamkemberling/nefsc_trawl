# WARMEM Trawl Survey Analyses Navigation

This index page serves as a navigation guide to the various html reports created for the WARMEM project that utilize long-term fisheries independent data. These reports are hosted via github pages. Links to the reports with brief descriptions are provided below.

---

## Survey Data Cleanup Code

[SURVDAT Prep Functions](https://adamkemberling.github.io/nefsc_trawl/01_Survdat_Standard_Cleanup.html)

## Presentations

[WARMEM 2021 Update](https://adamkemberling.github.io/nefsc_trawl/presentations/Northeast_Trawl_Size_Spectrum.html)


## Weight at Length and Survdat Cleanup QA/QC

The following markdown documents were made to verify various aspects of the research workflow. These were made to validate integral steps like the validity of length-weight relationships, the application of area-stratification, and the investigation of patterns that stand out within these reports.


 1. [NEFSC size spectrum build EDA](https://adamkemberling.github.io/nefsc_trawl/R/qaqc_reports/nefsc_eda.html): A look at the overall data trends following "normal" data prep for the NEFSC bottom trawl data. Includes overall abundance and biomass as well as the catch composition and stratified area abundances.
 
 2. [Aggregation Level Validation & LW Relationship Fits](https://adamkemberling.github.io/nefsc_trawl/R/qaqc_reports/stratification_validation.html): A Report confirming the validity of getting area stratified abundance and biomass at an individual or size level, and then aggregating by year or season later. 


 3. [Identifying Possible Issues with Fall Haddock Numbers](https://adamkemberling.github.io/nefsc_trawl/R/qaqc_reports/Haddock_check.html): When digging into the data for Haddock it appeared that the biomass for haddock in the Henry Bigelow Era, as derived from length-weight relationships was distinctly higher than for in the spring as well as during the Albatross years. This report dug into whether the relationship was tied to a specific survdat pull or not. The results from this markdown indicated that there was indeed an incorrect application of catch adjustment for some species, which was then corrected resulting in a new survdat pull.


 4. [SURVDAT Pull Diagnostics](https://adamkemberling.github.io/nefsc_trawl/R/qaqc_reports/survdat_pull_check.html): Because on the investigation into problems seen in the haddock data it was decided to go back and do comparisons across different data pulls. This side-by-side comparison of all the survdat pulls we have in inventory was done to identify whether or not the abundance and biomass conversions accounting for the survey vessel switch were being applied as intended and across different data sets.
 

 5. [Comparing Albatross to Bigelow Conversions using Haddock Data](https://adamkemberling.github.io/nefsc_trawl/R/qaqc_reports/albatross_bigelow_conversions.html): This report seeks to check the accuracy of the spring and fall conversions of haddock biomass from the nefsc groundfish survey data. This really identifies the inconsistent application of conversion factors for fall haddock.
 
 
---
 
# Size Spectrum Analyses


[Size Spectrum Methods Comparison](https://adamkemberling.github.io/nefsc_trawl/R/qaqc_reports/sizespectra_methods_comparison.html)

This document was done to directly compare the results of calculating size-spectrum characteristics using different methods. 

[Biomass Allocation Across Body  Size for the Study Regions](https://adamkemberling.github.io/nefsc_trawl/R/nmfs_size_spectra/bodymass_allocation.html)

This report digs into what the outcomes of the size spectrum estimations look like, and how their results break out into different stanzas. There's also timelines of each region showing relative prevalence of differently sized individuals.


[Species Composition within Size Bins & Bin-Size Gap Handling](https://adamkemberling.github.io/nefsc_trawl/R/nmfs_size_spectra/spectra_composition_suppl.html)

[Patterns Across Entire Survey Area](https://adamkemberling.github.io/nefsc_trawl/R/nmfs_size_spectra/size_spectrum_story.html)

[Size Spectrum Surfaces](https://adamkemberling.github.io/nefsc_trawl/R/nmfs_size_spectra/size_spectra_surfaces.html)

### Manuscript:

[WARMEM Size Spectrum Manuscript](https://adamkemberling.github.io/nefsc_trawl/R/nmfs_size_spectra/warmem_size_spectrum_manu.html)

---

# Size at Age Analyses

[Von-Bert Size at Age relationships](https://adamkemberling.github.io/nefsc_trawl/R/size_at_age/size_at_age_exploration.html)

[Temperature and Depth Preferences EDA using Observable.js](https://adamkemberling.github.io/nefsc_trawl/R/size_at_age/species_env_pref.html)

[Center of Biomass Changes Under Two Thermal Regimes](https://adamkemberling.github.io/nefsc_trawl/R/size_at_age/encounter_temperatures.html)

### Manuscript:

[Decadal Variability in Growth Manuscript](https://adamkemberling.github.io/nefsc_trawl/R/size_at_age/size_at_age_regimes.html)

---

##  Summary of Retrospective Issues Across Different SURVDAT Resources

These two reports seeks to identify for which species (if any), of which we have published literature, where there are notable differences in either the abundance and biomass between different pulls of the "survdat" data. The most recent data from NMFS has corrected the issues we've flagged, which largely stemmed from the absence or mis-specification of abundance and/or biomass corrections due to the survey vessel change.

[2019 Abundance and Biomass Retrospective Check](https://adamkemberling.github.io/nefsc_trawl/R/qaqc_reports/abundance_biomass_check_2019.html)

&

[2020 Abundance and Biomass Retrospective Check](https://adamkemberling.github.io/nefsc_trawl/R/qaqc_reports/abundance_biomass_check_2020.html) 
 
 

---

## Preliminary Exploration - Outdated

The following reports are no longer relevant as either the input data has changed or their premise is no longer relevant.

[Size spectrum slope sensitivity testing](https://adamkemberling.github.io/nefsc_trawl/R/qaqc_reports/ss_sensitivity.html)

When calculating size spectrum slopes the question arose of whether to truncate the data by minimum body size to account for size-selectivity in the survey gear. This report looked into the sensitivity of results to adjustments to the size cutoff threshold.
