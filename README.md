# InformativeCensoring

[![Build Status](https://travis-ci.org/scientific-computing-solutions/InformativeCensoring.svg?branch=master)](https://travis-ci.org/scientific-computing-solutions/InformativeCensoring)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/InformativeCensoring)](http://cran.r-project.org/package=InformativeCensoring)


Multiple Imputation for Informative Censoring

This R package implement two methods for multiple imputation of survival data.
* Gamma imputation from Jackson et al. [1]
* Risk score imputation from Hsu et al. [2]

## Contributor (alphabetical order)
Barlett, Jonathan; Burkoff, Nikolas; Jackson, Dan; Jones, Edmund; Law, Martin; 
Metcalfe, Paul; Ruau, David;

## Gamma imputation (Jackson 2014)
The Gamma imputation method implementation was developed in collaboration 
between AstraZeneca, the MRC Biostatistics Unit and the University of Cambridge.

This implementation was validated to the best of our effort following good coding
practice and thorough user testing.

## Risk Score Imputation (Hsu 2009)
We implemented the method described in Chiu-Hsieh Hsu and Jeremy Taylor (2009)
following the publication.

This implementation was validated to the best of our effort following good coding
practice and thorough user testing.

## Installation

To install the development version from github
```R
# install.packages("devtools")
devtools::install_github("hadley/ggplot2")
```

[1] Dan Jackson, Ian White, Shaun Seaman, Hannah Evans, Kathy Baisley, and James Carpenter. Relaxing the independent censoring assumption in the Cox proportional hazards model using multiple imputation. Statistics in Medicine, 33(27):4681–4694, 2014.

[2] Chiu-Hsieh Hsu and Jeremy MG Taylor. Nonparametric comparison of two survival functions with dependent censoring via nonparametric multiple imputation. Statistics in Medicine, 28(3):462–475, 2009.
