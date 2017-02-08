#This file contains the roxygen comments for the package and
#internal data frames


##' Perform methods of multiple imputation for
##' time to event data
##'
##' See Nonparametric comparison of two survival functions with
##' dependent censoring via nonparametric multiple imputation. Hsu and Taylor
##' Statistics in Medicine (2009) 28:462-475 for Hsu's method
##'
##' See Relaxing the independent censoring assumption in the Cox proportional
##' hazards model using multiple imputation. Jackson et al., Statistics in Medicine
##' (2014) 33:4681-4694 for Jackson's method
##'
##' @name InformativeCensoring-package
##' @aliases InformativeCensoring
##' @docType package
##' @title Perform methods of multiple imputation for
##' time to event data
##' @author \email{David.Ruau@@astrazeneca.com}
NULL


##' Simulated time to event data with 5 time independent covariates
##'
##' This dataset is inspired by the simulation described in Hsu and Taylor,
##' Statistics in Medicine (2009) 28:462-475 with an additional DCO.time column
##'
##' @name ScoreInd
##' @docType data
##' @format A data.frame containing a row per subject with eleven columns:
##' @field  Id subject identifier
##' @field arm factor for treatment group control=0, active=1
##' @field Z1 binary time independent covariate
##' @field Z2 continuous time independent covariate
##' @field Z3 binary time independent covariate
##' @field Z4 continuous time independent covariate
##' @field Z5 binary time independent covariate
##' @field event event indicator (1 yes, 0 no)
##' @field time subject censoring/event time (in years)
##' @field to.impute logical, should an event time be imputed for this subject?
##' (this is ignored if subject has event time)
##' @field DCO.time The time the subject would have been censored if they had not
##' had an event or been censored before the data cut off date
NULL

##' Simulated time dependent variables for time to event data
##'
##' This data set contains time dependent covariates for the
##' \code{\link{ScoreInd}} time to event data.
##'
##' @name ScoreTimeDep
##' @docType data
##' @format A data.frame containing 1 row per subject-visit
##' @field Id The Subject Id
##' @field start The covariate given in each row are for a given subject from time 'start'...
##' @field end ... until time end
##' @field W1 The value of a (binary) time dependent variable
##' for the subject with the given 'Id' between times 'start' and 'end'
##' @field W2 The value of a (continuous) time dependent variable
##' for the subject with the given 'Id' between times 'start' and 'end'
NULL
