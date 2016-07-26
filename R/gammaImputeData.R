#This file contains the code and object documentation associated with
#generating the imputed data sets usign gamma imputation see validationGamma.R
#for validation routines, and internalGamma.R for subfunctions

##' \code{GammaImputedData} object
##' 
##' An object which contains 
##' 
##' @slot data A data frame containing the time to event data
##' with 3 new columns impute.time and impute.event, the imputed event/censoring times and event indicators
##' (for subjects whose data is not imputed these columns contain the unchanged event/censoring time and
##' event indicator) and internal_gamma_val which is the value of gamma used for each subject in this data set
##' @slot default.formula The default model formula which will be used when fitting the imputed data  
##' @name GammaImputedData.object
NULL


##' \code{GammaImputedSet} object
##' 
##' An object which contains the set of gamma imputed data frames.
##' Use the \code{ExtractSingle} function to extract a single 
##' \code{GammaImputedData} objects. Use the ImputeStat function to fit models
##' to the entire set of imputed data frames 
##' 
##' It contains the following:
##' @slot data A data frame containing the unimputed time to event data (along with a column internal_gamma_val 
##' which is the value of gamma used for each subject in this data set)
##' @slot m The number of imputed data sets
##' @slot gamma.factor The value of gamma.factor used with the imputation
##' @slot impute.time A matrix (1 column per imputed data set) containing the imputed times 
##' @slot impute.event A matrix (1 column per imputed data set) containing the imputed event indicators
##' @slot default.formula The default model formula which will be used when fitting the imputed data   
##' @name GammaImputedSet.object
##' @seealso \code{\link{GammaImputedData.object}}
NULL

##' Perform gamma-Imputation for a given data set
##' 
##' This function performs the Imputation described in 
##' Relaxing the independent censoring assumptions in the Cox proportional 
##' hazards model using multiple imputation. (2014) D. Jackson et al. Statist. Med. (33)
##' 4681-4694
##' 
##' @details See the Gamma Imputation vignette for further details 
##' 
##' @param formula The model formula to be used when fitting the models to calculate
##' the cumulative hazard. A formula for coxph can include strata terms but not 
##' cluster or tt and only right-censored \code{Surv} objects can be used.
##' Note the function does not allow multiple strata to be written as \code{strata(W1)+strata(W2)}, 
##' use \code{strata(W1,W2)} instead  
##' @param data A time to event data set for which event times are to be imputed 
##' @param m The number of imputations to be created
##' @param gamma Either column name containing the value of gamma or a vector of values giving the subject specific 
##' size of the step change in the log hazard at censoring. If a subject has NA in this column then no imputation is performed 
##' for this subject (i.e. the subject's censored time remains unchanged after imputation). If a subject has already had an 
##' event then the value of gamma is ignored. If \code{gamma.factor} is also used then the subject specific gamma
##' are all multipled by \code{gamma.factor}. At least one of \code{gamma} and \code{gamma.factor} must be included.
##' @param gamma.factor If used, a single numeric value. If no \code{gamma} then the step change in log 
##' hazard for all subjects at censoring is given by \code{gamma.factor}. If \code{gamma} is used  
##' then for each subject, the step change in log hazard is given by \code{gamma.factor} multiplied by the subject specific gamma.
##' At least one of \code{gamma} and \code{gamma.factor} must be included.    
##' @param bootstrap.strata The strata argument for stratified bootstrap sampling, see argument \code{strata}
##' for the function \code{boot::boot} for further details. If argument is not used then standard sampling with
##' replacement will be used
##' @param DCO.time Either column name containing the subject's data cutoff time or a vector
##' of DCO.times for the subjects or a single number to be used as the DCO.time for all subjects 
##' (if imputed events are > this DCO.time then subjects are censored at DCO.time in imputed data sets)
##' @param ... Additional parameters to be passed into the model fit function
##' @param parallel The type of parallel operation to be used (if any).
##' @param ncpus integer: number of processes to be used in parallel operation: typically one would chose this to be
##' the number of available CPUs
##' @param cl An optional parallel or snow cluster for use if \code{parallel="snow"}. If not supplied, a
##' cluster on the local machine is created for the duration of the call.
##' @return A \code{GammaImputedSet.object} containing the imputed data sets 
##' 
##' @seealso \code{\link{GammaImputedSet.object}} \code{\link{GammaImputedData.object}} 
##' @examples 
##' data(nwtco)
##' nwtco <- nwtco[1:500,]
##'  
##' #creating 2 imputed data sets (m=2) for speed, would normally create more
##' ans <- gammaImpute(formula=Surv(edrel,rel)~histol + instit,
##'                    data = nwtco, m=2, gamma.factor=1, DCO.time=6209)
##' 
##' #subject specific gamma (multiplied by gamma.factor to give the jump)
##' #NA for subjects that are not to be imputed
##' jumps <- c(rep(NA,10),rep(1,490))                                      
##' DCO.values <- rep(6209,500)
##'                                                                               
##' ans.2 <- gammaImpute(formula=Surv(edrel,rel)~histol + instit + strata(stage),
##'                    data = nwtco, m=2, bootstrap.strata=strata(nwtco$stage),
##'                    gamma=jumps, gamma.factor=1, DCO.time=DCO.values)                    
##' 
##' #can also use column names
##' nwtco$gamma <- jumps
##' nwtco$DCO.time <- DCO.values
##' ans.3 <- gammaImpute(formula=Surv(edrel,rel)~histol + instit + strata(stage),
##'                    data = nwtco, m=2, bootstrap.strata=strata(nwtco$stage),
##'                    gamma="gamma", DCO.time="DCO.time")                                        
##'                    
##' @export
gammaImpute <- function(formula, data, m, gamma, gamma.factor,
                        bootstrap.strata=rep(1,nrow(data)), DCO.time,...,
                        parallel = c("no", "multicore", "snow")[1], ncpus = 1L, cl = NULL){
  
  #First sort out missing arguments for gamma
  if(missing(gamma) && missing(gamma.factor)){
    stop("At least one of gamma and gamma.factor must be included")
  }
  if(missing(gamma.factor)) gamma.factor <- 1
  if(!missing(gamma) && is.null(gamma)) stop("gamma cannot be NULL")
  if(missing(gamma)) gamma <- NULL
  
  #validate RHS of formula
  .validRHSFormula(getCovarFormula(formula))
  
  #Get the Surv object
  surv.object <- model.frame(formula,data)[,1]
  surv.times <- as.matrix(surv.object)
  if(class(surv.object) != "Surv" ||  attr(surv.object,"type")!="right" || ncol(surv.times)!=2){
    stop("Can only use a right censored Surv object on the left hand side of the formula")
  }
  
  #perform main validation of input arguments
  validate.Gamma.arguments(data,surv.times,m,gamma,bootstrap.strata,
                                   gamma.factor,DCO.time,match.call(expand.dots = TRUE))
  validate.parallel.arguments(parallel,ncpus,cl)
  
  #setup DCO.time column if needed
  if(!is.character(DCO.time)){
    data[,"internalDCO.time"] <- if(length(DCO.time)>1) DCO.time else rep(DCO.time,nrow(data)) 
    DCO.time <- "internalDCO.time"
  }
  
  #create the list of bootstrapped data indices
  boot.indices <- boot::boot(data=data,R=m,statistic=function(data,indices){return(indices)},strata=bootstrap.strata)
  boot.indices <- lapply(seq_len(m),function(x){boot.indices$t[x,]})
  
  #set up gamma
  data$internal_gamma_val <- getGamma(data,gamma,gamma.factor)
  
  #perform imputation
  if(parallel=="no"){
    imputes <- lapply(boot.indices,.singleGammaImpute,data=data,
                      surv.times=surv.times,DCO.time=DCO.time,formula=formula,...)
  }
  else{
    imputes <- parallelRun(parallel,ncpus,cl,lapply.list=boot.indices,
                           FUN=.singleGammaImpute,data,surv.times,DCO.time,formula,...)
  }
  
  #package up results
  retVal <- .createImputedSet(m,col.control=NULL,data,imputes,"Gamma",getCovarFormula(formula))
  retVal$gamma.factor <- gamma.factor
  retVal
}

##' @name ExtractSingle
##' @export
ExtractSingle.GammaImputedSet <- function(x,index){
  retVal <- .internalExtract(x,index,fit=FALSE)
  class(retVal) <- "GammaImputedData"
  retVal
}

