#This fle contains the functions used to validate the gammaimpute arguments
#(see validation.R for more some subfunctions (e.g. is.valid.call)) and
#parallel.R for validation of arguments associated with parallelism


validate.Gamma.arguments <- function(data,surv.times,m,gamma,strata,gamma.factor,DCO.time,Call){
  
  is.valid.call(Call)
  
  if(nrow(data)==0){
    stop("Empty data frame!")
  }
  
  if(any(c("impute.time","impute.event","internal_gamma_val","internalDCO.time") %in% colnames(data))){
    stop("Cannot use a data frame with columns impute.time, impute.event",
          " internalDCO.time or internal_gamma_val")
  } 
  
  if(!.internal.is.finite.number(m) ||!.internal.is.wholenumber(m) || m < 2){
    stop("m must be an integer and at least 2")
  }
  
  if(!is.numeric(gamma.factor) || is.na(gamma.factor) || length(gamma.factor)>1){
    stop("gamma.factor must be a single number that is multiplied by the values",
         " of the gamma argument in order to create subject specific jumps in the hazard rate.",
         " see help(gammaImpute) for examples")
  }
  
  #other errors on bootstrap.strata argument will be caught by boot
  if(length(strata) != nrow(data)){
    stop("Invalid strata argument it must be a vector the same length as the number of rows in the dataset")
  }
  
  validateDCO.time(DCO.time=DCO.time,data=data,times=surv.times[,1])
 
  if(!is.null(gamma)){
    validateGammaVal(gamma,data)
  }  
  
  #check time is positive
  if(any(surv.times[,1]<= 0)){
    stop("Time on study must be positive")
  }
  
}


validateDCO.time <- function(DCO.time,data,times){
  
  #attempt to get DCO.time to be a vector of DCO.times 
  #even if single number or column name is input
  if(is.character(DCO.time)){
    if(length(DCO.time) != 1 || !DCO.time %in% colnames(data)){
      stop("Invalid DCO.time argument")
    }
    DCO.time <- data[,DCO.time]
  }
  else if(length(DCO.time)==1){
    DCO.time <- rep(DCO.time,nrow(data)) 
  }
  else if(length(DCO.time)!=nrow(data)){
    stop("Invalid DCO.time length")
  }

  #Now validate the DCO.time
  if(any(!is.numeric(DCO.time) | is.infinite(DCO.time))){
    stop("DCO times must be numeric and finite")
  }

  #DCO.time is <= time
  if(any(DCO.time < times)){
    stop("DCO.time must be >= time for all subjects ")
  }
}

validateGammaVal <- function(gamma,data){
  #first attempt to get gamma a vector of gamma values
  #whatever format gamma was input as
  if(is.character(gamma)){
    if(length(gamma)!= 1 || !gamma %in% colnames(data)){
      stop("Invalid gamma column name")
    }
    gamma <- data[,gamma]
  }

  #then check gamma is the right length and numeric
  if(length(gamma)!= nrow(data)){
    stop("Invalid length of gamma its length must be the number of subjects in the",
         " the dataset. use the gamma.factor argument if you want to use a single number",
         " for each subject's jump in hazard. See help(gammaImpute) for futher details")
  }

  if(any(!is.na(gamma) & !is.numeric(gamma))){
    stop("gamma must be numeric")
  }
}