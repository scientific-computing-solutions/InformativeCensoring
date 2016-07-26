#The code contains the functions called by gammaImpute
#(not including the validation routines)

# Main function which performs gamma imputation to generate imputed times
# for a single Imputation
# @inheritParams gammaImpute
# @param indices The indicies of the rows of the data set used to get the bootstrapped data for fitting the Cox model 
# @param surv.times, matrix representation of the Surv object of the original event times
# and indicator
# @param DCO.time the column of data cutoff times of the data frame
# @return A data frame with the imputed times and events for this imputed dataset
# (subjects who are not imputed return unchanged values)
# @details This function performs the following steps:
# 1) fit a cox model to the bootstrapped data frame
# 2) calculates the baseline hazards by calling .SetupStrataBaseHaz
# 3) calculates the linear predictors (lp)
# 4) For each subject calculate an imputed time and event indicator (using .CoxMImpute)
# 5) wrap up the event times into a dataframe
.singleGammaImpute <- function(indices,data,surv.times,DCO.time,formula,...){
  
  boot.data <- data[indices,]
  model <- coxph(formula=formula,data=boot.data,model=TRUE,na.action=na.fail,...)
  
  stratas <- untangle.specials(model$terms,"strata")$vars 
  if(length(stratas) > 1){
    stop("Cannot have multiple strata arguments. Use strata(x,y) instead of strata(x)+strata(y)")
  }
  
  basehaz.details <- .SetupStrataBaseHaz(data,formula,stratas,basehaz(model))
  lp <- predict(model, data, type="lp",reference="sample")+data$internal_gamma_val
  
  #perform imputation for one subject at a time
  imputes <- lapply(seq_len(nrow(data)),function(x){
    
    my.time <- surv.times[x,1]
    
    #if have event or not imputing then do not impute
    if(surv.times[x,2]==1 || is.na(data[x,"internal_gamma_val"])){
      return(list(event=surv.times[x,2],time=my.time))
    }
    
    #(basehaz.details$index is a vector of which basehaz.details$base.hazard data frame
    #to be used for each subject)
    imputed.time <- .CoxMImpute(my.time,lp[x],
                                basehaz.details$base.hazard[[basehaz.details$index[x]]])
    DCO.time <- data[x,DCO.time]
    
    return(list(event=(imputed.time < DCO.time),
                time=min(DCO.time,imputed.time)))  
  })
  
  data.frame(impute.event=vapply(imputes,"[[","event",FUN.VALUE = numeric(1)),
             impute.time=vapply(imputes,"[[","time",FUN.VALUE = numeric(1)))  
}

#function which creates the baseline hazard data frames, one per strata
#and which data frame is required by each data point
# @param data The initial data frame for which imputation is required
# @param formula The formula used for the model fit
# @param stratas The stratas of the model (untangle.specials(terms,"strata"))
# @param base.hazard The output of basehaz(model)
# @return a list with two elements, base.hazard: a list of stratified baseline hazard
# data frames, named by the their strata values and index: A nrow(data) length
# vector where index[i] is the strata value required for subject in row i. If the
# baseline hazard function is not stratified then base.hazard is a list containing one data frame
# named "all" and index is a vector with every element = "all"
.SetupStrataBaseHaz <- function(data,formula,stratas,base.hazard){
  if(!"strata" %in% colnames(base.hazard)){
   return(list(index=rep("all",nrow(data)),
        base.hazard=list(all=base.hazard)))
  }
  
  
  base.hazard <- split(base.hazard,f=base.hazard$strata)
  
  #which basehaz data frame should each subject use?
  #cannot find a nicer way to do this...
  mf <- model.frame(formula = formula,data = data)
  
  index <- vapply(seq_len(nrow(data)),
                  function(i) paste(unlist(lapply(stratas,
                      function(x){as.character(mf[i,x])})),collapse=", "),
                  FUN.VALUE = character(1))
  
  if(any(!index %in% names(base.hazard))){
    stop("Cannot find baseline hazard function for a given strata ",
         "check your bootstrap.strata argument as you may be producing invalid ",
         "stratified bootstrapped samples")
  }
  
  list(index=index,
       base.hazard=base.hazard)
}
                      

# Function which performs the imputation for a single
# subject (we already assume a subject is to have time imputed)
# @param my.time The current time of censoring 
# @param beta The lp prediction from the model fit including 
# the step change gamma (= \hat{beta_j}Z_i+ gamma_i) see equation (4) in paper
# @param basehaz A data frame containing the baseline hazard function (time, hazard, [possibly strata])
# for the strata of this subject's strata
# @return A proposed imputed event time (higher up this could be replaced
# by censoring at DCO.time) 
.CoxMImpute <- function(my.time,beta,basehaz){
  
  basehaz$strata <- NULL
  basehaz$time <- basehaz$time-my.time
  
  if(any(basehaz$time<=0)){
    basehaz$hazard <- basehaz$hazard - max(basehaz$hazard[basehaz$time<=0])
  }
  
  #The 0,0 is inserted because all cumulative hazards pass through the origin. 
  #This also has the effect of making gamma=Inf give instant failure
  basehaz <- rbind(data.frame(time=0,hazard=0),
                   basehaz[basehaz$time>0,],
                   data.frame(time=Inf,hazard=Inf))
  
  quant <- -log(runif(1))*exp(-beta)
  my.time + min(basehaz$time[basehaz$hazard>= quant])
  
}

#simple input function which calculates subject specific jumps in 
#hazard given (already validated arguments)
getGamma <- function(data,gamma,gamma.factor){
  if(is.null(gamma)){
    return(rep(gamma.factor,nrow(data)))
  }
  if(is.character(gamma)){
    return(data[,gamma]*gamma.factor)
  }
  gamma*gamma.factor
}   


#from nlme, output the rhs of a formula as a formula object
getCovarFormula <- function(form){
  if (!(inherits(form, "formula"))) {
    stop("invalid formula")
  }
  form <- form[[length(form)]]
  if (length(form) == 3 && form[[1]] == as.name("|")) {
    form <- form[[2]]
  }
  eval(substitute(~form))
}