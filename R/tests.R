#The test statistics for fitting models to the imputed data
#they should return a list with 4 arguments which will
#be set to be a ScoreStat.object or GammaStat.object at a higher level (and the
#test statistic(s) will be calculated as estimate/sqrt(var) if needed)

.Testlogrank <- function(object,formula,...){
  .Testsurvdiff(object,formula,rho=0,method="logrank",...)
}

#Note this is the Peto & Peto modification of the Gehan-Wilcoxon test
.Testwilcoxon <- function(object,formula,...){
  .Testsurvdiff(object,formula,rho=1,method="Wilcoxon",...)
}


.Testsurvdiff <- function(object,formula,rho,method,...){
  
  model <- survdiff(formula=formula,data=object$data,rho=rho,...)
  
  list(model=model,
       method=paste(method,"(estimator for O-E)"),
       estimate=if(class(model$obs)=="matrix") sum(model$obs[2,]-model$exp[2,]) else (model$obs[2]-model$exp[2]),
       var=model$var[2,2])
}


.TestWeibull <- function(object,formula,...){
  .Testsurv("weibull",survreg,object,formula,dist="weibull",...)
}

.TestExponential <- function(object,formula,...){
  .Testsurv("exponential",survreg,object,formula,dist="exponential",...)
}

.Testcox <- function(object,formula,...){
  .Testsurv("Cox",coxph,object,formula,...)
}

.Testsurv <- function(method,model.function,object,formula,...){
  model <- model.function(formula,data=object$data,model=TRUE,...)
  list(model=model,
       method=method,
       estimate=model$coefficients,
       var=diag(vcov(model))[1:length(model$coefficients)]) 
  #we need to use [1:length(...)] in line above as the variances of the scale variable(s)
  #are included in the vcov matrix if the Weibull method is used
}


#if no formula argument is given then we are using the default from the object also checking if a formula
#argument is given then it is valid see .validRHSformula function for details
.getFormula <- function(formula,arm,method){
  
  #check the formula is valid
  .validRHSFormula(formula,arm)
  
  if(method %in% c("logrank","Wilcoxon")){
    #checking all but treatment arm are strata if using logrank/Wilcoxon
    tms<-terms(formula,specials=c("strata"))
    if(length(untangle.specials(tms,special = "strata")$vars) != length(attr(terms(formula),"term.labels")) - 1){
      stop("Cannot include non-stratified covariates in logrank or Wilcoxon method")  
    }
  }
  
  #Add the appropriate left hand side of the formula
  update(formula,paste("Surv(impute.time,impute.event)~."))
}


     