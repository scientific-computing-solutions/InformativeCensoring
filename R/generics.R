#This file contains the functions which are common to both imputation
#methods - these are mainly S3 generics and their internals 

##' @import survival
NULL

##' @importFrom boot boot
NULL

##' @importFrom dplyr inner_join
NULL

.internal.is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol

.internal.is.finite.number <- function(x){
  if(!is.numeric(x) || is.na(x) || length(x)>1 || is.infinite(x)){
    return(FALSE)
  }
  return(TRUE)
}

##' Extract a single risk score/gamma imputed data set/model fit
##' @param x The multiple imputed object
##' @param index Integer, which imputed data set/model fit should be returned
##' @return The individual data set/model fit 
##' @export
ExtractSingle <- function(x,index){
  UseMethod("ExtractSingle")  
}

##' @export
ExtractSingle.default <- function(x,index){
  stop("No method ExtractImputedData for this object")
}

#performs the ExtractSingle generic for various x
#@param x the S3 object for which ExtractSingle was used
#@param index - which individual object should be output
#@param fit - logical, if TRUE outputting a *Stat object, if not
#@ a *ImputedData object
.internalExtract <- function(x,index,fit){
  
  if(!.internal.is.finite.number(index) || !.internal.is.wholenumber(index) || index <= 0 || index > x$m){
    stop("invalid index")
  }
  
  if(fit){
    return(x$fits[[index]])
  }
  
  retVal <- list(data=x$data,
                 col.control=x$col.control,
                 default.formula=x$default.formula)
  
  retVal$data$impute.time <- x$impute.time[,index]
  retVal$data$impute.event <- x$impute.event[,index]
  retVal
}

##' S3 generic to fit model(s) to risk score/gamma Imputed objects
##'  
##' @param object A \code{ScoreImputedData}, \code{ScoreImputedSet}, \code{GammaImputedData} or \code{GammaImputedSet} object 
##' to fit the model to   
##' @param method The type of statistical model to fit. There are three methods which can be performed when using
##' Risk Score imputation   \cr
##' "logrank": a logrank test using \code{survival::survdiff} \cr
##' "Wilcoxon": Peto & Peto modification of the Gehan-Wilcoxon test using \code{survival::survdiff}
##' with \code{rho=1} \cr
##' "Cox": Fit a cox model using \code{survival::coxph}  \cr 
##' 
##' For gamma imputation the model can be "Cox" (using \code{survival::coxph}), 
##' "weibull" or "exponential" both using \code{survival::coxph}
##'  
##'  
##' @param formula The model formula to fit.
##' If no formula argument is used, then object$default.formula will be used.
##' For risk score imputation this is \code{~ treatment.group} and for gamma imputation
##' this is the formula used when fitting the Cox model
##' 
##' For \code{method="Cox"}, additional covariates can be included by explictily giving a 
##' formula argument. For logrank/Wilcoxon only additional strata terms can be
##' included. 
##'  
##' In all cases only the right hand side of the formula is required
##' The survival object on the left hand side is created automatically
##' E.g. for a Cox model could use formula=~arm + covar1. The cluster and tt options cannot be used
##' See the vignettes for further details
##' @param ... Additional arguments which are passed into the model fit function
##' @seealso \code{\link{ScoreStat.object}} \code{\link{ScoreImputedData.object}} 
##' @export
ImputeStat <- function(object,method=c("logrank","Wilcoxon","Cox","weibull","exponential")[1],formula,...){
  UseMethod("ImputeStat")  
}

##' @export
ImputeStat.default <- function(object,method=c("logrank","Wilcoxon","Cox","weibull","exponential")[1],formula,...){
  stop("No method ImputeStat for this object")
}

#Internal function for fitting the imputed data sets
#see ImputeStat function documentation
.imputeStat.internal <- function(object,method,formula,...){
  formula <- .getFormula(if(is.null(formula)) object$default.formula else formula,
                         object$col.control$arm,method)
  switch(method,
         logrank=.Testlogrank ,
         Wilcoxon=.Testwilcoxon ,
         Cox=.Testcox,
         exponential=.TestExponential,
         weibull=.TestWeibull,
         function(object,formula,...){stop("Unknown test type")})(object,formula,...)
}


#extract the imputed data and perform the fit
.internalImputeStatset <- function(object,method,formula,...,parallel,ncpus,cl){
  validate.parallel.arguments(parallel,ncpus,cl)
  
  extract_stat <- function(x,method,formula){
    ImputeStat(ExtractSingle(object,x),method=method,formula=formula,...)}
  
  if(parallel=="no"){
    fits <- lapply(seq_len(object$m),extract_stat,method,formula)
  }
  else{
    fits <- parallelRun(parallel,ncpus,cl,lapply.list=seq_len(object$m),
                        FUN=extract_stat,method,formula,...)
  }
}


#constructor function for the Score/Gamma ImputedSet objects
.createImputedSet <- function(m,col.control,data,imputes,classprefix,default.formula){
  retVal <- list(m=m,
               col.control=col.control,
               default.formula=default.formula,
               data=data,
               impute.time=vapply(imputes,"[[","impute.time",FUN.VALUE=numeric(nrow(data))),
               impute.event=vapply(imputes,"[[","impute.event",FUN.VALUE=numeric(nrow(data))))
  class(retVal) <- paste(classprefix,"ImputedSet",sep="")
  retVal
}

##' Test Cox proportional hazards assumption
##' 
##' See cox.zph function in the survival package
##' @inheritParams survival::cox.zph 
##' @param ... Additional arguments to cox.zph, for example \code{index} if
##' fit is a \code{GammaStatList}  object 
##' @seealso \code{\link[survival]{cox.zph}}
##' @export
cox.zph <- function(fit,transform="km",global=TRUE,...){
  UseMethod("cox.zph")
}

##' @export
cox.zph.default <- function(fit,transform="km",global=TRUE,...){
  survival::cox.zph(fit,transform=transform,global=global)
}

