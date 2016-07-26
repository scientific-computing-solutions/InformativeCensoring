#This file contains the code used to fit an analyze gamma imputed
#data sets. Most of these functions call functions inside generic.R
#and the Roxygen comments for the generics are found there

##' @name ImputeStat
##' @export
ImputeStat.GammaImputedData <- function(object,method=c("Cox","weibull","exponential")[1],formula=NULL,...){
  if(!method %in% c("Cox","weibull","exponential")){
    stop("Invalid method for gamma imputation must be one of Cox, weibull or exponential")
  }  
  
  test.stats <- .imputeStat.internal(object,method,formula,...)
  class(test.stats) <- "GammaStat"
  test.stats
}


##' @name ImputeStat
##' @export
ImputeStat.GammaImputedSet <- function(object,method=c("Cox","weibull","exponential")[1],formula=NULL,...,
                                       parallel = c("no", "multicore", "snow")[1], ncpus = 1L, cl = NULL){
  
  fits <- .internalImputeStatset(object,method,formula,...,parallel=parallel,ncpus=ncpus,cl=cl)
  
  #function to extract the estimates and variances from the fit
  .setupStats <- function(val){
    retVal <- do.call(rbind,lapply(fits,function(x){matrix(x[[val]],nrow=1)}))
    colnames(retVal) <- names(fits[[1]]$estimate) 
    retVal
  }
  
  retVal <- list(fits=fits,
                 statistics=list(estimates=.setupStats("estimate"),
                                 vars=.setupStats("var") ),
                 m=object$m)
  class(retVal) <- "GammaStatList"
  return(retVal)
}

##' @export
summary.GammaStatList <- function(object,...){
  
  M <- object$m
  estimates <- colMeans(object$statistics$estimates)
  mean.variances <- colMeans(object$statistics$vars)
  var.of.estimates <- apply(object$statistics$estimates,2,var)
  
  se <- sqrt(mean.variances + (1 + (1/M)) * var.of.estimates) 
  t <- estimates/se
  df <- (M-1)*(1+(mean.variances/((1+(1/M))*var.of.estimates)))^2
  
  t.val <- qt(1-(1-0.95)/2,df=df)
  
  retVal <- data.frame(est=estimates,
                       se=se,
                       t=t,
                       df=df,
                       "Pr(>|t|)"=2*(1-pt(abs(t),df=df)),
                       "lo 95"=estimates-se*t.val,
                       "hi 95"=estimates+se*t.val,
                       check.names=FALSE) 
  as.matrix(retVal,ncol=ncol(retVal))
}


##' \code{GammaStat} object
##' 
##' An S3 object which contains the point estimate 
##' and test statistic after fitting a model to 
##' a \code{GammaImputedData} object.
##' 
##' The function \code{print.GammaStat} has been implemented
##' 
##' The object contains the following:
##' @slot model The model used to create the fit
##' @slot method The model used for the fit
##' @slot estimate A point estimate of the test statistic
##' @slot var The estimate for the variance of the test statistic
##' @name GammaStat.object
NULL

##' \code{GammaStatList} object
##' 
##' The object containing the results of fitting models to
##' a \code{GammaImputedSet} object.
##' 
##' A \code{summary.GammaStatList} has been implemented which performs
##' Rubin's multiple imputation rules. 
##' 
##' The object contains the following
##' @slot fits A list of \code{GammaStat} objects containing the model fits for
##' the imputed data sets
##' @slot statistics A list with two elements: estimates and vars which contain the coefficient
##' estimates and their variances one column per covariate one row per imputed data set  
##' @slot m The number of model fits
##' @name GammaStatList.object
NULL

##' @export
print.GammaStat <- function(x,...){
  print(x$model)  
}
 
##' @name ExtractSingle
##' @export
ExtractSingle.GammaStatList <- function(x,index){
  .internalExtract(x,index,fit=TRUE)
}

##' @export 
cox.zph.GammaStat <- function(fit,transform="km",global=TRUE,...){
  if(!inherits(fit$model, "coxph"))  stop("The model fit is not Cox!")
  cox.zph(fit$model,transform=transform,global=global)
}

##' @export
cox.zph.GammaStatList <- function(fit,transform="km",global=TRUE,index,...){
  cox.zph(ExtractSingle(fit,index),transform=transform,global=global,...)
}

