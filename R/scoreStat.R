#This file contains the functions associated with a model fit
#of a single risk score imputed data set

##' ScoreStat object
##' 
##' An S3 object which contains the point estimate 
##' and test statistic after fitting a model to 
##' a \code{ScoreImputedData} object.
##' 
##' The functions \code{print.ScoreStat} and \code{as.vector.ScoreStat}
##' have been included 
##' 
##' The object contains the following:
##' @slot model The model used to create the fit
##' @slot method The method used for the fit
##' @slot estimate A point estimate of the test statistic
##' @slot var The estimate for the variance of the test statistic
##' @slot statistic The test statistic given by \code{estimate/sqrt(var)} 
##' 
##' @details The test statistic should be normally distributed and hence for
##' the logrank test Z = (O_2 - E_2)/sqrt(V_2), i.e. the square root of the standard
##' Chi squared statistic (with the appropriate sign)
##' 
##' @name ScoreStat.object
NULL

##' @export
ImputeStat.ScoreImputedData <- function(object,method=c("logrank","Wilcoxon","Cox")[1],formula=NULL,...){
  if(!method %in% c("logrank","Wilcoxon","Cox")){
    stop("Invalid method for risk score imputation it must be logrank, Wilcoxon or Cox")
  }  
  
  test.stats <- .imputeStat.internal(object,method,formula,...)
  #only take the first covariate: the treatment arm
  test.stats$estimate <- test.stats$estimate[1]
  test.stats$var <- test.stats$var[1]
  test.stats$statistic <- test.stats$estimate/sqrt(test.stats$var)
  
  class(test.stats) <- "ScoreStat"
  test.stats
}


##' @method  as.vector ScoreStat
##' @export
as.vector.ScoreStat <- function(x,mode="any"){
  ans <- c(x$estimate,x$var,x$statistic)
  names(ans) <- c("estimate","var","statistic")
  ans
}


##' @export
print.ScoreStat <- function(x,...){
  cat("Method used:",x$method,fill = TRUE)
  cat("Point Estimate:",x$estimate,fill=TRUE)
  cat("Variance estimate:",x$var,fill=TRUE)
  cat("Z stastistic:",x$statistic,fill=TRUE)
  cat("Use x$model to view the model fit")
} 