#This file contains code and documentation
#associated with the ScoreImputedSet object (a set of risk score imputed
#datasets) 

##' \code{ScoreImputedSet} object
##' 
##' An object which contains the set of score imputed data frames.
##' Use the \code{ExtractSingle} function to extract a single 
##' \code{ScoreImputedData} object. Use the \code{ScoreStat} function to fit models
##' to the entire set of imputed data frames 
##' 
##' It contains the following:
##' @slot data A data frame containing the unimputed time to event data
##' @slot col.control The list of column names the score imputation method requires see \code{\link{col.headings}}
##' for further details
##' @slot m The number of imputed data sets
##' @slot impute.time A matrix (1 column per imputed data set) containing the imputed times 
##' @slot impute.event A matrix (1 column per imputed data set) containing the imputed event indicators
##' @slot default.formula The default model formula which will be used when fitting the imputed data using a Cox model 
##' @name ScoreImputedSet.object
##' @seealso \code{\link{ScoreImputedData.object}}
NULL

 
##' @name ExtractSingle
##' @export
ExtractSingle.ScoreImputedSet <- function(x,index){
  retVal <- .internalExtract(x,index,fit=FALSE)
  class(retVal) <- "ScoreImputedData"
  retVal
}



##' @param parallel The type of parallel operation to be used (if any), can be used for \code{GammaImputedSet} and \code{ScoreImputedSet} 
##' @param ncpus integer: number of processes to be used in parallel operation: typically one would chose this to be
##' the number of available CPUs, can be used for \code{GammaImputedSet} and \code{ScoreImputedSet}.
##' @param cl An optional parallel or snow cluster for use if \code{parallel="snow"}. If not supplied, a
##' cluster on the local machine is created for the duration of the call, can be used for \code{GammaImputedSet} and \code{ScoreImputedSet}.
##' @name ImputeStat
##' @export
ImputeStat.ScoreImputedSet <- function(object,method=c("logrank","Wilcoxon","Cox")[1],formula=NULL,...,
                                       parallel = c("no", "multicore", "snow")[1], ncpus = 1L, cl = NULL){
  
  fits <- .internalImputeStatset(object,method,formula,...,parallel=parallel,ncpus=ncpus,cl=cl)
  statistics <- ScoreStatSet(t(vapply(fits,function(x){as.vector(x)},FUN.VALUE=numeric(3))))

  retVal <- list(fits=fits,
                 statistics=statistics,
                 m=object$m)
  
  class(retVal) <- "ScoreStatList"
  return(retVal)
}

