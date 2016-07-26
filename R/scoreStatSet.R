#This file contains the functions associated with
#the ScoreStatList object - the statistics from a set of 
#analyzed risk score imputed data sets. The summary.ScoreStatList object
#code and object documentation are found here (it
#implements the meth1, meth2 averaging methods in the Hsu and Taylor paper) 

##' ScoreStatList
##' 
##' The object containing the results of fitting models to
##' a \code{ScoreImputedSet} object.
##' 
##' A \code{summary.ScoreStatList} has been implemented.
##' 
##' The object contains the following
##' @slot fits A list of \code{ScoreStat} objects containing the model fits for
##' the imputed data sets
##' @slot statistics A \code{ScoreStatSet} object containing the statistics
##' @slot m The number of model fits  
##' @name ScoreStatList.object
##' @seealso \code{\link{ScoreStatSet.object}} \code{\link{ScoreStat.object}}
NULL

##' @export
summary.ScoreStatList <- function(object,...){
  summary(object$statistics)
}

##' @name ExtractSingle
##' @export
ExtractSingle.ScoreStatList <- function(x,index){
  .internalExtract(x,index,fit=TRUE)
}

##' An object which contains the test statistic and estimators for
##' a set of model fits to imputed data using risk score imputation
##' 
##' The object is a Mx3 matrix, one row per imputed data set
##' and columns: estimate (the point estimates), var (their variances)
##' and Z (the test statistic). M must be > 4
##' 
##' Note the Z should be ~ standard normal (so we do not use the chi_squared
##' test statistic see \code{\link{ScoreStat.object}} for further details)
##' 
##' The summary.ScoreStatSet function will apply the MI averaging procedures
##' and estimates of the test statistic and p-value    
##' @seealso \code{\link{summary.ScoreStatSet}}            
##' @name ScoreStatSet.object 
NULL

##' S3 generic to create a \code{ScoreStatSet} object
##' @param x The object to convert into a \code{ScoreStatSet} object
##' @return A ScoreStatSet object
##' @seealso \code{\link{ScoreStatSet.object}} 
##' @export
ScoreStatSet <- function(x){
  UseMethod("ScoreStatSet")  
}

##' @export
ScoreStatSet.default <- function(x){
  stop("No method ScoreStatSet for this object")
}


##' @export
ScoreStatSet.matrix <- function(x){
  if(ncol(x)!=3){
    stop("Matrix must have 3 columns")
  }
  
  if(any(x[,2]<0)){
    stop("The second column of the matrix is a column of variances, they must be positive")
  }
  
  if(nrow(x)<5){
    stop("Results for at least five data sets are required for the MI averaging calculation")
  }
  
  if(any(x[,1]/sqrt(x[,2])!= x[,3])){
    stop("The Z statistic [third col] should be the estimate [first col]/ sqrt(var [2nd col] )")
  }
  
  colnames(x) <- c("estimate","var","Z")
  class(x) <- "ScoreStatSet"
  x  
}


##' Summary object of \code{ScoreStatSet} object
##' 
##' This object contains the multiple imputed 
##' averages/p-values of a set of estimates from 
##' risk score imputed data sets.
##' 
##' A \code{print.summary.ScoreStatSet} function has been implemented
##' 
##' This object contains three lists meth1 and meth2 and methRubin
##' meth1 averages the point estimates to produce an F test statistic,
##' meth2 averages the test statistics and prodcues a t test statistic
##' and methRubin follows Rubin's standard rules and is used for calculating
##' confidence intervals
##' 
##' See the vignette for further details.
##'  
##' 
##' meth1, meth2 and methRubin are lists with the following elements:
##' estimate: average estimator for meth1, NOTE: for meth2 this is the average test statistic, \cr
##' var: estimate of variance of "estimate" field \cr
##' test.stat: test statistic \cr
##  df: degrees of freedom of test statistic \cr
##' distribution: distribution of statistical test (i.e. F or t) \cr
##' p.value: p-value of test \cr  
##' 
##' 
##' @name summary.ScoreStatSet
##' @aliases summary.ScoreStatSet.object
NULL

##' @export
summary.ScoreStatSet <- function(object,...){
  
  retVal <- list()
  
  #the same names are used as in Hsu's paper
  M <- nrow(object)
  
  theta.bar <- mean(object[,"estimate"])
  U1 <- mean(object[,"var"])
  B1 <- var(object[,"estimate"])
  V1 <- U1 + (1+1/M)*B1
  D <- theta.bar^2/V1
  t <- M-1
  r <- (1+1/M)*B1/U1
  v1 <- 4+(t-4)*(1+(1-2/t)/r)^2
  
  if(M==5 && r==0){
    warning("Cannot calculate degrees of freedom for meth1 when r=0 and M=5")
  }
  
  retVal[["meth1"]] <- list(estimate=theta.bar,
                            var=V1,
                            test.stat=D,
                            df=c(1,v1),
                            distribution="F",
                            p.value=1-pf(D,df1=1,df2=v1))
 
  Z.bar <- mean(object[,"Z"])
  B2 <- var(object[,"Z"])
  V2 <- 1 + (1+1/M)*B2
  v2 <- (M-1)*(1 + (M/(M+1))*(1/B2) )^2
  
  retVal[["meth2"]] <-list(estimate=Z.bar,
                           var=V2,
                           test.stat=Z.bar/sqrt(V2),
                           df=v2,
                           distribution="t",
                           p.value=2*(1-pt(abs(Z.bar/sqrt(V2)),v2)))
  
  rubin.df <- (M-1)*(1+ U1/((1+1/M)*B1))^2
  
  retVal[["methRubin"]] <- list(estimate=theta.bar,
                                var=V1,
                                test.stat=theta.bar/sqrt(V1),
                                df=rubin.df,
                                distribution="t",
                                p.value=2*(1-pt(abs(theta.bar/sqrt(V1)),rubin.df)))
  
  class(retVal) <- "summary.ScoreStatSet"
  retVal
} 


##' @export
print.summary.ScoreStatSet <- function(x,...){
  cat("Summary statistics of Imputed Data Sets\n")
  cat("Method 1 (averaging of point estimates):\n")
  cat("Estimator:",x$meth1$estimate,"\n")
  cat("Var of Estimator:",x$meth1$var,"\n")
  cat("Test statistic:",x$meth1$test.stat,"\n")
  cat("Distribution:",x$meth1$distribution,"\n")
  cat("with ",x$meth1$df[1],", ",x$meth1$df[2]," degrees of freedom\n",sep="")
  cat("giving a p-value of",x$meth1$p.value,"\n\n")

  cat("Method 2 (averaging of test statistics):\n")
  cat("Test statistic:",x$meth2$test.stat,"\n")
  cat("Distribution:",x$meth2$distribution,"\n")
  cat("with",x$meth2$df,"degrees of freedom\n")
  cat("giving a p-value of",x$meth2$p.value,"\n\n")
}

##' @export
confint.summary.ScoreStatSet <- function(object, parm, level = 0.95, ...){
  if(!.internal.is.finite.number(level) || level <= 0 || level >= 1){
    stop("Invalid level argument")
  }
  
  standard.error <- sqrt(object$methRubin$var)
  t_val <- qt(1-(1-level)/2,df=object$methRubin$df)
  retVal <- object$methRubin$estimate + c(-1,1)*t_val*standard.error
  names(retVal) <- c(paste(100*(1-level)/2,"%",sep=""),paste(100*(1 - (1-level)/2) ,"%",sep=""))
  retVal
}
