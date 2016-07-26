#This file contains the main function to perform the risk
#score imputation and the roxygen comments for the objects
#associated with generating the risk score imputation
#See internalScore.R for subfunctions and validation.R for 
#many validation routines

##' \code{ScoreImputedData} object
##' 
##' An object which contains 
##' 
##' @slot data A data frame containing the time to event data
##' with 2 new columns impute.time and impute.event, the imputed event/censoring times and event indicators
##' (for subjects whose data is not imputed these columns contain the unchanged event/censoring time and
##' event indicator )
##' @slot col.control The list of column names the risk score imputation method requires see \code{\link{col.headings}}
##' for further details. If censor.type was not used then \code{col.control$censor.type="using_has.event_col"}
##' @slot default.formula The default model formula which will be used when fitting the imputed data using a Cox model 
##' @name ScoreImputedData.object
NULL


##' Perform risk score multiple imputation method
##' @param data The data set for which imputation is required
##' @param event.model The right hand side of the formula to be used for fitting the Cox model for calculating the time to 
##' event score e.g. ~Z1+Z2+Z3.
##' @param censor.model The right hand side of the formula to be used for fitting the Cox model for calculating the time to 
##' censoring score if not included then \code{event.model} will be used
##' @param col.control A list of the columns names of \code{data} which are used by the imputation algorithm
##' See example below and for further details of these columns and their purpose see \code{\link{col.headings}}
##' @param NN.control Parameters which control the nearest neighbour algorithm. See \code{\link{NN.options}}
##' @param time.dep A ScoreTD object, to be included if the time dependent score imputation method is to be used, otherwise it should be NULL 
##' @param m The number of data sets to impute
##' @param bootstrap.strata When performing the bootstrap procedure for fitting the models,
##' how should the data be stratified (see strata argument to \code{boot::boot}). if argument
##' is not used then standard sampling with replacement is used to generate the bootstrap data
##' @param ... Additional arguments passed into the Cox model Note the subset and na.action arguments should not be used 
##' (na.fail will be used when fitting the Cox model)
##' @param parallel The type of parallel operation to be used (if any).
##' @param ncpus integer: number of processes to be used in parallel operation: typically one would chose this to be
##' the number of available CPUs
##' @param cl An optional parallel or snow cluster for use if \code{parallel="snow"}. If not supplied, a
##' cluster on the local machine is created for the duration of the call.
##' @return A \code{ScoreImputedSet} object
##' @seealso \code{\link{ScoreImputedSet.object}}
##' @details  Note that coxph may fail to converge and the following output
##' Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :  
##' Ran out of iterations and did not converge
##' 
##' It is possible to use ridge regression by including a ridge term in the model formula
##' (e.g. \code{~Z1+ridge(Z2,theta=1)}). See \code{\link[survival]{ridge}} for further details
##' @examples 
##' 
##' data(ScoreInd)
##' 
##' col.control <- col.headings(has.event="event", time="time",
##'                             Id="Id",arm="arm",
##'                             DCO.time="DCO.time", 
##'                             to.impute="to.impute")
##' 
##' ans <- ScoreImpute(data=ScoreInd,event.model=~Z1+Z2+Z3+Z4+Z5,
##'                    col.control=col.control, m=5,
##'                    bootstrap.strata=ScoreInd$arm,
##'                    NN.control=NN.options(NN=5,w.censoring = 0.2))
##' 
##' @export
ScoreImpute <- function(data,event.model,censor.model=event.model,
                        col.control, NN.control=NN.options(),time.dep=NULL,m,
                        bootstrap.strata=rep(1,nrow(data)),...,
                        parallel = c("no", "multicore", "snow")[1], ncpus = 1L, cl = NULL){
  
  #validate arguments
  validate.Score.Arguments(data,col.control,NN.control,time.dep,match.call(expand.dots = TRUE),m)
  validate.parallel.arguments(parallel,ncpus,cl)
  
  #setup censor model event type if not given
  if(is.null(col.control$censor.type)){
    col.control$censor.type <- "using_has.event_col"
    data[,col.control$censor.type] <- 1-data[,col.control$has.event]
  }
  
  #indices for bootstrap sample 
  boot.indices <- boot::boot(data=data,R=m,statistic=function(data,indices){return(indices)},strata=bootstrap.strata)
  boot.indices <- lapply(seq_len(m),function(x){boot.indices$t[x,]})

  #perform imputation
  if(parallel=="no"){
    imputes <- lapply(boot.indices,Sfn,data,event.model,
                      censor.model,col.control,
                      NN.control,time.dep,...)
  }
  else{
    imputes <- parallelRun(parallel,ncpus,cl,lapply.list=boot.indices,
                           FUN=Sfn,data,event.model,censor.model,col.control,
                           NN.control,time.dep,...)
  }
    
  #remove censor.type column if created earlier  
  if(col.control$censor.type=="using_has.event_col"){
    data[,col.control$censor.type] <- NULL
  }
  
  #package up result
  .createImputedSet(m,col.control,data,imputes,"Score",formula(paste("~",col.control$arm)))
}



