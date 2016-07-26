#This file containts functions associated with col.control and NN.options arguments
#for srisk score imputation

##' Specify the columns of the data frame required by score imputation method
##' @param arm column name which will contain the subject's treatment group
##' @param has.event column name which will contain whether 
##' the subject has an event (1) or not(0)
##' @param time column name of censoring/event time 
##' @param Id column name of subject Id
##' @param DCO.time column name of the time at which the subject would have been
##' censored had they not had an event before data cut off
##' @param to.impute column name of the logical column as to whether events should
##' be imputed
##' @param censor.type column name of the column containing the reason for censoring,
##' 0=had event, 1=non-administrative censoring 2=administrative censoring -- only subjects
##' with 1 in this column count as having an `event' in the Cox model for censoring 
##' (optionally used -- if not used then all subjects who are censored are used)
##' @return A list contain the given arguments 
##' @export
col.headings <- function(arm,has.event,time,Id,DCO.time,to.impute,censor.type=NULL){
  
  retVal <- list(arm=arm,
                 has.event=has.event,
                 time=time,
                 Id=Id,
                 DCO.time=DCO.time,
                 to.impute=to.impute,
                 censor.type=censor.type)
  
  .validcoloption(retVal)
  retVal
}

.validcoloption <- function(retVal){

  if(is.null(retVal) || !is.list(retVal)){
    stop("Invalid column control")
  }
  
  which.cols <- c("arm","has.event","time","Id","DCO.time","to.impute","censor.type")
  
  if(any(!which.cols %in% names(retVal))){
    stop("Invalid column control")
  }
  
  lapply(which.cols,function(name){
    if(name!="censor.type" && is.null(retVal[[name]])){
      stop("Invalid argument to col.headings only censor.type can be NULL")
    }
    if(!is.null(retVal[[name]]) && 
         (length(retVal[[name]])!=1|| !is.character(retVal[[name]]) || retVal[[name]]=="")){
        stop("Invalid argument to col.headings")  
    }
  })
}


.dummy.col.headings <- function(){
  col.headings(arm=" ",has.event=" ",time=" ",Id=" ",DCO.time=" ",to.impute=" ",censor.type = " ")
}
  

##' Create a list of options which control the nearest neighbour algorithm
##' for risk score imputation
##' @param NN The (maximum) number of subjects to be included in the
##' risk set
##' @param w.censoring The weighting on the censoring risk score when 
##' calculating distances for the nearest neighbour calculation
##' A weighting of \code{(1-w.censoring)} is used for the event risk score  
##' @param min.subjects If using time dependent score imputation include at least
##' this number of subjects when fitting the Cox model (i.e. include some subjects who were censored/had event
##' earlier than the cenosred observation if neccessary)
##' @return A list of options used within the ScoreImputedData
##' function
##' @export
NN.options <- function(NN=5,w.censoring=0.2,min.subjects=20){
  retVal <- list(NN=NN,
       w.censoring=w.censoring,
       min.subjects=min.subjects)
  
  .validNNoption(retVal)
  retVal
}

.validNNoption <- function(val){
  if(is.null(val) || !is.list(val)){
    stop("Invalid NNoption")
  }
  
  if(any(! c("NN","w.censoring","min.subjects") %in% names(val))){
    stop("NNoption must contain the entries NN, w.censoring and min.subjects")
  }
  
  if(!.internal.is.finite.number(val$w.censoring) || val$w.censoring < 0 || val$w.censoring >1){
    stop("Invalid argument w.censoring must be in [0,1]")
  }

  if(!.internal.is.finite.number(val$NN) || !.internal.is.wholenumber(val$NN) || val$NN <=0){
    stop("Invalid argument NN must be positive integer")
  }
  
  
  if(!.internal.is.finite.number(val$min.subjects) || !.internal.is.wholenumber(val$min.subjects) || val$min.subjects <=0){
    stop("Invalid argument min.subjects must be positive integer")
  }
}
