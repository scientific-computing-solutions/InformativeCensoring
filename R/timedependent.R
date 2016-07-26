#This file contains code associated with the time dependent risk score imputation
#method, the ScoreTD object (data frame containing the time dependent covariates)
#

##' A \code{ScoreTD} object
##' 
##' This data frame holds time dependent covariates for
##' use with risk score imputation
##' 
##' The data frame contains the following columns:
##' 'Id' for subject ID \cr
##' 'time.start' and 'time.end' the range of time for which 
##' the covariate values are valid - i.e. [time.start,time.end] \cr 
##' Additional columns are the time dependent covariates 
##'  
##' All data for a single subject should be stored in consecutive rows, sorted 
##' by time and the starting time of a row should match the ending time of the previous row  
##'     
##' @name ScoreTD.object
##' @seealso \code{\link{MakeTimeDepScore}}
NULL


##' Create a valid \code{ScoreTD} object
##' @param data A data frame of time dependent covariates
##' @param Id The column name of the subject Id 
##' @param time.start The covariates are valid for the time [time.start,time.end] where
##' time.start is the column name of time.start
##' @param time.end The covariates are valid for the time [time.start,time.end] where
##' time.end is the column name of time.end
##' @return A \code{ScoreTD} object
##' @export
MakeTimeDepScore <- function(data,Id,time.start,time.end){
  if(any(!c(Id,time.start,time.end) %in% colnames(data))){
    stop("Invalid arguments, Id, time.start and time.end should be column names in the data frame ")
  }
  
  sort.cols <- function(data,new.name,given.name){
    if(given.name != new.name){
      if(new.name %in% colnames(data)) stop("The column",new.name,"already exists!")
      data[,new.name] <- data[,given.name]
      data[,given.name] <- NULL
    }
    data
  }
  
  data <- sort.cols(data,"Id",Id)
  data <- sort.cols(data,"time.start",time.start)
  data <- sort.cols(data,"time.end",time.end)
 
  data$Id <- factor(data$Id)
  data$time.start <- as.numeric(data$time.start)
  data$time.end <- as.numeric(data$time.end)
  
  checkContiguous(data$Id)
  by(data,data$Id,checkPanelling)
  
  class(data) <- c("ScoreTD",class(data))
  data
}


#For a give data frame data and ScoreTD object, time.dep, 
#output a merged data frame containing data and the time dependent covariate
#values at time my.time if my.time is NULL then output the time dependent covariates
#at the time of subject having an event/being censored
.getTimeDepDataSet <- function(data,time.dep,dataID,datatime,my.time){
  
  #if have no subjects (as current subject has max time, then output empty data frame)
  #we only fit a model if the number of rows of data frame > NN, so 
  #this edge case behaves as expected
  if(nrow(data)==0){ 
    return(data)
  }
  
  #when merging and matching need to be careful about factor IDs not being matched correctly
  time.dep$Id <- as.character(time.dep$Id)
  data$internal_char_IDXX <-as.character(data[,dataID])
  
  
  time.dep <- time.dep[time.dep$Id %in% data$internal_char_IDXX,] 
  time.dep$internal.Score.time <- if(is.null(my.time)) vapply(time.dep$Id,function(x){
                        data[which(data$internal_char_IDXX==x),datatime]  
                        },numeric(1)) else my.time 
  
  retVal <- time.dep[time.dep$time.start < time.dep$internal.Score.time 
                     & time.dep$time.end >= time.dep$internal.Score.time,]
  
  s <- sort(unique(retVal$Id))
  s2 <- sort(unique(data$internal_char_IDXX)) 
  
  if(length(s) != length(s2) || any(s != s2)){
    msg <- if(is.null(my.time)) " the value of a covariate for a subject at the time of event/censoring" 
           else paste(" the value of a covariate for an at risk subject at time",my.time) 
    stop("Error in creating data set with time dependent covariates. Cannot determine", msg) 
  } 
  retVal$internal.Score.time <- NULL
  retVal$time.start <- NULL
  retVal$time.end <- NULL
  
  
  data <- inner_join(data,retVal,by=c("internal_char_IDXX"="Id"))
  data$internal_char_IDXX  <- NULL
  data
}

