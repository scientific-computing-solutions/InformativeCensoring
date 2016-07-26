#This file contains the code to validate the riskscore imputation arguments and
#the formula arguments and various subfunctions

#check the formula is valid (specifically there is no left hand side, should not be tt or cluster
#and if arm is not NULL then the 'arm' is the first term on the right hand side and there are no
#interactions between arm and other covariates)
.validRHSFormula <- function(formula,arm=NULL){
  
  if(class(formula)!="formula") stop("Invalid formula is not of type formula")
  
  if(length(.getResponse(formula))!=0){
    stop("formula cannot have any variables on the left hand side.")
  }
  
  tms<-terms(formula,specials=c("strata","cluster","tt"))
  
  if(length(untangle.specials(tms,"tt")$vars)>0 ||
     length(untangle.specials(tms,"cluster")$vars)>0){
    stop("Cannot use tt or cluster in formula for gamma/score imputation")
  }
  
  
  if(is.null(arm)) return(NULL)
  
  #validate the right hand side of the formula
  k <- attr(terms(formula),"term.labels")
  if(length(k)==0){
    stop("Empty formula!")
  }
  
  first_term <- k[[1]]
  if(first_term!=arm){
    stop("The first term of the formula must be the treatment group")
  }
  
  if(length(k)>1){
    covariates <- tail(k,-1)
    
    ans <- unlist(lapply(covariates,function(x){
      first_term %in% unlist(strsplit(x,split=c(":")))
    }))
    
    if(any(ans)){
      stop(error=paste("Model formula cannot include interactions between",arm, "and covariates.",
                       "The", arm, "must be the first model variable"))
    }
  }
  
}


#from http://stackoverflow.com/questions/13217322/how-to-reliably-get-dependent-variable-name-from-formula-object
.getResponse <- function(formula) {
  tt <- terms(formula)
  vars <- as.character(attr(tt, "variables"))[-1] ## [1] is the list call
  response <- attr(tt, "response") # index of response var
  vars[response] 
}

#check that the columns [time.start,time.end] for a single
#subject are ordered correctly, the initial time.start = 0, 
#time.end > time.start and there are no gaps (i.e. [0,1], [2,10])
#@param data a data frame with the time dep. covariates for a
#single subject
checkPanelling <- function(data){
  if(data$time.start[1] != 0){
    stop("The time.start for the first row for subject",data$Id[1],"must be 0")
  }
  if(nrow(data)>1 && !all.equal(data$time.start[2:nrow(data)],data$time.end[1:(nrow(data)-1)])){
    stop("There cannot be gaps in the [time.start,time.end) intervals, check subject",data$Id[1])      
  }
  if(any(data$time.end <= data$time.start)){
    stop("The time.end cannot be <= time.start, check subject",data$Id[1])
  }
}


# Check duplicated elements of a vector are contiguous
# 
# for example c(1,1,1,2,2,4,4,5,5) is contiguous
# but c(1,1,2,2,1,4,5,4,5) is not
# @param vec a vector of values
checkContiguous <- function(vec){
  if(length(vec)<2) return(TRUE)
  known.vals <- vec[1]
  for(i in 2:length(vec)){
    if(vec[i] != vec[i-1]){
      if(!vec[i] %in% known.vals){
        known.vals <- c(known.vals,vec[i])
      }
      else{
        stop(paste("The order of the data frame is incorrect. All rows with the same",
                   "Subject ID must be contiguous"))
      }
    }
  }
}

is.valid.call <- function(Call){
  
  indx <- match(c("subset", "na.action"), names(Call),nomatch=0L)
  if(any(indx[1:2]!=0)){
    stop("Cannot use na.action or subset argument with this imputation ",
         "function (na.fail is used when fitting the Cox models)")
  }
}


.additionalScore.validate <- function(data,col.control,Call){
  
  is.valid.call(Call)
  
  if(any(c("impute.time","impute.event") %in% colnames(data))){
    stop("Cannot use a data frame with columns impute.time or impute.event")
  } 
  
  #validate col.control (note this is needed here in case people do not use col.headings function) 
  .validcoloption(col.control)
  
  #validate col.control with data
  lapply(col.control,function(x){
    if(!is.null(x) && !x %in% colnames(data))stop("Invalid column name '",x,"'not found in data frame")})
  
  if(nrow(data)==0){
    stop("Empty data frame!")
  }
  
  if(any(!is.numeric(data[,col.control$time]))){
    stop("time on study must be numeric")
  }
  
  if(any(!is.numeric(data[,col.control$DCO.time]))){
    stop("DCO time must be numeric")
  }
  
  #check time is positive
  if(any(data[,col.control$time]<= 0)){
    stop("Time on study must be positive")
  }
  
  #DCO.time is <= time
  if(any(data[,col.control$DCO.time] < data[,col.control$time])){
    stop("DCO.time must be >= time for all subjects ")
  }
  
  #unique IDs
  if(nrow(data)!=length(unique(data[,col.control$Id]))){
    stop("Subject Ids must be unique")
  }
  
  if(any(!data[,col.control$has.event] %in% c(0,1))){
    stop("Event indicator column can only contain 0s and 1s")
  }
  
  if(any(!is.logical(data[,col.control$to.impute]))){
    stop("To impute column can only contain TRUE and FALSE")
  }
  
}


validate.Score.Arguments <- function(data,col.control,NN.control,time.dep,Call,m){
  
  #validation routines additional
  .additionalScore.validate(data,col.control,Call=Call)

  #check m
  if(!.internal.is.finite.number(m) ||!.internal.is.wholenumber(m) || m < 5){
    stop("m must be at least 5")
  }
  
  #check have event.model
  indx <- match(c("event.model"), names(Call),nomatch=0L)
  if(indx[1]==0){
    stop("Missing event.model argument")
  }
  
  #event.model and censor.model are validated later
  #validate NN.control  (note this is needed here in case people do not use NN.options function)
  .validNNoption(NN.control)

  #check arm is two level factor
  if(!is.factor(data[,col.control$arm])|| length(levels(data[,col.control$arm]))!=2){
    stop("Treatment group must be a two level factor variable")
  }
  
  #censor.type validation
  if("using_has.event_col" %in% colnames(data)){
    stop("Cannot use a data frame with column name using_has.event_col")
  }
  
  if(!is.null(col.control$censor.type)){
    if(any(!data[,col.control$censor.type] %in% 0:2)){
      stop("Censor type must be 0, 1 or 2")
    }
    if(any(data[data[,col.control$has.event]==1,col.control$censor.type]!=0)){
      stop("Censor type for subjects who have event must=0")
    }
    if(any(data[data[,col.control$has.event]==0,col.control$censor.type]==0)){
      stop("Censor type for subjects who do not have event cannot be=0")
    }
  }
  
  
  if(!is.null(time.dep)){
    #timedep class
    if(!"ScoreTD" %in% class(time.dep)){
      stop("time.dep must be of type ScoreTD")  
    }
    
    ans <- intersect(colnames(time.dep),colnames(data))
    if(length(ans)>1 || (length(ans)==1 && ans[1] != "Id") ){
      stop("Cannot have columns with the same name in data and time.dep data frames")
    }
    
  } 

}  

