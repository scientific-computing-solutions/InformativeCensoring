#This file contains the functions to produce a single imputed 
#data frame using risk score imputation - the main entry point is
#function Sfn below which is called from ScoreImpute in scoreImputedData.R

#The Sfn function splits the data frame by arm and each arm is imputed separately using
#the .imputeTimes function. This is where the control flow becomes more complex. See
#comments in the @details section of .imputeTimes for details.


#Function which performs the imputation for a single imputed data set
#@inheritparams ScoreImpute
#@param indices The indicies of the rows of the data set used to get the bootstrapped data for fitting the Cox model
#@return A list of imputation and event times
Sfn <- function(indices,data,event.model,censor.model,col.control,NN.control,time.dep,...){
  
  #get imputation for both arms separately
  imputes <- by(data,data[,col.control$arm],function(x){
    .imputeTimes(x,data[indices,] ,event.model,censor.model,
                 col.control,NN.control,time.dep=time.dep,...)},simplify=FALSE)
  
  imputed.df <- do.call(rbind,imputes)
  retVal <- inner_join(data,imputed.df,by=setNames(col.control$Id,col.control$Id))
  
  list(impute.time=retVal$impute.time,
       impute.event=retVal$impute.event)
}



# Function which performs the imputation for a single arm
# @inheritParams ScoreImpute
# @param data All the subjects for a given arm from the ScoreImpute data frame
# @param boot.data The bootstrapped data frame (both arms )
# At this point all subjects in the data frame are in the same treatment group
# @return A data frame with 3 columns, subject ID, impute.time,impute.event
# which will be merged into the input data frame in ScoreImputedData
# @details The control flow of this function is reasonably complex (due to the complex
# nature of the algorithm) so the following is performed:
#
# 1) subset the bootstrapped data set ot only include the treatment group we are considering
##########################
# For time independent case
# 2) use .calcscores to fit Cox models and return a data frame of risk scores
# and a list of the event and censor model fits using the bootstrapped data (this list
# is called modelfits)
# 3) set df = data, it is to be the data frame from which the risk set is chosen 
#(in time dep case df is not data)
###########################
# For time dependent case
# 2') use .getTimeDepDataSet to attach the time dependent covariates onto the data frame
# (not boot.data) taking the values at which 
# 3') calculate the maximum time we are going to use to subset the data when fitting the models
# to ensure there are at least NN.control$min.subjects subjects in boot.data frame when fitting
# the Cox models this  uses the .getLatestTimeCutoff function
# ##########################
# The for each subject (i.e. row in 'data' data frame) we generate an imputed event time
# and indicator in the following way
# 4) If already had event or not imputing then leave as is
# 5) In timedep case then subset the boot.data frame to subjects who are still on the trial
# beyond this subject's censor time (or cur at maximum time to enusre NN.control$min.subjects is satisfied)
# and as in "2)" use .calcscores to return a data frame of risk scores and model fits (this list is called modelfits)
# and a list of the event and censor model fits using the bootstrapped data
# 6) In timedep case set df=the data frame used in 5). Now we have a df and modelfits whether we used time
# dependent covariates or not
# 7) Use .calcDistances to calcuate the weighted distance from our given subject we are trying to impute for to all
# subjects in data frame df 
# 8) Use .getRiskSet to get the risk set for our subject
# 9) Use .kmi to perform the Kaplan-Meier imputation
################################
# 10) package up the outputs from 9 into a data frame, adding subject Id so that the data can be merged
# with the other treatment arm in function Sfn above 
.imputeTimes <- function(data,boot.data,event.model,censor.model,col.control,NN.control,time.dep,...){
  
  #get only the boot data from the arm we are imputing for
  boot.data <- boot.data[boot.data[,col.control$arm]==data[1,col.control$arm],]
  
  #if no time dependent then perform the model fits once here
  if(is.null(time.dep)){
    modelfits <- .calcscores(list(event.model,censor.model),boot.data,
                            time=col.control$time,has.event=col.control$has.event,
                            has.censoring=col.control$censor.type,
                            NN=NN.control$NN,...)
    
    #df is the data frame from which the risk sets are extracted from 
    df <- boot.data
    
  }
  else{ #otherwise add the values of the time dep. covariates at time of event/censoring to the data frame
    data <- .getTimeDepDataSet(data,time.dep=time.dep,col.control$Id,col.control$time,my.time=NULL)
    #when fitting the Cox models, here we calculate a time such that, if a subject has a censored
    #time > latest.time.cutoff we subset the data at latest.time.cutoff instead to ensure
    #there are at least min.subjects for the model fit.
    latest.time.cutoff <- .getLatestTimeCutoff(boot.data[,col.control$time],NN.control$min.subjects)
  }
 
  #now for each subject we create imputed data times
  kmimputes <- lapply(seq_len(nrow(data)),function(x){
    
    #current time
    my.time <- data[x,col.control$time]
    
    #if have event or not imputing then do not impute
    if(data[x,col.control$has.event]==1 || !data[x,col.control$to.impute]){
      return(list(event=data[x,col.control$has.event],time=my.time))
    }
   
    #time dependent case here - we need to fit individual model
    if(!is.null(time.dep)){
      #data frame for model fit and used for selecting risk set  
      #note we subset at min(my.time,latest.time.cutoff)
      df <- .getTimeDepDataSet(boot.data[boot.data[,col.control$time] > min(my.time,latest.time.cutoff),],time.dep,
                               col.control$Id,col.control$time,min(my.time,latest.time.cutoff))
      
      modelfits <- .calcscores(list(event.model,censor.model),df,
                              time=col.control$time,has.event=col.control$has.event,
                              has.censoring=col.control$censor.type,
                              NN=NN.control$NN,...)
    }
    
    #calculate the distances and hence risk set from current subject (data[x,]) to the subjects
    #used to fit the model (scores are in modelfits$raw.scores)
    my.dist <- .calcDistances(modelfits,NN.control$w.censoring,data[x,])
    risk.set.df <- .getRiskSet(my.dist,df[,col.control$time],my.time,NN.control$NN,df[,col.control$has.event])
    
    #apply Kaplan-Meier imputation
    .kmi(risk.set.df,
         my.time=my.time,
         my.DCO.time=data[x,col.control$DCO.time])
  })
  
  kmimputes <- data.frame(impute.event=vapply(kmimputes,"[[","event",FUN.VALUE = numeric(1)),
                          impute.time=vapply(kmimputes,"[[","time",FUN.VALUE = numeric(1)))  
  kmimputes[,col.control$Id] <- data[,col.control$Id]
  
  return(kmimputes)
  
}

# Perform the model fit and score calculation for data from a single arm
#@param models this is list(event.model,censor.model)
#@param data the data for the model fit, in the time dependent case this will be
#subjects in the same arm with still at risk at a given subject's censor time
# @param time, has.event, has.censoring - the column names of the data frame for the Surv object
# @param NN the number of subjects to be included in risk sets
# @param ... additional arguments to be passed into coxph
# @return a list with 2 objects: scores a data frame of risk scores
# and model a list of the event and censor model fits
.calcscores <- function(models,data,time,has.event,has.censoring,NN,...){
  
  model.fits <-  if(nrow(data) > NN) #only need model if more than NN subjects
                   mapply(.fitPHmodel,models,c(TRUE,FALSE), 
                          MoreArgs=c(list(data=data,time=time,has.event=has.event,has.censoring=has.censoring),list(...)),
                          SIMPLIFY=FALSE)
                  else list(NULL,NULL) 
  
  lin.comb <- function(model.fit){
    if(is.null(model.fit)) return(rep(0,nrow(data)))
    #Note do not need to include the constant term as we are normalizing
    ans <- predict(model.fit,type="lp")#+attr(predict(model.fit,type="terms"),"constant")
    if(length(ans)!=nrow(data)){
      stop("predict.coxph is not handling duplicate rows correctly please change implementation")
    }
    ans
  }
  
  return(list(raw.scores=data.frame(Rs.f=lin.comb(model.fits[[1]]),
                                    Rs.c=lin.comb(model.fits[[2]])),
              model=model.fits))
}

# Fit a model for calculating the risk scores using
# @param formula as the model for coxph
# @param event a logical flag for whether this is event (TRUE) or censoring 
# model (FALSE)
# @param data the data frame for Coxph model
# @param time, has.event,has.censoring - the column names of the data frame for the Surv object
# @param ... additional arguments to be passed into coxph (will not be na.action or subset they have been caught higher up)
# @return a Cox model fit or NULL if no model could be fitted
# a vector of raw scores
.fitPHmodel <- function(formula,event,data,time,has.event,has.censoring,...){
  .validRHSFormula(formula)
  formula <- update(formula,paste("Surv(",time,",",if(event) has.event else paste(has.censoring,"==1"),") ~."))
  
  #if no events/censoring
  if((event && all(data[,has.event]==0)) || (!event && all(data[,has.censoring]!=1))){
    return(NULL)
  }
  return(coxph(formula,droplevels(data),model=TRUE,na.action=na.fail,...))
}



#given a data frame of scores model.fits$scores (column names Rs.f and Rs.c) and the 
#weighting for the censored score (Rs.c) output a distance matrix for subject given in the row of
#a data frame, data.row (using model.fits$model)
.calcDistances <- function(model.fits,w.censoring,data.row){
  
  #get scores for subject whose row is data.row
  my.scores <- vapply(model.fits$model,function(x){
    if(is.null(x)) return(0)
    predict(x,type="lp",newdata = data.row)#+attr(predict(x,type="terms"),"constant")
  },FUN.VALUE = numeric(1))
  
  my <- data.frame(Rs.f=my.scores[1],Rs.c=my.scores[2])
  #add my scores to the scores calculated from the model fit data
  .internalDistances(rbind(model.fits$raw.scores,my),w.censoring)
}  

#raw.scores contains the data frame of unnormalized scores,
#the final row of which contains the scores of the subject whos 
#risk set we are calculating
.internalDistances <- function(raw.scores,w.censoring){
  scores <- normalize.scores(raw.scores)
  
  scores$Rs.f <- scores$Rs.f*sqrt(1-w.censoring)
  scores$Rs.c <- scores$Rs.c*sqrt(w.censoring)
  
  #The last entry is the scores for the data point of interest
  index <- nrow(scores)
  
  distances <- sqrt((scores$Rs.f - scores$Rs.f[index])^2 +
         (scores$Rs.c - scores$Rs.c[index])^2)
  
  return(distances[1:(index-1)]) #do not need distances[index] (=0 as point of interest)
}  


normalize.scores <- function(raw.scores){
  .normfunc <- function(x){
    y <- x[1:(length(x)-1)] #not including data point of interest in
                          #calculation of mean and variance for normalization
    if(is.na(var(y)) ||var(y)==0) return(rep(0,length(x)))
    (x-mean(y))/sqrt(var(y))   
  }
  data.frame(Rs.f=.normfunc(raw.scores$Rs.f),Rs.c=.normfunc(raw.scores$Rs.c))
}

#Output the indices of the risk set
#@param distances A vector of distances one for each subject in consideration
#@param times A vector of times  one for each subject in consideration
#@param my.time The time of censoring 
#@param NN Number of subjects to be included in the risk set
#@param has.event A vector of event indicators one for each subject in consideration
#@return A data frame containing the time and event indicator for a given risk set
.getRiskSet <- function(distances,times,my.time,NN,has.event){
  possible.set <- which(times > my.time)
  risk.set.index <- if(length(possible.set)<=NN) possible.set
                    else intersect(which(distances <=  sort(distances[possible.set],partial=NN)[NN]),possible.set) 
  data.frame(times=times[risk.set.index],
             has.event=has.event[risk.set.index])
}

#perform Kaplan-Meier imputation for a single subject S
#it is already assumed at this point that S was censored and their event time is to
#be imputed
#@param risk.set.df A data frame with columns `times' and `has.event' giving the times
#and event indicators fo the risk set
#@param my.time, subject's current censored time
#@param my.DCO.time The time at which S is to be censored if imputed time is large
#@return A list with two elements, event and time: the imputed event indicator and
#time for S
.kmi <- function(risk.set.df,my.time,my.DCO.time){
  
  if(nrow(risk.set.df)==0){
    return(list(event=0,time=my.time))
  }
  
  km <- survfit(Surv(times,has.event)~1,data=risk.set.df)
  impute.time <- quantile(km,probs=runif(n=1))$quantile
  names(impute.time) <- NULL
  event <- as.numeric(!is.na(impute.time) && impute.time < my.DCO.time)
  
  if(is.na(impute.time)){
    impute.time <- max(risk.set.df$times)  
  }
  
  return(list(event=event,
              time=min(impute.time,my.DCO.time)))
}


#Given a vector of times and an integer variable min.subjects
#find the min.subjects largest time, remove all values of time >= this value
#and return the largest time remaining. 
.getLatestTimeCutoff <- function(times,min.subjects){
  N <- 1+max(0,length(times)-min.subjects)
  if(N==1)return(0)
  max(times[times<sort(times,partial=N)[N]])  
}
