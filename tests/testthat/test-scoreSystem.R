context("ScoreSystem")

test_that("arguments_to_ScoreImpute",{
  data(ScoreInd)
  
  #first testing the arm and gamma arguments are valid in col.control for score imputation
  
  #invalid arm
  expect_error(ScoreImpute(data=ScoreInd,event.model=~Z1+Z2+Z3+Z4+Z5,
                              col.control=col.headings(has.event="event",
                                                       time="time",
                                                       Id="Id",
                                                       arm=23,
                                                       DCO.time="DCO.time",
                                                       to.impute="to.impute"),
                              NN.control=NN.options(NN=5,w.censoring = 0.2),
                              bootstrap.strata=ScoreInd$arm,m=5))
  
  #using gamma
  expect_error(ScoreImpute(data=ScoreInd,event.model=~Z1+Z2+Z3+Z4+Z5,
                              col.control=col.headings(has.event="event",
                                                       time="time",
                                                       Id="Id",
                                                       arm="arm",
                                                       gamma="arm",
                                                       DCO.time="DCO.time",
                                                       to.impute="to.impute"),
                              NN.control=NN.options(NN=5,w.censoring = 0.2),
                              bootstrap.strata=ScoreInd$arm,m=5))
  
  #missing arm
  expect_error(ScoreImpute(data=ScoreInd,event.model=~Z1+Z2+Z3+Z4+Z5,
                              col.control=col.headings(has.event="event",
                                                       time="time",
                                                       Id="Id",
                                                       DCO.time="DCO.time",
                                                       to.impute="to.impute"),
                              NN.control=NN.options(NN=5,w.censoring = 0.2),
                              bootstrap.strata=ScoreInd$arm,m=5))
  
  #Also testing:
  
  col.control <- col.headings(has.event="event", time="time",
                           Id="Id",arm="arm", DCO.time="DCO.time",
                           to.impute="to.impute")
  
  #use subset
  expect_error(ScoreImpute(data=ScoreInd,event.model=~Z1+Z2+Z3+Z4+Z5,
                           col.control=col.control,NN.control=NN.options(NN=5,w.censoring = 0.2),
                           subset=Z1==1,bootstrap.strata=ScoreInd$arm,m=5))
  
  #use NA.action
  expect_error(ScoreImpute(data=ScoreInd,event.model=~Z1+Z2+Z3+Z4+Z5,
                           col.control=col.control,NN.control=NN.options(NN=5,w.censoring = 0.2),
                           bootstrap.strata=ScoreInd$arm,m=5,na.action=na.exclude))
  
  #missing event model
  expect_error(ScoreImpute(data=ScoreInd,col.control=col.control,
                           NN.control=NN.options(NN=5,w.censoring = 0.2),
                           bootstrap.strata=ScoreInd$arm,m=5))

})

test_that("ScoreImputedDataOutput",{
  data(ScoreInd)
  
  ScoreInd$to.impute[1:40] <- FALSE
  set.seed(25)
  
  ans <- ScoreImpute(data=ScoreInd,event.model=~Z1+Z2+Z3+Z4+Z5,
                     col.control=col.headings(has.event="event", time="time",
                                  Id="Id",arm="arm", DCO.time="DCO.time",to.impute="to.impute"),
                     NN.control=NN.options(NN=5,w.censoring = 0.2),m=5,bootstrap.strata=ScoreInd$arm)
  
  ans <- ExtractSingle(ans,index=4)
  
  expect_equal("ScoreImputedData",class(ans))
  expect_equal(c("data","col.control","default.formula"),names(ans))
  
  expect_equal(col.headings(has.event="event",
                            time="time",
                            Id="Id",
                            arm="arm",
                            DCO.time="DCO.time",
                            to.impute="to.impute",
                            censor.type="using_has.event_col"),ans$col.control)
  
  expect_equal(c(colnames(ScoreInd),"impute.time","impute.event"),colnames(ans$data))
  expect_equal(nrow(ScoreInd),nrow(ans$data)) 
  
  df <- ans$data
  df$impute.event <- NULL
  df$impute.time <- NULL
  
  expect_equal(ScoreInd,df)
  
  df <- ans$data
  expect_true(all(df$impute.event[df$event==1]==1))
  expect_true(all(df$impute.event[!df$to.impute]==df$event[!df$to.impute]))
  expect_true(all(df$time[!df$to.impute]==df$impute.time[!df$to.impute]))
  
  expect_true(all(df$time <= df$impute.time & df$impute.time <= df$DCO.time))
  #if impute to DCO.time then don't have event
  expect_true(all(df$impute.event[df$time!=df$DCO.time & df$impute.time==df$DCO.time]==0  ))
})


test_that("algorithm_is_stochastic",{
  data(ScoreInd)
  set.seed(25)
  
  col.control <-col.headings(has.event="event",
                  time="time",Id="Id",
                  arm="arm", DCO.time="DCO.time", to.impute="to.impute")
  
  ans <- ScoreImpute(data=ScoreInd,event.model=~Z1+Z2+Z3+Z4+Z5,
                     col.control = col.control,m=5,bootstrap.strata = ScoreInd$arm,
                     NN.control=NN.options(NN=5,w.censoring = 0.2))
  
  set.seed(26)
  ans2 <- ScoreImpute(data=ScoreInd,event.model=~Z1+Z2+Z3+Z4+Z5,
                     col.control = col.control,m=5,bootstrap.strata = ScoreInd$arm,
                     NN.control=NN.options(NN=5,w.censoring = 0.2))
  
  expect_false(all(ExtractSingle(ans,index=1)$data$impute.time==ExtractSingle(ans2,index=1)$data$impute.time))
  
})


test_that("factor_numeric_character_Id_same_answer",{
  data(ScoreInd)
  set.seed(25)
  
  col.control <-col.headings(has.event="event",
                             time="time",Id="Id",
                             arm="arm", DCO.time="DCO.time", to.impute="to.impute")
  
  ans <- ScoreImpute(data=ScoreInd,event.model=~Z1+Z2+Z3+Z4+Z5,
                        col.control=col.control,m=5,bootstrap.strata=ScoreInd$arm,
                        NN.control=NN.options(NN=5,w.censoring = 0.2))
  
  set.seed(25)
  ScoreInd$Id <- as.character(ScoreInd$Id)
  ScoreInd$Id[1] <- "a"
  
  ans2 <- ScoreImpute(data=ScoreInd,event.model=~Z1+Z2+Z3+Z4+Z5,
                      col.control=col.control,m=5,bootstrap.strata=ScoreInd$arm,
                      NN.control=NN.options(NN=5,w.censoring = 0.2))
  
  expect_equal(ExtractSingle(ans,3)$data$impute.time,ExtractSingle(ans2,3)$data$impute.time)
  
  
  set.seed(25)
  ScoreInd$Id <- factor(ScoreInd$Id)
  ans3 <- ScoreImpute(data=ScoreInd,event.model=~Z1+Z2+Z3+Z4+Z5,
                      col.control=col.control,m=5,bootstrap.strata=ScoreInd$arm,
                      NN.control=NN.options(NN=5,w.censoring = 0.2))
  
  expect_equal(ExtractSingle(ans,3)$data$impute.time,ExtractSingle(ans3,3)$data$impute.time)
 
})

test_that("using_censor_type_gives_same_result_as_not_if_no_administrative_censoring",{
  data(ScoreInd)
  set.seed(125)
  
  ScoreInd$ctype <- 1 - ScoreInd$event 
  
  col.control <- col.headings(has.event="event",
                    time="time", Id="Id",
                    arm="arm", DCO.time="DCO.time",
                    to.impute="to.impute")
  
  ans <- ScoreImpute(data=ScoreInd,event.model=~Z1+Z2+Z3+Z4+Z5,m=5,
                          bootstrap.strata=ScoreInd$arm,
                          col.control=col.control,
                          NN.control=NN.options(NN=5,w.censoring = 0.2))
  
  set.seed(125)
  col.control$censor.type <- "ctype"
  ans2 <- ScoreImpute(data=ScoreInd,event.model=~Z1+Z2+Z3+Z4+Z5,m=5,
                     bootstrap.strata=ScoreInd$arm,
                     col.control=col.control,
                     NN.control=NN.options(NN=5,w.censoring = 0.2))
  
  #only difference should be col.control$censor.type
  expect_equal("ctype",ans2$col.control$censor.type)
  expect_equal("using_has.event_col",ans$col.control$censor.type)
  ans2$col.control$censor.type <- ans$col.control$censor.type
  expect_equal(ans,ans2)
  
})


test_that("Sfn_time_dep",{
  data(ScoreInd)
  data(ScoreTimeDep)
  set.seed(25)
  
  time.dep <- MakeTimeDepScore(ScoreTimeDep,Id="Id",
                             time.start="start",
                             time.end="end")
  
  #ok if time dep not used 
  #note do not get same answer without timedep as still use separate
  #model fits for each censored observation if timedep is not NULL
  expect_warning(ans <- ScoreImpute(data=ScoreInd,event.model=~Z1+Z3+Z5,
                        col.control=col.headings(has.event="event",
                                                 time="time",
                                                 Id="Id",
                                                 arm="arm",
                                                 DCO.time="DCO.time",
                                                 to.impute="to.impute"),
                        NN.control=NN.options(NN=5,w.censoring = 0.2),
                        time.dep = time.dep,m=5,
                        bootstrap.strata=ScoreInd$arm))
  
  #same answer if ID is factor (matching time.dep$ID which is always a factor)
  set.seed(25)
  
  ScoreInd$Id <- factor(ScoreInd$Id)
 
  expect_warning(ans2 <- ScoreImpute(data=ScoreInd,event.model=~Z1+Z3+Z5,
                                       col.control=col.headings(has.event="event",
                                                                time="time",
                                                                Id="Id",
                                                                arm="arm",
                                                                DCO.time="DCO.time",
                                                                to.impute="to.impute"),
                                       NN.control=NN.options(NN=5,w.censoring = 0.2),
                                       time.dep = time.dep,m=5,
                                     bootstrap.strata=ScoreInd$arm))
  expect_equal(ans$data$impute.time,ans2$data$impute.time)
  expect_equal(ans$data$impute.event,ans2$data$impute.event)
})

test_that("ordering_timeindep",{
  
  data(ScoreInd)
  set.seed(250)
  
  col.control <- col.headings(has.event="event",
                              time="time", Id="Id",
                              arm="arm", DCO.time="DCO.time",
                              to.impute="to.impute")
  
  #random reordering
  ScoreInd <- ScoreInd[sample(nrow(ScoreInd)),]
  
  #randomly set DCOs
  ScoreInd$DCO.time <- ScoreInd$DCO.time + runif(n=nrow(ScoreInd),min = 0,max=10)
  
  ans <- ScoreImpute(data=ScoreInd,event.model=~Z1+Z2+Z3+Z4+Z5,m=5,
                     bootstrap.strata=ScoreInd$arm,
                     col.control=col.control,
                     NN.control=NN.options(NN=5,w.censoring = 0.2))
  
  ans <- ExtractSingle(ans,index=2)
  
  expect_equal(ans$data$Id,ScoreInd$Id)
  expect_equal(ans$data$arm,ScoreInd$arm)
  
  expect_true(all(ans$data$impute.time >= ans$data$time))
  expect_true(all(ans$data$DCO.time>= ans$data$impute.time))
  
})


