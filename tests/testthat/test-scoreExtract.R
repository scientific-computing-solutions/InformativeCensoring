context("ScoreExtract")


test_that("ScoreImputeSet_extract",{
  set.seed(25)
  data(ScoreInd)
  
  ans <- ScoreImpute(data=ScoreInd,event.model=~Z1+Z2+Z3+Z4+Z5,
                     col.control=col.headings(has.event="event",
                                              time="time",
                                              Id="Id",
                                              arm="arm",
                                              DCO.time="DCO.time",
                                              to.impute="to.impute"),
                     NN.control=NN.options(NN=5,w.censoring = 0.2),
                     m=6,bootstrap.strata=ScoreInd$arm)
 
  expect_equal("ScoreImputedSet",class(ans))

  expect_equal(6,ans$m)
  
  expect_equal(col.headings(has.event="event",
               time="time",
               Id="Id",
               arm="arm",
               DCO.time="DCO.time",
               to.impute="to.impute",
               censor.type="using_has.event_col"),ans$col.control)
  
  expect_equal(ScoreInd,ans$data)
  
  expect_equal("matrix",class(ans$impute.time))
  expect_equal("matrix",class(ans$impute.event))
  expect_equal(400,nrow(ans$impute.time))
  expect_equal(6,ncol(ans$impute.event))
  
  
  expect_error(ExtractSingle(ScoreInd,index=5))
  expect_error(ExtractSingle(ans,index=4.5))
  expect_error(ExtractSingle(ans,index=7))
  expect_error(ExtractSingle(ans,index=0))
  
  my.data <- ExtractSingle(ans,index=5)
  expect_equal("ScoreImputedData",class(my.data))
  expect_equal(my.data$data$impute.time,ans$impute.time[,5])
  expect_equal(my.data$data$impute.event,ans$impute.event[,5])
  
  my.data <- ExtractSingle(ans,index=3)
  expect_equal(my.data$data$impute.time,ans$impute.time[,3])
  expect_equal(my.data$data$impute.event,ans$impute.event[,3])
})


test_that("ExtractSingle.Stat",{
  set.seed(25)
  data(ScoreInd)
  
  ans <- ScoreImpute(data=ScoreInd,
                     event.model=~Z1+Z2+Z3+Z4+Z5,
                     col.control=col.headings(has.event="event",
                                              time="time",
                                              Id="Id",
                                              arm="arm",
                                              DCO.time="DCO.time",
                                              to.impute="to.impute"),
                     NN.control=NN.options(NN=5,w.censoring = 0.2),
                     m=5,bootstrap.strata=ScoreInd$arm)
  fits <- ImputeStat(ans,method="Cox")
  
  onefit <- ExtractSingle(fits,index=3)
  
  onedata <- ExtractSingle(ans,index=3)
  expect_equal(onefit,ImputeStat(onedata,method="Cox"))
  
   
}) 