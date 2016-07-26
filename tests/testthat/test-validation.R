context("validation")

test_that("checkContiguous",{
  expect_error(checkContiguous(c(1,2,1)))
  expect_error(checkContiguous(c(2,2,1,1,1,3,3,3,2)))
  expect_error(checkContiguous(c("A","","A")))
  expect_error(checkContiguous(c(1,1,1,2,2,2,3,3,3,2,3,3,3)))
  expect_that(checkContiguous(c(1,2,3,4)), not(throws_error()))
  expect_that(checkContiguous(c(41)), not(throws_error()))
  expect_that(checkContiguous(c("A","A","C","C","C","B")), not(throws_error()))
})

test_that("checkPanelling_invalid",{
  make.df <- function(s,e){
    data.frame(time.start=s,time.end=e)
  }
  
  #negative start
  expect_error(checkPanelling(make.df(c(-5),c(4))))
  
  #not start at zero
  expect_error(checkPanelling(make.df(c(1),c(4))))
  
  #incorrect order
  expect_error(checkPanelling(make.df(c(0,2,1),c(1,3,2))))
  
  #more than one subject
  expect_error(checkPanelling(make.df(c(0,0),c(2,5))))
  
  #non-contiguous
  expect_error(checkPanelling(make.df(c(0,2,5),c(2,4,6))))

  #start and end mismatch
  expect_error(checkPanelling( make.df(c(0,10,20),c(10,15,30))))
    
  #invalid last end
  expect_error(checkPanelling(make.df(c(0,10,20,30),c(10,20,30,25))))
  
  #invalid interval of 0
  expect_error(checkPanelling(make.df(c(0,10,20,20),c(10,20,20,25))))
  
})

test_that("checkPanelling_valid",{
  make.df <- function(s,e){
    data.frame(time.start=s,time.end=e)
  }
  
  expect_that(checkPanelling(make.df(c(0),c(4))), not(throws_error()))
  expect_that(checkPanelling(make.df(c(0,2.5,5,7),c(2.5,5,7,10))), not(throws_error()))
  expect_that(checkPanelling(make.df(c(0,2.5,5,7),c(2.5,5,7,7.01))), not(throws_error()))
})


test_that(".getResponse",{
  expect_equal("y",.getResponse(formula(y~x)))
  expect_equal("Surv(x, y)",.getResponse(formula(Surv(x,y)~x+y+r*y)))
  expect_equal("w + y * z",.getResponse(formula(w+y*z~x)))
  expect_equal(0,length(.getResponse(formula(~x))))
})

test_that(".validRHSFormula",{
  #first without arm argument LHS must be empty
  expect_error(.validRHSFormula(formula(y~x)))  
  expect_that(.validRHSFormula(formula(~a+b+c)),not(throws_error()) )
  expect_error(.validRHSFormula(formula(~x+cluster(y))))
  expect_error(.validRHSFormula(formula(~tt(x)+y)))
  expect_that(.validRHSFormula(formula(~a+strata(b)+c)),not(throws_error()) )
  #if do have arm argument then it must be the first on the rhs
  #and no interaction terms with it
  expect_error(.validRHSFormula(formula(y~a),arm="a"))
  expect_that(.validRHSFormula(formula(~a+b+c),arm="a"),not(throws_error()) )
  expect_error(.validRHSFormula(formula(~b+arm),arm="arm"))
  expect_error(.validRHSFormula(formula(~b),arm="arm"))
  expect_error(.validRHSFormula(formula(~arm+b*arm),arm="arm"))
  expect_error(.validRHSFormula(formula(~b+b*arm),arm="b"))
})


test_that(".additionalScore.validate_control",{
  data(ScoreInd)
  
  col.control <- col.headings(has.event="event",
                              time="time",
                              Id="Id",
                              arm="arm",
                              DCO.time="DCO.time",
                              to.impute="to.impute")
  
  
  Call <- call("mycall",event.model=~Z1+Z2)
  
  #col.control not matching columns in data frame
  col.control$time <- "my.time"
  expect_error(.additionalScore.validate(ScoreInd,col.control=col.control,Call))
  
  col.control$time <- "time"
  col.control$gamma <- "BOO"
  expect_error(.additionalScore.validate(ScoreInd,col.control=col.control,Call)) 
  
  col.control$to.impute <- NULL
  expect_error(.additionalScore.validate(ScoreInd,col.control=col.control,Call))  

})

test_that(".additionalScore.validate",{
  col.control <- col.headings(has.event="event",
                              time="time",
                              Id="Id",
                              arm="arm",
                              DCO.time="DCO.time",
                              to.impute="to.impute")
 
 
  Call <- call("mycall",event.model=~Z1+Z2)
  #non-unique Id
  df <- data.frame(Id=c(1,6,9,1),event=c(0,0,1,1),time=c(4,5,6,7),
                   arm=factor(c(0,1,0,1)),DCO.time=c(5,6,7,8),to.impute=c(TRUE,TRUE,FALSE,FALSE))
  
  expect_error(.additionalScore.validate(df,col.control=col.control,Call))
  
  #negative time
  df <- data.frame(Id=c(1,6,9,21),event=c(0,0,1,1),time=c(-4,5,6,7),
                   arm=factor(c(0,1,0,1)),DCO.time=c(5,6,7,8),to.impute=c(TRUE,TRUE,FALSE,FALSE))
  
  expect_error(.additionalScore.validate(df,col.control=col.control,Call))
  
  #zero time
  df <- data.frame(Id=c(1,6,9,21),event=c(0,0,1,1),time=c(4,5,0,7),
                   arm=factor(c(0,1,0,1)),DCO.time=c(5,6,7,8),to.impute=c(TRUE,TRUE,FALSE,FALSE))
  
  expect_error(.additionalScore.validate(df,col.control=col.control,Call))
  
  #DCO.time < time
  df <- data.frame(Id=c(1,6,9,21),event=c(0,0,1,1),time=c(4,5,6,7),
                   arm=factor(c("A","B","A","B")),DCO.time=c(5,3,7,8),to.impute=c(TRUE,TRUE,FALSE,FALSE))
  
  expect_error(.additionalScore.validate(df,col.control=col.control,Call))
  
  #ok if DCO.time = time
  df <- data.frame(Id=c(1,6,9,21),event=c(0,0,1,1),time=c(4,5,6,7),
                   arm=factor(c("A","B","A","B")),DCO.time=c(4,5,6,7),to.impute=c(TRUE,TRUE,FALSE,FALSE))
  
  expect_that(.additionalScore.validate(df,col.control=col.control,Call),not(throws_error()))
  
  
  df$arm <- factor(c("A","B","A","B"))             
 
  #call has subset
  expect_error(.additionalScore.validate(df,col.control=col.control,
                                         Call=call("my.func",event.model=~Z1,subset="a")))
  
  #toimpute invalid
  df <- data.frame(Id=c(1,6,9,21),event=c(0,0,1,1),time=c(4,5,6,7),
                   arm=factor(c(1,0,0,1)),DCO.time=c(4,5,6,7),to.impute=c(7,TRUE,FALSE,FALSE))
  expect_error(.additionalScore.validate(df,col.control=col.control,Call))
  
  #event indicator incorrect
  df <- data.frame(Id=c(1,6,9,21),event=c(0,5,1,1),time=c(4,5,6,7),
                   arm=factor(c(1,0,0,1)),DCO.time=c(4,5,6,7),to.impute=c(TRUE,TRUE,FALSE,FALSE))
  expect_error(.additionalScore.validate(df,col.control=col.control,Call))
})


test_that("validate_Score_arguments_control",{
  data(ScoreInd)
  
  col.control <- col.headings(has.event="event",
                              time="time",
                              Id="Id",
                              arm="arm",
                              DCO.time="DCO.time",
                              to.impute="to.impute")
  
  Call <- call("mycall",event.model=~Z1+Z2)
  
  expect_error(validate.Score.Arguments(ScoreInd,col.control=col.control,NN.control=NULL,NULL,Call,m=5))
  
  NN.control <- c(10,20)
  expect_error(validate.Score.Arguments(ScoreInd,col.control=col.control,NN.control=NN.control,NULL,Call,m=5))
  NN.control <- list(NN=10)
  expect_error(validate.Score.Arguments(ScoreInd,col.control=col.control,NN.control=NN.control,NULL,Call,m=5))
  
  NN.control <- NN.options()
  expect_that(validate.Score.Arguments(ScoreInd,col.control=col.control,NN.control=NN.control,NULL,Call,m=5),not(throws_error()))
  
})


test_that("validate_Score_arguments_data_and_Call",{
  Call <- call("mycall",event.model=~Z1+Z2)
  col.control <- col.headings(has.event="event",
                              time="time",
                              Id="Id",
                              arm="arm",
                              DCO.time="DCO.time",
                              to.impute="to.impute")
  NN.control <- NN.options()
  
  
  df <- data.frame(Id=c(1,6,9,21),event=c(0,0,1,1),time=c(4,5,6,7),
                   arm=factor(c("A","B","A","B")),DCO.time=c(4,5,6,7),to.impute=c(TRUE,TRUE,FALSE,FALSE))
  
  expect_that(validate.Score.Arguments(df,col.control=col.control,NN.control=NN.control,NULL,
                                       Call=Call,m=5),not(throws_error()))
  
  #invalid m
  expect_error(validate.Score.Arguments(df,col.control=col.control,NN.control=NN.control,NULL,
                                        Call=Call,m=-5))
  
  expect_error(validate.Score.Arguments(df,col.control=col.control,NN.control=NN.control,NULL,
                                        Call=Call,m=c(3,4,5)))
  
  #m must be > 4
  expect_error(validate.Score.Arguments(df,col.control=col.control,NN.control=NN.control,NULL,
                                        Call=Call,m=4))
               
               
  
  #no event model in call
  expect_error(validate.Score.Arguments(df,col.control=col.control,NN.control=NN.control,NULL,
                                      Call=call("my.func"),m=5))
  
  #arm not a factor
  df <- data.frame(Id=c(1,6,9,21),event=c(0,0,1,1),time=c(4,5,6,7),
                   arm=c(1,0,1,1),DCO.time=c(4,5,"a6",7),to.impute=c(TRUE,TRUE,FALSE,FALSE))
  expect_error(validate.Score.Arguments(df,col.control=col.control,NN.control=NN.control,NULL,Call,m=5))
  
  #arm not two level factor
  df <- data.frame(Id=c(1,6,9,21),event=c(0,0,1,1),time=c(4,5,6,7),
                   arm=factor(c(1,1,1,1)),DCO.time=c(4,5,6,7),to.impute=c(TRUE,TRUE,FALSE,FALSE))
  expect_error(validate.Score.Arguments(df,col.control=col.control,NN.control=NN.control,NULL,Call,m=5))
  df <- data.frame(Id=c(1,6,9,21),event=c(0,0,1,1),time=c(4,5,6,7),
                   arm=factor(c(1,4,7,1)),DCO.time=c(4,5,6,7),to.impute=c(TRUE,TRUE,FALSE,FALSE))
  expect_error(validate.Score.Arguments(df,col.control=col.control,NN.control=NN.control,NULL,Call,m=5))
})


test_that("Invalid_censortype",{
  Call <- call("mycall",event.model=~Z1+Z2)
  col.control <- col.headings(has.event="event",
                              time="time",
                              Id="Id",
                              arm="arm",
                              DCO.time="DCO.time",
                              to.impute="to.impute",
                              censor.type="ctype")
  NN.control <- NN.options()
  
  
  df <- data.frame(Id=c(1,6,9,21),event=c(0,0,1,1),time=c(4,5,6,7),
                   arm=factor(c("A","B","A","B")),DCO.time=c(4,5,6,7),to.impute=c(TRUE,TRUE,FALSE,FALSE))
  
  
  #invalid column name
  df2 <- df
  df2$using_has.event_col <- rep(1,nrow(df))
  expect_error(validate.Score.Arguments(df2,col.control=col.control,NN.control=NN.control,NULL,Call=Call,m=5) ) 
  
  #column contains something other than 0,1 or 2
  df$ctype <- c(4,1,0,1)
  expect_error(validate.Score.Arguments(df,col.control=col.control,NN.control=NN.control,NULL,Call=Call,m=5) ) 
  
  df$ctype <- c(-1,"hello",1,1)
  expect_error(validate.Score.Arguments(df,col.control=col.control,NN.control=NN.control,NULL,Call=Call,m=5) ) 
  
  #error if have event and censor type != 0
  df$ctype <- c(1,0,2,1)
  expect_error(validate.Score.Arguments(df,col.control=col.control,NN.control=NN.control,NULL,Call=Call,m=5) ) 
  
  #and error if do not have event and censor type = 0
  df$ctype <- c(0,0,0,1)
  expect_error(validate.Score.Arguments(df,col.control=col.control,NN.control=NN.control,NULL,Call=Call,m=5) ) 
  
})

test_that("validate_Score_arguments_timedep",{
  Call <- call("mycall",event.model=~Z1+Z2)
  col.control <- col.headings(has.event="event",
                              time="time",
                              Id="Id",
                              arm="arm",
                              DCO.time="DCO.time",
                              to.impute="to.impute")
  NN.control <- NN.options()
  
  df <- data.frame(Id=c(1,6,9,21),event=c(0,0,1,1),time=c(14,15,16,17),
                   arm=factor(c("A","B","A","B")),DCO.time=c(24,25,26,27),to.impute=c(TRUE,TRUE,FALSE,FALSE))
  
  
  time.dep.df <- data.frame(Id=c(1,1,6,6,6,21,21,9),
                            time=c(0,5,0,2,8,0,10,0),
                            end=c(5,14,2,8,15,10,17,16),
                            W1=c(1,2,3,4,1,2,2,1))
  
  time.dep <- MakeTimeDepScore(time.dep.df,Id="Id",time.start="time",time.end="end")
  
  expect_error(validate.Score.Arguments(df,col.control=col.control,NN.control=NN.control,time.dep=time.dep.df,Call,m=5))
  
  #ok if both columns have Id
  expect_that(validate.Score.Arguments(df,col.control=col.control,NN.control=NN.control,time.dep=time.dep,Call,m=5),
                             not(throws_error()))
  
  #invalid if both have same column names 
  df$W1 <- c(2,3,4,5)
  expect_error(validate.Score.Arguments(df,col.control=col.control,NN.control=NN.control,time.dep=time.dep,Call,m=5))
  
  df$W1 <- NULL
  time.dep$arm <- c(1,1,1,1,0,0,0,0)
  expect_error(validate.Score.Arguments(df,col.control=col.control,NN.control=NN.control,time.dep=time.dep,Call,m=5))
  
})

