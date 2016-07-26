context("internalScore")

test_that("getLatestTimeCutoff",{
  
  expect_equal(0,.getLatestTimeCutoff(c(9,1.5,7,10,12),5))
  expect_equal(10,.getLatestTimeCutoff(c(9,1.5,7,10,12),1))
  expect_equal(7,.getLatestTimeCutoff(c(9,9,9,11,1.5,7,10,12),4))
  expect_equal(7,.getLatestTimeCutoff(c(9,9,9,11,1.5,7,10,12),5))
  expect_equal(7,.getLatestTimeCutoff(c(9,9,9,11,1.5,7,10,12),6))
  expect_equal(1.5,.getLatestTimeCutoff(c(9,9,9,11,1.5,7,10,12),7))
})

test_that(".kmi_deterministiccases",{
  #empty risk set
  df <- data.frame(times=numeric(0),has.event=numeric(0))
  
  ans <- .kmi(df,10,20)
  expect_equal(list(event=0,time=10),ans)
  
  #risk set with no events
  df <- data.frame(times=c(20,30,10),has.event=c(0,0,0))
  ans <- .kmi(df,5,40)
  expect_equal(list(event=0,time=30),ans)
  ans <- .kmi(df,5,15)
  expect_equal(list(event=0,time=15),ans)
  
  #risk set with 1 subject with event
  df <- data.frame(times=20,has.event=1)
  ans <- .kmi(df,5,15)
  expect_equal(list(event=0,time=15),ans)
  ans <- .kmi(df,5,30)
  expect_equal(list(event=1,time=20),ans)
  
  #risk set with 1 subject with no events
  df <- data.frame(times=40,has.event=0)
  ans <- .kmi(df,5,15)
  expect_equal(list(event=0,time=15),ans)
  ans <- .kmi(df,5,50)
  expect_equal(list(event=0,time=40),ans)
  
  #risk set 1 event, last one
  df <- data.frame(times=c(40,20,30),has.event=c(1,0,0))
  ans <- .kmi(df,5,15)
  expect_equal(list(event=0,time=15),ans)
  ans <- .kmi(df,5,50)
  expect_equal(list(event=1,time=40),ans)
})


test_that("kmi_stochastic",{
  df <- data.frame(times=c(10,20,30,40,50),has.event=c(1,1,0,1,1))
  
  set.seed(10)
  #get the 20 random U which will be used below
  Us <- runif(n=20) 
  expect_vals <- vapply(Us,function(x){
    if(x<0.2) return(10)
    if(x<0.4) return(20)
    if(x<0.7) return(40)
    return(50)
  },FUN.VALUE = numeric(1))
  
  #reset the seed back
  set.seed(10)
  ans <- vapply(1:20,function(x){.kmi(df,5,60)$time},FUN.VALUE = numeric(1))
  
  expect_equal(expect_vals,ans)
  
})


test_that("getriskset",{
  
  distances <- c(10,5,12,1.5,5,10,7)
  times <- c(10,20,30,40,50,60,5)
  has.event <- c(1,0,0,1,1,0,1)
  
  #test if NN is large to include everything
  rs <- .getRiskSet(distances,times,my.time=1,NN=10,has.event)
  expect_equal(data.frame(times=times,has.event=has.event),rs)
  
  #Only include subjects who have time strictly > my.time
  rs <- .getRiskSet(distances,times,my.time=5,NN=10,has.event)
  expect_equal(data.frame(times=c(10,20,30,40,50,60),
                          has.event=c(1,0,0,1,1,0)),rs)
  #Empty data set
  rs <- .getRiskSet(distances,times,my.time=60,NN=10,has.event)
  expect_equal(data.frame(times=numeric(0),
                          has.event=numeric(0)),rs)
  
  #NN are output if possible
  rs <- .getRiskSet(distances,times,my.time=5,NN=5,has.event)
  expect_equal(data.frame(times=c(10,20,40,50,60),
                          has.event=c(1,0,1,1,0)),rs)
 
  #handling ties as described in vignette
  rs <- .getRiskSet(distances,times,my.time=5,NN=4,has.event)
  expect_equal(data.frame(times=c(10,20,40,50,60),
                          has.event=c(1,0,1,1,0)),rs)
  
  #extreme NN (=1)
  rs <- .getRiskSet(distances,times,my.time=40,NN=1,has.event)
  expect_equal(data.frame(times=c(50),
                          has.event=c(1)),rs)
})

test_that(".fitPHmodel",{
  my.data <- data.frame(
    my.time=c(10,20,80,40,50,50,70,90),
    my.event=c(1,1,1,1,0,0,0,0),
    my.censoring=c(0,0,0,0,1,1,1,1),
    c1=c(0,1,1,0,1,1,0,0)
  )
  
  #error if any NA
  my.data$c1[1] <- NA
  expect_error(.fitPHmodel(formula(~c1),event=FALSE,my.data,time="my.time",has.event="my.event",has.censoring="my.censoring"))
  
  my.data$c1[1] <- 0
  #return NULL if cannot fit event
  expect_true(is.null(.fitPHmodel(formula(~c1),event=FALSE,my.data[1:4,],time="my.time",has.event="my.event",has.censoring="my.censoring")))
  #or censoring
  expect_true(is.null(.fitPHmodel(formula(~c1),event=TRUE,my.data[5:8,],time="my.time",has.event="my.event",has.censoring="my.censoring")))
  #but ok in reverse case
  expect_false(is.null(.fitPHmodel(formula(~c1),event=TRUE,my.data[1:4,],time="my.time",has.event="my.event",has.censoring="my.censoring")))
  
  samePH <- function(x,y){
    expect_equal(x$coefficients,y$coefficients)
    expect_equal(x$var,y$var)
    expect_equal(x$residuals,y$residuals)
  }
  
  #check model fit works
  my.mod <- .fitPHmodel(formula(~c1),event=TRUE,my.data,time="my.time",has.event="my.event",has.censoring="my.censoring")
  surv.mod <- coxph(Surv(my.time,my.event)~c1,data=my.data)
  samePH(surv.mod,my.mod)

  #and that it works for censored model
  my.mod <- .fitPHmodel(formula(~c1),event=FALSE,my.data,time="my.time",has.event="my.event",has.censoring="my.censoring",ties="breslow")
  surv.mod <- coxph(Surv(my.time,my.censoring==1)~c1,data=my.data,ties="breslow")
  samePH(surv.mod,my.mod)
  expect_equal("breslow",my.mod$method)
  
  #and censored model works with my.censoring including administrative censoring
  my.data$my.censoring[8] <- 2
  my.mod <- .fitPHmodel(formula(~c1),event=FALSE,my.data,time="my.time",has.event="my.event",has.censoring="my.censoring")
  surv.mod <- coxph(Surv(my.time,my.censoring==1)~c1,data=my.data)
  samePH(surv.mod,my.mod)
  
})

#note there will always be at least one row in the raw.scores data frame
test_that("normalize.scores",{
  
  my.df <- data.frame(Rs.f=c(1,2,3,4,5),
                      Rs.c=c(5,5,5,5,5))
  
  norm.df <- normalize.scores(my.df)
  expect_equal(rep(0,5),norm.df$Rs.c)
  expect_equal(((1:5)-2.5)/sqrt(5/3),norm.df$Rs.f)
  
  norm.df <- normalize.scores(my.df[1,])
  expect_equal(0,norm.df$Rs.f)
  expect_equal(0,norm.df$Rs.c)
})

test_that(".internalDistances",{
  
  my.df <- data.frame(Rs.f=c(1,2,3,4,5),
                      Rs.c=c(5,5,5,5,5))
  
  s <- 1/sqrt(5/3)
  
  expect_equal(rep(0,4),.internalDistances(my.df,w.censoring = 1))
  expect_equal((4:1)*s,.internalDistances(my.df,w.censoring = 0))
  expect_equal((4:1)*s*sqrt(0.5),.internalDistances(my.df,w.censoring = 0.5))
 
  my.df$Rs.c <- my.df$Rs.f
  expect_equal((4:1)*s,.internalDistances(my.df,w.censoring = 0.5))
  
  my.df$Rs.c <- 5:1
  expect_equal((4:1)*s,.internalDistances(my.df,w.censoring = 1))
  
})
