#These tests don't work on R cmd check but do in devtools::test so they are
#not run gnerally (they are testing the parallelism)
context("parallel")

test_that("valid_parallel_args",{
  
  expect_error(validate.parallel.arguments(parallel="no",ncpus=-5,cl=NULL))
  expect_error(validate.parallel.arguments(parallel="no",ncpus=34.5,cl=NULL))
  expect_error(validate.parallel.arguments(parallel="no",ncpus=c(2,6,10),cl=NULL))
  expect_error(validate.parallel.arguments(parallel="no",ncpus=5,cl="hello"))
  
  expect_error(validate.parallel.arguments(parallel="no",ncpus=5,cl=NULL))
  expect_error(validate.parallel.arguments(parallel="mue",ncpus=5,cl=NULL))
  
})

test_that("parallel",{
  
  f <- function(x,y){x*y}
  
  RNGkind("L'Ecuyer-CMRG")
  
  ans.1 <- parallelRun(parallel="multicore",ncpus=2,cl=NULL,lapply.list=1:10,FUN=f,y=5)
  expect_equal(unlist(ans.1),seq(5,50,5))
  
  ans.2 <- parallelRun(parallel="snow",ncpus=2,cl=NULL,lapply.list=1:10,FUN=f,y=5)
  expect_equal(unlist(ans.2),seq(5,50,5))
  
  RNGkind("Mersenne-Twister")
  expect_warning(parallelRun(parallel="multicore",ncpus=2,cl=NULL,lapply.list=1:10,FUN=f,y=5))
  
})



test_that("parallel_score",{
  RNGkind("L'Ecuyer-CMRG")
  set.seed(104)
  data(ScoreInd)
  
  col.control <- col.headings(has.event="event", time="time",Id="Id",arm="arm",
                              DCO.time="DCO.time", to.impute="to.impute")
  
  d1 <-  ScoreImpute(data=ScoreInd,event.model=~Z1+Z2+Z3+Z4+Z5,
                     col.control=col.control, m=5,
                     bootstrap.strata=ScoreInd$arm,
                     NN.control=NN.options(NN=5,w.censoring = 0.2),
                     parallel="snow",ncpus=2)
  
  expect_that(ImputeStat(d1,parallel="multicore",ncpus=2),not(throws_error()))
})


test_that("parallel_gamma",{
  RNGkind("L'Ecuyer-CMRG")
  set.seed(10)
  load("gamma_test.rda")
  d1 <- gammaImpute(formula=Surv(Yi,delta)~Z+strata(W1),
                    data = gamma.dataset,
                    m=5, DCO.time="DCO.time",
                    bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                    gamma.factor=1,gamma="gamma",
                    parallel="multicore",ncpus=2)
  
  expect_that(ImputeStat(d1,method="weibull",parallel="snow",ncpus=2),not(throws_error()))
  
  
})

test_that("parallel_seed",{
  
  RNGkind("L'Ecuyer-CMRG")
  
  load("gamma_test.rda")
  #snow is reproducible
  set.seed(10)
  d1 <- gammaImpute(formula=Surv(Yi,delta)~Z+strata(W1),
                    data = gamma.dataset,
                    m=9, DCO.time="DCO.time",
                    bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                    gamma.factor=1,gamma="gamma",
                    parallel="snow",ncpus=2)
  
  set.seed(10)
  d2 <- gammaImpute(formula=Surv(Yi,delta)~Z+strata(W1),
                    data = gamma.dataset,
                    m=9, DCO.time="DCO.time",
                    bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                    gamma.factor=1,gamma="gamma",
                    parallel="snow",ncpus=2)
  
  expect_equal(d1,d2)
  
  #multicore is reproducible
  set.seed(10)
  d1 <- gammaImpute(formula=Surv(Yi,delta)~Z+strata(W1),
                    data = gamma.dataset,
                    m=5, DCO.time="DCO.time",
                    bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                    gamma.factor=1,gamma="gamma",
                    parallel="multicore",ncpus=2)
  
  set.seed(10)
  d2 <- gammaImpute(formula=Surv(Yi,delta)~Z+strata(W1),
                    data = gamma.dataset,
                    m=5, DCO.time="DCO.time",
                    bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                    gamma.factor=1,gamma="gamma",
                    parallel="multicore",ncpus=2)
  
  expect_equal(d1,d2)
  
  
})