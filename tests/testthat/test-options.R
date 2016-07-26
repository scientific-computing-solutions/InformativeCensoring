context("options")

test_that(".internal.is.finite.number",{
  expect_true(.internal.is.finite.number(10))
  expect_true(.internal.is.finite.number(0))
  expect_true(.internal.is.finite.number(-8))
  expect_true(.internal.is.finite.number(4.6764))
  expect_false(.internal.is.finite.number("w"))
  expect_false(.internal.is.finite.number(c(4,5,6)))
  expect_false(.internal.is.finite.number(TRUE))
  expect_false(.internal.is.finite.number(Inf))
  expect_false(.internal.is.finite.number(NA))
  expect_false(.internal.is.finite.number(as.numeric(NA)))
})

test_that("NN.options_invalid",{
  expect_error(NN.options(NN=0))
  expect_error(NN.options(NN=-1))
  expect_error(NN.options(NN=3.5))
  expect_error(NN.options(NN=c(5,6)))
  expect_error(NN.options(NN="fr"))
  expect_error(NN.options(NN=5,w.censoring = -1.2))
  expect_error(NN.options(NN=5,w.censoring = 1.4))
  expect_error(NN.options(NN=5,w.censoring = c(0.5,0.3)))
  expect_error(NN.options(NN=5,w.censoring = FALSE))
  expect_error(NN.options(NN=5,w.censoring = 0.2,min.subjects=c(2,10)))
  expect_error(NN.options(NN=5,w.censoring = 0.2,min.subjects="hello"))
  expect_error(NN.options(NN=5,w.censoring = 0.2,min.subjects=0))
  expect_error(NN.options(NN=5,w.censoring = 0.2,min.subjects=4.5))
})

test_that("NN.options",{
  #default
  n <- NN.options()
  expect_equal(5,n$NN)
  expect_equal(0.2,n$w.censoring)
  expect_equal(20,n$min.subjects)
  expect_equal(c("NN","w.censoring","min.subjects"),names(n))
  
  #non default 
  n <- NN.options(NN=6,w.censoring = 0.5,min.subjects=15)
  expect_equal(6,n$NN)
  expect_equal(0.5,n$w.censoring)
  expect_equal(15,n$min.subjects)
})

test_that("dummy.col.headings",{
  expect_that(.dummy.col.headings(),not(throws_error()))
})


test_that("col.headings",{
  #note most tests require a data set with col.headings list
  #and so are tested elsewhere
  expect_error(col.headings(arm="a",has.event="",time="t",Id="I",DCO.time="D",to.impute ="ti"))
  expect_error(col.headings(arm="a",has.event="h",time="t",Id=c("I","J"),DCO.time="D",to.impute ="ti"))
  col.control <- col.headings(arm="a",has.event="h",time="t",Id="I",DCO.time="D",to.impute ="ti")
  expect_equal(c("arm","has.event","time","Id","DCO.time","to.impute","censor.type"),names(col.control))
  expect_equal("a",col.control$arm)
  expect_equal("h",col.control$has.event)
  expect_equal("t",col.control$time)
  expect_equal("I",col.control$Id)
  expect_equal("D",col.control$DCO.time)
  expect_equal("ti",col.control$to.impute)
  
  #No arm
  expect_error(col.headings(has.event="event",
                            time="time", Id="Id",
                            DCO.time="DCO.time",
                            to.impute="to.impute"))
  
})


