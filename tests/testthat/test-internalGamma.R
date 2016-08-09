context("internalGamma")


test_that(".CoxMImpute",{
  #Here we check we get the same result as the example code given by Jackson

  Jackson.snippet <- function(mytime,lp,hazard_boot,U){
    hazard_boot$time=hazard_boot$time-mytime
    haz_to_sub=subset(hazard_boot, time<=0)$hazard
    hazard_boot=subset(hazard_boot, time>0)
    if(length(haz_to_sub)>0) hazard_boot$hazard=hazard_boot$hazard-max(haz_to_sub)
    hazard_boot=rbind(c(0,0), hazard_boot, c(Inf, Inf))
    hazard_boot=data.frame(hazard=hazard_boot[,1], time=hazard_boot[,2])
    quant=-log(U)*exp(-lp)
    answer=min(subset(hazard_boot, hazard>=quant)$time)
    mytime+answer
  }

  check.same <- function(my.time,beta,basehaz,seed){
    set.seed(seed)
    Jans <- Jackson.snippet(my.time,beta,basehaz,runif(1))

    set.seed(seed)
    myans <- .CoxMImpute(my.time,beta,basehaz)

    expect_equal(Jans,myans)
  }


   basehaz <- data.frame(hazard=c(0.1,0.5,1.0,1.5,2.5,4),
                         time=c(1,2,3,4,5,6))

   check.same(my.time=3,beta=0.4,basehaz,seed=10)
   check.same(my.time=8,beta=0.8,basehaz,seed=23)
   check.same(my.time=5,beta=-0.8,basehaz,seed=123)
   check.same(my.time=5.3,beta=1,basehaz,seed=1323)
   check.same(my.time=0.5,beta=1,basehaz,seed=1323)


})

test_that(".SetupStrataBaseHaz.no.strata",{

  data <- data.frame(id=1:5,
                     time=1:5,
                     event=c(0,1,1,0,1),
                     x=c(1,0,1,1,0),
                     y=c(1,1,1,0,0))

  formula <- formula(Surv(time,event)~x)

  basehaz <- data.frame(hazard=c(1,2.2,4.5,8),
                        time=c(1:4))

  stratas <- character(0)
  ans <- .SetupStrataBaseHaz(data,formula,stratas,basehaz)

  expect_equal(rep("all",5),ans$index)
  expect_equal("all",names(ans$base.hazard))

})


test_that(".SetupStrataBaseHaz.strata.error",{

  data <- data.frame(id=1:5,
                     time=1:5,
                     event=c(0,1,1,0,1),
                     x=c(1,0,1,1,0),
                     y=c(1,1,1,0,0))

  formula <- formula(Surv(time,event)~strata(x))

  basehaz <- data.frame(hazard=c(1,2.2,4.5,2,4,7,11),
                        time=c(1:3,1:4),
                        strata=c(rep("x=1",3),rep("x=2",4)))

  stratas <- "strata(x)"

  expect_error(.SetupStrataBaseHaz(data,formula,stratas,basehaz))

})


test_that(".SetupStrataBaseHaz.one.strata",{

  data <- data.frame(id=1:5,
                     time=1:5,
                     event=c(0,1,1,0,1),
                     x=c(1,0,1,1,0),
                     y=c(1,1,1,0,0))

  formula <- formula(Surv(time,event)~strata(x))
  basehaz <- data.frame(hazard=c(1,2.2,4.5,2,4,7,11),
                        time=c(1:3,1:4),
                        strata=c(rep("x=0",3),rep("x=1",4)))

  stratas <- "strata(x)"

  ans <- .SetupStrataBaseHaz(data,formula,stratas,basehaz)
  expect_equal(c("x=1","x=0","x=1","x=1","x=0"),ans$index)
  expect_equal(c("x=0","x=1"),names(ans$base.hazard))

  expect_equal(ans$base.hazard$`x=0`,basehaz[basehaz$strata=="x=0",])


  #Now redo with x as factor
  data$x <- factor(data$x)
  basehaz$strata <- c(rep(0,3),rep(1,4))
  ans <- .SetupStrataBaseHaz(data,formula,stratas,basehaz)
  expect_equal(c("1","0","1","1","0"),ans$index)
  expect_equal(c("0","1"),names(ans$base.hazard))

  expect_equal(ans$base.hazard$`1`,basehaz[basehaz$strata=="1",])

  #Now check strata(x*y) works
  data$x <- c(1,0,1,1,0)
  formula <- formula(Surv(time,event)~strata(x*y))
  stratas <- "strata(x * y)"
  basehaz$strata <- c(rep("x * y=0",3),rep("x * y=1",4))
  ans <- .SetupStrataBaseHaz(data,formula,stratas,basehaz)
  expect_equal(c("x * y=1","x * y=0","x * y=1","x * y=0","x * y=0"),ans$index)
  expect_equal(c("x * y=0","x * y=1"),names(ans$base.hazard))
})

test_that(".SetupStrataBaseHaz.two.strata",{
  data <- data.frame(id=1:5,
                     time=1:5,
                     event=c(0,1,1,0,1),
                     x=c(1,0,1,1,0),
                     y=c(1,1,1,0,0))

  data$x <- factor(data$x)
  formula <- formula(Surv(time,event)~strata(x)+strata(y))
  stratas <- c("strata(x)","strata(y)")

  basehaz <- data.frame(hazard=c(1,3,2,1),
                        time=1:4,
                        strata=c("0, y=0","0, y=1","1, y=0","1, y=1"))

  ans <- .SetupStrataBaseHaz(data,formula,stratas,basehaz)
  expect_equal(c("1, y=1","0, y=1","1, y=1","1, y=0","0, y=0"),ans$index)
  expect_equal(c("0, y=0","0, y=1","1, y=0","1, y=1"),names(ans$base.hazard))
  expect_equal(ans$base.hazard$`0, y=1`,basehaz[basehaz$strata=="0, y=1",])
})

test_that("untangle.specials.behaves.as.expected",{

  data <- data.frame(id=1:9,
                     time=1:9,
                     event=c(0,1,1,0,1,1,0,1,1),
                     x=c(1,0,1,1,0,1,1,0,1),
                     y=c(1,1,1,0,0,1,1,0,1))

  formula <- formula(Surv(time,event)~x)
  m <- coxph(formula,data)
  expect_equal(character(0),untangle.specials(m$terms,"strata")$var)


  formula <- formula(Surv(time,event)~x+strata(y))
  expect_warning(m <- coxph(formula,data)) #not interested in lack of convergence
  expect_equal("strata(y)",untangle.specials(m$terms,"strata")$var)

  formula <- formula(Surv(time,event)~x*strata(y))
  expect_warning(m <- coxph(formula,data)) #not interested in lack of convergence
  expect_equal("strata(y)",untangle.specials(m$terms,"strata")$var)

  formula <- formula(Surv(time,event)~strata(x)*strata(y))
  suppressWarnings(m <- coxph(formula,data)) #not interested in lack of convergence
  expect_equal(c("strata(x)","strata(y)"),untangle.specials(m$terms,"strata")$var)

  formula <- formula(Surv(time,event)~strata(x*y))
  m <- coxph(formula,data)
  expect_equal("strata(x * y)",untangle.specials(m$terms,"strata")$var)

  formula <- formula(Surv(time,event)~strata(x, y))
  m <- coxph(formula,data)
  expect_equal("strata(x, y)",untangle.specials(m$terms,"strata")$var)

  formula <- formula(Surv(time,event)~strata(y) + strata( x))
  m <- coxph(formula,data)
  expect_equal(c("strata(y)","strata(x)"),untangle.specials(m$terms,"strata")$var)

})
