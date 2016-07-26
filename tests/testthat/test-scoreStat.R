context("ScoreStat")

mockImputed <- function(){
  data(ScoreInd)
  
  df <- ScoreInd
  
  df$impute.event <- df$event
  df$impute.time <- df$time
  
  df$my.arm <- df$arm
  df$arm <- NULL
  
  col.control <- list(has.event="event",
                      time="time",
                      Id="Id",
                      arm="my.arm",
                      DCO.time="DCO.time",
                      to.impute="to.impute")
  
  retVal <- list(data=df,
    col.control=col.control,
    default.formula=formula("~ my.arm"))
  
  class(retVal) <- "ScoreImputedData"
  retVal
}


test_that("creation_invalid",{
  
  scoreData <- mockImputed()
  expect_error(ImputeStat(scoreData$data,method="logrank"))
  expect_error(ImputeStat(scoreData,method="log"))
  expect_error(ImputeStat(scoreData,method="exponential"))
  expect_error(ImputeStat(scoreData,method="weibull"))
  
  #No formula apart from stratified for Wilcoxon/logrank
  expect_error(ImputeStat(scoreData,method="logrank",formula=~my.arm+Z1))
  expect_error(ImputeStat(scoreData,method="Wilcoxon",formula=~my.arm+Z1))
  expect_error(ImputeStat(scoreData,method="Wilcoxon",formula=~my.arm+strata(Z1)+Z2))
  expect_error(ImputeStat(scoreData,method="Wilcoxon",formula=~strata(Z1)+my.arm))
  
  #Invalid Cox formulae
  expect_error(ImputeStat(scoreData,method="Cox",formula=~Z1))
  expect_error(ImputeStat(scoreData,method="Cox",formula=~my.arm+Z1*my.arm))
  expect_error(ImputeStat(scoreData,method="Cox",formula=Surv(impute.time,impute.event)~my.arm))
})

test_that("logrank",{
  scoreData <- mockImputed()
  
  a <- ImputeStat(scoreData,method="logrank",formula=~my.arm)
  expect_equal("ScoreStat",class(a))
  expect_equal(c("model","method","estimate","var","statistic"),names(a))
  
  expect_equal("logrank (estimator for O-E)",a$method)
  
  my.mod <- survdiff(Surv(impute.time,impute.event)~my.arm,data=scoreData$data)
  
  #hack call to then show model is same
  my.mod$call <- a$model$call
  expect_equal(my.mod,a$model)
  
  #Z^2 is correct
  expect_equal(my.mod$chisq,a$statistic^2)
  
  expect_equal(a$estimate/sqrt(a$var),a$statistic )
  expect_equal(a$var,my.mod$var[2,2])
  expect_equal(a$estimate,my.mod$obs[2]-my.mod$exp[2])
  
  #default test is logrank and default formula is ~my.arm
  a1 <- ImputeStat(scoreData)
  expect_equal(a,a1)
  
})

test_that("asvector",{
  
  #mock object
  a <- list(model=NULL,
            estimate=6,
            var=9,
            statistic=2,
            test=NULL)
  
  class(a) <- "ScoreStat"
  ans <- as.vector(a)
  
  expects <- c(6,9,2)
  names(expects) <- c("estimate","var","statistic")
  
  expect_equal(expects,ans)
  
})

test_that("Wilcoxon",{
  scoreData <- mockImputed()
  
  #remove time and event cols as they should not be being used
  scoreData$data$time <- NULL
  scoreData$data$event <- NULL
  
  a <- ImputeStat(scoreData,method="Wilcoxon")
  expect_equal("ScoreStat",class(a))
  expect_equal(c("model","method","estimate","var","statistic"),names(a))
  
  expect_equal("Wilcoxon (estimator for O-E)",a$method)
  
  my.mod <- survdiff(Surv(impute.time,impute.event)~my.arm,data=scoreData$data,rho = 1)
  
  #hack call to then show model is same
  my.mod$call <- a$model$call
  expect_equal(my.mod,a$model)
  
  #Z^2 is correct
  expect_equal(my.mod$chisq,a$statistic^2)
  
  expect_equal(a$estimate/sqrt(a$var),a$statistic )
  expect_equal(a$var,my.mod$var[2,2])
  expect_equal(a$estimate,my.mod$obs[2]-my.mod$exp[2])
})

samePH <- function(x,y){
  expect_equal(x$coefficients,y$coefficients)
  expect_equal(x$var,y$var)
  expect_equal(x$residuals,y$residuals)
}

test_that("Cox_defaultformula",{
  scoreData <- mockImputed()
  a <- ImputeStat(scoreData,method="Cox",ties="breslow")
  expect_equal("ScoreStat",class(a))
  
  expect_equal("Cox",a$method)
  
  my.mod <- coxph(Surv(impute.time,impute.event)~my.arm,data=scoreData$data,ties="breslow")
  
  expect_equal("breslow",a$model$method)
  
  samePH(my.mod,a$model)
  expect_equal(a$statistic,a$estimate/sqrt(a$var))
  expect_equal(a$estimate,my.mod$coefficients["my.arm1"])
  expect_true(abs(a$var-my.mod$var[1,1])<1e-12)
})


test_that("Cox_usersuppliedformula",{
  scoreData <- mockImputed()
  a <- ImputeStat(scoreData,method="Cox",formula=~my.arm+Z1+Z4)
  expect_equal("ScoreStat",class(a))
  
  expect_equal("Cox",a$method)
  
  my.mod <- coxph(Surv(impute.time,impute.event)~my.arm+Z1+Z4,data=scoreData$data)
  expect_equal("efron",a$model$method)
  
  samePH(my.mod,a$model)
  expect_equal(a$statistic,a$estimate/sqrt(a$var))
  expect_equal(a$estimate,my.mod$coefficients["my.arm1"])
  expect_true(abs(a$var-my.mod$var[1,1])<1e-12)
})


test_that("Logrank_Wilcoxon_usersupplied_formula",{
  scoreData <- mockImputed()
  a <- ImputeStat(scoreData,method="logrank",formula=~my.arm+strata(Z1))
  expect_equal(a$statistic,sum(a$model$obs[2,]-a$model$exp[2,])/sqrt(a$model$var[2,2]))
  
  scoreData <- mockImputed()
  a <- ImputeStat(scoreData,method="logrank",formula=~my.arm+strata(Z1,Z3))
  a2 <- ImputeStat(scoreData,method="logrank",formula=~my.arm+strata(Z3) + strata(Z1))
  expect_equal(as.vector(a),as.vector(a2))
  
})
