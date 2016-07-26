context("gammaStat")

test_that("GammaStat",{
  load("gamma_test.rda")
  
  set.seed(10)
  imputed <- gammaImpute(formula=Surv(Yi,delta)~Z+W1,
                         gamma.dataset,m=2,DCO.time="DCO.time",
                         bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                         gamma.factor=1,gamma="gamma")
  
  df <- ExtractSingle(imputed,1)
  
  expect_error(ImputeStat(df,method="Wilcoxon"))
  expect_error(ImputeStat(df,method="logrank"))
  expect_that(ImputeStat(df,method="Cox"),not(throws_error())) #no error uses default formula
  
  expect_error(ImputeStat(df,method="Cox",formula=~cluster(W1)))
  expect_error(ImputeStat(df,method="Cox",formula=~Z+tt(W1)))
  expect_error(ImputeStat(df,method="Cox",formula=Surv(impute.time,impute.event)~Z))
  
  
  ans <- ImputeStat(df,method="Cox",formula=~Z+W1)
  expect_equal("GammaStat",class(ans))
  expect_equal("Cox",ans$method)
  
  mod <- coxph(Surv(impute.time,impute.event)~Z+W1,data=df$data,model=TRUE)
  
  #change Call so that mod should equal ans$model
  mod$call <- ans$model$call
  expect_equal(mod,ans$model)
  
  vars <- ans$var
  expect_equal(c("Z1","Z2","W1"),names(vars))
  names(vars) <- NULL
  expect_true(abs(vcov(mod)[1,1]-vars[1])<1e-8)
  expect_true(abs(vcov(mod)[3,3]-vars[3])<1e-8)
  
  expect_equal(coefficients(mod),ans$estimate)

})


test_that("weibull_and_exponential",{
  load("gamma_test.rda")

  set.seed(10)
  imputed <- gammaImpute(formula=Surv(Yi,delta)~Z+W1,
                         gamma.dataset,m=2,gamma="gamma",
                         bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                         gamma.factor=1,DCO.time="DCO.time")
                         
  df <- ExtractSingle(imputed,1)
  ans <- ImputeStat(df,method="weibull")

  expect_equal("weibull",ans$method)
  mod <- survreg(Surv(impute.time,impute.event)~Z+W1,data=df$data,dist="weibull")
  expect_equal(mod$coefficients,ans$estimate)
  
  #Do not pull out scales at least for now
  expect_equal(ans$var,diag(vcov(mod))[1:4])
  
  #Test strata
  ans <- ImputeStat(df,method="weibull",formula=~Z+strata(W1))
  mod <- survreg(Surv(impute.time,impute.event)~Z+strata(W1),data=df$data,dist="weibull")
  expect_equal(mod$coefficients,ans$estimate)
  expect_equal(ans$var,diag(vcov(mod))[1:3])
  
  #Test exponential
  ans <- ImputeStat(df,method="exponential")
  expect_equal("exponential",ans$method)
  mod <- survreg(Surv(impute.time,impute.event)~Z+W1,data=df$data,dist="exponential")
  expect_equal(mod$coefficients,ans$estimate)
  expect_equal(ans$var,diag(vcov(mod)))
})


test_that("default_formula_behaves_as_expected",{
  load("gamma_test.rda")
 
  set.seed(10)
  imputed <- gammaImpute(formula=Surv(Yi,delta)~Z,
                         gamma.dataset,m=2,gamma="gamma",
                         bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                         gamma.factor=1,DCO.time="DCO.time")
  
  df <- ExtractSingle(imputed,1)
  ans <- ImputeStat(df,method="Cox",formula=~Z)
 
  expect_equal(formula("Surv(impute.time,impute.event)~Z"),ans$model$formula)
 
  ans2 <- ImputeStat(df,method="Cox",formula=~Z+W1)
  expect_equal(formula("Surv(impute.time,impute.event)~Z+W1"),ans2$model$formula)
})


test_that("GammaStatSet",{
  load("gamma_test.rda")
  set.seed(10)
  imputed <- gammaImpute(formula=Surv(Yi,delta)~Z+W1,
                         gamma.dataset,m=2,gamma="gamma",
                         bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                         gamma.factor=1,DCO.time="DCO.time")
  
  ans <- ImputeStat(imputed,method="Cox",formula=~W1+strata(Z))

  expect_equal("GammaStatList",class(ans))
  expect_equal(2,ans$m)
  
  df <- ExtractSingle(imputed,1)
  expect_equal(ExtractSingle(ans,1),ImputeStat(df,method="Cox",formula=~W1+strata(Z)))
  
  expect_equal(2,length(ans$fits))
  expect_equal(c("estimates","vars"),names(ans$statistics))
  expect_equal("matrix",class(ans$statistics$estimates))
  expect_equal("matrix",class(ans$statistics$vars))
  expect_equal(2,nrow(ans$statistics$estimates))
  expect_equal(1,ncol(ans$statistics$vars))
  expect_equal("W1",colnames(ans$statistics$vars))
  
  expect_equal(ans$fits[[1]]$estimate,ans$statistics$estimates[1,1])
  expect_equal(ans$fits[[2]]$estimate,ans$statistics$estimates[2,1])
  expect_equal(ans$fits[[1]]$var,ans$statistics$vars[1,1])
  expect_equal(ans$fits[[2]]$var,ans$statistics$vars[2,1])
  
})

test_that("GammaStatSet_multiplecovar",{
  load("gamma_test.rda")
  set.seed(10)
  imputed <- gammaImpute(formula=Surv(Yi,delta)~Z+W1,
                         gamma.dataset,m=3,gamma="gamma",
                         bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                         gamma.factor=1,DCO.time="DCO.time")
  
  
  ans <- ImputeStat(imputed,method="Cox",formula=~Z+strata(W1))
  
  expect_equal(3,length(ans$fits))
  expect_equal("matrix",class(ans$statistics$estimates))
  expect_equal("matrix",class(ans$statistics$vars))
  expect_equal(3,nrow(ans$statistics$estimates))
  expect_equal(2,ncol(ans$statistics$vars))
  expect_equal(c("Z1","Z2"),colnames(ans$statistics$vars))
  
  expect_equal(ans$fits[[1]]$estimate[1],ans$statistics$estimates[1,1])
  expect_equal(ans$fits[[2]]$estimate[1],ans$statistics$estimates[2,1])
  expect_equal(ans$fits[[3]]$estimate[2],ans$statistics$estimates[3,2])
  expect_equal(ans$fits[[1]]$var[2],ans$statistics$vars[1,2])
  expect_equal(ans$fits[[2]]$var[1],ans$statistics$vars[2,1])
  expect_equal(ans$fits[[3]]$var[2],ans$statistics$vars[3,2])
  
})

test_that("summary_1_covar",{
  #A mock GammaStatList object
  
  estimates <- matrix(c(2,4),ncol=1)
  vars <- matrix(c(0.5,0.8),ncol=1)
  
  colnames(estimates) <- "W1"
  colnames(vars) <- "W1"
  
  
  fits <- list(m=2,
               statistics=list(estimates=estimates,
                               vars=vars))
  class(fits) <- "GammaStatList"
  
  ans <- summary(fits)
  expect_equal("matrix",class(ans))

  expect_equal("W1",rownames(ans))
  expect_equal(c("est","se","t","df","Pr(>|t|)","lo 95","hi 95"),colnames(ans))
  
  expect_equal(3,ans[1,1])
  expect_equal(sqrt(3.65),ans[1,2])
  expect_equal(3/sqrt(3.65),ans[1,3])
  expect_equal((1+0.65/3)^2,ans[1,4])
  expect_equal(2*(1-pt(3/sqrt(3.65),df=(1+0.65/3)^2)),ans[1,5]) 
  tval <- qt(0.975,df=(1+0.65/3)^2)
  expect_equal(3-tval*sqrt(3.65),ans[1,6])
  expect_equal(3+tval*sqrt(3.65),ans[1,7])

})

test_that("summary_3_covar",{
  
  #A mock GammaStatList object
  
  estimates <- matrix(1:9,ncol=3)
  vars <- matrix(seq(0.1,0.9,0.1),ncol=3)
  
  colnames(estimates) <- c("W1","R4","Q")
  colnames(vars) <- c("W1","R4","Q")
  fits <- list(m=3,
               statistics=list(estimates=estimates,
                               vars=vars))
  class(fits) <- "GammaStatList"
  
  ans <- summary(fits)
  expect_equal(c("W1","R4","Q"),rownames(ans))
  expect_equal(5,ans[2,1])
  expect_equal(sqrt(0.8+4/3),ans[3,2])
  expect_equal(2/sqrt(0.2+4/3),ans[1,3])
  expect_equal(2*(1+0.5*3/4)^2,ans[2,4])
})