context("gammaImputeData")

#see validate_gamma_arguments in test-validation for
#input parameter validation


test_that("GammaImputeSet_object",{
  
  set.seed(25)
  load("gamma_test.rda")
  
  imputed.data.sets <- gammaImpute(formula=Surv(Yi,delta)~Z+strata(W1),
                                   data = gamma.dataset,
                                   m=2, DCO.time="DCO.time",
                                   bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                                   gamma.factor=1,gamma="gamma")
  
  expect_equal("GammaImputedSet",class(imputed.data.sets))
  expect_equal(2,imputed.data.sets$m)
  
  ans <- imputed.data.sets$data
  ans$internal_gamma_val <- NULL
  
  expect_equal(gamma.dataset,ans)
  expect_equal("matrix",class(imputed.data.sets$impute.time))
  expect_equal("matrix",class(imputed.data.sets$impute.event))
  
  expect_equal(2,ncol(imputed.data.sets$impute.time))
  expect_equal(nrow(gamma.dataset),nrow(imputed.data.sets$impute.event))
  
  expect_equal(1,imputed.data.sets$gamma.factor)
  expect_equal(formula(~Z+strata(W1)),imputed.data.sets$default.formula)
  
})

test_that("gamma.factor",{
  load("gamma_test.rda")
  
  set.seed(410)
  gamma.dataset$gamma <- runif(nrow(gamma.dataset))
  
  set.seed(10)
  imputed.data.sets <- gammaImpute(formula=Surv(Yi,delta)~Z+strata(W1),
                                   data = gamma.dataset,
                                   m=2,gamma="gamma",
                                   bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                                   gamma.factor=-1.765,DCO.time="DCO.time")
  
  expect_equal(imputed.data.sets$data$internal_gamma_val,-1.765*gamma.dataset$gamma)
  
  gamma.dataset$gamma <- gamma.dataset$gamma*-1.765
  set.seed(10)
  imputed.data.sets2 <-gammaImpute(formula=Surv(Yi,delta)~Z+strata(W1),
                                   data = gamma.dataset,
                                   m=2,gamma="gamma",
                                   bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                                   gamma.factor=1,DCO.time="DCO.time")
  
  
  expect_equal(imputed.data.sets$impute.time,imputed.data.sets2$impute.time)
  expect_equal(imputed.data.sets$impute.event,imputed.data.sets2$impute.event)
})

test_that("inf.gamma.factor",{
  load("gamma_test.rda")
  
  set.seed(10)
  inf.gamma <- gammaImpute(formula=Surv(Yi,delta)~Z+strata(W1),
                           data = gamma.dataset,
                           m=2,gamma="gamma",DCO.time="DCO.time",
                           bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                           gamma.factor=Inf)
  set.seed(10)
  neg.inf.gamma <- gammaImpute(formula=Surv(Yi,delta)~Z+strata(W1),
                               data = gamma.dataset,
                               m=2,gamma="gamma",DCO.time="DCO.time",
                               bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                               gamma.factor=-Inf)
  
  inf.gamma.df <- ExtractSingle(inf.gamma,1)
  neg.inf.gamma.df <- ExtractSingle(neg.inf.gamma,2)
  
  index <- which(inf.gamma.df$data$to.impute==1)
  
  #+Inf gamma implies have event at censoring
  expect_equal(inf.gamma.df$data$impute.time[index],inf.gamma.df$data$Yi[index])
  
  #-Inf gamma implies have no event
  expect_equal(neg.inf.gamma.df$data$impute.time[index],neg.inf.gamma.df$data$DCO.time[index])
})

test_that("impute.times.non_stochastic_tests",{
  load("gamma_test.rda")
  
  df <- gamma.dataset
  df$gamma[1:49] <- as.numeric(NA)
  df$gamma[50] <- NA
  df$gamma[50:100] <- 1
  
  #check data set is suitable for tests below
  expect_false(all(is.na(df$gamma[df$delta])))
  expect_true(any(is.na(df$gamma[df$delta])))
  expect_false(all(is.na(df$gamma[!df$delta])))
  expect_true(any(is.na(df$gamma[!df$delta])))
  
  df$DCO.time <- df$Yi + runif(nrow(df),1,5) 
  
  ans <- ExtractSingle(gammaImpute(formula=Surv(Yi,delta)~Z+strata(W1),df,
              m=2,gamma="gamma",DCO.time="DCO.time",
              bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
              gamma.factor=1),index=1)
  
  expect_equal("GammaImputedData",class(ans))
  expect_equal(formula("~Z+strata(W1)"),ans$default.formula)
  
  #gamma = NA implies no imputation
  index <- which(is.na(df$gamma))
  expect_equal(ans$data$Yi[index],ans$data$impute.time[index])
  expect_equal(as.numeric(ans$data$delta[index]),ans$data$impute.event[index])
  
  expect_true(all(ans$data$DCO.time>=ans$data$impute.time))
  expect_true(all(ans$data$Yi<=ans$data$impute.time))
  
  #All DCO.time are censored
  expect_true(all(!ans$data$impute.event[ans$data$DCO.time==ans$data$impute.time]))
  
  #gamma=NA, with already having event remains unchanged
  index <- which(ans$data$delta)
  expect_true(all(as.logical(ans$data$impute.event[index])))
  expect_true(all(ans$data$Yi[index]==ans$data$impute.time[index]))
})

test_that("internal_gamma_val_same_with_colname_or_vect",{
  load("gamma_test.rda")
  
  set.seed(2323)
  gamma.dataset$G <- runif(n=nrow(gamma.dataset))
  
  set.seed(2323)
  use.gamma.col <- gammaImpute(formula=Surv(Yi,delta)~Z+strata(W1),
                               data = gamma.dataset,
                               m=2,DCO.time="DCO.time",gamma="G",
                               bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1))
  
  gamma.vec <- gamma.dataset$G
  
  set.seed(2323)
  use.gamma.vec <- gammaImpute(formula=Surv(Yi,delta)~Z+strata(W1),
                               data = gamma.dataset,
                               m=2,DCO.time="DCO.time",gamma=gamma.vec,
                               bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1))
  
  expect_equal(use.gamma.col,use.gamma.vec)
  expect_equal(use.gamma.col$data$internal_gamma_val,gamma.vec)
  
})


test_that("gamma_and_gamma_factor",{
  load("gamma_test.rda")
  
  gamma.dataset$G <- gamma.dataset$gamma
  
  set.seed(10)
  use.gamma.factor <- gammaImpute(formula=Surv(Yi,delta)~Z+strata(W1),
                           data = gamma.dataset,
                           m=2,DCO.time="DCO.time",
                           bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                           gamma.factor=5)

  set.seed(10)
  gamma.dataset$G <- 5
  use.gamma <- gammaImpute(formula=Surv(Yi,delta)~Z+strata(W1),
                           data = gamma.dataset,
                           m=2,DCO.time="DCO.time",gamma="G",
                           bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1))
  
  gvec <- gamma.dataset$G
  set.seed(10)
  use.gamma.vector <-  gammaImpute(formula=Surv(Yi,delta)~Z+strata(W1),
                                   data = gamma.dataset,
                                   m=2,DCO.time="DCO.time",gamma=gvec,
                                   bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1))
  
  expect_equal(use.gamma.vector,use.gamma)
  
  #everything except gamma.factor and gamma column are the same
  use.gamma$gamma.factor <- NULL
  use.gamma.factor$gamma.factor <- NULL
  use.gamma$data$G <- use.gamma.factor$data$G
  expect_equal(use.gamma,use.gamma.factor)
  
})


test_that("DCO.time_vector",{
  load("gamma_test.rda")
  
  set.seed(10)
  
  dco.time <- 3+runif(n = nrow(gamma.dataset))
  
  DCO.ans <- gammaImpute(formula=Surv(Yi,delta)~Z+strata(W1),
                         data = gamma.dataset,
                         m=2,DCO.time=dco.time,
                         bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                         gamma.factor=-Inf)
  
  expect_equal(DCO.ans$data$internalDCO.time,dco.time)
  dco.one <- ExtractSingle(DCO.ans,1)
  expect_equal(dco.time[gamma.dataset$to.impute],
               dco.one$data$impute.time[gamma.dataset$to.impute])
})



test_that("DCO.time_different_inputs_agree",{
  load("gamma_test.rda")
  
  set.seed(10)
  DCO.col <- gammaImpute(formula=Surv(Yi,delta)~Z+strata(W1),
                                  data = gamma.dataset,
                                  m=2,DCO.time="DCO.time",
                                  bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                                  gamma.factor=1)
  set.seed(10)
  DCO.single <- gammaImpute(formula=Surv(Yi,delta)~Z+strata(W1),
                            data = gamma.dataset,
                            m=2,DCO.time=3,
                            bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                            gamma.factor=1)
  
  set.seed(10)
  DCO.all <- gammaImpute(formula=Surv(Yi,delta)~Z+strata(W1),
                         data = gamma.dataset,
                         m=2,DCO.time=rep(3,nrow(gamma.dataset)),
                         bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                         gamma.factor=1)
  
  expect_equal(DCO.single,DCO.all)
  #remove internal_DCO from DCO.single should equal DCO.col
  DCO.single$data$internalDCO.time <- NULL
  expect_equal(DCO.col,DCO.single)
})

