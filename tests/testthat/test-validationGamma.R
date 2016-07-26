context("validationGamma")

test_that("invalid_rhs_formula",{
  load("gamma_test.rda")
  
  
  expect_error(gammaImpute(formula=Surv(Yi,delta)~Z+strata(W1)+strata(W2),
                           gamma.dataset,m=2,gamma="gamma",DCO.time="DCO.time",
                           bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                           gamma.factor=1))
  
  expect_error(gammaImpute(formula=Surv(Yi,delta)~Z+cluster(W1)+strata(W2),
                           gamma.dataset,m=2,gamma="gamma",DCO.time="DCO.time",
                           bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                           gamma.factor=1))
  
  expect_error(gammaImpute(formula=Surv(Yi,delta)~~Z+tt(W1)+strata(W2),
                           gamma.dataset,m=2,gamma="gamma",DCO.time="DCO.time",
                           bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                           gamma.factor=1))
  
  
  expect_that(gammaImpute(formula=Surv(Yi,delta)~Z+strata(W1,W2),
                          gamma.dataset,DCO.time="DCO.time",m=2,gamma="gamma",
                          bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                          gamma.factor=1),not(throws_error()))
  
  
})


test_that("validate_gamma_arguments_lhs_formula",{
  load("gamma_test.rda")

  #ok to have expression in Surv
  expect_that(gammaImpute(formula=Surv(Yi,delta==W2)~Z+strata(W1),
                          gamma.dataset,DCO.time="DCO.time",m=2,gamma="gamma",
                          bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                          gamma.factor=1),not(throws_error()))
  
  #no lhs
  expect_error(gammaImpute(formula=~Z+strata(W1),
                          gamma.dataset,DCO.time="DCO.time",m=2,gamma="gamma",
                          bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                          gamma.factor=1))
  
  #lhs not right censored
  expect_error(gammaImpute(formula=Surv(Yi-5,Yi,delta)~Z+strata(W1),
                           gamma.dataset,DCO.time="DCO.time",m=2,gamma="gamma",
                           bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                           gamma.factor=1))
  #lhs not Surv
  expect_error(gammaImpute(formula=Yi~Z+strata(W1),
                           gamma.dataset,DCO.time="DCO.time",m=2,gamma="gamma",
                           bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                           gamma.factor=1))
})

test_that("validate_gamma_arguments_m_gammafactor_and_strata_call",{
  load("gamma_test.rda")
  
  surv.times <- as.matrix(model.frame(formula(Surv(Yi,delta)~1),data=gamma.dataset))
  
  Call <- call("mycall",event.model=~Z1+Z2)
  
  
  expect_that(validate.Gamma.arguments(data=gamma.dataset,surv.times=surv.times,
                                       m=2,gamma="gamma",
                                       strata=rep(1,nrow(gamma.dataset)),
                                       gamma.factor=1,DCO.time="DCO.time",Call=Call), not(throws_error()))
  
  
  #m
  expect_error(validate.Gamma.arguments(data=gamma.dataset,surv.times=surv.times,
                                        m=-4,gamma="gamma",
                                        strata=rep(1,nrow(gamma.dataset)),
                                        gamma.factor=1,DCO.time="DCO.time",Call=Call))
  
  expect_error(validate.Gamma.arguments(data=gamma.dataset,surv.times=surv.times,
                                        m=1,gamma="gamma",
                                        strata=rep(1,nrow(gamma.dataset)),
                                        gamma.factor=1,DCO.time="DCO.time",Call=Call))
  
  expect_error(validate.Gamma.arguments(data=gamma.dataset,surv.times=surv.times,
                                        m=1.6,gamma="gamma",
                                        strata=rep(1,nrow(gamma.dataset)),
                                        gamma.factor=1,DCO.time="DCO.time",Call=Call))
  
  #gamma.factor
  expect_error(validate.Gamma.arguments(data=gamma.dataset,col.control=col.control,
                                        m=4,gamma="gamma",
                                        strata=rep(1,nrow(gamma.dataset)),
                                        gamma.factor=c(3,4,5),Call=Call))
  
  expect_error(validate.Gamma.arguments(data=gamma.dataset,col.control=col.control,
                                        m=4,gamma="gamma",
                                        strata=rep(1,nrow(gamma.dataset)),
                                        gamma.factor="hello",Call=Call))
  
  
  #Call
  Call <- call("mycall",event.model=~Z1+Z2,subset="x==8")
  expect_error(validate.Gamma.arguments(data=gamma.dataset,surv.times=surv.times,
                                       m=2,gamma="gamma",
                                       strata=rep(1,nrow(gamma.dataset)),
                                       gamma.factor=1,DCO.time="DCO.time",Call=Call))
  
  Call <- call("mycall",event.model=~Z1+Z2,na.action="boo")
  expect_error(validate.Gamma.arguments(data=gamma.dataset,surv.times=surv.times,
                                        m=2,gamma="gamma",
                                        strata=rep(1,nrow(gamma.dataset)),
                                        gamma.factor=1,DCO.time="DCO.time",Call=Call))
  
  
  #strata (errors will be caught by boot if not caught here)
  df <- gamma.dataset
  expect_error(validate.Gamma.arguments(data=df,surv.times=surv.times,
                                        m=2,gamma="gamma",
                                        strata=rep(1,nrow(gamma.dataset)+6),
                                        gamma.factor=1,DCO.time="DCO.time",Call=Call))

})


test_that("negative_time",{
  load("gamma_test.rda")
  
  surv.times <- as.matrix(model.frame(formula(Surv(Yi,delta)~1),data=gamma.dataset))
  Call <- call("mycall",event.model=~Z1+Z2)
  
  surv.times[1,1] <- -8
  
  df <- gamma.dataset
  df$internal_gamma_val <- rep(1,nrow(df))
  expect_error(validate.Gamma.arguments(data=df,surv.times=surv.times,
                                        m=2,gamma="gamma",
                                        strata=rep(1,nrow(gamma.dataset)),
                                        gamma.factor=1,DCO.time="DCO.time",Call=Call))
})

test_that("validate_gamma_arguments_data",{
  load("gamma_test.rda")
  
  surv.times <- as.matrix(model.frame(formula(Surv(Yi,delta)~1),data=gamma.dataset))
  
  Call <- call("mycall",event.model=~Z1+Z2)
  
  #data 
  expect_error(validate.Gamma.arguments(data=gamma.dataset[numeric(0),],surv.times=surv.times,
                                        m=2,gamma="gamma",
                                        strata=rep(1,nrow(gamma.dataset)),
                                        gamma.factor=1,DCO.time="DCO.time",Call=Call))
  
  
  df <- gamma.dataset
  df$internal_gamma_val <- rep(1,nrow(df))
  expect_error(validate.Gamma.arguments(data=df,surv.times=surv.times,
                                        m=2,gamma="gamma",
                                        strata=rep(1,nrow(gamma.dataset)),
                                        gamma.factor=1,DCO.time="DCO.time",Call=Call))
  
  df <- gamma.dataset
  df$impute.time <- rep(1,nrow(df))
  expect_error(validate.Gamma.arguments(data=df,surv.times=surv.times,
                                        m=2,gamma="gamma",
                                        strata=rep(1,nrow(gamma.dataset)),
                                        gamma.factor=1,DCO.time="DCO.time",Call=Call))
  
  
  df <- gamma.dataset
  df$gamma <- rep("HELLO",nrow(df))
  expect_error(validate.Gamma.arguments(data=df,surv.times=surv.times,
                                        m=2,gamma="gamma",
                                        strata=rep(1,nrow(gamma.dataset)),
                                        gamma.factor=1,DCO.time="DCO.time",Call=Call))
  
  
})

test_that("validate_DCO.time",{
  load("gamma_test.rda")
  
  surv.times <- as.matrix(model.frame(formula(Surv(Yi,delta)~1),data=gamma.dataset))
  
  Call <- call("mycall",event.model=~Z1+Z2)
  
  #No DCO time column
  df <- gamma.dataset
  expect_error(validate.Gamma.arguments(data=df,surv.times=surv.times,
                                        m=2,gamma="gamma",
                                        strata=rep(1,nrow(gamma.dataset)),
                                        gamma.factor=1,DCO.time="O.time",Call=Call))
  
  #Inf DCO time
  df$DCO.time[1] <- Inf
  expect_error(validate.Gamma.arguments(data=df,surv.times=surv.times,
                                        m=2,gamma="gamma",
                                        strata=rep(1,nrow(gamma.dataset)),
                                        gamma.factor=1,DCO.time="DCO.time",Call=Call))
  
  #DCO less than time
  df <- gamma.dataset
  df$DCO.time[5] <- 0.5*df$Yi[5] 
  expect_error(validate.Gamma.arguments(data=df,surv.times=surv.times,
                                        m=2,gamma="gamma",
                                        strata=rep(1,nrow(gamma.dataset)),
                                        gamma.factor=1,DCO.time="DCO.time",Call=Call))
  
  #invalid length of DCO.time  
  expect_error(validate.Gamma.arguments(data=df,surv.times=surv.times,
                                        m=2,gamma="gamma",
                                        strata=rep(1,nrow(gamma.dataset)),
                                        gamma.factor=1,DCO.time=c(1,2),Call=Call))
  
  
})


test_that("Validate_gamma",{
  
  load("gamma_test.rda")
  
  #gamma length incorrect for character string
  expect_error(gammaImpute(formula=Surv(Yi,delta)~Z+strata(W1),
                          gamma.dataset,DCO.time="DCO.time",m=2,gamma=rep("gamma",nrow(gamma.dataset)),
                          bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                          gamma.factor=1))
  
  #gamma length incorrect for vector
  expect_error(gammaImpute(formula=Surv(Yi,delta)~Z+strata(W1),
                           gamma.dataset,DCO.time="DCO.time",m=2,gamma=1,
                           bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                           gamma.factor=1))
  
  #cannot be single number
  expect_error(gammaImpute(formula=Surv(Yi,delta)~Z+strata(W1),
                           gamma.dataset,DCO.time="DCO.time",m=2,gamma=1,
                           bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                           gamma.factor=1))
  
  expect_error(gammaImpute(formula=Surv(Yi,delta)~Z+strata(W1),
                           gamma.dataset,DCO.time="DCO.time",m=2,gamma=c(1,6),
                           bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                           gamma.factor=1))
  
  #missing col
  expect_error(gammaImpute(formula=Surv(Yi,delta)~Z+strata(W1),
                           gamma.dataset,DCO.time="DCO.time",m=2,gamma="missing",
                           bootstrap.strata=strata(gamma.dataset$Z,gamma.dataset$W1),
                           gamma.factor=1))
  
})
