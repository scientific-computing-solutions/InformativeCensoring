context("ScoreStatSet")

test_that("from_matrix_invalid",{
  #wrong number of columns
  expect_error(ScoreStatSet(matrix(rep(1,20),ncol = 2)))
  expect_error(ScoreStatSet(matrix(rep(1,20),ncol = 4)))

  #fewer than 5 rows
  expect_error(ScoreStatSet(matrix(rep(1,12),ncol = 3)))
  
  #negative variances
  expect_error(ScoreStatSet(matrix(c(rep(1,9),-1,rep(1,5)),ncol = 3)))
  
  #Z must equal estimate/sqrt(var)
  m <- matrix(rep(1,15),ncol=3)
  m[1,3] <- 2
  expect_error(ScoreStatSet(m))
})

.get.test.matrix <- function(){
  matrix(c(0.5,-0.6,0.7,0.8,0.9,
        0.1^2,0.11^2,0.12^2,0.13^2,0.14^2,
        0.5/0.1,-0.6/0.11,0.7/0.12,0.8/0.13,0.9/0.14),ncol = 3)}


test_that("from_matrix",{
  m <- .get.test.matrix()
  sss <- ScoreStatSet(m)
  
  expect_equal("ScoreStatSet",class(sss))
  expect_equal(3,ncol(sss))
  expect_equal(c(0.5,-0.6,0.7,0.8,0.9),sss[,"estimate"])
  expect_equal(m[,2],sss[,"var"])
  expect_equal(m[,3],sss[,"Z"])
})


#check the summary.ScoreStatSet creation works
test_that("summary_creation",{
  m <- .get.test.matrix()
  s <- summary(ScoreStatSet(m))
  expect_equal("summary.ScoreStatSet",class(s))
  expect_equal(3,length(s))
  
  expect_equal(c("meth1","meth2","methRubin"),names(s))
  
  meth1 <- s$meth1
  meth2 <- s$meth2
  methRubin <- s$methRubin
  
  expect_equal(c("estimate","var","test.stat","df","distribution","p.value"),names(meth1))
  expect_equal(names(meth1),names(meth2))
  expect_equal(names(meth1),names(methRubin))
  
  expect_equal("F",meth1$distribution)
  expect_equal("t",meth2$distribution)
  expect_equal("t",methRubin$distribution)
})

test_that("summary_all_rows_same",{
  m <- matrix(rep(c(1,0.25,2),5),ncol=3,byrow = TRUE)
  #get a warning here as v1 is undefined (yes, this is caught...)
  expect_warning(s5 <- summary(ScoreStatSet(m)))
  
  m <- matrix(rep(c(1,0.25,2),6),ncol=3,byrow = TRUE)
  s6 <- summary(ScoreStatSet(m))
  m <- matrix(rep(c(1,0.25,2),10),ncol=3,byrow = TRUE)
  s10 <- summary(ScoreStatSet(m))
  
  expect_equal(s6,s10) #these should be the same 
  expect_equal(s5$meth2,s6$meth2) #s5[[1]] has NAN df -> this is correct 0/0 = NAN
  
  expect_equal(1,s6$meth1$estimate)
  expect_equal(0.25,s6$meth1$var)
  expect_equal(4,s6$meth1$test.stat)
  expect_equal(c(1,Inf),s6$meth1$df)
  expect_equal(2*(1-pnorm(abs(2))),s6$meth1$p.value)
  
  
  expect_equal(2,s6$meth2$estimate)
  expect_equal(1,s6$meth2$var)
  expect_equal(2,s6$meth2$test.stat)
  expect_true(is.infinite(s6$meth2$df))
  expect_equal(2*(1-pnorm(abs(2))),s6$meth2$p.value)
  
  expect_equal(1,s6$methRubin$estimate)
  expect_equal(0.25,s6$methRubin$var)
  expect_equal(2,s6$methRubin$test.stat)
  expect_true(is.infinite(s6$methRubin$df))
  expect_equal(2*(1-pnorm(abs(2))),s6$methRubin$p.value)
  
  
})

test_that("summary_system_test",{
  m <- .get.test.matrix()
  s <- summary(ScoreStatSet(m))
  
  expect_equal(0.46,s$meth1$estimate)
  #0.0146 + 6*var(c(0.5,-0.6,0.7,0.8,0.9))/5
  expect_equal(0.4622,s$meth1$var)

  #0.46*0.46/0.4622
  expect_equal(0.457810471657,s$meth1$test.stat)
  #as (t-4) = 0
  expect_equal(c(1,4),s$meth1$df)
  expect_equal(1-pf(0.457810471657,1,4),s$meth1$p.value)
  
  me <- mean(c(0.5/0.1,-0.6/0.11,0.7/0.12,0.8/0.13,0.9/0.14))
  expect_equal(me,s$meth2$estimate)
  
  va <- var(c(0.5/0.1,-0.6/0.11,0.7/0.12,0.8/0.13,0.9/0.14))
  expect_equal(1+(6*va/5),s$meth2$var)
  
  expect_equal( me/sqrt(1+(6*va/5)),s$meth2$test.stat)
  
  expect_equal( 4*(1+5/(va*6))^2,s$meth2$df)
  expect_equal( 2*(1-pt(me/sqrt(1+(6*va/5)),4*(1+5/(va*6))^2)),s$meth2$p.value  )
})

test_that("methRubin",{
  m <- .get.test.matrix()
  s <- summary(ScoreStatSet(m))
  
  expect_equal(s$methRubin$estimate,s$meth1$estimate)
  expect_equal(s$methRubin$var,s$meth1$var)
  expect_equal(s$methRubin$test.stat^2,s$meth1$test.stat)
  
  #4*(1+0.0146/(6*0.373/5))^2
  expect_true(abs(4.265203-s$methRubin$df)<1e-6)
})

test_that("meth1.df",{
  
  m <- matrix(c(5:12,seq(0.1,0.8,0.1),rep(1,8)),ncol=3)
  m[,3] <- m[,1]/sqrt(m[,2])
  s <- summary(ScoreStatSet(m))
  
  r <- (9*6)/(0.45*8)
  
  expect_equal(4+3*(1+5/(7*r))^2, s$meth1$df[2])
  
})

test_that("confint",{
  m <- matrix(c(5:12,seq(0.1,0.8,0.1),rep(1,8)),ncol=3)
  m[,3] <- m[,1]/sqrt(m[,2])
  s <- summary(ScoreStatSet(m))
  
  expect_error(confint(s,level=Inf))
  expect_error(confint(s,level=0))
  expect_error(confint(s,level=1))
  expect_error(confint(s,level="Hello"))
  
  ans <- confint(s,level=0.001)
  expect_true(s$methRubin$estimate > ans[1])
  expect_true(s$methRubin$estimate < ans[2])
  
  ans.95 <- confint(s)
  expect_equal(c("2.5%","97.5%"),names(ans.95))
  #symmetric
  expect_true(all.equal(s$methRubin$estimate,mean(ans.95)))
  
})