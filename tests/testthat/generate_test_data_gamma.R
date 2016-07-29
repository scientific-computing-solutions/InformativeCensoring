#code for generating gamma imputation test set, this is not exported
#for users but is used in the testthat tests

set.seed(6110)

N <- 1000

Z <- sample(x=0:2,size = N,replace = TRUE,prob = c(0.5,0.3,0.2))

Ci <- rexp(n = N,rate = 0.3)
Ti <- vapply(Z,function(x){
  rate <- 0.03
  if(x == 1) rate <- 0.05
  if(x == 2) rate <- 0.09
  rexp(1,rate=rate)
} ,FUN.VALUE = numeric(1))

Yi <- pmin(Ti,Ci,3)
delta <- (Ti < Ci & Ti <3)

W1 <- rbinom(n = N,size=1,p=0.5)
W2 <- rbinom(n = N,size=1,p=0.5)

gamma.dataset <- data.frame(Id=1:N,
                    Yi=Yi,
                    delta=delta,
                    Z=factor(Z),
                    W1=W1,
                    W2=W2,
                    to.impute=(delta==0 & Yi <3),
                    DCO.time=rep(3,N),
                    gamma=rep(1,N))

#save(gamma.dataset,file="gamma_test.rda")
#must move file into testthat directory
