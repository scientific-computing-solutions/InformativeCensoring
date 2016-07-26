#Code used to generate time independent score imputation data set
#code for time dependent dataset is below


Gen.Cova <- function(N){
  data.frame(Z1=rbinom(n=N,size=1,prob=0.5),
             Z2=runif(n=N, min = 0,max=1),
             Z3=rbinom(n=N,size=1,prob=0.5),
             Z4=runif(n=N, min = 0,max=1),
             Z5=rbinom(n=N,size=1,prob=0.5))
}

Get.times <- function(degree,values){
  ((degree+1)*(-log(runif(length(values)))/values))^(1/(degree+1))
}

Gen.Event <- function(df){
  TRT <- 1-df$arm
  val <- exp(0.75*TRT-2*df$Z1+0.5*df$Z2-2*df$Z3+2*df$Z4+2*df$Z5)
  Get.times(4,val)
}

Gen.Censor <- function(df){
  TRT <- 1-df$arm
  val <- exp(-3*(TRT+0.1)*df$Z1+0.5*df$Z2-2*(TRT+0.1)*df$Z3+1.5*df$Z4+ 2*(TRT+0.1)*df$Z5 )
  Get.times(3,val)
}



set.seed(1120)

N.placebo <- 200
N.active <- 200
N <- N.placebo+N.active
df <- data.frame(Id=1:N,
                 arm=c(rep(0,N.placebo),rep(1,N.active)))

df <- cbind(df,Gen.Cova(N))

event.times <- Gen.Event(df)
censor.times <- Gen.Censor(df)

df$event <- ifelse(event.times < censor.times,1,0)
df$time <- pmin(event.times,censor.times)
df$arm <- as.factor(df$arm)
df$to.impute <- df$event==0
df$DCO.time <- df$time*runif(N,min = 1,max=2)

ScoreInd <- df

#####################################################

#Generate time dependent covariate data frame

getW1 <- function(num){
  U <- rpois(n=1,lambda = 3)
  if(U>num){
    return(rep(0,num))
  }
  if(U==0){
    return(rep(1,num))
  }
  c(rep(0,U-1),rep(1,num-U+1))
  
}

getW2 <- function(num){
  start <- 1+rgamma(n = 1,shape=0.5,scale=1)
  start*(1:num)^1.2
}



visit <- 0.2  #in years

ans <- lapply(1:N,function(x){
  num <- 1 + floor(df[x,"time"]/visit)
  start <- visit*(0:(num-1))
  end <- start + visit
  end[num] <- df[x,"time"]
  
  my.df <- data.frame(Id=rep(x,num),
                      start=start,
                      end=end,
                      W1=getW1(num),
                      W2=getW2(num))
})

ScoreTimeDep <-  do.call("rbind", ans) 
rownames(ScoreTimeDep) <- NULL

devtools::use_data(ScoreInd,overwrite = TRUE)
devtools::use_data(ScoreTimeDep,overwrite = TRUE)
