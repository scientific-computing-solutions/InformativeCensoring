context("timedependent")

#not testing that data is contiguous or 
#that the pannelling is valid - 
#these are tested in the test-validation.R
#file
test_that("MakeTimeDepScore_invalid",{
  df <- data.frame(Subject=c(1,1,2,2),
                   c1=c(10,20,20,30),
                   my.start=c(0,10,0,10),
                   my.end=c(10,12,10,15))
  
  #invalid Id
  expect_error(MakeTimeDepScore(data=df,Id="Id",time.start="my.start",time.end="my.end"))
  
  #invalid time.start
  expect_error(MakeTimeDepScore(data=df,Id="Subject",time.start="c",time.end="my.end"))
 
  df$Id <- df$my.start
  
  #cannot have "Id" column and use it for non "Id" column
  expect_error(MakeTimeDepScore(data=df,Id="Subject",time.start="my.start",time.end="my.end"))
  
  df$Id <- df$Subject
  df$time.start <- df$my.start
  
  #cannot have "time.start" column if not used as time.start 
  expect_error(MakeTimeDepScore(data=df,Id="Id",time.start="my.start",time.end="my.end"))
  
  #test if not numeric time.start/time.end  
  df$time.start <- NULL
  
  df$my.start[1] <- "r"
  expect_warning(expect_error(MakeTimeDepScore(data=df,Id="Id",time.start="my.start",time.end="my.end")))
})

test_that("MakeTimeDepScore",{
  df <- data.frame(Subject=c(1,1,2,2),
                   c1=c(10,20,20,30),
                   my.start=c(0,10,0,10),
                   my.end=c(10,12,10,15),
                   W1=c(0,1,1,0))
  
  df$W1 <- factor(df$W1)
  
  td <- MakeTimeDepScore(data=df,Id="Subject",time.start="my.start",time.end="my.end")
  expect_equal(c("ScoreTD","data.frame"),class(td))
  
  expect_equal(c("c1","W1","Id","time.start","time.end"),colnames(td))
  expect_equal(factor(c(1,1,2,2)),td$Id) #ID has become a factor if it wasn't already
  expect_equal(c(0,10,0,10),td$time.start)
  expect_equal(c(10,20,20,30),td$c1)
  expect_equal(factor(c(0,1,1,0)),td$W1)
  
  df <- data.frame(Id=c(1,1,2,2),
                   c1=c(10,20,20,30),
                   time.start=c(0,10,0,10),
                   my.end=c(10,12,10,15),
                   W1=c(0,1,1,0))
  
  df$W1 <- factor(df$W1)
  td <- MakeTimeDepScore(data=df,Id="Id",time.start="time.start",time.end="my.end")
  
  expect_equal(c("Id","c1","time.start","W1","time.end"),colnames(td))
  expect_equal(df$my.end,td$time.end)
  expect_equal(df$Id,as.numeric(td$Id))
  
})


test_that(".getTimeDepDataSet_invalid",{
  #validation such as ensuring data and time.dep do not have same
  #colnames etc. has already taken place before .getTimeDepDataSet is 
  #called. my.time should also be valid for the given data set
  
  df <- data.frame(Id=c(1,1,2,2),
                   c1=c(10,20,20,30),
                   time.start=c(0,10,0,10),
                   my.end=c(10,12,10,15),
                   W1=c(0,1,1,0))
  
  baseline.df <- data.frame(Sub=c(1,2,3),
                            Z1=c(0,1,1),
                            time=c(12,15,20),
                            DCO.time=c(20,25,40),
                            has.event=c(1,0,0),
                            arm=c(1,1,1))
  td <- MakeTimeDepScore(data=df,Id="Id",time.start="time.start",time.end="my.end")
  
  #missing Id
  expect_error(.getTimeDepDataSet(baseline.df,td,"Sub","time",my.time=NULL))
  
  expect_that(.getTimeDepDataSet(baseline.df[1:2,],td,"Sub","time",my.time=NULL),not(throws_error()))
  
  #missing someone
  expect_error(.getTimeDepDataSet(baseline.df[1:2,],td,"Sub","time",my.time=14)) 
  
  
})

test_that(".getTimeDepDataSet_nullmy.time",{
  df <- data.frame(Id=c(1,1,2,2),
                   c1=c(10,20,20,30),
                   time.start=c(0,10,0,10),
                   my.end=c(10,12,10,15),
                   W1=c(0,1,1,0))
  df$W1 <- factor(df$W1)
  
  td <- MakeTimeDepScore(data=df,Id="Id",time.start="time.start",time.end="my.end")
  
  baseline.df <- data.frame(Sub=c(1,2),
                            Z1=c(0,1),
                            time=c(12,15),
                            DCO.time=c(20,25),
                            has.event=c(1,0),
                            arm=c(1,1)) 
  
  baseline.df$Sub <- factor(baseline.df$Sub) 
  ans <- .getTimeDepDataSet(baseline.df,td,"Sub","time",my.time=NULL)
  
  expect_equal(2,nrow(ans))
  expect_equal(c("Sub","Z1","time","DCO.time", "has.event","arm","c1","W1"),colnames(ans))
               
  expect_equal(c(20,30),ans$c1)             
  expect_equal(factor(c(1,0)),ans$W1)
  expect_equal(c(0,1),ans$Z1)
  
  #now reorder Sub the merge does not sort the data set
  baseline.df <- baseline.df[2:1,]
  ans <- .getTimeDepDataSet(baseline.df,td,"Sub","time",my.time=NULL)
  
  expect_equal(factor(2:1),ans$Sub)
  
})

test_that(".getTimeDepDataSet_my.time",{
  df <- data.frame(Id=c(1,1,2,2,3,3,4,4,4),
                   c1=c(10,20,20,30,45,10,1,4,3),
                   time.start=c(0,10,0,10,0,5,0,7,12),
                   my.end=c(10,12,10,15,5,15,7,12,15),
                   W1=c(0,1,1,0,1,0,0,1,1))
  
  baseline.df <- data.frame(Sub=c(1,2,3,4),
                            Z1=c(0,1,1,0),
                            time=c(12,15,15,15),
                            DCO.time=c(20,25,25,25),
                            has.event=c(1,0,0,0),
                            arm=c(1,1,1,1)) 
  
  td <- MakeTimeDepScore(data=df,Id="Id",time.start="time.start",time.end="my.end")
  ans <- .getTimeDepDataSet(baseline.df,td,"Sub","time",my.time=4)
  
  expect_equal(c(10,20,45,1),ans$c1)
  expect_equal(c(0,1,1,0),ans$W1)
  
  #also OK with Id in both columns
  baseline.df$Id <- baseline.df$Sub
  baseline.df$Sub <- NULL
  ans <- .getTimeDepDataSet(baseline.df,td,"Id","time",my.time=12)
  expect_equal(1:4,ans$Id)
  
  expect_equal(8,ncol(ans))
  
  expect_equal(c(20,30,10,4),ans$c1)
  expect_equal(c(1,0,0,1),ans$W1)
  
})

test_that("bootstrap_df_OK_with_timedep",{
  df <- data.frame(Id=c(1,1,2,2,3,3,4,4,4),
                   c1=c(10,20,20,30,45,10,1,4,3),
                   time.start=c(0,10,0,10,0,5,0,7,12),
                   my.end=c(10,12,10,15,5,15,7,12,15),
                   W1=c(0,1,1,0,1,0,0,1,1))
  
  baseline.df <- data.frame(Sub=c(1,2,3,4),
                            Z1=c(0,1,1,0),
                            time=c(12,15,15,15),
                            DCO.time=c(20,25,25,25),
                            has.event=c(1,0,0,0),
                            arm=c(1,1,1,1))
  
  baseline.df <- baseline.df[rep(1:4,2),]
  
  td <- MakeTimeDepScore(data=df,Id="Id",time.start="time.start",time.end="my.end")
  ans <- .getTimeDepDataSet(baseline.df,td,"Sub","time",my.time=4)
  
  expect_equal(rep(1:4,2),ans$Sub)
  expect_equal(c(0,1,1,0,0,1,1,0),ans$Z1)
  expect_equal(c(10,20,45,1,10,20,45,1),ans$c1)
  expect_equal(c(0,1,1,0,0,1,1,0),ans$W1)
})
