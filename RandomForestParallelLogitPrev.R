#' Running RandomForest in parallel for max prevalence (if pathogen spreads to another individual, what
#' determines severity of outbreak?)
#' November 21, 2017
#' @author Lauren White, whit1951@umn.edu

#Load libraries
library(parallel)
library(doParallel)
library(foreach)
library(randomForest)

#Detect cores
ncores <- parallel::detectCores()
doParallel::registerDoParallel(ncores)

# setwd("C:/Users/law2y/OneDrive/R Code/IBM/August 2017")
# merged_data<-read.csv("IBMsummary.csv")

#Load data
#source('/home/forester/whit1951/IBMsummary.csv')
merged_data<-read.csv("IBMsummary.csv")
merged_data<-merged_data[-which(merged_data$beta2==5),] #Take out saturating beta2
merged_data<-merged_data[which(merged_data$max_I>1),]

# x <- matrix(runif(500), 100)
# y <- gl(2, 50)

x<-merged_data[,c(3:11)] #covariates
#y<-merged_data[,12] #response variable, duration
y<-merged_data[,14] #response variable, max_prev


tic=Sys.time()
set.seed(123)
rf <- foreach(ntree=rep(500, 20), .combine=combine, .multicombine=TRUE,
              .packages='randomForest') %dopar% {
                randomForest(x, y, ntree=ntree, importance=TRUE)
              }
print(difftime(Sys.time(),tic,units="mins"))
save(rf, file = "RF_logit_prev.RData")

#varImpPlot(rf, type=2, main= "Duration")

imp_maxprev<-rf$importance
write.csv(imp_maxprev, "logit_imp_maxprev.csv")


rm(list=ls(all=T)) #clear workspace