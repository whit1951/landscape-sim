#Run IBM_loop_function in parallel via clusters
#August 24, 2017
#Updated with \beta_3\ values

library(parallel)
library(foreach)
library(doParallel)
library(SDMTools)

#load functions
#source('C:/Users/law2y/OneDrive/R Code/IBM/May 2017/MSI Code/sim_parallel_functions.R')
source('/home/forester/whit1951/sim_parallel_functions_2017_08_22.R')

maxtime<-10000
nsim<-100
k<-8 #dimensions of landscape= 2^k+1
density<-0.25
inf_prob<-0.2
rec_rate<-0.1
percep<- 1
p<-0.25 #0.5, 0.75
H<-c(0.1, 0.5, 0.9)
beta1<-c(0, 3, 6)
beta2<-c(0, 5)
beta3<-c(0, -1)

params<-expand.grid(maxtime=maxtime, nsim=nsim, k=k, density=density, inf_prob=inf_prob, rec_rate=rec_rate, p=p, H=H, beta1=beta1, beta2=beta2, beta3=beta3, percep=percep)
params<-params[-which(params$beta2==5 & params$beta3==0),]
#summary_data<-list(NA)

## Apply the declared function in parallel

ncores <- parallel::detectCores()
doParallel::registerDoParallel(ncores)
tic=Sys.time()
foreach(n = 1:nrow(params)) %dopar% IBM_loop_function(params[n,]$maxtime, params[n,]$nsim, params[n,]$k, params[n,]$density, params[n,]$inf_prob, params[n,]$rec_rate, params[n,]$p, params[n,]$H, params[n,]$beta1, params[n,]$beta2, params[n,]$beta3, params[n,]$percep)
print(difftime(Sys.time(),tic,units="mins"))

#save(summary_data, file="move_sim1.RData")
#save(params, file="parameters_sim1.RData")

rm(list=ls(all=T)) #clear workspace