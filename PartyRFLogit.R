#Load libraries
library(party)

#Load data
merged_data<-read.csv("IBMsummary.csv")
merged_data<-merged_data[-which(merged_data$beta2==5),] #Take out saturating beta2

#Is an outbreak successful or not?
merged_data$outbreak<-ifelse(merged_data$max_I>1, 1,0)

set.seed(213)

tic=Sys.time()
fit.cf1 <- cforest(outbreak ~ p + H+ beta1 +beta2 +beta3+ density + rec_rate + percep, data=merged_data, controls=cforest_unbiased(ntree=1000))
v1 <- varimp(fit.cf1, conditional= FALSE)
v1<-v1[order(v1)]
write.csv(v1, "partyRF_logit1000.csv")

print(difftime(Sys.time(),tic,units="mins"))


rm(list=ls(all=T)) #clear workspace