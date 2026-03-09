#model flood risk time trends as a function of wetland flood value loss / difference with bank
library(tidyverse)
library(lfe)
library(modelsummary)
library(fixest)

#first assemble bank statistics dataset
loss=readRDS("data/Summary data/loss_summary.rds")
colnames(loss)[1:25]=paste0(colnames(loss)[1:25],"_loss")

bank=readRDS("data/Summary data/bank_summary.rds")
colnames(bank)[2:28]=paste0(colnames(bank)[2:28],"_bank")

wetlanddata=merge(loss,bank,by.x="sa_id",by.y="bank_name")

#read in coefficients shorter time period and coverage control
coefs=read.csv(file="data/hurdle_model_2009_2021_censusFEs_policycontrol_weighted.csv")
#merge in wetland data
coefs=merge(coefs,wetlanddata,by.x="banks_upstream",by.y="sa_id",all=FALSE)
coefs$weight=1/sqrt(coefs$std.error)

#cut spurious tail of weights
coefs$weight[which(coefs$weight>20)]=NA

#merge in clustering variable based on banks serving same service area
pooling=readRDS("data/pool_map.rds")
coefs=merge(coefs,pooling[,c(1,3)],by.x="banks_upstream",by.y="bank_name",all.x=TRUE,all.y=FALSE)
coefs$pool_id[which(is.na(coefs$pool_id))]=coefs$banks_upstream[which(is.na(coefs$pool_id))]

#greater flood damages over time, controling for NFIP coverage, is associated with estimated upstream lost wetland value
mod_zero_simple=feols(fml=estimate~I(log(sum_housing_value_loss)),data=coefs[which(coefs$stage=="zero"&coefs$sum_housing_value_loss>0),],weights = coefs$weight[which(coefs$stage=="zero"&coefs$sum_housing_value_loss>0)],cluster=coefs$pool_id[which(coefs$stage=="zero"&coefs$sum_housing_value_loss>0)])
mod_pos_simple=feols(fml=estimate~I(log(sum_housing_value_loss)),data=coefs[which(coefs$stage=="positive"&coefs$sum_housing_value_loss>0&is.finite(coefs$weight)),],weights = coefs$weight[which(coefs$stage=="positive"&coefs$sum_housing_value_loss>0&is.finite(coefs$weight))],cluster=coefs$pool_id[which(coefs$stage=="positive"&coefs$sum_housing_value_loss>0&is.finite(coefs$weight))])

#now look at any additional effect of differential flood valuation in lost wetlands vs bank
#the larger this variable, the larger trends in flood damages we would expect
coefs$diff_val=coefs$median_housing_value_bank-coefs$median_housing_value_loss
#normalize 
coefs$diff_val=(coefs$diff_val-mean(coefs$diff_val,na.rm=T))/sd(coefs$diff_val,na.rm=T)

mod_zero=feols(estimate~I(log(sum_housing_value_loss))*diff_val,data=coefs%>%filter(stage=="zero"&sum_housing_value_loss>0),weights = coefs$weight[which(coefs$stage=="zero"&coefs$sum_housing_value_loss>0)],cluster=coefs$pool_id[which(coefs$stage=="zero"&coefs$sum_housing_value_loss>0)])
mod_pos=feols(estimate~I(log(sum_housing_value_loss))*diff_val,data=coefs%>%filter(stage=="positive"&sum_housing_value_loss>0),weights = coefs$weight[which(coefs$stage=="positive"&coefs$sum_housing_value_loss>0)],cluster=coefs$pool_id[which(coefs$stage=="positive"&coefs$sum_housing_value_loss>0)])

#output summary table

mods=list("Simple Model, Binary"=mod_zero_simple,"Simple Model, Flood Damage"=mod_pos_simple, "Full Model, Binary"=mod_zero,"Full Model, Flood Damage"=mod_pos)

msummary(mods,
         output="summarytable.docx",
         stars=TRUE,
         coef_omit=1,
         coef_rename=c("Log Lost Flood Protection/n(Housing Value)","Bank Median-Loss Median","Log Lost Flood Protection*/n(Bank Median-Loss Median)"))
