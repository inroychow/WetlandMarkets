#model flood risk time trends as a function of wetland flood value loss / difference with bank

#first assemble bank statistics dataset
loss=readRDS("data/Summary data/loss_summary.rds")
colnames(loss)[1:25]=paste0(colnames(loss)[1:25],"_loss")

bank=readRDS("data/Summary data/bank_summary.rds")
colnames(bank)[2:28]=paste0(colnames(bank)[2:28],"_bank")

wetlanddata=merge(loss,bank,by.x="sa_id",by.y="bank_name")

#read in coefficients shorter time period and coverage control
coefs=read.csv(file="data/hurdle_model_2009_2021_censusFEs_policycontrol.csv")
#merge in wetland data
coefs=merge(coefs,wetlanddata,by.x="banks_upstream",by.y="sa_id",all=FALSE)
coefs$weight=1/sqrt(coefs$std.error)

#greater flood damages over time, controling for NFIP coverage, is associated with estimated upstream lost wetland value
mod_zero=lm(estimate~I(log(sum_housing_value_loss)),data=coefs%>%filter(stage=="zero"&sum_housing_value_loss>0),weights = coefs$weight[which(coefs$stage=="zero"&coefs$sum_housing_value_loss>0)])
mod_pos=lm(estimate~I(log(sum_housing_value_loss)),data=coefs%>%filter(stage=="positive"&sum_housing_value_loss>0),weights = coefs$weight[which(coefs$stage=="positive"&coefs$sum_housing_value_loss>0)])

#now look at any additional effect of differential flood valuation in lost wetlands vs bank
#the larger this variable, the larger trends in flood damages we would expect
coefs$diff_val=coefs$sum_housing_value_bank-coefs$sum_housing_value_loss
#normalize 
coefs$diff_val=(coefs$diff_val-mean(coefs$diff_val,na.rm=T))/sd(coefs$diff_val,na.rm=T)

mod_zero=lm(estimate~I(log(sum_housing_value_loss))+diff_val,data=coefs%>%filter(stage=="zero"&sum_housing_value_loss>0),weights = coefs$weight[which(coefs$stage=="zero"&coefs$sum_housing_value_loss>0)])
mod_pos=lm(estimate~I(log(sum_housing_value_loss))+diff_val,data=coefs%>%filter(stage=="positive"&sum_housing_value_loss>0),weights = coefs$weight[which(coefs$stage=="positive"&coefs$sum_housing_value_loss>0)])

summary(mod_zero)
summary(mod_pos)
#could add a control for total remaining flood value in service area?
#look at combined effect from hurdle model
         