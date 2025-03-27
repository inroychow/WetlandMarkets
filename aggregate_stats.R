library(future)
library(future.apply)
library(ggplot2)
library(tidyverse)

datadir="C:\\Users\\fmoore\\Box\\ECDF Functions"

#get unique bank names
names=list.files(paste0(datadir,"\\SA ECDF Functions"))
names = sub("^ecdf_fn_(.*)\\.rds$", "\\1", names)

#read in an example - bank, SA and SA loss
testname=names[754] #randomly chosen bank

bankfile=paste0(datadir,"\\Bank ECDF Functions\\bank_ecdf_fn_",testname,".rds")
banktest=readRDS(bankfile)

safile=paste0(datadir,"\\SA ECDF Functions\\ecdf_fn_",testname,".rds")
satest=readRDS(safile)

lossfile=paste0(datadir,"\\SA ECDF Functions - Loss\\ecdf_fn_",testname,".rds")
losstest=readRDS(lossfile)

plot(losstest$housing_units)
plot(banktest$housing_units,add=T,col="red")
plot(satest$housing_units,add=T,col="blue")
legend("bottomright",legend=c("Bank","Service Area","Wetland Loss"),col=c("red","blue","black"),bty="n",lwd=1)

#mean values
bankmean=mean(environment(banktest$housing_units)$x)
lossmean=mean(environment(losstest$housing_units)$x)
samean=mean(environment(satest$housing_units)$x)

#calculate differences in means, medians and area between ECDFs for:
#1. lost wetlands vs service areas
#2. lost wetlands vs banks
#3. service areas vs banks

#following written with assistance from ChatGPT
# Functions
get_ecdf_data <- function(ecdf_obj) {
  dat=environment(ecdf_obj)$x
  if(length(dat)>1e6) dat=sample(dat,1e6,replace=FALSE)
  return(dat)
}

signed_area_diff <- function(ecdf2, ecdf1, grid_size = 1000) {
  x1 <- get_ecdf_data(ecdf1)
  x2 <- get_ecdf_data(ecdf2)
  x_range <- range(c(x1, x2))
  grid <- seq(x_range[1], x_range[2], length.out = grid_size)
  delta_x <- diff(range(grid)) / (grid_size - 1)
  sum(ecdf1(grid) - ecdf2(grid)) * delta_x
}

quantile_diff <- function(x1, x2, probs = c(0.25, 0.5, 0.75,0.95,0.99)) {
  q1 <- quantile(x1, probs = probs)
  q2 <- quantile(x2, probs = probs)
  setNames(q1 - q2, paste0("qdiff_", probs*100))
}

make_row <- function(unit_name, label, ecdf1, ecdf2) {
  if (!is.function(ecdf1) || !is.function(ecdf2)) {
    return(data.frame(
      unit = unit_name,
      comparison = label,
      mean_diff = NA,
      area_diff = NA,
      qdiff_25 = NA,
      qdiff_50 = NA,
      qdiff_75 = NA,
      qdiff_95 = NA,
      qdiff_99 = NA
    ))
  }
  
  x1 <- get_ecdf_data(ecdf1)
  x2 <- get_ecdf_data(ecdf2)
  qdiffs <- quantile_diff(x1, x2)
  
  data.frame(
    unit = unit_name,
    comparison = label,
    mean_diff = mean(x1) - mean(x2),
    area_diff = signed_area_diff(ecdf1, ecdf2),
    qdiff_25 = qdiffs["qdiff_25"],
    qdiff_50 = qdiffs["qdiff_50"],
    qdiff_75 = qdiffs["qdiff_75"],
    qdiff_95 = qdiffs["qdiff_95"],
    qdiff_99 = qdiffs["qdiff_99"]
  )
}

bankfiles=list.files(paste0(datadir,"\\Bank ECDF Functions\\"),full.names=TRUE)
safiles=list.files(paste0(datadir,"\\SA ECDF Functions\\"),full.names=TRUE)
lossfiles=list.files(paste0(datadir,"\\SA ECDF Functions - Loss\\"),full.names=TRUE)

process_unit <- function(unit_name, datadir, bankfiles=bankfiles, safiles=safiles, lossfiles=lossfiles) {
  bankfile <- paste0(datadir, "\\Bank ECDF Functions\\bank_ecdf_fn_", unit_name, ".rds")
  safile   <- paste0(datadir, "\\SA ECDF Functions\\ecdf_fn_", unit_name, ".rds")
  lossfile <- paste0(datadir, "\\SA ECDF Functions - Loss\\ecdf_fn_", unit_name, ".rds")
  
  bank <- if (bankfile %in% bankfiles) readRDS(bankfile)$housing_units else NA
  sa   <- if (safile %in% safiles) readRDS(safile)$housing_units else NA
  loss <- if (lossfile %in% lossfiles) readRDS(lossfile)$housing_units else NA
  
  do.call(rbind, list(
    make_row(unit_name, "loss_vs_sa",   loss, sa),
    make_row(unit_name, "loss_vs_bank", loss, bank),
    make_row(unit_name, "sa_vs_bank",   sa, bank)
  ))
}
plan(multisession)

results_df_long <- future_lapply(
  names,
  FUN = function(n) process_unit(n, datadir, bankfiles, safiles, lossfiles)
)
results_df_long <- do.call(rbind, results_df_long)

#plots
a=ggplot(results_df_long,aes(x=mean_diff))+geom_histogram()+facet_wrap(~comparison)

df_quantiles=results_df_long%>%
  select(unit,comparison,starts_with("qdiff"))%>%
  pivot_longer(cols=starts_with("qdiff"),names_to="quantile",names_prefix = "qdiff_",values_to="value")

b=ggplot(df_quantiles,aes(y=value,x=quantile,col=quantile,fill=quantile))+geom_violin()+facet_wrap(~comparison)+theme_bw()+
  stat_summary(fun="mean",geom="point",color="black")

#next compare: mean bank values (i.e. reasonable mitigation) against different quantiles of loss distribution
#look at ratios - i.e. what is ratio of flood benefit at different parts of the wetland loss distribution?
  