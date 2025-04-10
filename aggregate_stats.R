library(future)
library(future.apply)
library(ggplot2)
library(tidyverse)
library(tibble)
library(progressr)
library(data.table)
library(dplyr)

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
plan(sequential)

#plots
a=ggplot(results_df_long,aes(x=mean_diff))+geom_histogram()+facet_wrap(~comparison)

df_quantiles=results_df_long%>%
  select(unit,comparison,starts_with("qdiff"))%>%
  pivot_longer(cols=starts_with("qdiff"),names_to="quantile",names_prefix = "qdiff_",values_to="value")

b=ggplot(df_quantiles,aes(y=value,x=quantile,col=quantile,fill=quantile))+geom_violin()+facet_wrap(~comparison)+theme_bw()+
  stat_summary(fun="mean",geom="point",color="black")

x_grid <- 10^seq(-2, 3.5, length.out = 750)
#next compare: mean bank values (i.e. reasonable mitigation) against different quantiles of loss distribution
#look at ratios - i.e. what is ratio of flood benefit at different parts of the wetland loss distribution?
process_unit <- function(unit_name, datadir=datadir, bankfiles=bankfiles, lossfiles=lossfiles, x_grid=x_grid) {
  # File paths
  bankfile <- paste0(datadir, "\\Bank ECDF Functions\\bank_ecdf_fn_", unit_name, ".rds")
  lossfile <- paste0(datadir, "\\SA ECDF Functions - Loss\\ecdf_fn_", unit_name, ".rds")
  
  # Load data if files exist
  bank_data <- if (bankfile %in% bankfiles) readRDS(bankfile) else NA
  loss_data <- if (lossfile %in% lossfiles) readRDS(lossfile) else NA
  
  # Validate required objects
  if (!is.list(bank_data) || !is.list(loss_data)) return(NULL)
  if (!is.function(bank_data$housing_units) || !is.function(loss_data$housing_units)) return(NULL)
  if (!is.function(bank_data$housing_value) || !is.function(loss_data$housing_value)) return(NULL)
  
  # Extract raw values
  bank_vals  <- environment(bank_data$housing_units)$x
  bank_hval  <- environment(bank_data$housing_value)$x
  loss_vals  <- environment(loss_data$housing_units)$x
  loss_hval  <- environment(loss_data$housing_value)$x
  
  # Compute means
  bank_mean_vals <- mean(bank_vals, na.rm = TRUE)
  bank_mean_hval <- mean(bank_hval, na.rm = TRUE)
  
  # Normalize accordingly
  loss_vals_norm <- loss_vals / bank_mean_vals
  loss_hval_norm <- loss_hval / bank_mean_hval
  
  n_obs <- length(loss_vals_norm)
  
  # Optional downsampling
  if (n_obs > 1e5) {
    sampled_indices <- sample(seq_len(n_obs), 1e5)
    loss_vals_norm <- loss_vals_norm[sampled_indices]
    loss_hval_norm <- loss_hval_norm[sampled_indices]
  }
  
  # Evaluate ECDFs
  ecdf_vals  <- ecdf(loss_vals_norm)
  ecdf_hvals <- ecdf(loss_hval_norm)
  
  # Output both sets as tibbles
  result_units <- tibble(
    unit = unit_name,
    comparison = "normalized_loss_units",
    x = x_grid,
    y = ecdf_vals(x_grid),
    n_obs = n_obs
  )
  
  result_values <- tibble(
    unit = unit_name,
    comparison = "normalized_loss_values",
    x = x_grid,
    y = ecdf_hvals(x_grid),
    n_obs = n_obs
  )
  
  # Return combined result
  bind_rows(result_units, result_values)
}
plan(multisession)  # or sequential/multicore/etc.
handlers("txtprogressbar")  # or "progress" for shiny/fancier

with_progress({
  p <- progressor(along = names)
  
  results_long <- future_lapply(names, function(n) {
    p(message = n)
    process_unit(n, datadir, bankfiles, lossfiles, x_grid)
  })
})

results_long <- Filter(Negate(is.null), results_long)
results_long_normalized <- bind_rows(results_long)

fwrite(results_long_normalized,file=paste0(datadir,"\\normalized_ecdfs_loss.csv"))

ecdf_summary <- results_long_normalized %>%
  group_by(x,comparison) %>%
  summarize(
    mean_y = weighted.mean(x=y,w=n_obs)
  )

newlabs=c("Flood Protection: Housing Units","Flood Protection: Housing Values")
names(newlabs)=c("normalized_loss_units","normalized_loss_values")

a=ggplot(
  results_long_normalized,
  aes(x = x, y = y,group=unit,lwd=n_obs,col=comparison)) + geom_line(alpha=0.2)+
  geom_vline(xintercept = 1, linetype = "dashed",col="grey10") + geom_line(data=ecdf_summary,aes(x=x,y=mean_y),col="black",inherit.aes=FALSE)+
  scale_x_log10() + facet_wrap(~comparison, labeller=labeller(comparison=newlabs))+labs(y="Cumulative Density",x="Value Relative to Bank Mean")+
  theme_bw()+scale_linewidth_continuous(guide="none")+scale_color_discrete(guide="none")+theme(strip.background =element_rect(fill="white"))+
  geom_vline(xintercept = 1.5, linetype = "dashed",col="grey50")+geom_vline(xintercept = 2, linetype = "dashed",col="grey80")

