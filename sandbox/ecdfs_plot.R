library(future)
library(future.apply)
library(ggplot2)
library(tidyverse)
library(tibble)
library(progressr)
library(data.table)
library(dplyr)

#datadir = "your_directory"
datadir="C:\\Users\\indumati\\Box\\ECDF Functions"

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
  dplyr::select(unit,comparison,starts_with("qdiff"))%>%
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
  bank_pop   <- environment(bank_data$population)$x
  loss_vals  <- environment(loss_data$housing_units)$x
  loss_hval  <- environment(loss_data$housing_value)$x
  loss_pop   <- environment(loss_data$population)$x
  
  # Compute means
  bank_mean_vals <- mean(bank_vals, na.rm = TRUE)
  bank_mean_hval <- mean(bank_hval, na.rm = TRUE)
  bank_mean_pop <- mean(bank_pop, na.rm = TRUE)
  
  # Normalize accordingly
  loss_vals_norm <- loss_vals / bank_mean_vals
  loss_hval_norm <- loss_hval / bank_mean_hval
  loss_pop_norm <- loss_pop / bank_mean_pop
  
  n_obs <- length(loss_vals_norm)
  
  # Optional downsampling
  if (n_obs > 1e5) {
    sampled_indices <- sample(seq_len(n_obs), 1e5)
    loss_vals_norm <- loss_vals_norm[sampled_indices]
    loss_hval_norm <- loss_hval_norm[sampled_indices]
    loss_pop_norm <- loss_pop_norm[sampled_indices]
  }
  
  # Evaluate ECDFs
  ecdf_vals  <- ecdf(loss_vals_norm)
  ecdf_hvals <- ecdf(loss_hval_norm)
  ecdf_pop <- ecdf(loss_pop_norm)
  
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
  
  result_pop <- tibble(unit = unit_name,
                       comparison = "normalized_pop_values",
                       x = x_grid,
                       y = ecdf_pop(x_grid),
                       n_obs = n_obs
  )
    
  
  # Return combined result
  bind_rows(result_units, result_values, result_pop)
}
plan(multisession, workers = 3)  # or sequential/multicore/etc.
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

fwrite(results_long_normalized,file=paste0(datadir,"\\normalized_ecdfs_loss_5.22.csv"))

ecdf_summary <- results_long_normalized %>%
  group_by(x,comparison) %>%
  summarize(
    mean_y = weighted.mean(x=y,w=n_obs)
  )

newlabs=c("Flood Protection: Housing Units","Flood Protection: Housing Values", "Flood Protection: Population")
names(newlabs)=c("normalized_loss_units","normalized_loss_values", "normalized_pop_values")

a=ggplot(
  results_long_normalized,
  aes(x = x, y = y,group=unit,lwd=n_obs,col=comparison)) + geom_line(alpha=0.2)+
  geom_vline(xintercept = 1, linetype = "dashed") + geom_line(data=ecdf_summary,aes(x=x,y=mean_y),col="black",inherit.aes=FALSE)+
  scale_x_log10() + facet_wrap(~comparison, labeller=labeller(comparison=newlabs))+labs(y="Cumulative Density",x="Value Relative to Bank Mean")+
  theme_bw()+scale_linewidth_continuous(guide="none")+scale_color_discrete(guide="none")+theme(strip.background =element_rect(fill="white"))
a


# ------- Indu code begins here --------

results_long_normalized = read.csv(paste0(datadir,"\\normalized_ecdfs_loss.csv"))

# Highest number of loss cells:
loss_cell_counts <- results_long_normalized %>%
  filter(comparison == "normalized_loss_units") %>%
  group_by(unit) %>%
  slice(1) %>%             # take the first row per unit
  ungroup() %>%
  dplyr::select(unit, n_obs) %>%
  arrange(desc(n_obs))
loss_cell_counts
######## cutoff thresholds  ------------------

cut_pts <- c(`gt1` = 1, `gte1.5` = 1.5, `gte2` = 2, `gt10` = 10)

# helper function: 1 - ECDF(threshold)
tail_share <- function(df, threshold) {
  1 - approx(x = df$x, y = df$y, xout = threshold, rule = 2)$y
}

# compute tail shares per unit
shares_per_unit <- results_long_normalized %>%
  filter(comparison == "normalized_loss_units") %>%
  group_by(unit) %>%
  group_modify(~{
    unit_n <- max(.x$n_obs, na.rm = TRUE)
    shares <- map_dbl(cut_pts, function(thresh) tail_share(.x, thresh))
    tibble(n_obs = unit_n, !!!setNames(shares * 100, names(cut_pts)))
  }) %>%
  ungroup()
# overall weighted mean
overall_wtd <- shares_per_unit %>%
  select(-unit) %>%
  summarise(across(-n_obs, ~ weighted.mean(.x, w = n_obs, na.rm = TRUE)))



# ------------------------


# ------------------------

ecdf_summary <- results_long_normalized %>%
  group_by(x,comparison) %>%
  summarize(
    mean_y = weighted.mean(x=y,w=n_obs)
  )
newlabs=c("Flood Protection: Housing Units","Flood Protection: Housing Values", "Flood Protection: Population")
names(newlabs)=c("normalized_loss_units","normalized_loss_values", "normalized_loss_population")

# 1. Define the reference x-values
ref_values <- c(1, 1.5, 2)

crossings <- ecdf_summary %>%
 
  filter(x > 0) %>%
  group_by(comparison) %>%
  summarize(
    # approx(..., rule=2) will clamp outside x-range rather than return NA
    below_1   = approx(x, mean_y, xout = 1,   rule = 2)$y,
    below_1_5 = approx(x, mean_y, xout = 1.5, rule = 2)$y,
    below_2   = approx(x, mean_y, xout = 2,   rule = 2)$y,
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = starts_with("below_"),
    names_to = "key",
    values_to = "cdf_value"
  ) %>%
  mutate(
    ref_x = case_when(
      key == "below_1"   ~ 1,
      key == "below_1_5" ~ 1.5,
      key == "below_2"   ~ 2
    ),
    percent_below = round(100 * cdf_value, 1),
    percent_above = round(100 * (1 - cdf_value), 1),
    label = paste0(percent_below, "% below\n", percent_above, "% above")
  )

#  Add horizontal lines & labels to existing plot ‘a’
#    anchor the left side at x = 0.1 for readability in log-scale plots:
xmin_val <- min(results_long_normalized$x[results_long_normalized$x > 0], na.rm = TRUE)

# ---- Same plot but only housing units -----

# Filter only housing unit rows
results_hu <- results_long_normalized %>%
  filter(comparison == "normalized_loss_units")

# Recalculate ecdf_summary for housing units only
ecdf_summary_hu <- results_hu %>%
  group_by(x, comparison) %>%
  summarize(mean_y = weighted.mean(y, w = n_obs), .groups = "drop")
big_cypress_line <- results_hu %>% filter(unit == "Big_Cypress_MB_Phase_I-V")

big_cypress_label_point <- big_cypress_line %>%
  filter(y >= 0.4 & y <= 0.6) %>%
  slice_min(x)

# Assign colors to the reference lines
ref_line_colors <- c("gray20", "gray50", "gray80")
names(ref_line_colors) <- c("below_1", "below_1_5", "below_2")

# Merge colors into crossings
crossings <- crossings %>%
  mutate(line_color = ref_line_colors[key])
crossings_hu <- crossings %>% filter(comparison == "normalized_loss_units")

a_hu_improved <- ggplot(
  results_hu,
  aes(x = x, y = y, group = unit, lwd = n_obs)
) +
  geom_line(alpha = 0.2, aes(color = "Big Cypress Mitigation Bank")) +  # label the raw lines
  geom_line(data = ecdf_summary_hu, aes(x = x, y = mean_y, color = "Mean ECDF"), inherit.aes = FALSE, linewidth = 1) +
  
  # Vertical dashed lines from xmin to crossing point only
  geom_segment(
    data = crossings_hu,
    aes(x = ref_x, xend = ref_x, y = 0, yend = cdf_value, color = key),
    linetype = "dashed",
    linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  
  # Horizontal dashed lines from xmin to vertical line, color-coded
  geom_segment(
    data = crossings_hu,
    aes(x = xmin_val, xend = ref_x, y = cdf_value, yend = cdf_value, color = key),
    linetype = "dashed",
    linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  annotate("text",
           x = big_cypress_label_point$x + 50,
           y = big_cypress_label_point$y + 0.01,
           label = "Big Cypress\nMitigation Bank",
           color = "black",
           fontface = "bold",
           size = 4,
           hjust = 0.3) +
  annotate("segment",
           x = big_cypress_label_point$x + 26,  # moved from +20 to +23
           xend = big_cypress_label_point$x + 6,  # was just big_cypress_label_point$x
           y = big_cypress_label_point$y + 0.01,
           yend = big_cypress_label_point$y,
           color = "black",
           arrow = arrow(length = unit(0.25, "cm")))+
  
  # Color mapping for reference lines
  scale_color_manual(
    values = c(
      "below_1" = "gray20",
      "below_1_5" = "gray50",
      "below_2" = "gray80",
      "Mean ECDF" = "black",
      "Big Cypress Mitigation Bank" = "#F8766D"
    ),
    guide = "none"  # remove legend
  )+
  
  labs(
    y = "Cumulative Density",
    x = "Value Relative to Bank Mean"
  ) +
  scale_x_log10(
    breaks = c(0.1, 1, 1.5, 2, 10, 100, 1000),
    labels = c("0.1", "1", "1.5", "2", "10", "100", "1000"),
    expand = c(0, 0)
  )+
  scale_y_continuous(expand = c(0, 0))+
  theme_bw(base_size = 14) +  # make all text larger
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_blank(),  # remove facet title
    strip.background = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_linewidth_continuous(guide = "none")

a_hu_improved 

# # Identify the thick line that has 10x higher value in SA than in bank ----
# 
# # Get both x values at y=0.1 and y=0.9, then compute the steepness range
# find_steep_midrange <- function(df) {
#   if (n_distinct(df$y) < 2 || all(is.na(df$y))) return(tibble(x10 = NA, x90 = NA, range = NA))
#   x10 <- approx(df$y, df$x, xout = 0.1, rule = 2)$y
#   x90 <- approx(df$y, df$x, xout = 0.9, rule = 2)$y
#   tibble(x10 = x10, x90 = x90, range = x90 - x10)
# }
# 
# candidates <- results_long_normalized %>%
#   filter(comparison == "normalized_loss_units") %>%
#   group_by(unit) %>%
#   group_modify(~ {
#     x_vals <- find_steep_midrange(.x)
#     x_vals$n_obs <- unique(.x$n_obs)
#     x_vals
#   }) %>%
#   filter(!is.na(x10), x10 > 5, x90 < 40) %>%
#   arrange(desc(n_obs), range)
# 
# print(candidates)

highlight_unit <- "Big_Cypress_MB_Phase_I-V"

df_highlight <- results_long_normalized %>%
  filter(unit == highlight_unit, comparison == "normalized_loss_units")

a_hu_improved +
  geom_line(data = df_highlight, aes(x = x, y = y), color = "black", size = 1.5) +
  ggtitle("ECDF Plot with Big_Cypress_MB_Phase_I-V Highlighted")

# Check median??

median(df_highlight$x)





# ---------- Same plot but HV and POP ------
results_hv_pop <- results_long_normalized %>%
  filter(comparison %in% c("normalized_loss_values", "normalized_pop_values"))

# Recalculate ecdf_summary for hv and pop 
ecdf_summary_hv_pop <- results_hv_pop %>%
  group_by(x, comparison) %>%
  summarize(mean_y = weighted.mean(y, w = n_obs), .groups = "drop")

a_hv_pop_improved <- ggplot(
  results_hv_pop, 
  aes(x = x, y = y, group = unit, lwd = n_obs)
) + 
  geom_line(aes(color = comparison), alpha = 0.2) + 
  geom_line(
    data = ecdf_summary_hv_pop, 
    aes(x = x, y = mean_y, color = "Mean ECDF"), 
    inherit.aes = FALSE,
    linewidth = 1
  ) + 
  
  # Vertical dashed lines at fixed values
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray20", linewidth = 0.5) + 
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "gray50", linewidth = 0.5) + 
  geom_vline(xintercept = 2, linetype = "dashed", color = "gray80", linewidth = 0.5) + 
  
  # Horizontal dashed segments
  geom_segment(
    data = crossings %>% filter(comparison %in% c("normalized_loss_values", "normalized_pop_values")),
    aes(x = xmin_val, xend = ref_x, y = cdf_value, yend = cdf_value, color = comparison),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  facet_wrap(~comparison, labeller = labeller(comparison = newlabs)) + 
  
  labs(y = "Cumulative Density", x = "Value Relative to Bank Mean") + 
  
  scale_x_log10(
    breaks = c(0.1, 1, 1.5, 2, 10, 100, 1000),
    labels = c("0.1", "1", "1.5", "2", "10", "100", "1000"),
    expand = c(0, 0)
  ) + 
  scale_y_continuous(expand = c(0, 0)) +
  
  theme_bw(base_size = 14) + 
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14),
    strip.background = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  
  scale_color_manual(
    values = c(
      "normalized_pop_values" = "#00BA38",     # green
      "normalized_loss_values" = "#619CFF",    # blue
      "Mean ECDF" = "black"
    ),
    guide = "none"
  ) +
  
  scale_linewidth_continuous(guide = "none")

a_hv_pop_improved

a_hu_improved



## AREA WEIGHTED DIFFERENCE AND MAX VALUES ##

# ----------  Build one light-weight summary row per mitigation bank ----------
library(dplyr)
library(glue)
library(purrr)
library(tibble)

ratio_df <- map_dfr(names, function(nm) {
  bank_pth <- file.path(datadir, "Bank ECDF Functions",      paste0("bank_ecdf_fn_", nm, ".rds"))
  loss_pth <- file.path(datadir, "SA ECDF Functions - Loss", paste0("ecdf_fn_",    nm, ".rds"))
  if (!file.exists(bank_pth) || !file.exists(loss_pth)) return(NULL)
  
  bank_vals <- environment(readRDS(bank_pth)$housing_units)$x
  loss_vals <- environment(readRDS(loss_pth)$housing_units)$x
  
  tibble(
    unit        = nm,
    mean_bank   = mean(bank_vals,  na.rm = TRUE),
    mean_loss   = mean(loss_vals,  na.rm = TRUE),
    n_pixels    = length(loss_vals)
  )
}) %>%
  ## drop pathological or empty cases up-front
  filter(is.finite(mean_bank), mean_bank > 0,
         is.finite(mean_loss),  mean_loss > 0) %>%
  mutate(ratio_mean = mean_loss / mean_bank)     # safe now

# ---------- Aggregate across banks with area (loss) weighting ----------------------
area_weighted_avg <- with(ratio_df, weighted.mean(ratio_mean, w = n_pixels, na.rm = TRUE))
max_ratio         <- max(ratio_df$ratio_mean,  na.rm = TRUE)          # most extreme bank-level case
p95_ratio         <- quantile(ratio_df$ratio_mean, 0.95, na.rm = TRUE)  # “typical high end” 

cat(glue(
  "Wetlands lost to development tend to be near previously developed areas and ",
  "therefore deliver higher downstream flood-protection value than the wetlands ",
  "created to compensate for them – on average about ",
  "{round(area_weighted_avg, 1)}× as much when weighted by the area of loss, ",
  "with some banks showing ratios as high as {round(max_ratio, 1)}×."
))


