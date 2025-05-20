library(fixest)
library(tidyverse)
library(data.table)
library(broom)

dat=readRDS("data/nfipbankdata_unlist.rds")

policies=readRDS("data/policies.rds")

dat_policies_merged=merge(dat,policies%>%select(tract,policyeffective_year,total_coverage_2025),by.x=c("tract","claim_year"),by.y=c("tract","policyeffective_year"),all.x=TRUE,all.y=FALSE)

# Function to fit hurdle model for one group
fit_hurdle_model <- function(df_group,policy_control=FALSE) {
  # Complete all tract-year combinations
  df_complete <- df_group %>%
    complete(tract, claim_year,
             fill = list(
               total_paid_2025 = 0,
               n_claims = 0
             )) %>%
    select(tract,claim_year,total_paid_2025)%>%
    mutate(has_claims = total_paid_2025 > 0)%>%
    mutate(tract=as.factor(tract))%>%
    filter(claim_year%in%(1985:2021))
  
  if(policy_control==TRUE){
    #merge in relevant data from policies
    df_complete=merge(df_complete,policies%>%filter(tract%in%unique(df_complete$tract))%>%
                      select(tract,policyeffective_year,total_coverage_2025),
                      by.x=c("tract","claim_year"),by.y=c("tract","policyeffective_year"),all.x=TRUE,all.y=FALSE)
  }
  
  if(policy_control==FALSE){
    # First stage: binary model
    model_zero <- tryCatch({
      feglm(has_claims ~ claim_year | tract, family = binomial(), data = df_complete)
    }, error = function(e) NULL)
    
    # Second stage: log-linear model for positive values
    model_positive <- tryCatch({
      feols(log(total_paid_2025) ~ claim_year | tract, data = df_complete %>% filter(total_paid_2025 > 0))
    }, error = function(e) NULL)
  }
  # Tidy model outputs if both succeeded
  if (!is.null(model_zero) && !is.null(model_positive)) {
    zero_tidy <- tidy(model_zero) %>%
      mutate(stage = "zero")
    
    pos_tidy <- tidy(model_positive) %>%
      mutate(stage = "positive")
    
    bind_rows(zero_tidy, pos_tidy)
  } else {
    tibble(term = NA, estimate = NA, std.error = NA, p.value = NA, stage = NA)
  }
  
  if(policy_control==TRUE){
    # First stage: binary model
    model_zero <- tryCatch({
      feglm(has_claims ~ claim_year+log(total_coverage_2025) | tract, family = binomial(), data = df_complete)
    }, error = function(e) NULL)
    
    # Second stage: log-linear model for positive values
    model_positive <- tryCatch({
      feols(log(total_paid_2025) ~ claim_year+log(total_coverage_2025) | tract, data = df_complete %>% filter(total_paid_2025 > 0))
    }, error = function(e) NULL)
    
    # Tidy model outputs if both succeeded
    if (!is.null(model_zero) && !is.null(model_positive)) {
      zero_tidy <- tidy(model_zero) %>%
        mutate(stage = "zero")%>%
        filter(term=="claim_year")
      
      pos_tidy <- tidy(model_positive) %>%
        mutate(stage = "positive")%>%
        filter(term=="claim_year")
      
      bind_rows(zero_tidy, pos_tidy)
    } else {
      tibble(term = NA, estimate = NA, std.error = NA, p.value = NA, stage = NA)
    }

  }
}

# Apply to each banks_upstream group
model_results <- dat %>%
  # filter(banks_upstream%in%c("Sharp_Mountain_Creek_-_GWTF",
  #                            "Prater_Island_Single_Client_Bank",
  #                            "Conasauga_Bend",
  #                            "Conasauga_River",
  #                            "Old_Creek_Place_Mitigation_Bank",
  #                            "Upper_Coosa_Mitigation_Bank"  ))%>%
  group_by(banks_upstream) %>%
  group_split() %>%
  map_dfr(
    ~ fit_hurdle_model(.x,policy_control=FALSE) %>%
      mutate(banks_upstream = unique(.x$banks_upstream)),
    .id = "group_id"
  ) %>%
  select(banks_upstream, stage, term, estimate, std.error, p.value)

write.csv(model_results,file="data/hurdle_model_1985_2021_censusFEs.csv")

a=ggplot(model_results%>%filter(!is.na(estimate)),aes(x=estimate))+geom_histogram()+facet_wrap(~stage)

#version controling for policy coverage - limited to post 2008 time trends

model_results <- dat %>%
  group_by(banks_upstream) %>%
  group_split() %>%
  map_dfr(
    ~ fit_hurdle_model(.x,policy_control=TRUE) %>%
      mutate(banks_upstream = unique(.x$banks_upstream)),
    .id = "group_id"
  ) %>%
  select(banks_upstream, stage, term, estimate, std.error, p.value)

write.csv(model_results,file="data/hurdle_model_2009_2021_censusFEs_policycontrol.csv")

a=ggplot(model_results%>%filter(!is.na(estimate)),aes(x=estimate))+geom_histogram()+facet_wrap(~stage)
