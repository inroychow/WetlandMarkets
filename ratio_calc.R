library(readr)      # parse_number() is handy if any numeric columns came in as character
library(dplyr) 

#setwd()

ledger = read.csv("C:\\Users\\indumati\\Box\\Paper2\\ledgerdata\\ledger.csv")
data = readRDS("C:\\Users\\indumati\\Box\\Paper2\\Examples\\CONUS example\\full_summary_2021.rds")
names = unique(data$Name)

ledger <- ledger %>%
  mutate(Bank.Name.Clean = gsub("[^A-Za-z0-9]+", "_", Bank.Name)) %>%
  mutate(Bank.Name.Clean = gsub("_+", "_", Bank.Name.Clean),
         Bank.Name.Clean = gsub("_$", "", Bank.Name.Clean))  # remove trailing _ if any

# Then match:
ledger_matched <- ledger %>% 
  filter(Bank.Name.Clean %in% names)
head(ledger_matched)


ledger_ratios <- ledger_matched %>% 
  ## 1.  Make sure all the credit / acre columns are numeric
  mutate(across(
    c(Initiation.Credits, Available.Credits, Released.Credits,
      Withdrawn.Credits, Initiation.Acres),
    ~ parse_number(as.character(.x))
  )) %>% 
  
  ## 2.  One row per bank 
  group_by(Bank.Name.Clean) %>%            
  summarise(
    Initiation_Credits = dplyr::last(Initiation.Credits),
    Initiation_Acres   = dplyr::last(Initiation.Acres),
    Available_Credits  = dplyr::last(Available.Credits),
    Released_Credits   = dplyr::last(Released.Credits),
    Withdrawn_Credits  = dplyr::last(Withdrawn.Credits),
    
    ## 3.  Core ratios -------------------------------------------
    
    ## Credit density of the bank (credits generated per acre restored / preserved)
    Acres_per_credit      = Initiation_Acres/ Initiation_Credits,
    
    ## How many of the originally approved credits have been released?
    Share_Released       = Released_Credits   / Initiation_Credits,
    
    ## How many credits are still on the ledger?
    Share_Available      = Available_Credits  / Initiation_Credits,
    
    ## How many credits have been used (withdrawn) so far?
    Share_Withdrawn      = Withdrawn_Credits  / Initiation_Credits
  ) %>% 
  ungroup()%>% 
  select(Bank.Name.Clean, Acres_per_credit, Share_Released, Share_Available, Share_Withdrawn)

saveRDS(ledger_ratios, "C:\\Users\\indumati\\Box\\Paper2_final\\Bank Footprints\\ledger_ratios.rds")
