library(tidyverse)
library(readr)

wellness_df <- 
  read_csv("C:/Users/lheidenreich/OneDrive - LA TROBE UNIVERSITY/Study 1_Acute/Meta Data (REDCap)/SPPS1DailyWellnessMo_DATA_2024-08-23_1317.csv")

test_dates <- c("27/07/2024", "28/07/2024", "03/08/2024", "04/08/2024",
                "05/08/2024", "10/08/2024", "11/08/2024")

test_dates1 <- c("27/07/2024", "03/08/2024", "04/08/2024", "10/08/2024")

wellness_df %>% 
  group_by(initials) %>% 
  summarise(n = n())

wellness_df_test <- wellness_df %>% 
  filter(well_date %in% test_dates1)

wellness_df_test %>% 
  group_by(well_date) %>% 
  summarise(n = n(), 
            initials = paste(unique(initials), collapse = ", "),
            .groups = "drop")

wellness_df_test %>% 
  group_by(well_date) %>% 
  summarise(n = n(), 
            initials = paste(unique(initials), collapse = ", "),
            .groups = "drop")
