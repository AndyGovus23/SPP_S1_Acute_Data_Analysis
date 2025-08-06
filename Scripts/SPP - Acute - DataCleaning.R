###La Trobe University - The Swimmer's Phenomics Project###

#Authors: Laine Heidenreich (La Trobe University), Nathan Lawler (Murdoch University), Lachlan Mitchell (VIS), Louise Cato (VIS), David Pyne (UC), Andrew Govus (La Trobe University)
#Purpose: To map the time course of metabolite responses after a high-intensity swimming session (7 x 200 m) in highly-trained male and female swimmers

#####DATA CLEANING SCRIPT#####

#To Do:

#1. Reshape data
#2. Extract and clean metabolite names
#3. Determine data quality

####PART 1 - Load Libraries####

#Read Data
library(readxl)

#Data Tables
library(DT)

#Missing Value Imputation & Outlier detection

#Using Random Forests
library(missForest)
library(outForest)

#Data Visualisation
library(tidyverse)
library(ggpubr)

#Data Cleaning#
library(janitor)
library(broom)
library(stringr)
library(forcats)

#Descriptive stats#
library(skimr)

####PART 2 - Import Data####

rp_pos_data <- readxl::read_xlsx("Data/SPP_S1 - UntargetedWorkbook.xlsx", sheet="RPpos") #Reverse phase - positive mode
rp_neg_data <- readxl::read_xlsx("Data/SPP_S1 - UntargetedWorkbook.xlsx", sheet="RPneg") #Reverse phase - negative mode
hil_pos_data <- readxl::read_xlsx("Data/SPP_S1 - UntargetedWorkbook.xlsx", sheet="HILpos") #HELIC - positive mode
hil_neg_data <- readxl::read_xlsx("Data/SPP_S1 - UntargetedWorkbook.xlsx", sheet="HILneg") #HELIC - negative mode
meta_data <- readxl::read_xlsx("Data/SPP_S1 - metaData.xlsx", sheet="Data") #Athlete meta data

####PART 3 - Clean Meta Data####

#Import - metadata#
meta_data <- readxl::read_xlsx("Data/SPP_S1 - metaData.xlsx", sheet="Data")

#Coerce columns to named factors with correct labels for each level

#Create factors
meta_data <- meta_data %>%
  dplyr::mutate(across(c(squad,
                         initial,
                         sex,
                         session), as.factor))

#Create factors with specific levels and labels
meta_data$time <- factor(meta_data$time, 
          levels=c(0,1, 2,3),
          labels=c("Pre", "Post", "8h post", "22h post"))

#Rename inital column for later join
meta_data <- meta_data%>%
  dplyr::rename(subject_id=initial)

#Coerce subject_id column to lowercase
meta_data$subject_id <- tolower(meta_data$subject_id)

#Remove unnecessary columns
meta_data_1 <- meta_data%>%
  dplyr::select(-c(time,
                   session,
                   date,
                   tube_label,
                   source_name))%>%
  dplyr::group_by(subject_id, sex, squad)%>%
  dplyr::summarise_at(vars(age, body_mass), mean)

#Coerce to tibble
meta_data_1 <- as_tibble(meta_data_1)

####PART 4 - Function - Clean Metabolite Data####

#Create column name cleaning function - extract participant ID, run order (1:96), and run (1 or 2)
extract_col_name <- function(headers) {
  library(stringr)
  
  # Extract "HPPp01" followed by digits (the prefix)
  run_prefix <- str_extract(headers, "HPPp01\\d+")
  
  # Extract everything after that prefix and underscore using capturing group
  extracted_suffix <- str_replace(headers, ".*HPPp01\\d+_", "")
  
  # Combine prefix and suffix
  extracted_text <- ifelse(!is.na(run_prefix) & extracted_suffix != headers, 
                           paste(run_prefix, extracted_suffix, sep = "_"), 
                           NA)
  
  return(extracted_text)
}

####PART 5 - Clean Metabolite RP_POS_DATA####

#Step 1: Apply the function to rename selected columns in each data frame to extract Unique HPP
rp_pos_data <- rp_pos_data %>%
  dplyr::rename_with(~extract_col_name(.), .cols = "performance_spp1_PLA_MS-AI-RPPOS@fNMR_HPPr12_HPPp014_PQC_2":"performance_spp1_PLA_MS-AI-RPPOS@fNMR_HPPr12_HPPp015_Blank2_58")

#Step 2: Clean column names
rp_pos_data <- janitor::clean_names(rp_pos_data)

#Step 3: Pivot data and relocate data to desired format - RP_POS
rp_pos_data <- rp_pos_data %>%
  tidyr::pivot_longer(
    cols =  hp_pp014_pqc_2:hp_pp015_blank2_58,  # Only pivot these measurement columns
    names_to = "unique_id",
    values_to = "Value"
  ) %>%
  dplyr::relocate(idx:qc_rsd_percent, .after = Value)  # Move metadata columns after Value

#Step 4: Add idx to unique id column
rp_pos_data <- rp_pos_data %>%
  dplyr::mutate(
    run_no = word(unique_id, 2, sep = "_"),
    run_order = word(unique_id, -1, sep = "_"),  # Extracts the last number
    sample = case_when(
      stringr::str_detect(word(unique_id, 3, sep = "_"), "^hpp") ~ "sample",
      word(unique_id, 3, sep = "_") == "pqc" ~ "pqc",
      word(unique_id, 3, sep = "_") == "blank" ~ "blank",
      TRUE ~ NA_character_
    ),
    subject_id = case_when(
      !sample %in% c("pqc", "blank") & str_detect(word(unique_id, 4, sep = "_"), "[a-zA-Z]") ~ str_sub(word(unique_id, 4, sep = "_"), 1, 2),
      TRUE ~ ""
    ),
    week = case_when(
      !str_detect(unique_id, "pqc|blank") ~ str_sub(word(unique_id, 4, sep = "_"), 3, 3),
      TRUE ~ ""
    ),
    sample_time = case_when(
      !str_detect(unique_id, "pqc|blank") ~ str_sub(word(unique_id, 4, sep = "_"), 4, 4),
      TRUE ~ ""
    )
  )

#Step 6: Add metadata using 'subject_id' as the join key
rp_pos_data <- rp_pos_data %>%
  dplyr::left_join(meta_data_1, by = "subject_id")%>%
  dplyr::rename(metabolite = name)%>%
  dplyr::rename(value = Value)%>%
  dplyr::mutate(mode="rp_pos")

#Step 7: Relocate columns
rp_pos_data <- rp_pos_data%>%
  dplyr::relocate(run_no:body_mass, .after=unique_id)

#Step 8: Check unique cases in metabolite column
dplyr::n_distinct(rp_pos_data$metabolite)
#n = 69 - unique metabolites - annotated already

#Step 9: Export rp_pos_data to .csv file
write.csv(rp_pos_data, 
          "I:/My Drive/Research/Projects - Statistical Analyses/SPP - Acute_Data_Analysis/Data/SPP_S1 - rp_pos_data.csv",
          row.names = FALSE)

####PART 6 - Clean Metabolite RP_NEG_DATA####

#Step 1: Apply the function to rename selected columns in each data frame to extract Unique HPP
rp_neg_data <- rp_neg_data %>%
  dplyr::rename_with(~extract_col_name(.), .cols = "performance_spp1_PLA_MS-AI-RPNEG@fNMR_HPPr12_HPPp014_HPP00515_lk20-1-capillary-03-08-2024_4":"performance_spp1_PLA_MS-AI-RPNEG@fNMR_HPPr12_HPPp015_Blank2_58")

#Step 2: Clean column names
rp_neg_data <- janitor::clean_names(rp_neg_data)

#Step 3: Pivot data and relocate data to desired format - rp_neg
rp_neg_data <- rp_neg_data %>%
  tidyr::pivot_longer(
    cols =  hp_pp014_hpp00515_lk20_1_capillary_03_08_2024_4:hp_pp015_blank2_58,  # Only pivot these measurement columns
    names_to = "unique_id",
    values_to = "Value"
  ) %>%
  dplyr::relocate(idx:qc_rsd_percent, .after = Value)  # Move metadata columns after Value

#Step 4: Add idx to unique id column
rp_neg_data <- rp_neg_data %>%
  dplyr::mutate(
    run_no = word(unique_id, 2, sep = "_"),
    run_order = word(unique_id, -1, sep = "_"),  # Extracts the last number
    sample = case_when(
      stringr::str_detect(word(unique_id, 3, sep = "_"), "^hpp") ~ "sample",
      word(unique_id, 3, sep = "_") == "pqc" ~ "pqc",
      word(unique_id, 3, sep = "_") == "blank" ~ "blank",
      TRUE ~ NA_character_
    ),
    subject_id = case_when(
      !sample %in% c("pqc", "blank") & str_detect(word(unique_id, 4, sep = "_"), "[a-zA-Z]") ~ str_sub(word(unique_id, 4, sep = "_"), 1, 2),
      TRUE ~ ""
    ),
    week = case_when(
      !str_detect(unique_id, "pqc|blank") ~ str_sub(word(unique_id, 4, sep = "_"), 3, 3),
      TRUE ~ ""
    ),
    sample_time = case_when(
      !str_detect(unique_id, "pqc|blank") ~ str_sub(word(unique_id, 4, sep = "_"), 4, 4),
      TRUE ~ ""
    )
  )

#Step 6: Add metadata using 'subject_id' as the join key
rp_neg_data <- rp_neg_data %>%
  dplyr::left_join(meta_data_1, by = "subject_id")%>%
  dplyr::rename(metabolite = name)%>%
  dplyr::rename(value = Value)%>%
  dplyr::mutate(mode="rp_neg")

#Step 7: Relocate columns
rp_neg_data <- rp_neg_data%>%
  dplyr::relocate(run_no:body_mass, .after=unique_id)

#Step 8: Check unique cases in metabolite column
dplyr::n_distinct(rp_neg_data$metabolite)
#n = 132 - unique metabolites - annotated already

#Step 9: Export rp_neg_data to .csv file
write.csv(rp_neg_data, 
          "I:/My Drive/Research/Projects - Statistical Analyses/SPP - Acute_Data_Analysis/Data/SPP_S1 - rp_neg_data.csv",
          row.names = FALSE)

####PART 7 - Clean Metabolite HIL_POS_DATA####

#Step 1: Apply the function to rename selected columns in each data frame to extract Unique HPP
hil_pos_data <- hil_pos_data %>%
  dplyr::rename_with(~extract_col_name(.), .cols = "performance_spp1_PLA_MS-AI-HILPOS@fNMR_HPPr12_HPPp014_PQC_2":"performance_spp1_PLA_MS-AI-HILPOS@fNMR_HPPr12_HPPp015_Blank2_58")

#Step 2: Clean column names
hil_pos_data <- janitor::clean_names(hil_pos_data)

#Step 3: Pivot data and relocate data to desired format - hil_pos
hil_pos_data <- hil_pos_data %>%
  tidyr::pivot_longer(
    cols =  hp_pp014_pqc_2:hp_pp015_blank2_58,  # Only pivot these measurement columns
    names_to = "unique_id",
    values_to = "Value"
  ) %>%
  dplyr::relocate(idx:qc_rsd_percent, .after = Value)  # Move metadata columns after Value

#Step 4: Add idx to unique id column
hil_pos_data <- hil_pos_data %>%
  dplyr::mutate(
    run_no = word(unique_id, 2, sep = "_"),
    run_order = word(unique_id, -1, sep = "_"),  # Extracts the last number
    sample = case_when(
      stringr::str_detect(word(unique_id, 3, sep = "_"), "^hpp") ~ "sample",
      word(unique_id, 3, sep = "_") == "pqc" ~ "pqc",
      word(unique_id, 3, sep = "_") == "blank" ~ "blank",
      TRUE ~ NA_character_
    ),
    subject_id = case_when(
      !sample %in% c("pqc", "blank") & str_detect(word(unique_id, 4, sep = "_"), "[a-zA-Z]") ~ str_sub(word(unique_id, 4, sep = "_"), 1, 2),
      TRUE ~ ""
    ),
    week = case_when(
      !str_detect(unique_id, "pqc|blank") ~ str_sub(word(unique_id, 4, sep = "_"), 3, 3),
      TRUE ~ ""
    ),
    sample_time = case_when(
      !str_detect(unique_id, "pqc|blank") ~ str_sub(word(unique_id, 4, sep = "_"), 4, 4),
      TRUE ~ ""
    )
  )

#Step 6: Add metadata using 'subject_id' as the join key
hil_pos_data <- hil_pos_data %>%
  dplyr::left_join(meta_data_1, by = "subject_id")%>%
  dplyr::rename(metabolite = name)%>%
  dplyr::rename(value = Value)%>%
  dplyr::mutate(mode="hil_pos")

#Step 7: Relocate columns
hil_pos_data <- hil_pos_data%>%
  dplyr::relocate(run_no:body_mass, .after=unique_id)

#Step 8: Check unique cases in metabolite column
dplyr::n_distinct(hil_pos_data$metabolite)
#n = 152 - unique metabolites - annotated already

#Step 9: Export hil_pos_data to .csv file
write.csv(hil_pos_data, 
          "I:/My Drive/Research/Projects - Statistical Analyses/SPP - Acute_Data_Analysis/Data/SPP_S1 - hil_pos_data.csv",
          row.names = FALSE)

####PART 8 - Clean Metabolite HIL_NEG_DATA####

#Step 1: Apply the function to rename selected columns in each data frame to extract Unique HPP
hil_neg_data <- hil_neg_data %>%
  dplyr::rename_with(~extract_col_name(.), .cols = "performance_spp1_PLA_MS-AI-HILNEG@fNMR_HPPr13_HPPp014_PQC_2":"performance_spp1_PLA_MS-AI-HILNEG@fNMR_HPPr13_HPPp015_Blank2_58")

#Step 2: Clean column names
hil_neg_data <- janitor::clean_names(hil_neg_data)

#Step 3: Pivot data and relocate data to desired format - hil_neg
hil_neg_data <- hil_neg_data %>%
  tidyr::pivot_longer(
    cols =  hp_pp014_pqc_2:hp_pp015_blank2_58,  # Only pivot these measurement columns
    names_to = "unique_id",
    values_to = "Value"
  ) %>%
  dplyr::relocate(idx:qc_rsd_percent, .after = Value)  # Move metadata columns after Value

#Step 4: Add idx to unique id column
hil_neg_data <- hil_neg_data %>%
  dplyr::mutate(
    run_no = word(unique_id, 2, sep = "_"),
    run_order = word(unique_id, -1, sep = "_"),  # Extracts the last number
    sample = case_when(
      stringr::str_detect(word(unique_id, 3, sep = "_"), "^hpp") ~ "sample",
      word(unique_id, 3, sep = "_") == "pqc" ~ "pqc",
      word(unique_id, 3, sep = "_") == "blank" ~ "blank",
      TRUE ~ NA_character_
    ),
    subject_id = case_when(
      !sample %in% c("pqc", "blank") & str_detect(word(unique_id, 4, sep = "_"), "[a-zA-Z]") ~ str_sub(word(unique_id, 4, sep = "_"), 1, 2),
      TRUE ~ ""
    ),
    week = case_when(
      !str_detect(unique_id, "pqc|blank") ~ str_sub(word(unique_id, 4, sep = "_"), 3, 3),
      TRUE ~ ""
    ),
    sample_time = case_when(
      !str_detect(unique_id, "pqc|blank") ~ str_sub(word(unique_id, 4, sep = "_"), 4, 4),
      TRUE ~ ""
    )
  )

#Step 6: Add metadata using 'subject_id' as the join key
hil_neg_data <- hil_neg_data %>%
  dplyr::left_join(meta_data_1, by = "subject_id")%>%
  dplyr::rename(metabolite = name)%>%
  dplyr::rename(value = Value)%>%
  dplyr::mutate(mode="hil_neg")

#Step 7: Relocate columns
hil_neg_data <- hil_neg_data%>%
  dplyr::relocate(run_no:body_mass, .after=unique_id)

#Step 8: Check unique cases in metabolite column
dplyr::n_distinct(hil_neg_data$metabolite)
#n = 88 - unique metabolites - annotated already

#Step 9: Export hil_neg_data to .csv file
write.csv(hil_neg_data, 
          "I:/My Drive/Research/Projects - Statistical Analyses/SPP - Acute_Data_Analysis/Data/SPP_S1 - hil_neg_data.csv",
          row.names = FALSE)

###PART 9 - Merge Data####

#Import Cleaned Data
rp_pos_data_clean <- readr::read_csv("Data/SPP_S1 - rp_pos_data.csv",) #Reverse phase - positive mode
rp_neg_data_clean <- readr::read_csv("Data/SPP_S1 - rp_neg_data.csv") #Reverse phase - negative mode
hil_pos_data_clean <- readr::read_csv("Data/SPP_S1 - hil_pos_data.csv") #HELIC - positive mode
hil_neg_data_clean <- readr::read_csv("Data/SPP_S1 - hil_neg_data.csv") #HELIC - negative mode

#Check metabolites present in all datasets
common_metabolites <- Reduce(intersect, list(rp_pos_data_clean$metabolite, 
                                             rp_neg_data_clean$metabolite,
                                             hil_pos_data_clean$metabolite,
                                             hil_neg_data_clean$metabolite))

#Bind three analysis modes together
long_data_raw <- rbind(rp_pos_data_clean,
                       rp_neg_data_clean,
                       hil_pos_data_clean,
                       hil_neg_data_clean)

#Check - Distinct metabolites
dplyr::n_distinct(long_data_raw$metabolite)
#n = 366 unique metabolites

#rp_pos = 69 annotated
#rp_neg = 132 annotated
#hil_pos = 152 annotated
#hil_neg = 88 annotated
#Sum = 441

#366 metabolites - less than the sum of the three datasets due to overlaps
#Many metabolites unannotated that require annotation

###PART 10 - Export long_data_raw to .csv file####
write.csv(long_data_raw, 
          "Data/SPP_S1 - long_data_raw.csv",
          row.names = FALSE)