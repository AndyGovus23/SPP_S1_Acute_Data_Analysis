###The Swimmer's Phenomics Project - Study 1###

#Data analysis for the Swimmer's Phenomics Project - Study 1 - Acute Post-Exercise Tiem Course Study
#Purpose: Data clean to create DF of untargeted data

#Analysis plan:

#1. Remove features with >30% QC RSD
#2. Remove features with excessive zeros
#3. Remove features that are flagged as include == FALSE
#4. Perform PCA to check QC vs. samples
#5. Perform outlier detection and RF imputation once untargeted data is combined with lipid data

#Data Frames:

#Cleaned data frame with features removed: long_data_cleaned

#Data Analysis:

#PCA of QCs vs. Samples in XXXX script.

####PART X - Load Libraries####

#Data Tables
library(DT)

#Missing Value Imputation & Outlier detection
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

#Multivariate statistics
library(mixOmics)

#Import Data
rp_pos_data <- readr::read_csv("Data/SPP_S1 - rp_pos_data.csv",) #Reverse phase - positive mode
rp_neg_data <- readr::read_csv("Data/SPP_S1 - rp_neg_data.csv") #Reverse phase - negative mode
hil_pos_data <- readr::read_csv("Data/SPP_S1 - hil_pos_data.csv") #HELIC - positive mode
hil_neg_data <- readr::read_csv("Data/SPP_S1 - hil_neg_data.csv") #HELIC - negative mode
long_data_raw <- readr::read_csv("Data/SPP_S1 - long_data_raw.csv") #Combined dataset containing metabolites detected in rp_pos, rp_neg and hil_pos modes

####PART X - Handling of QC_RSD values####

#Create cleaner dataset
# long_data_clean <- long_data_raw%>%
#   dplyr::select(unique_id,
#                 idx,
#                 subject_id,
#                 run_no,
#                 run_order,
#                 squad,
#                 age,
#                 sex, 
#                 sample, 
#                 week,
#                 sample_time,
#                 mode,
#                 metabolite,
#                 value,
#                 qc_rsd_percent)

#Create idx_mode column that combines idx and mode for better filtering
long_data_clean <- long_data_raw %>%
  dplyr::mutate(idx_mode = paste(idx, mode, sep = "_"))

#Add pcq or blank to subject_id column
long_data_clean <- long_data_clean %>%
  dplyr::mutate(subject_id = case_when(
    sample %in% c("pqc", "blank") ~ paste(sample, run_order, mode,sep = "_"),  # Append idx_mode
    TRUE ~ subject_id  # Keep original subject_id otherwise
  ))

# long_data_clean <- long_data_clean %>%
#   dplyr::mutate(idx_mode = paste(idx, mode, sep = "_"))

#Check - Distinct features
dplyr::n_distinct(long_data_clean$idx_mode)
#10885 features

#Calculate descriptive statistics and display QC_RSD%
long_data_clean_desc <- long_data_clean %>%
  dplyr::group_by(idx_mode)%>%
  dplyr::summarise(
    na = sum(is.na(value)),                          # Count of NAs
    zero = sum(value == 0, na.rm = TRUE),            # Count of zeros
    mean = mean(value, na.rm = TRUE),                # Mean
    sd = sd(value, na.rm = TRUE),                    # Standard deviation
    n = sum(!is.na(value)),                          # Count of non-NA values
    zero_percent = (zero/n)*100,
    na_percent = (na/n)*100,
    p25 = quantile(value, 0.25, na.rm = TRUE),       # 25th percentile
    p50 = quantile(value, 0.50, na.rm = TRUE),       # Median (50th percentile)
    p75 = quantile(value, 0.75, na.rm = TRUE),       # 75th percentile
    p100 = max(value, na.rm = TRUE),                 # Maximum (100th percentile)
    upper_95_ci = mean + 1.96 * (sd / sqrt(n)),      # Upper 95% confidence interval
    lower_95_ci = mean - 1.96 * (sd / sqrt(n)),      # Lower 95% confidence interval
    qc_rsd_percent = mean(qc_rsd_percent)
  )

#Check - Distinct features
dplyr::n_distinct(long_data_clean_desc$idx_mode)
#10885 unique features - match

#Filter - Metabolites with QC_RSD >30%
long_data_clean_desc_qc_rsd30 <- long_data_clean_desc%>%
  dplyr::filter(qc_rsd_percent>30)

#Plot - Metabolites with QC_RSD >30% - rp_pos
long_data_clean_desc_qc_rsd30 %>%
  dplyr::filter(str_detect(idx_mode, "rp_pos") & qc_rsd_percent > 30) %>%
  ggplot(aes(y = qc_rsd_percent, x = reorder(idx_mode, qc_rsd_percent))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  xlab("") +
  scale_y_continuous("QC RSD > 30%") +
  geom_hline(yintercept = 30, color = "red")

#Plot - Metabolites with QC_RSD >30% - rp_neg
long_data_clean_desc_qc_rsd30 %>%
  dplyr::filter(mode=="rp_neg" & qc_rsd_percent>30)%>%
  ggplot(aes(y=qc_rsd_percent, x=reorder(idx, qc_rsd_percent)))+
  geom_bar(stat="identity")+
  coord_flip()+
  xlab("")+
  scale_y_continuous("QC RSD > 30%")+
  geom_hline(yintercept = 30, color = "red")

#Plot - Metabolites with QC_RSD >30% - hil_pos
long_data_clean_desc_qc_rsd30 %>%
  dplyr::filter(mode=="hil_pos" & qc_rsd_percent>30)%>%
  ggplot(aes(y=qc_rsd_percent, x=reorder(idx, qc_rsd_percent)))+
  geom_bar(stat="identity")+
  coord_flip()+
  xlab("")+
  scale_y_continuous("QC RSD > 30%")+
  geom_hline(yintercept = 30, color = "red")

#Get names of metabolites with RSD >30%
long_data_clean_exclude <- long_data_clean_desc%>%
  dplyr::filter(qc_rsd_percent>30)

#Check unique features to be excluded
dplyr::n_distinct(long_data_clean_exclude$idx_mode)
#n = 5206 features to be removed
#n = 10885 unique features in long_data_clean
#Remaining = 5679

#Remove features with QC_RSD >30%
long_data_clean_2 <- long_data_clean%>%
  dplyr::filter(!idx_mode %in% long_data_clean_exclude$idx_mode)

#Check remaining unique features
dplyr::n_distinct(long_data_clean_2$idx_mode)
#5679 = matches

####PART X - Handling of Zero Values####

#Omit features where zero values >10%

#Calculate descriptive statistics
long_data_clean_desc_2 <- long_data_clean_2 %>%
  dplyr::group_by(idx_mode)%>%
  dplyr::summarise(
    na = sum(is.na(value)),                          # Count of NAs
    zero = sum(value == 0, na.rm = TRUE),            # Count of zeros
    mean = mean(value, na.rm = TRUE),                # Mean
    sd = sd(value, na.rm = TRUE),                    # Standard deviation
    n = sum(!is.na(value)),                          # Count of non-NA values
    zero_percent = (zero/n)*100,
    p25 = quantile(value, 0.25, na.rm = TRUE),       # 25th percentile
    p50 = quantile(value, 0.50, na.rm = TRUE),       # Median (50th percentile)
    p75 = quantile(value, 0.75, na.rm = TRUE),       # 75th percentile
    p100 = max(value, na.rm = TRUE),                 # Maximum (100th percentile)
    upper_95_ci = mean + 1.96 * (sd / sqrt(n)),      # Upper 95% confidence interval
    lower_95_ci = mean - 1.96 * (sd / sqrt(n)),      # Lower 95% confidence interval
    qc_rsd_percent = mean(qc_rsd_percent)
  )

#Filter - Metabolites with zeros
long_data_clean_desc_zero <- long_data_clean_desc_2%>%
  dplyr::filter(zero_percent>10)

#Plot - Zeros - rp_neg
long_data_clean_desc_zero %>%
  dplyr::filter(str_detect(idx_mode, "rp_neg") & zero_percent>10)%>%
  ggplot(aes(y=zero_percent, x=reorder(idx_mode, zero_percent)))+
  geom_bar(stat="identity")+
  coord_flip()+
  xlab("")+
  ylab("% zeros")+
  geom_hline(yintercept = 10, color = "red")

#Plot - Zeros - rp_pos
long_data_clean_desc_zero %>%
  dplyr::filter(str_detect(idx_mode, "rp_pos") & zero_percent>10)%>%
  ggplot(aes(y=zero_percent, x=reorder(idx_mode, zero_percent)))+
  geom_bar(stat="identity")+
  coord_flip()+
  xlab("")+
  ylab("% zeros")+
  geom_hline(yintercept = 10, color = "red")

#Plot - Zeros - hil_pos
long_data_clean_desc_zero %>%
  dplyr::filter(str_detect(idx_mode, "hil_pos") & zero_percent>10)%>%
  ggplot(aes(y=zero_percent, x=reorder(idx_mode, zero_percent)))+
  geom_bar(stat="identity")+
  coord_flip()+
  xlab("")+
  ylab("% zeros")+
  geom_hline(yintercept = 10, color = "red")

#Plot - Zeros - hil_neg
long_data_clean_desc_zero %>%
  dplyr::filter(str_detect(idx_mode, "hil_neg") & zero_percent>10)%>%
  ggplot(aes(y=zero_percent, x=reorder(idx_mode, zero_percent)))+
  geom_bar(stat="identity")+
  coord_flip()+
  xlab("")+
  ylab("% zeros")+
  geom_hline(yintercept = 10, color = "red")

#Exclude features with >10% zeros
long_data_clean_exclude_zero <- long_data_clean_desc_2%>%
  dplyr::filter(zero_percent>10)
n_distinct(long_data_clean_exclude_zero)
#n = 188 features excludes due to >10% zeros

#Check unique cases in indx_mode column
dplyr::n_distinct(long_data_clean_2$idx_mode)
#n = 5679 features
#Expect after zeros excluded n = 5491

#Filter metabolites with RSD >30% or zeros>5 from long_data_clean
long_data_clean_2 <- long_data_clean_2%>%
  dplyr::filter(!idx_mode %in% long_data_clean_exclude_zero$idx_mode)

#Check - Distinct features
dplyr::n_distinct(long_data_clean_2$idx_mode)
#n = 5491 - matches expected

####PART X - Remove Flagged Features####

#Remove features where include=="FALSE"
long_data_clean_3 <- long_data_clean_2%>%
  dplyr::filter(!include=="FALSE")

#Check - Distinct features
n_distinct(long_data_clean_3$idx_mode)
#n = 3743
#n = 2471 after flagged exclude tagged features are removed

#Remove flag and include column
long_data_clean_3 <- long_data_clean_3 %>%
  dplyr::select(-include, -flags)

#Create unique metabolite name
long_data_clean_3 <- long_data_clean_3 %>%
  dplyr::mutate(
    metabolite = case_when(
      is.na(metabolite) | metabolite == "" ~ idx_mode,  # If metabolite is blank, set to idx_mode
      !str_detect(metabolite, paste0("_", idx_mode)) ~ paste0(metabolite, "_", idx_mode),  # Append idx_mode if not already present
      TRUE ~ metabolite  # Otherwise, keep it unchanged
    )
  )

#Unique features
n_distinct(long_data_clean_3$metabolite)
#n=3743

#Function - Clean metabolite names in columns
clean_column_values <- function(data, col) {
  data %>%
    mutate({{ col }} := {{ col }} %>%
             str_to_lower() %>%
             str_replace_all("[^a-z0-9]+", "_") %>%
             str_replace_all("_+", "_") %>%
             str_replace_all("^_|_$", ""))
}

#Clean - metabolite names
long_data_cleaned <- clean_column_values(long_data_clean_3, metabolite)


#Export .csv file
write.csv(long_data_cleaned, 
          "Data/SPP_S1 - long_data_cleaned.csv",
          row.names = FALSE)

####PART X - Check QC RSD for each metabolite####


#Descriptives
long_data_cleaned_desc <- long_data_cleaned %>%
  dplyr::group_by(metabolite)%>%
  dplyr::summarise(
    qc_rsd_percent = mean(qc_rsd_percent),
    mz = mean(m_z_meas),
    rt = mean(rt_min),
    ccs = mean(ccs_a2),
    na = sum(is.na(value)),                          # Count of NAs
    zero = sum(value == 0, na.rm = TRUE),            # Count of zeros
    mean = mean(value, na.rm = TRUE),                # Mean
    sd = sd(value, na.rm = TRUE),                    # Standard deviation
    n = sum(!is.na(value)),                          # Count of non-NA values
    zero_percent = (zero/n)*100,
    p25 = quantile(value, 0.25, na.rm = TRUE),       # 25th percentile
    p50 = quantile(value, 0.50, na.rm = TRUE),       # Median (50th percentile)
    p75 = quantile(value, 0.75, na.rm = TRUE),       # 75th percentile
    p100 = max(value, na.rm = TRUE),                 # Maximum (100th percentile)
    upper_95_ci = mean + 1.96 * (sd / sqrt(n)),      # Upper 95% confidence interval
    lower_95_ci = mean - 1.96 * (sd / sqrt(n)),      # Lower 95% confidence interval
    )


#Filter data to exclude unnamed metabolities
long_data_cleaned_named <- long_data_cleaned %>%
  dplyr::filter(!str_detect(metabolite, "^\\d{1,4}_(rp_pos|rp_neg|hil_pos|hil_neg)$"))

#Clean data to remove metabolite
long_data_cleaned_named <- long_data_cleaned_named %>%
  dplyr::mutate(metabolite = str_replace(metabolite, "_[^_]+_[^_]+_[^_]+$", ""))

#Descriptive statistics - long_data_cleaned_named
long_data_cleaned_named_desc <- long_data_cleaned_named %>%
  dplyr::group_by(metabolite, mode)%>%
  dplyr::summarise(
    qc_rsd_percent = mean(qc_rsd_percent),
    mz = mean(m_z_meas),
    rt = mean(rt_min),
    ccs = mean(ccs_a2),
    na = sum(is.na(value)),                          # Count of NAs
    zero = sum(value == 0, na.rm = TRUE),            # Count of zeros
    mean = mean(value, na.rm = TRUE),                # Mean
    sd = sd(value, na.rm = TRUE),                    # Standard deviation
    n = sum(!is.na(value)),                          # Count of non-NA values
    zero_percent = (zero/n)*100,
    p25 = quantile(value, 0.25, na.rm = TRUE),       # 25th percentile
    p50 = quantile(value, 0.50, na.rm = TRUE),       # Median (50th percentile)
    p75 = quantile(value, 0.75, na.rm = TRUE),       # 75th percentile
    p100 = max(value, na.rm = TRUE),                 # Maximum (100th percentile)
    upper_95_ci = mean + 1.96 * (sd / sqrt(n)),      # Upper 95% confidence interval
    lower_95_ci = mean - 1.96 * (sd / sqrt(n)),      # Lower 95% confidence interval
  )

# Distinct metabolites
n_distinct(long_data_cleaned_named_desc$metabolite)

# Filter - Metabolites detected in multiple ionisation modes by lowest qc_rsd_percent
df_filtered <- long_data_cleaned_named_desc %>%
  dplyr::group_by(metabolite) %>%
  dplyr::slice_min(qc_rsd_percent, with_ties = FALSE) %>%
  dplyr::ungroup()

# Distinct metabolites
n_distinct(df_filtered) #249 unique named metabolites

# Filter - Long data based on distinct metabolites
filtered_dataset1 <- long_data_cleaned_named %>%
  dplyr::semi_join(df_filtered, by = c("metabolite", "mode"))

# Descriptive statistics - Filtered data
filtered_dataset1_desc <- filtered_dataset1 %>%
  dplyr::group_by(metabolite, mode)%>%
  dplyr::summarise(
    qc_rsd_percent = mean(qc_rsd_percent),
    mz = mean(m_z_meas),
    rt = mean(rt_min),
    ccs = mean(ccs_a2),
    na = sum(is.na(value)),                          # Count of NAs
    zero = sum(value == 0, na.rm = TRUE),            # Count of zeros
    mean = mean(value, na.rm = TRUE),                # Mean
    sd = sd(value, na.rm = TRUE),                    # Standard deviation
    n = sum(!is.na(value)),                          # Count of non-NA values
    zero_percent = (zero/n)*100,
    p25 = quantile(value, 0.25, na.rm = TRUE),       # 25th percentile
    p50 = quantile(value, 0.50, na.rm = TRUE),       # Median (50th percentile)
    p75 = quantile(value, 0.75, na.rm = TRUE),       # 75th percentile
    p100 = max(value, na.rm = TRUE),                 # Maximum (100th percentile)
    upper_95_ci = mean + 1.96 * (sd / sqrt(n)),      # Upper 95% confidence interval
    lower_95_ci = mean - 1.96 * (sd / sqrt(n)),      # Lower 95% confidence interval
  )

#Export .csv of named metabolites for 1 of the ionisation modes (lowest qc_rsd_percent)
write.csv(filtered_dataset1, "Data/SPP_S1 - NamedMetabs_long.csv",
          row.names = FALSE)
