###The Swimmer's Phenomics Project - Study 1###

#Data analysis for the Swimmer's Phenomics Project - Study 1 - Acute Time Course Study.
#Purpose: Check QC vs. samples for combined untargeted data set using a PCA

#Load Libraries

#Data manipulation
library(dplyr)
library(tidyr)

#Multivariate statistics
library(mixOmics)
library(ropls)

#Data imputation
library(missForest)

#Import Data
long_data_cleaned <- readr::read_csv("Data/SPP_S1 - long_data_cleaned.csv")

####PART 2 - Create Wide Data Frame####

#Pivot long_data_cleaned to wide format
wide_data_clean_1 <- long_data_cleaned %>%
  group_by(sample, subject_id, sample_time, week, squad, age, sex, metabolite) %>%
  summarise(value = mean(value), na.rm=FALSE, .groups = "drop") %>%  # Aggregate values
  pivot_wider(names_from = metabolite, values_from = value) %>%
  ungroup()%>%
  dplyr::select(-na.rm)

#Clean column names
wide_data_clean_1 <-janitor::clean_names(wide_data_clean_1)

#Count metabolites in wide df to ensure no data loss
metab_no_wide <- wide_data_clean_1%>%
  dplyr::select(-c(sample:sex))%>%
  ncol()
#n = 3743 = long_data_cleaned - No data loss

#Count number of NAs and zeros in each metabolite column
na_zero_counts <- wide_data_clean_1 %>%
  dplyr::group_by(subject_id,sample_time, week)%>%
  dplyr::summarise(
    across(8:3747,  # Proper range selection - should  be column 8:2478
           list(
             na_count = ~sum(is.na(.)), 
             zero_count = ~sum(. == 0, na.rm = TRUE)
           ), 
           .names = "{.col}_{.fn}"
    )
  )

####PART 2 - Impute zero values using Random Forest imputation####

#Create metabolite data matrix
wide_data_clean_metab <-  wide_data_clean_1 %>% 
  dplyr::select(8:3747) %>% 
  dplyr::mutate_all(as.numeric) %>% 
  as.matrix()

#NB: 170*3740 matrix

# Convert zeros to NA
wide_data_clean_metab[wide_data_clean_metab == 0] <- NA

# Random forest non-parametric imputation to replace NAs
# wide_data_clean_imp <- missForest::missForest(wide_data_clean_metab, 
#                                               xtrue = wide_data_clean_metab, 
#                                               verbose = T)
# 
# wide_data_clean_imp <- data.frame(wide_data_clean_imp[[1]], check.names = T)


# Random forest imputation using missRanger
wide_data_clean_metab <- as.data.frame(wide_data_clean_metab)

wide_data_clean_imp <- missRanger::missRanger(
  data = wide_data_clean_metab,
  num.trees = 500,  # Default number of trees
  verbose = TRUE
)

#Reappend meta data
meta_data_wide <- wide_data_clean_1%>%
  dplyr::select(sample:sex)

#Bind meta data with data
wide_data_imp <- cbind(meta_data_wide,
                       wide_data_clean_imp)

#Export - RF imputed dataset
readr::write_csv(wide_data_imp, "Data/SPP_S1 - WideData_RF_imputed.csv")

###PART 2 - PCA of PCQ vs. Samples###

#Analysis plan

#1. Split data into ionisation mode
#2. Check PQC clustering in each ionisation mode separately

#Import - Imputed data frame (can be imported without having to impute again)
wide_data_imp <- readr::read_csv("Data/SPP_S1 - WideData_RF_imputed.csv")

#Create copy of wide_data_imp
wide_data_pca <- wide_data_imp

#Get wide metabolite column names
col_names <- wide_data_pca%>%
  dplyr::select(8:3747)%>%
  colnames()

# #Append "m" in front of metabolites to give metabolite numbers
# wide_data_pca <- wide_data_pca %>%
#   rename_with(~ str_c("m", seq_along(.)+1, "_", .), .cols = all_of(col_names))
# 
# #Get wide metabolite column names
# col_names_m <- colnames(wide_data_pca[8:110])
# 
# #Replace metabolite names with metabolite number
# wide_data_pca <- wide_data_pca %>%
#   rename_with(~ sub("_.*", "", .x), .cols= all_of(col_names_m))

#Remove blanks
wide_data_pca <- wide_data_pca%>%
  dplyr::filter(sample !="blank")

#Coerce relevant meta data to factors
wide_data_pca <- wide_data_pca %>%
  dplyr::mutate(across(.cols=c(sample:sex), .fns=as.factor))

#Filter - metabolites only#
x_metabs <- wide_data_pca %>%
  dplyr::select(8:3747)

#Select samples only (not PQC)
y_labels <- wide_data_pca %>%
  dplyr::select(sample)

#Create ID matrix
m_id <- wide_data_pca %>%
  dplyr::select(subject_id)

#Model 1 - PCA - All metabolites - 5 component model
pca_metabs <- mixOmics::pca(x_metabs,
                            scale=TRUE,
                            center = TRUE,
                            ncomp=5)
#Check - Proportion of variance explained by each component
pca_metabs$prop_expl_var
#PC1: 14%, PC2: 8%

#Check - Cumulative Proportion of variance explained by each component
pca_metabs$cum.var
#2 PCs: ~14% variance explained

#Plot - Scree Plot
plot(pca_metabs)
#2PCs required

#Model 2 - PCA - All metabolites - 2 component model
pca_metabs_2 <- mixOmics::pca(x_metabs,
                              scale=TRUE,
                              center = TRUE,
                              ncomp=2)

#Plot - Individual responses
plotIndiv(pca_metabs_2,
          group=y_labels$sample,
          legend=TRUE,
          ellipse = TRUE,
          title="PCA on Metabolites")