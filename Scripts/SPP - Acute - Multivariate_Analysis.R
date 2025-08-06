###The Swimmer's Phenomics Project - Study 1###

#Data analysis for the Swimmer's Phenomics Project - Study 1 - Acute Post-Exercise Tiem Course Study
#Purpose: Data clean to create DF of untargeted data

#Analysis plan:

#1. Calculate CV for each metabolite for resting values
#2. ASCA: Time (Pre, Post, 8 h post, 22 h post) x Condition (Control, S1, S2)
#3. Consider sex effects later
#4. Consider squad effects

####PART X - Load Libraries####
library(tidyverse)
library(tidyplots)
library(ALASCA)
library(readr)

####PART X - Import Data#####

#Import - Imputed data frame
wide_data_imp <- readr::read_csv("Data/SPP_S1 - WideData_RF_imputed.csv")

####PART X - Clean Data####

#Rename week to session and filter out PQCs
wide_data_imp_1 <- wide_data_imp %>%
  dplyr::rename(session=week)%>%
  dplyr::filter(sample=="sample")%>%
  dplyr::select(-sample)

#Pivot wide to long format
long_data_imp_1 <- wide_data_imp_1 %>% 
  tidyr::pivot_longer(!subject_id:sex, names_to = "variable", values_to = "value")

#Create factors - Sample time & session
long_data_imp_1$sample_time <- factor(long_data_imp_1$sample_time,
                                      levels=c(0, 1, 2, 3),
                                      labels=c("Pre", "Post", "8 h post", "22 h post"))

long_data_imp_1$session <- factor(long_data_imp_1$session,
                                      levels=c(0, 1, 2),
                                      labels=c("Control", "Session 1", "Session 2"))

#####PART X - ALASCA Model - Time*Session - All, no filtering#####

#Covert DF to DT
long_data_imp_asca <- setDT(long_data_imp_1)

#ASCA Model
alasca_time_session <- ALASCA(formula = value ~ sample_time*session + (1|subject_id),
                          df = long_data_imp_asca,
                          participant_column="subject_id",
                          scale_function = "sdt1",
                          reduce_dimensions.limit=TRUE,
                          n_validation_runs=5,
                          max_PC=5,
                          validate = TRUE,
                          use_Rfast= TRUE,
                          save=TRUE,
                          filename = "SPP_S1_sampleTime_session")

#Plot - Effects - Model Scores and Loadings
plot(alasca_time_session, 
     effect=1, 
     component=c(1,2), 
     type="effect")

#Plot - Marginal means from regression
plot(alasca_time_session, 
     component = 1,
     effect = 1,
     type = 'prediction',
     n_limit = 20,
     variable = c())

#Import ASCA Results  - Loadings
alasca_time_session_loadings <- 
  read_csv("ALASCA/SampleTime_Session/loadings_effect_1.csv")

# Rename "covars" column
alasca_time_session_loadings <- alasca_time_session_loadings %>%
  dplyr::rename(metabolite=covars)

#Import ASCA Results - Predictions
alasca_time_session_pred <- 
  read_csv("ALASCA/SampleTime_Session/model_prediction.csv")

# Create factors - Sample time & session
alasca_time_session_pred$sample_time <- factor(alasca_time_session_pred$sample_time,
                                                 levels=c("Pre", "Post", "8 h post", "22 h post"),
                                                 labels=c("Pre", "Post", "8 h post", "22 h post"))

alasca_time_session_pred$session <- factor(alasca_time_session_pred$session,
                                             levels=c("Control", "Session 1", "Session 2"),
                                             labels=c("Control", "Session 1", "Session 2"))

alasca_time_session_pred_1 <- alasca_time_session_pred %>%
  dplyr::filter(str_detect(variable, "_acid"))

# Plot - Individual metabolites
dodge_width <- 0.3 #position dodge

plot_pred <- alasca_time_session_pred %>%
  ggplot(aes(x = sample_time, y = pred, color = session, group = session, shape = session)) +
  geom_point(size = 3, position = position_dodge(width = dodge_width)) +
  geom_line(aes(linetype = session), linewidth = 1, position = position_dodge(width = dodge_width)) +
  geom_errorbar(aes(ymin = low, ymax = high),
                width = 0.1, linewidth = 0.8,
                position = position_dodge(width = dodge_width)) +
  scale_color_manual(values = c("Control" = "#E69F00", "Session 1" = "#56B4E9", "Session 2" = "black")) + # Colorblind-safe
  scale_shape_manual(values = c(16, 17, 18)) +
  scale_linetype_manual(values = c("Control" = "solid", "Session 1" = "dashed", "Session 2" = "dotdash")) +
  ylab("Std. Value") +
  xlab("") +
  facet_wrap(~variable, scales = "free_y") +
  theme_minimal(base_size = 14) +
  theme(
    strip.background = element_rect(color = "black", fill = "white", linewidth = 1),
    strip.text = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)),
         shape = guide_legend(override.aes = list(size = 4)))


# Import ASCA Results - Scores
alasca_time_session_scores <- 
  read_csv("ALASCA/SampleTime_Session/scores_effect_1.csv")

# Create factors - Sample time & session
alasca_time_session_scores$sample_time <- factor(alasca_time_session_scores$sample_time,
                                      levels=c("Pre", "Post", "8 h post", "22 h post"),
                                      labels=c("Pre", "Post", "8 h post", "22 h post"))

alasca_time_session_scores$session <- factor(alasca_time_session_scores$session,
                                  levels=c("Control", "Session 1", "Session 2"),
                                  labels=c("Control", "Session 1", "Session 2"))

## PC1
alasca_time_session_scores <- alasca_time_session_scores %>% 
  filter(PC == "1")

# Plot - Scores - tidyplots
alasca_time_session_scores |>
  tidyplot(x = sample_time, y = score, color = session) |>
  add_mean_line() |>
  add_mean_dot() |>
  add_sem_ribbon() |>
  adjust_legend_title("Session") |>
  adjust_x_axis_title("Sample Time") |>
  adjust_y_axis_title("Score")

#Plot - Effects - Model Scores and Loadings
plot(alasca_time_session, 
     effect=1, 
     component=1, 
     type="effect")

#Plot - Model regression residuals
plot(alasca_time_session, 
     effect=1, component=c(1,2), 
     type="residuals", 
     variable="sample_time")

#Plot - Marginal means from regression
plot(alasca_time_session, 
     component = 1,
     effect = 1,
     type = 'prediction',
     n_limit = 20,
     variable = c())

####PART X - Named Metabolites Only####

#Import - Imputed data frame
named_metabs_long <- readr::read_csv("Data/SPP_S1 - NamedMetabs_long.csv")

####PART X - Clean Data####

#Rename week to session and filter out PQCs
named_metabs_long_1 <- named_metabs_long %>%
  dplyr::filter(sample=="sample") %>%
  dplyr::relocate(metabolite, .before = value) %>%
  dplyr::select(-c(unique_id:sample)) %>%
  dplyr::select(-c(idx:idx_mode)) %>%
  dplyr::rename(session=week)%>%
  dplyr::rename(variable=metabolite)

#Create factors - Sample time & session
named_metabs_long_1$sample_time <- factor(named_metabs_long_1$sample_time,
                                      levels=c(0, 1, 2, 3),
                                      labels=c("Pre", "Post", "8 h post", "22 h post"))

named_metabs_long_1$session <- factor(named_metabs_long_1$session,
                                  levels=c(0, 1, 2),
                                  labels=c("Control", "Session 1", "Session 2"))

#####PART X - ALASCA Model - Time*Session - All, no filtering#####
long_named_data_asca <- data.table::setDT(named_metabs_long_1)

#ASCA Model
alasca_named_time_session <- ALASCA(formula = value ~ sample_time*session + (1|subject_id),
                              df = long_named_data_asca,
                              participant_column="subject_id",
                              scale_function = "sdt1",
                              reduce_dimensions.limit=TRUE,
                              n_validation_runs=5,
                              max_PC=5,
                              validate = TRUE,
                              use_Rfast= TRUE,
                              save=TRUE,
                              filename = "SPP_S1_sampleTime_session_Named")

#Plot - Effects - Model Scores and Loadings
plot(alasca_named_time_session, 
     effect=1, 
     component=c(1,2), 
     type="effect")

#Plot - Model regression residuals
plot(alasca_named_time_session, 
     effect=1, component=c(1,2), 
     type="residuals", 
     variable="sample_time")

#Plot - Marginal means from regression - PC1
plot(alasca_named_time_session, 
     component = 1,
     effect = 1,
     type = 'prediction',
     n_limit = 20,
     variable = c())

#Plot - Marginal means from regression - PC2
plot(alasca_named_time_session, 
     component = 2,
     effect = 1,
     type = 'prediction',
     n_limit = 20,
     variable = c())
