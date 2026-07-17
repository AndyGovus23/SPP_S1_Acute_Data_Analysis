###The Swimmer's Phenomics Project - Study 1###

#Data analysis for the Swimmer's Phenomics Project - Study 1 - Acute Post-Exercise Time Course Study
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
library(data.table)

####PART X - Import Data#####

# #Import - Imputed data frame
# wide_data_imp <- readr::read_csv("Data/SPP_S1 - WideData_RF_imputed.csv")
# 
# ####PART X - Clean Data####
# 
# #Rename week to session and filter out PQCs
# wide_data_imp_1 <- wide_data_imp %>%
#   dplyr::rename(session=week)%>%
#   dplyr::filter(sample=="sample")%>%
#   dplyr::select(-sample)
# 
# #Pivot wide to long format
# long_data_imp_1 <- wide_data_imp_1 %>% 
#   tidyr::pivot_longer(!subject_id:sex, names_to = "variable", values_to = "value")
# 
# #Create factors - Sample time & session
# long_data_imp_1$sample_time <- factor(long_data_imp_1$sample_time,
#                                       levels=c(0, 1, 2, 3),
#                                       labels=c("Pre", "Post", "8 h post", "22 h post"))
# 
# long_data_imp_1$session <- factor(long_data_imp_1$session,
#                                       levels=c(0, 1, 2),
#                                       labels=c("Control", "Session 1", "Session 2"))
# 
# #####PART X - ALASCA Model - Time*Session - All, no filtering#####
# 
# #Covert DF to DT
# long_data_imp_asca <- setDT(long_data_imp_1)

gc()

#Read in all long data
long_data_imp_asca <- readr::read_csv("Data/SPP_S1 - long_data_imp_asca.csv") %>% 
  select(-`...1`)

set.seed(123)

#ASCA Model
alasca_time_session <- ALASCA(formula = value ~ sample_time*session + (1|subject_id),
                              df = long_data_imp_asca,
                              participant_column = "subject_id",
                              scale_function = "sdt1",
                              reduce_dimensions = TRUE,
                              reduce_dimensions.limit = 0.95,
                              pca_function = "irlba", 
                              n_validation_runs = 1000,
                              max_PC = 5,
                              validate = TRUE,
                              use_Rfast = TRUE,
                              save = TRUE,
                              filename = "SPP_S1_sampleTime_session")


# #Plot - Effects - Model Scores and Loadings
# plot(alasca_time_session, 
#      effect=1, 
#      component=c(1,2), 
#      type="effect")
# 
# #Plot - Marginal means from regression
# plot(alasca_time_session, 
#      component = 1,
#      effect = 1,
#      type = 'prediction',
#      n_limit = 20,
#      variable = c())


### LOADINGS ###

# Import ASCA Results - Loadings
alasca_time_session_loadings <-
  read_csv("ALASCA/sampleTime_session/loadings_effect_1.csv")

# Rename "covars" column
alasca_time_session_loadings <- alasca_time_session_loadings %>%
  dplyr::rename(metabolite = covars)

## PC1 ##

# Filter for PC1
alasca_time_session_loadings_pc1 <- alasca_time_session_loadings %>% 
  filter(PC == "1")
# Order data
alasca_time_session_loadings_pc1 <- 
  alasca_time_session_loadings_pc1[order(-alasca_time_session_loadings_pc1$loading), ]
# Top 20 variables (10 positive, 10 negatives)
alasca_time_session_loadings_pc1_top20 <- 
  alasca_time_session_loadings_pc1[c(1:10,3731:3740), ]
# Order from negative to positive
alasca_time_session_loadings_pc1_top20 <- 
  alasca_time_session_loadings_pc1_top20[order(alasca_time_session_loadings_pc1_top20$loading), ]
# Rename metabolites 
alasca_time_session_loadings_pc1_top20 <- alasca_time_session_loadings_pc1_top20 %>% 
  mutate(metabolite = metabolite %>% gsub("_[0-9]+_(hil|rp)_(pos|neg)$", "", .))
# Order of the names 
alasca_time_session_loadings_pc1_order <- 
  alasca_time_session_loadings_pc1_top20$metabolite

# Plot - loadings
alasca_time_session_loadings_pc1_plot <- 
  ggplot(alasca_time_session_loadings_pc1_top20, aes(x = loading, y = factor(metabolite, levels = alasca_time_session_loadings_pc1_order))) + 
  geom_point() + 
  geom_errorbar(aes(xmin = low, xmax = high), width = .1) + 
  scale_x_continuous(breaks = seq(-0.15, 0.27, by = 0.03)) + 
  labs(x = "Loading PC1 (47.6%)", y = "Metabolite") + 
  geom_vline(xintercept = 0.00, linetype = "dotted") + 
  theme(axis.line = element_line(colour = "black"), 
        panel.background = element_blank(), 
        axis.title.x = element_text(margin = ggplot2::margin(t = 10)), 
        axis.title.y = element_text(margin = ggplot2::margin(r = 10)))
  

## PC2 ##

# Filter for PC2
alasca_time_session_loadings_pc2 <- alasca_time_session_loadings %>% 
  filter(PC == "2")
# Order data
alasca_time_session_loadings_pc2 <- 
  alasca_time_session_loadings_pc2[order(-alasca_time_session_loadings_pc2$loading), ]
# Top 20 variables (10 positive, 10 negatives)
alasca_time_session_loadings_pc2_top20 <- 
  alasca_time_session_loadings_pc2[c(1:10,3731:3740), ]
# Order from negative to positive
alasca_time_session_loadings_pc2_top20 <- 
  alasca_time_session_loadings_pc2_top20[order(alasca_time_session_loadings_pc2_top20$loading), ]
# Rename metabolites 
alasca_time_session_loadings_pc2_top20 <- alasca_time_session_loadings_pc2_top20 %>% 
  mutate(metabolite = metabolite %>% gsub("_[0-9]+_(hil|rp)_(pos|neg)$", "", .))
# Order of the names 
alasca_time_session_loadings_pc2_order <- 
  alasca_time_session_loadings_pc2_top20$metabolite

# Plot - loadings
alasca_time_session_loadings_pc2_plot <- 
  ggplot(alasca_time_session_loadings_pc2_top20, aes(x = loading, y = factor(metabolite, levels = alasca_time_session_loadings_pc2_order))) + 
  geom_point() + 
  geom_errorbar(aes(xmin = low, xmax = high), width = .1) + 
  scale_x_continuous(breaks = seq(-0.15, 0.45, by = 0.03)) + 
  labs(x = "Loading PC2 (27.3%)", y = "Metabolite") + 
  geom_vline(xintercept = 0.00, linetype = "dotted") + 
  theme(axis.line = element_line(colour = "black"), 
        panel.background = element_blank(), 
        axis.title.x = element_text(margin = ggplot2::margin(t = 10)), 
        axis.title.y = element_text(margin = ggplot2::margin(r = 10)))


### SCORES ###

# Import ASCA Result - scores
alasca_time_session_scores <-
  read_csv("ALASCA/sampleTime_session/scores_effect_1.csv")

# Set dodge width
dodge <- position_dodge(width = 0.3)

## PC1 ##
alasca_time_session_scores_pc1 <- alasca_time_session_scores %>% 
  filter(PC == "1") %>% 
  mutate(sample_time = factor(sample_time, levels = c("Pre", "Post", "8 h post", "22 h post")))

# Plot - scores
alasca_time_session_scores_pc1_plot <- 
  ggplot(alasca_time_session_scores_pc1, 
         aes(x = sample_time, y = score, group = session, colour = session, 
             shape = session, linetype = session)) + 
  geom_point(position = dodge, size = 3) + geom_line(position = dodge, size = 0.8) + 
  geom_errorbar(aes(ymin = low, ymax = high), width = 0, position = dodge, linewidth = 0.8) +
  # geom_ribbon(aes(ymax = high, ymin = low, fill = session), alpha = 0.2, colour = NA, position = dodge, inherit.aes = T) +
  labs(x = "Sample Time", y = "Score PC1 (47.6%)") + 
  scale_colour_manual(name = "Session", values = c("Control" = "#E69F00", "Session 1" = "#56B4E9", "Session 2" = "black")) +
  scale_shape_manual(name = "Session", values = c(16, 17, 18)) +
  scale_linetype_manual(name = "Session", values = c("Control" = "solid", "Session 1" = "dashed", "Session 2" = "dotdash")) +
  theme(panel.background = element_blank(), 
        axis.line = element_line(), 
        axis.title = element_text(margin = ggplot2::margin(t = 10), size = 15),
        legend.position = "bottom", 
        legend.title = element_text(size = 13), 
        legend.text = element_text(size = 13))

## PC2 ##
alasca_time_session_scores_pc2 <- alasca_time_session_scores %>% 
  filter(PC == "2") %>% 
  mutate(sample_time = factor(sample_time, levels = c("Pre", "Post", "8 h post", "22 h post")))

# Plot - scores
alasca_time_session_scores_pc2_plot <- 
  ggplot(alasca_time_session_scores_pc2, 
         aes(x = sample_time, y = score, group = session, colour = session, 
             shape = session, linetype = session)) + 
  geom_point(position = dodge, size = 3) + geom_line(position = dodge, size = 0.8) + 
  geom_errorbar(aes(ymin = low, ymax = high), width = 0, position = dodge, size = 0.8) +
  # geom_ribbon(aes(ymax = high, ymin = low, fill = session), alpha = 0.2, colour = NA, position = dodge, inherit.aes = T) +
  labs(x = "Sample Time", y = "Score PC2 (27.3%)") + 
  scale_colour_manual(name = "Session", values = c("Control" = "#E69F00", "Session 1" = "#56B4E9", "Session 2" = "black")) +
  scale_shape_manual(name = "Session", values = c(16, 17, 18)) +
  scale_linetype_manual(name = "Session", values = c("Control" = "solid", "Session 1" = "dashed", "Session 2" = "dotdash")) +
  theme(panel.background = element_blank(), 
        axis.line = element_line(), 
        axis.title = element_text(margin = ggplot2::margin(t = 10), size = 15),
        legend.position = "bottom", 
        legend.title = element_text(size = 13), 
        legend.text = element_text(size = 13))

### PREDICTIONS ###

# Import ASCA Results - Predictions
alasca_time_session_pred <- 
  read_csv("ALASCA/sampleTime_session/model_prediction.csv")

# Create factors - Sample time & session
alasca_time_session_pred$sample_time <- factor(alasca_time_session_pred$sample_time,
                                               levels=c("Pre", "Post", "8 h post", "22 h post"),
                                               labels=c("Pre", "Post", "8 h post", "22 h post"))

alasca_time_session_pred$session <- factor(alasca_time_session_pred$session,
                                           levels=c("Control", "Session 1", "Session 2"),
                                           labels=c("Control", "Session 1", "Session 2"))

alasca_time_session_pred_top20_pc1 <- alasca_time_session_pred %>%
  mutate(variable = variable %>% gsub("_[0-9]+_(hil|rp)_(pos|neg)$", "", .)) %>% 
  filter(variable %in% alasca_time_session_loadings_pc1_order)

alasca_time_session_pred_top20_pc2 <- alasca_time_session_pred %>% 
  mutate(variable = variable %>% gsub("_[0-9]+_(hil|rp)_(pos|neg)$", "", .)) %>% 
  filter(variable %in% alasca_time_session_loadings_pc2_order)

dodge_width <- 0.3 #position dodge

# Plot - Individual metabolites
alasca_time_session_pred_pc1_plot <- alasca_time_session_pred_top20_pc1 %>%
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
  theme(strip.background = element_rect(color = "black", fill = "white", linewidth = 1),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4)),
         shape = guide_legend(override.aes = list(size = 4)))

alasca_time_session_pred_pc2_plot <- alasca_time_session_pred_top20_pc2 %>%
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
  theme(strip.background = element_rect(color = "black", fill = "white", linewidth = 1),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4)),
         shape = guide_legend(override.aes = list(size = 4)))



# # #Import ASCA Results - Predictions
# # alasca_time_session_pred <- 
# #   read_csv("ALASCA/SampleTime_Session/model_prediction.csv")
# 
# #Import ASCA Results - Predictions
# alasca_time_session_pred <- 
#   read_csv("ALASCA/2026-07-12T213759/model_prediction.csv")
# 
# # Create factors - Sample time & session
# alasca_time_session_pred$sample_time <- factor(alasca_time_session_pred$sample_time,
#                                                  levels=c("Pre", "Post", "8 h post", "22 h post"),
#                                                  labels=c("Pre", "Post", "8 h post", "22 h post"))
# 
# alasca_time_session_pred$session <- factor(alasca_time_session_pred$session,
#                                              levels=c("Control", "Session 1", "Session 2"),
#                                              labels=c("Control", "Session 1", "Session 2"))
# 
# alasca_time_session_pred_1 <- alasca_time_session_pred %>%
#   dplyr::filter(str_detect(variable, "_acid"))
# 
# # Plot - Individual metabolites
# dodge_width <- 0.3 #position dodge
# 
# plot_pred <- alasca_time_session_pred %>%
#   ggplot(aes(x = sample_time, y = pred, color = session, group = session, shape = session)) +
#   geom_point(size = 3, position = position_dodge(width = dodge_width)) +
#   geom_line(aes(linetype = session), linewidth = 1, position = position_dodge(width = dodge_width)) +
#   geom_errorbar(aes(ymin = low, ymax = high),
#                 width = 0.1, linewidth = 0.8,
#                 position = position_dodge(width = dodge_width)) +
#   scale_color_manual(values = c("Control" = "#E69F00", "Session 1" = "#56B4E9", "Session 2" = "black")) + # Colorblind-safe
#   scale_shape_manual(values = c(16, 17, 18)) +
#   scale_linetype_manual(values = c("Control" = "solid", "Session 1" = "dashed", "Session 2" = "dotdash")) +
#   ylab("Std. Value") +
#   xlab("") +
#   facet_wrap(~variable, scales = "free_y") +
#   theme_minimal(base_size = 14) +
#   theme(
#     strip.background = element_rect(color = "black", fill = "white", linewidth = 1),
#     strip.text = element_text(face = "bold", size = 14),
#     legend.position = "bottom",
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor = element_blank()
#   ) +
#   guides(color = guide_legend(override.aes = list(size = 4)),
#          shape = guide_legend(override.aes = list(size = 4)))
# 
# 
# # Import ASCA Results - Scores
# alasca_time_session_scores <- 
#   read_csv("ALASCA/SampleTime_Session/scores_effect_1.csv")
# 
# # Create factors - Sample time & session
# alasca_time_session_scores$sample_time <- factor(alasca_time_session_scores$sample_time,
#                                       levels=c("Pre", "Post", "8 h post", "22 h post"),
#                                       labels=c("Pre", "Post", "8 h post", "22 h post"))
# 
# alasca_time_session_scores$session <- factor(alasca_time_session_scores$session,
#                                   levels=c("Control", "Session 1", "Session 2"),
#                                   labels=c("Control", "Session 1", "Session 2"))
# 
# ## PC1
# alasca_time_session_scores <- alasca_time_session_scores %>% 
#   filter(PC == "1")
# 
# # Plot - Scores - tidyplots
# alasca_time_session_scores |>
#   tidyplot(x = sample_time, y = score, color = session) |>
#   add_mean_line() |>
#   add_mean_dot() |>
#   add_sem_ribbon() |>
#   adjust_legend_title("Session") |>
#   adjust_x_axis_title("Sample Time") |>
#   adjust_y_axis_title("Score")
# 
# #Plot - Effects - Model Scores and Loadings
# plot(alasca_time_session, 
#      effect=1, 
#      component=1, 
#      type="effect")
# 
# #Plot - Model regression residuals
# plot(alasca_time_session, 
#      effect=1, component=c(1,2), 
#      type="residuals", 
#      variable="sample_time")
# 
# #Plot - Marginal means from regression
# plot(alasca_time_session, 
#      component = 1,
#      effect = 1,
#      type = 'prediction',
#      n_limit = 20,
#      variable = c())



####PART X - Named Metabolites Only####

# #Import - Imputed data frame
# named_metabs_long <- readr::read_csv("Data/SPP_S1 - NamedMetabs_long.csv")
# 
# ####PART X - Clean Data####
# 
# #Rename week to session and filter out PQCs
# named_metabs_long_1 <- named_metabs_long %>%
#   dplyr::filter(sample=="sample") %>%
#   dplyr::relocate(metabolite, .before = value) %>%
#   dplyr::select(-c(unique_id:sample)) %>%
#   dplyr::select(-c(idx:idx_mode)) %>%
#   dplyr::rename(session=week)%>%
#   dplyr::rename(variable=metabolite)
# 
# #Create factors - Sample time & session
# named_metabs_long_1$sample_time <- factor(named_metabs_long_1$sample_time,
#                                       levels=c(0, 1, 2, 3),
#                                       labels=c("Pre", "Post", "8 h post", "22 h post"))
# 
# named_metabs_long_1$session <- factor(named_metabs_long_1$session,
#                                   levels=c(0, 1, 2),
#                                   labels=c("Control", "Session 1", "Session 2"))
# 
# #####PART X - ALASCA Model - Time*Session - All, no filtering#####
# long_named_data_asca <- data.table::setDT(named_metabs_long_1)


#Read in named data
long_named_data_asca <- readr::read_csv("Data/SPP_S1 - long_named_data_asca.csv")

#ASCA Model
alasca_named_time_session <- ALASCA(formula = value ~ sample_time*session + (1|subject_id),
                              df = long_named_data_asca,
                              participant_column="subject_id",
                              scale_function = "sdt1",
                              reduce_dimensions.limit=TRUE,
                              n_validation_runs=1000,
                              max_PC=5,
                              validate = TRUE,
                              use_Rfast= TRUE,
                              save=TRUE,
                              filename = "SPP_S1_sampleTime_session_Named")


### LOADINGS ###

# Import ASCA Results - Loadings
alasca_time_session_named_loadings <-
  read_csv("ALASCA/sampleTime_session_Named/loadings_effect_1.csv")

# Rename "covars" column
alasca_time_session_named_loadings <- alasca_time_session_named_loadings %>%
  dplyr::rename(metabolite = covars)

## PC1 ##

# Filter for PC1
alasca_time_session_named_loadings_pc1 <- alasca_time_session_named_loadings %>% 
  filter(PC == "1")
# Order data
alasca_time_session_named_loadings_pc1 <- 
  alasca_time_session_named_loadings_pc1[order(-alasca_time_session_named_loadings_pc1$loading), ]
# Top 10 variables (5 positive, 5 negatives)
alasca_time_session_named_loadings_pc1_top10 <- 
  alasca_time_session_named_loadings_pc1[c(1:5,245:249), ]
# Order from negative to positive
alasca_time_session_named_loadings_pc1_top10 <- 
  alasca_time_session_named_loadings_pc1_top10[order(alasca_time_session_named_loadings_pc1_top10$loading), ]
# Order of the names 
alasca_time_session_named_loadings_pc1_order <- 
  alasca_time_session_named_loadings_pc1_top10$metabolite

# Plot - loadings
alasca_time_session_named_loadings_pc1_plot <- 
  ggplot(alasca_time_session_named_loadings_pc1_top10, 
         aes(x = loading, y = factor(metabolite, levels = alasca_time_session_named_loadings_pc1_order))) + 
  geom_point() + 
  geom_errorbar(aes(xmin = low, xmax = high), width = .1) + 
  scale_x_continuous(breaks = seq(-0.15, 0.33, by = 0.03)) + 
  labs(x = "Loading PC1 (63.9%)", y = "Metabolite") + 
  geom_vline(xintercept = 0.00, linetype = "dotted") + 
  theme(axis.line = element_line(colour = "black"), 
        panel.background = element_blank(), 
        axis.title.x = element_text(margin = ggplot2::margin(t = 10)), 
        axis.title.y = element_text(margin = ggplot2::margin(r = 10)))


## PC2 ##

# Filter for PC2
alasca_time_session_named_loadings_pc2 <- alasca_time_session_named_loadings %>% 
  filter(PC == "2")
# Order data
alasca_time_session_named_loadings_pc2 <- 
  alasca_time_session_named_loadings_pc2[order(-alasca_time_session_named_loadings_pc2$loading), ]
# Top 10 variables (5 positive, 5 negatives)
alasca_time_session_named_loadings_pc2_top10 <- 
  alasca_time_session_named_loadings_pc2[c(1:5,245:249), ]
# Order from negative to positive
alasca_time_session_named_loadings_pc2_top10 <- 
  alasca_time_session_named_loadings_pc2_top10[order(alasca_time_session_named_loadings_pc2_top10$loading), ]
# Order of the names 
alasca_time_session_named_loadings_pc2_order <- 
  alasca_time_session_named_loadings_pc2_top10$metabolite

# Plot - loadings
alasca_time_session_named_loadings_pc2_plot <- 
  ggplot(alasca_time_session_named_loadings_pc2_top10, 
         aes(x = loading, y = factor(metabolite, levels = alasca_time_session_named_loadings_pc2_order))) + 
  geom_point() + 
  geom_errorbar(aes(xmin = low, xmax = high), width = .1) + 
  scale_x_continuous(breaks = seq(-0.30, 0.33, by = 0.03)) + 
  labs(x = "Loading PC2 (19.7%)", y = "Metabolite") + 
  geom_vline(xintercept = 0.00, linetype = "dotted") + 
  theme(axis.line = element_line(colour = "black"), 
        panel.background = element_blank(), 
        axis.title.x = element_text(margin = ggplot2::margin(t = 10)), 
        axis.title.y = element_text(margin = ggplot2::margin(r = 10)))


### SCORES ###

# Import ASCA Result - scores
alasca_time_session_named_scores <-
  read_csv("ALASCA/sampleTime_session_Named/scores_effect_1.csv")

# Set dodge width
dodge <- position_dodge(width = 0.3)

## PC1 ##
alasca_time_session_named_scores_pc1 <- alasca_time_session_named_scores %>% 
  filter(PC == "1") %>% 
  mutate(sample_time = factor(sample_time, levels = c("Pre", "Post", "8 h post", "22 h post")))

# Plot - scores
alasca_time_session_named_scores_pc1_plot <- 
  ggplot(alasca_time_session_named_scores_pc1, 
         aes(x = sample_time, y = score, group = session, colour = session, 
             shape = session, linetype = session)) + 
  geom_point(position = dodge, size = 3) + geom_line(position = dodge, size = 0.8) + 
  geom_errorbar(aes(ymin = low, ymax = high), width = 0, position = dodge, linewidth = 0.8) +
  # geom_ribbon(aes(ymax = high, ymin = low, fill = session), alpha = 0.2, colour = NA, position = dodge, inherit.aes = T) +
  labs(x = "Sample Time", y = "Score PC1 (63.9%)") + 
  scale_colour_manual(name = "Session", values = c("Control" = "#E69F00", "Session 1" = "#56B4E9", "Session 2" = "black")) +
  scale_shape_manual(name = "Session", values = c(16, 17, 18)) +
  scale_linetype_manual(name = "Session", values = c("Control" = "solid", "Session 1" = "dashed", "Session 2" = "dotdash")) +
  theme(panel.background = element_blank(), 
        axis.line = element_line(), 
        axis.title = element_text(margin = ggplot2::margin(t = 10), size = 15),
        legend.position = "bottom", 
        legend.title = element_text(size = 13), 
        legend.text = element_text(size = 13))

## PC2 ##
alasca_time_session_named_scores_pc2 <- alasca_time_session_named_scores %>% 
  filter(PC == "2") %>% 
  mutate(sample_time = factor(sample_time, levels = c("Pre", "Post", "8 h post", "22 h post")))

# Plot - scores
alasca_time_session_named_scores_pc2_plot <- 
  ggplot(alasca_time_session_named_scores_pc2, 
         aes(x = sample_time, y = score, group = session, colour = session, 
             shape = session, linetype = session)) + 
  geom_point(position = dodge, size = 3) + geom_line(position = dodge, size = 0.8) + 
  geom_errorbar(aes(ymin = low, ymax = high), width = 0, position = dodge, linewidth = 0.8) +
  # geom_ribbon(aes(ymax = high, ymin = low, fill = session), alpha = 0.2, colour = NA, position = dodge, inherit.aes = T) +
  labs(x = "Sample Time", y = "Score PC2 (19.7%)") + 
  scale_colour_manual(name = "Session", values = c("Control" = "#E69F00", "Session 1" = "#56B4E9", "Session 2" = "black")) +
  scale_shape_manual(name = "Session", values = c(16, 17, 18)) +
  scale_linetype_manual(name = "Session", values = c("Control" = "solid", "Session 1" = "dashed", "Session 2" = "dotdash")) +
  theme(panel.background = element_blank(), 
        axis.line = element_line(), 
        axis.title = element_text(margin = ggplot2::margin(t = 10), size = 15),
        legend.position = "bottom", 
        legend.title = element_text(size = 13), 
        legend.text = element_text(size = 13))

### PREDICTIONS ###

# Import ASCA Results - Predictions
alasca_time_session_named_pred <- 
  read_csv("ALASCA/sampleTime_session_Named/model_prediction.csv")

# Create factors - Sample time & session
alasca_time_session_named_pred$sample_time <- factor(alasca_time_session_named_pred$sample_time,
                                               levels=c("Pre", "Post", "8 h post", "22 h post"),
                                               labels=c("Pre", "Post", "8 h post", "22 h post"))

alasca_time_session_named_pred$session <- factor(alasca_time_session_named_pred$session,
                                           levels=c("Control", "Session 1", "Session 2"),
                                           labels=c("Control", "Session 1", "Session 2"))

alasca_time_session_named_pred_top10_pc1 <- alasca_time_session_named_pred %>%
  filter(variable %in% alasca_time_session_named_loadings_pc1_order)

alasca_time_session_named_pred_top10_pc2 <- alasca_time_session_named_pred %>% 
  filter(variable %in% alasca_time_session_named_loadings_pc2_order)

dodge_width <- 0.3 #position dodge

# Plot - Individual metabolites
alasca_time_session_named_pred_pc1_plot <- alasca_time_session_named_pred_top10_pc1 %>%
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
  theme(strip.background = element_rect(color = "black", fill = "white", linewidth = 1),
        strip.text = element_text(face = "bold", size = 10),
        legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4)),
         shape = guide_legend(override.aes = list(size = 4)))

alasca_time_session_named_pred_pc2_plot <- alasca_time_session_named_pred_top10_pc2 %>%
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
  theme(strip.background = element_rect(color = "black", fill = "white", linewidth = 1),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4)),
         shape = guide_legend(override.aes = list(size = 4)))


# #Plot - Effects - Model Scores and Loadings
# plot(alasca_named_time_session, 
#      effect=1, 
#      component=c(1,2), 
#      type="effect")
# 
# #Plot - Model regression residuals
# plot(alasca_named_time_session, 
#      effect=1, component=c(1,2), 
#      type="residuals", 
#      variable="sample_time")
# 
# #Plot - Marginal means from regression - PC1
# plot(alasca_named_time_session, 
#      component = 1,
#      effect = 1,
#      type = 'prediction',
#      n_limit = 20,
#      variable = c())
# 
# #Plot - Marginal means from regression - PC2
# plot(alasca_named_time_session, 
#      component = 2,
#      effect = 1,
#      type = 'prediction',
#      n_limit = 20,
#      variable = c())
# 
# #Import ASCA Results  - Loadings
# alasca_time_session_named_loadings <-
#   read_csv("ALASCA/sampleTime_session_Named/loadings_effect_1.csv")
# 
# # Rename "covars" column
# alasca_time_session_named_loadings <- alasca_time_session_named_loadings %>%
#   dplyr::rename(metabolite=covars)
# 
# #Import ASCA Results - Predictions
# alasca_time_session_named_pred <- 
#   read_csv("ALASCA/sampleTime_session_Named/model_prediction.csv")
# 
# # Create factors - Sample time & session
# alasca_time_session_named_pred$sample_time <- factor(alasca_time_session_named_pred$sample_time,
#                                                      levels=c("Pre", "Post", "8 h post", "22 h post"),
#                                                      labels=c("Pre", "Post", "8 h post", "22 h post"))
# 
# alasca_time_session_named_pred$session <- factor(alasca_time_session_named_pred$session,
#                                                  levels=c("Control", "Session 1", "Session 2"),
#                                                  labels=c("Control", "Session 1", "Session 2"))
# 
# alasca_time_session_named_pred_1 <- alasca_time_session_named_pred %>%
#   dplyr::filter(str_detect(variable, "_acid"))
# 
# # Plot - Individual metabolites
# dodge_width <- 0.3 #position dodge
# 
# plot_named_pred <- alasca_time_session_named_pred_1 %>%
#   ggplot(aes(x = sample_time, y = pred, color = session, group = session, shape = session)) +
#   geom_point(size = 3, position = position_dodge(width = dodge_width)) +
#   geom_line(aes(linetype = session), linewidth = 1, position = position_dodge(width = dodge_width)) +
#   geom_errorbar(aes(ymin = low, ymax = high),
#                 width = 0.1, linewidth = 0.8,
#                 position = position_dodge(width = dodge_width)) +
#   scale_color_manual(values = c("Control" = "#E69F00", "Session 1" = "#56B4E9", "Session 2" = "black")) + # Colorblind-safe
#   scale_shape_manual(values = c(16, 17, 18)) +
#   scale_linetype_manual(values = c("Control" = "solid", "Session 1" = "dashed", "Session 2" = "dotdash")) +
#   ylab("Std. Value") +
#   xlab("") +
#   facet_wrap(~variable, scales = "free_y") +
#   theme_minimal(base_size = 14) +
#   theme(
#     strip.background = element_rect(color = "black", fill = "white", linewidth = 1),
#     strip.text = element_text(face = "bold", size = 14),
#     legend.position = "bottom",
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor = element_blank()
#   ) +
#   guides(color = guide_legend(override.aes = list(size = 4)),
#          shape = guide_legend(override.aes = list(size = 4)))
# 
# # Import ASCA Results - Scores
# alasca_time_session_named_scores <- 
#   read_csv("ALASCA/sampleTime_session_Named/scores_effect_1.csv")
# 
# # Create factors - Sample time & session
# alasca_time_session_named_scores$sample_time <- factor(alasca_time_session_named_scores$sample_time,
#                                                        levels=c("Pre", "Post", "8 h post", "22 h post"),
#                                                        labels=c("Pre", "Post", "8 h post", "22 h post"))
# 
# alasca_time_session_named_scores$session <- factor(alasca_time_session_named_scores$session,
#                                                    levels=c("Control", "Session 1", "Session 2"),
#                                                    labels=c("Control", "Session 1", "Session 2"))
# 
# ## PC1
# alasca_time_session_named_scores <- alasca_time_session_named_scores %>% 
#   filter(PC == "1")
# 
# # Plot - Scores - tidyplots
# alasca_time_session_named_scores |>
#   tidyplot(x = sample_time, y = score, color = session) |>
#   add_mean_line() |>
#   add_mean_dot() |>
#   add_sem_ribbon() |>
#   adjust_legend_title("Session") |>
#   adjust_x_axis_title("Sample Time") |>
#   adjust_y_axis_title("Score")



