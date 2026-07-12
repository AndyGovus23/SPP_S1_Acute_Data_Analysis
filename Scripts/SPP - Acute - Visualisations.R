#### Load Packages ####
library(tidyverse)
library(readr)
library(patchwork)

#### Import Results Data ####
step_test_df <- read_csv("Data/SPP_S1 - stepTestData.csv")

#### Clean data ####

# Create gender column
females <- c("Lily Koch", "Semra Olowoniyi", "Aleisha Clark", "Madison Brand", 
             "Makayla Larkin")

step_test_df <- step_test_df %>% 
  mutate(sex = if_else(Swimmer %in% females, "Female", "Male"))

# Convert date column
step_test_df$Date <- as.Date(step_test_df$Date, format = "%d/%m/%Y")

# Create summary data by gender
sex_avg <- step_test_df %>% 
  group_by(Session, sex, Set) %>% 
  summarise(avg_time = mean(`200m Time (s)`),
            se_time = sd(`200m Time (s)`, na.rm = TRUE)/sqrt(n()),
            avg_bla = mean(Lactate), 
            se_bla = sd(Lactate, na.rm = TRUE)/sqrt(n()), 
            .groups = "drop")

# Function to format time
format_mmss <- function(seconds) {
  sprintf("%d:%04.1f", floor(seconds / 60), seconds %% 60)
}

y_breaks <- seq(115, 193, by = 3)

#### Create Plot ####
time_plot <- sex_avg %>% 
  ggplot() + 
  geom_line(aes(x = Set, y = avg_time, group = sex, colour = sex), linewidth = 1.2) + 
  geom_point(aes(x = Set, y = avg_time, colour = sex), size = 3) +
  geom_errorbar(aes(x = Set, y = avg_time, ymin = avg_time - se_time, 
                    ymax = avg_time + se_time, colour = sex), 
                width = 0.2, linewidth = 0.8) +
  scale_y_continuous(breaks = y_breaks, labels = format_mmss(y_breaks)) +
  scale_x_continuous(breaks = 1:7) +
  scale_color_manual(values = c("Female" = "grey", "Male" = 'black')) + 
  labs(y = "200 m Time (mm:ss)", x = "Step", title = "A") +
  facet_wrap(~Session, labeller = labeller(Session = c("1" = "INC1", "2" = "INC2"))) +
  theme(plot.background = element_rect(fill = 'white', colour = NA), 
        panel.background = element_rect(fill = 'white', colour = NA), 
        plot.title = element_text(size = 24, face = "bold"), 
        axis.line = element_line(), 
        axis.title.x = element_text(size = 14, face = "italic", margin = margin(t = 8)), 
        axis.title.y = element_text(size = 14, face = "italic", margin = margin(r = 8)), 
        axis.text = element_text(size = 13, colour = 'black'), 
        legend.title = element_blank(), 
        legend.position = "none", 
        strip.text = element_text(size = 14, face = "bold", colour = "black"), 
        strip.background = element_rect(fill = "lightgrey", colour = "black", linewidth = 0.5))

lactate_plot <- sex_avg %>% 
  ggplot() +
  geom_line(aes(x = Set, y = avg_bla, group = sex, colour = sex), 
            linewidth = 1.2, position = position_dodge(width = 0.5)) + 
  geom_point(aes(x = Set, y = avg_bla, group = sex, colour = sex), 
             size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(x = Set, y = avg_bla, ymin = avg_bla - se_bla, 
                    ymax = avg_bla + se_bla, colour = sex), 
                width = 0.2, linewidth = 0.8, position = position_dodge(width = 0.5)) +
  scale_x_continuous(breaks = 1:7) +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, by = 1)) + 
  scale_color_manual(values = c("Female" = "grey", "Male" = 'black')) + 
  labs(y = "Blood Lactate Concentration (mmol/L)", x = "Step", title = "B") + 
  facet_wrap(~Session, labeller = labeller(Session = c("1" = "INC1", "2" = "INC2"))) + 
  theme(plot.background = element_rect(fill = 'white', colour = NA), 
        panel.background = element_rect(fill = 'white', colour = NA), 
        plot.title = element_text(size = 24, face = "bold"), 
        axis.line = element_line(), 
        axis.title.x = element_text(size = 14, face = "italic", margin = margin(t = 8)), 
        axis.title.y = element_text(size = 14, face = "italic", margin = margin(r = 8)), 
        axis.text = element_text(size = 12, colour = 'black'), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14), 
        legend.position = "right", 
        strip.text = element_text(size = 14, face = "bold", colour = "black"), 
        strip.background = element_rect(fill = "lightgrey", colour = "black", linewidth = 0.5))

combined_step_test_plot <- time_plot + lactate_plot

ggsave("Images/time_bla_results.png", plot = combined_step_test_plot, width = 10, 
       height = 6, dpi = 600, bg = "white")
  


