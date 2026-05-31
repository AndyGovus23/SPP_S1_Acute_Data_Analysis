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
  group_by(sex, Set) %>% 
  summarise(avg_time = mean(`200m Time (s)`),
            avg_bla = mean(Lactate), 
            .groups = "drop")

# Function to format time
format_mmss <- function(seconds) {
  sprintf("%d:%04.1f", floor(seconds / 60), seconds %% 60)
}

y_breaks <- seq(120, 178, by = 3)

#### Create Plot ####
time_plot <- sex_avg %>% 
  ggplot() + 
  geom_line(aes(x = Set, y = avg_time, group = sex, colour = sex), linewidth = 1.2) + 
  geom_point(aes(x = Set, y = avg_time, group = sex, colour = sex), size = 3) +
  scale_y_continuous(breaks = y_breaks, labels = format_mmss(y_breaks)) +
  scale_x_continuous(breaks = 1:7) +
  scale_color_manual(values = c("Female" = "grey", "Male" = 'black')) + 
  labs(y = "200 m Time (mm:ss)", title = "Average 200 m Time for Step 1-7") + 
  theme(plot.background = element_rect(fill = 'white', colour = NA), 
        panel.background = element_rect(fill = 'white', colour = NA), 
        axis.line = element_line(), 
        plot.title = element_text(face = "bold", size = 16, hjust = 0.05), 
        axis.title.x = element_text(size = 13, face = "italic", margin = margin(t = 8)), 
        axis.title.y = element_text(size = 13, face = "italic", margin = margin(r = 8)), 
        axis.text = element_text(size = 12, colour = 'black'), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12), 
        legend.position = "right")

lactate_plot <- sex_avg %>% 
  ggplot() +
  geom_line(aes(x = Set, y = avg_bla, group = sex, colour = sex), 
            linewidth = 1.2, position = position_dodge(width = 0.5)) + 
  geom_point(aes(x = Set, y = avg_bla, group = sex, colour = sex), 
             size = 3, position = position_dodge(width = 0.5)) +
  scale_x_continuous(breaks = 1:7) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 1)) + 
  scale_color_manual(values = c("Female" = "grey", "Male" = 'black')) + 
  labs(y = "Blood Lactate Concentration (mmol/L)", title = "Average Blood Lactate for Step 1-7") + 
  theme(plot.background = element_rect(fill = 'white', colour = NA), 
        panel.background = element_rect(fill = 'white', colour = NA), 
        axis.line = element_line(), 
        plot.title = element_text(face = "bold", size = 16, hjust = 0.05), 
        axis.title.x = element_text(size = 13, face = "italic", margin = margin(t = 8)), 
        axis.title.y = element_text(size = 13, face = "italic", margin = margin(r = 8)), 
        axis.text = element_text(size = 12, colour = 'black'), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12), 
        legend.position = "right")

combined_step_test_plot <- time_plot + lactate_plot

ggsave("Images/time_bla_results.png", plot = combined_step_test_plot, width = 10, 
       height = 6, dpi = 600, bg = "white")
  


