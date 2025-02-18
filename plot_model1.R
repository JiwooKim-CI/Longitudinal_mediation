# Load required libraries
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(jtools)
library(plyr)
library(tikzDevice)
library(distributional)
library(ggdist)

# Set theme
theme_set(theme_ggdist())

### 1. Load and Process Simulation Data ###

# Function to load and preprocess datasets
load_and_preprocess <- function(file, method_name) {
  df <- read_csv(file)
  df[,1] <- method_name
  colnames(df) <- c("method", "value", "diff", "n")
  df <- aggregate(n ~ method + value, data = df, FUN = sum)
  return(df)
}

# Load datasets
ancova_sim1_1  <- load_and_preprocess("ancova_sim1_1.csv", "ANCOVA")
change_sim1_1  <- load_and_preprocess("change_sim1_1.csv", "Change (F)")
change2_sim1_1 <- load_and_preprocess("change2_sim1_1.csv", "Change (A)")
naive_sim1_1   <- load_and_preprocess("naive_sim1_1.csv", "Naive")
change1_sim1_1 <- load_and_preprocess("change1_sim1_1.csv", "Change (D)")

# Combine datasets
graph1 <- rbind(ancova_sim1_1, change_sim1_1, change2_sim1_1, naive_sim1_1, change1_sim1_1) %>%
  group_by(method) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

### 2. Create and Save First Plot ###
new_sim1 <- ggplot(graph1, aes(x = method, y = value, fill = method, weight = percent)) +
  stat_halfeye(adjust = 0.5, justification = -0.2, .width = 0, point_colour = NA, show.legend = FALSE) +
  geom_boxplot(width = 0.12, outlier.color = NA, alpha = 0.5, show.legend = FALSE) +
  labs(x = "Method", y = "Raw Bias") +
  coord_flip() +
  ylim(-0.5, 0.5) +
  theme_apa()

# Save plots
tikz("new_sim1.tex", width = 5.5, height = 5, standAlone = FALSE)
new_sim1
dev.off()

pdf("new_sim1.pdf", height = 6, width = 6)
new_sim1
dev.off()

### 3. Load and Process Second Simulation Data ###
ancova_sim1_c  <- load_and_preprocess("ancova_sim1_c.csv", "ANCOVA")
change_sim1_c  <- load_and_preprocess("change_sim1_c.csv", "Change (F)")
change2_sim1_c <- load_and_preprocess("change2_sim1_c.csv", "Change (A)")
naive_sim1_c   <- load_and_preprocess("naive_sim1_c.csv", "Naive")
change1_sim1_c <- load_and_preprocess("change1_sim1_c.csv", "Change (D)")

# Combine datasets
graph2 <- rbind(ancova_sim1_c, change_sim1_c, change2_sim1_c, naive_sim1_c, change1_sim1_c) %>%
  group_by(method) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

### 4. Create and Save Second Plot ###
new_sim2 <- ggplot(graph2, aes(x = method, y = value, fill = method, weight = percent)) +
  stat_halfeye(adjust = 5, justification = -0.2, .width = 0, point_colour = NA, show.legend = FALSE) +
  geom_boxplot(width = 0.12, outlier.color = NA, alpha = 0.5, show.legend = FALSE) +
  labs(x = "Method", y = "Raw Bias") +
  coord_flip() +
  ylim(-0.5, 0.5) +
  theme_apa()

# Save plots
tikz("new_sim2.tex", width = 5.5, height = 5, standAlone = FALSE)
new_sim2
dev.off()

pdf("new_sim2.pdf", height = 6, width = 6)
new_sim2
dev.off()

### 5. Load and Process Third Simulation Data ###
ancova_sim1_3  <- load_and_preprocess("ancova_sim1_3.csv", "ANCOVA")
change_sim1_3  <- load_and_preprocess("change_sim1_3.csv", "Change (F)")
change2_sim1_3 <- load_and_preprocess("change2_sim1_3.csv", "Change (A)")
naive_sim1_3   <- load_and_preprocess("naive_sim1_3.csv", "Naive")
change1_sim1_3 <- load_and_preprocess("change1_sim1_3.csv", "Change (D)")

# Combine datasets
graph3 <- rbind(ancova_sim1_3, change_sim1_3, change2_sim1_3, naive_sim1_3, change1_sim1_3) %>%
  group_by(method) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

### 6. Create and Save Third Plot ###
new_sim3 <- ggplot(graph3, aes(x = method, y = value, fill = method, weight = percent)) +
  stat_halfeye(adjust = 0.5, justification = -0.2, .width = 0, point_colour = NA, show.legend = FALSE) +
  geom_boxplot(width = 0.12, outlier.color = NA, alpha = 0.5, show.legend = FALSE) +
  labs(x = "Method", y = "Raw Bias") +
  coord_flip() +
  ylim(-0.5, 0.5) +
  theme_apa()

# Save plots
tikz("new_sim3.tex", width = 5.5, height = 5, standAlone = FALSE)
new_sim3
dev.off()

pdf("new_sim3.pdf", height = 6, width = 6)
new_sim3
dev.off()

### 7. Violation of Ignorability Assumption Plot ###

# Load grid simulation data
ancova_grid <- read_csv("ancova_sim1_1.csv") %>%
  select(-1) %>%
  rename(value = X1) %>%
  mutate(method = "ANCOVA")

change_grid <- read_csv("change_sim1_1.csv") %>%
  select(-1) %>%
  rename(value = X1) %>%
  mutate(method = "Change (F)")

change2_grid <- read_csv("change2_sim1_1.csv") %>%
  select(-1) %>%
  rename(value = X1) %>%
  mutate(method = "Change (A)")

change1_grid <- read_csv("change1_sim1_1.csv") %>%
  select(-1) %>%
  rename(value = X1) %>%
  mutate(method = "Change (D)")

naive_grid <- read_csv("naive_sim1_1.csv") %>%
  select(-1, -4) %>%
  rename(value = X1) %>%
  mutate(method = "Naive")

# Combine data
sim1 <- rbind(ancova_grid, change_grid, change1_grid, change2_grid, naive_grid)

change_summary <- sim1 %>%
  filter(method %in% c("Change (F)", "Change (A)", "Change (D)")) %>%
  mutate(method = "Change (Combined)")  # Label for the combined line

other_methods <- sim1 %>%
  filter(!method %in% c("Change (F)", "Change (A)", "Change (D)"))
combined_data <- bind_rows(other_methods, change_summary)

# Create bias plot
model_1_ig <- ggplot(combined_data, aes(x = diff, y = abs(value), color = method)) +
  stat_smooth(method = "loess", aes(weight = n), se = FALSE) +
  facet_grid("method") +
  scale_color_manual(
    values = c(
      "ANCOVA" = "#F8766D",        
      "Naive" = "#C77CFF",     
      "Change (Combined)" = "black"  
    )) +
  theme_apa() +
  labs(x = "Violation of Ignorability Assumption", y = "Bias")

# Save plot
tikz("sim_1_ig.tex", height = 6, width = 6)
model_1_ig
dev.off()

pdf("sim_1_ig.pdf", height = 6, width = 6)
model_1_ig
dev.off()
