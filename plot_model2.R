# Load required libraries
library(readr)
library(dplyr)
library(ggplot2)
library(jtools)
library(ggdist)
library(tikzDevice)

theme_set(theme_ggdist())  # Set theme

### 1. Function to Load & Process Data ###
load_and_preprocess <- function(file_list, method_name) {
  df <- do.call(rbind, lapply(file_list, read_csv))
  df <- df[,-c(1,3)]  # Remove unnecessary columns
  colnames(df)[1] <- "value"
  df <- aggregate(n ~ value, df, sum)
  df$perc <- df$n / sum(df$n)
  df$method <- method_name
  return(df)
}

# Define file names
ancova_files <- c("ancova_2_1.csv", "ancova_2_2.csv", "ancova_2_3.csv")
change_files <- c("change_2_1.csv", "change_2_2.csv", "change_2_3.csv")
change2_files <- c("change2_2_1.csv", "change2_2_2.csv", "change2_2_3.csv")
change1_files <- c("change1_2_1.csv", "change1_2_2.csv", "change1_2_3.csv")
naive_files <- c("naive_2_1.csv", "naive_2_2.csv", "naive_2_3.csv")

# Load and preprocess data
ancova_grid <- load_and_preprocess(ancova_files, "ANCOVA")
change_grid <- load_and_preprocess(change_files, "Change(F)")
change2_grid <- load_and_preprocess(change2_files, "Change(A)")
change1_grid <- load_and_preprocess(change1_files, "Change(D)")
naive_grid <- load_and_preprocess(naive_files, "Naive")

# Combine all data
graph_sim2 <- rbind(ancova_grid, change_grid, change2_grid, change1_grid, naive_grid)

### 2. Generate Plot ###
new_sim2_1 <- ggplot(graph_sim2, aes(x = method, y = value, fill = method, weight = n)) +
  stat_halfeye(adjust = 0.5, justification = -0.2, .width = 0, point_colour = NA, show.legend = FALSE) +
  geom_boxplot(width = 0.12, outlier.color = NA, alpha = 0.5, show.legend = FALSE) +
  labs(x = "Method", y = "Raw Bias") +
  coord_flip() +
  ylim(-0.5, 0.5) +
  theme_apa()

# Save Plot
tikz("new_sim2_1.tex", width = 5.5, height = 5, standAlone = TRUE)
new_sim2_1
dev.off()

pdf("new_sim2_1.pdf", height = 6, width = 6)
new_sim2_1
dev.off()

### 3. Process Ignorability Assumption Data ###
# Define new files for ignorability assumption
ancova_files_i <- c("ancova_2_1_i.csv", "ancova_2_2_i.csv", "ancova_2_3_i.csv")
change_files_i <- c("change_2_1_i.csv", "change_2_2_i.csv", "change_2_3_i.csv")
change2_files_i <- c("change2_2_1_i.csv", "change2_2_2_i.csv", "change2_2_3_i.csv")
change1_files_i <- c("change1_2_1_i.csv", "change1_2_2_i.csv", "change1_2_3_i.csv")
naive_files_i <- c("naive_2_1_i.csv", "naive_2_2_i.csv", "naive_2_3_i.csv")

# Load and preprocess data
ancova_grid_i <- load_and_preprocess(ancova_files_i, "ANCOVA")
change_grid_i <- load_and_preprocess(change_files_i, "Change(F)")
change2_grid_i <- load_and_preprocess(change2_files_i, "Change(A)")
change1_grid_i <- load_and_preprocess(change1_files_i, "Change(D)")
naive_grid_i <- load_and_preprocess(naive_files_i, "Naive")

# Combine all data
graph_sim2_i <- rbind(ancova_grid_i, change_grid_i, change2_grid_i, change1_grid_i, naive_grid_i)

### 4. Generate Ignorability Plot ###
new_sim2_3 <- ggplot(graph_sim2_i, aes(x = method, y = value, fill = method, weight = n)) +
  stat_halfeye(adjust = 0.5, justification = -0.2, .width = 0, point_colour = NA, show.legend = FALSE) +
  geom_boxplot(width = 0.12, outlier.color = NA, alpha = 0.5, show.legend = FALSE) +
  labs(x = "Method", y = "Raw Bias") +
  coord_flip() +
  ylim(-0.5, 0.5) +
  theme_apa()

# Save Plot
tikz("new_sim2_3.tex", width = 5.5, height = 5, standAlone = TRUE)
new_sim2_3
dev.off()

pdf("new_sim2_3.pdf", height = 6, width = 6)
new_sim2_3
dev.off()

### 5. Process Data for Bias vs Ignorability Assumption ###
sim2 <- rbind(ancova_grid, change_grid, change1_grid, change2_grid, naive_grid)

# Combine "Change" methods into one
change_summary <- sim2 %>%
  filter(method %in% c("Change(F)", "Change(A)", "Change(D)")) %>%
  mutate(method = "Change (Combined)")

# Keep other methods
other_methods <- sim2 %>%
  filter(!method %in% c("Change(F)", "Change(A)", "Change(D)"))

# Bind data together
combined_data <- bind_rows(other_methods, change_summary)

### 6. Generate Bias vs Ignorability Assumption Plot ###
model_2_ig_2 <- ggplot(combined_data, aes(x = diff, y = abs(value), linetype = method, color = method)) +
  stat_smooth(method = "loess", aes(weight = n), se = FALSE) +
  scale_color_manual(
    values = c(
      "ANCOVA" = "#F8766D",
      "Naive" = "#C77CFF",
      "Change (Combined)" = "black"
    )
  ) +
  theme_apa() +
  labs(
    x = "Violation of Ignorability Assumption",
    y = "Bias",
    color = "Method",   # Legend title for colors
    linetype = "Method" # Legend title for line types
  ) +
  ylim(0, 0.1)

# Save Plot
tikz("sim_2_ig_2.tex", width = 5.5, height = 5, standAlone = TRUE)
model_2_ig_2
dev.off()

pdf("sim_2_ig_2.pdf", height = 6, width = 6)
model_2_ig_2
dev.off()
