# lipid_qc_all_data.R

# Libraries ---------------------------------------------------------------

library(scales)
library(ggsci)
library(ggrepel)
library(ggpubr)
library(writexl)
library(readxl)
library(tidyverse)
setwd("~/Documents/R/LJ ")

# Minor stuff -------------------------------------------------------------

# PIQ base color
piq_color = "#DE246B"
accent_color <- "#03CACE"

# Set the default ggplto theme to avoid re-typing
my_ggtheme <- theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(plot.subtitle = element_text(face = "italic"), plot.caption = element_text(face = "italic"))
theme_set(my_ggtheme)

# Function to save standard pdf and png files for this run of scripts
my_ggsave <- function(figin, width = 8, height = 6, units = "in") {
  pdffileout = paste0("RESULTS/", figin, ".pdf")
  ggsave(pdffileout, width = width, height = height, units = units)
  pngfileout = paste0("RESULTS/", figin, ".png")
  ggsave(pngfileout, width = width, height = height, units = units, dpi = "retina")
}

# Data --------------------------------------------------------------------

# Load and munge the data file provided by Yuya 1/19
df.raw_data <- read_xlsx("LN001_MET_POS.xlsx", range = "A3:CG199", col_names = FALSE, col_types = c("text",rep("numeric", 84)))

# Load the com
v.col_names <- read_xlsx("LN001_MET_POS.xlsx", range = "A1:CG1", col_names = FALSE, col_types = c("text")) %>%
  transpose() %>%
  unlist()
v.col_names[1] <- "lipid_id"
colnames(df.raw_data) <- v.col_names

# Convert to long form and remove the NA, tag QC spike vs endogenous lipids, ln transform
df.data <- df.raw_data %>%
  pivot_longer(cols = c(PQC_1:last_col()), names_to = "sample_id", values_to = "intensity") %>%
  filter(! is.na(intensity)) %>%
  mutate(l_intensity = log(intensity)) %>%
  select(-intensity) %>%
  mutate(lipid_group = ifelse(str_detect(lipid_id, "QPRESS"), "Spike", "Endogenous")) %>%
  extract(sample_id, c("sample_type","sample_index"), "(\\w{3})_(\\d{1,2})", remove = FALSE) %>%
  mutate(sample_index = as.numeric(sample_index) - 1) %>%
  mutate(plate_id = (sample_index %/% 5) + 1) %>%
  mutate(plate_id = factor(plate_id)) %>%
  select(plate_id, sample_id, sample_type, sample_index, lipid_id, lipid_group, l_intensity)

# Plot the distribution of the mean intensity values by plate, calling out the spikes
df.plate_means <- df.data %>%
  group_by(plate_id, sample_type, lipid_group, lipid_id) %>%
  summarize(mean_intensity = mean(l_intensity), count = n(), sd_intensity = sd(l_intensity), cv_intensity = sqrt(exp(sd_intensity^2)-1), .groups = "drop")
 
#median intensity 
df.plate_intensity_medians <- df.plate_means %>%
  group_by(plate_id, sample_type) %>%
  summarize(median_intensity = median(mean_intensity, na.rm = TRUE), .groups = "drop") %>%
  mutate(median_intensity_label = round(median_intensity,2))

# Median normalization ----------------------------------------------------
# Get sample normalization factors based on common features
df.norm_factors <- df.data %>%
  mutate(max_samples = n_distinct(sample_id)) %>%
  group_by(lipid_id) %>%
  mutate(num_detected = n_distinct(sample_id)) %>%
  ungroup() %>%
  filter(num_detected == max_samples) %>%
  group_by(sample_id) %>%
  summarize(median_sample_intensity = median(l_intensity), .groups = "drop") %>%
  mutate(mean_median = mean(median_sample_intensity)) %>%
  group_by(sample_id) %>%
  summarize(norm_factor = mean_median/median_sample_intensity, .groups = "drop")
# Apply the norm factors to the raw data
df.metabolite_mednorm <- df.data %>%
  left_join(df.norm_factors, by = "sample_id") %>%
  mutate(l_intensity_mednorm = l_intensity * norm_factor) 

# Plot the normalized distribution of the mean intensity values by plate, calling out the spikes
df.plate_norm_means <- df.metabolite_mednorm %>%
  group_by(plate_id, sample_type, lipid_group, lipid_id) %>%
  summarize(mean_intensity = mean(l_intensity_mednorm), count = n(), sd_intensity = sd(l_intensity_mednorm), cv_intensity = sqrt(exp(sd_intensity^2)-1), .groups = "drop")

#median norm intensity 
df.plate_norm_intensity_medians <- df.plate_norm_means %>%
  group_by(plate_id, sample_type) %>%
  summarize(median_intensity = median(mean_intensity, na.rm = TRUE), .groups = "drop") %>%
  mutate(median_intensity_label = round(median_intensity,2))

# Plot the intensity ranges
fig_intensities <- ggplot(df.plate_means, aes(plate_id, mean_intensity)) +
  scale_fill_npg() +
  scale_color_npg() +
  geom_boxplot(aes(color = sample_type), show.legend = FALSE, width = 0.5) +
  geom_line(data = filter(df.plate_means, lipid_group == "Spike"), aes(group = lipid_id), color = "lightgray", size = 0.25) +
  geom_point(data = filter(df.plate_means, lipid_group == "Spike"), aes(fill = sample_type), shape = 21, alpha = 0.3, show.legend = FALSE) +
  geom_label(data = df.plate_intensity_medians, aes(y = median_intensity, label = median_intensity_label), size = 2.5) +
  facet_grid(sample_type ~.) +
  labs(x = "Plate ID",
       y = "Mean ln Intensity",
       title = "Distribution of Metabolite Intensities",
       subtitle = "Mean of intensities by plate are plotted",
       caption = "Spike metabolites are highlighted as circles and connected by lines")

# Save the plot
my_ggsave("fig_intensities")

# Plot the normalized intensity ranges
fig_norm_intensities <- ggplot(df.plate_norm_means, aes(plate_id, mean_intensity)) +
  scale_fill_npg() +
  scale_color_npg() +
  geom_boxplot(aes(color = sample_type), show.legend = FALSE, width = 0.5) +
  geom_line(data = filter(df.plate_norm_means, lipid_group == "Spike"), aes(group = lipid_id), color = "lightgray", size = 0.25) +
  geom_point(data = filter(df.plate_norm_means, lipid_group == "Spike"), aes(fill = sample_type), shape = 21, alpha = 0.3, show.legend = FALSE) +
  geom_label(data = df.plate_norm_intensity_medians, aes(y = median_intensity, label = median_intensity_label), size = 2.5) +
  facet_grid(sample_type ~.) +
  labs(x = "Plate ID",
       y = "Mean ln Intensity",
       title = "Distribution of Metabolite Intensities",
       subtitle = "Mean of intensities by plate are plotted",
       caption = "Spike metabolites are highlighted as circles and connected by lines")

# Save the plot
my_ggsave("fig_norm_intensities")

# Get the median CV's for plotting
df.plate_means_medians <- df.plate_means %>%
  group_by(plate_id, sample_type) %>%
  summarize(median_cv = median(cv_intensity, na.rm = TRUE), .groups = "drop") %>%
  mutate(median_label = percent(median_cv, accuracy = 0.1))

# Get the normalized median CV's for plotting
df.plate_norm_means_medians <- df.plate_norm_means %>%
  group_by(plate_id, sample_type) %>%
  summarize(median_cv = median(cv_intensity, na.rm = TRUE), .groups = "drop") %>%
  mutate(median_label = percent(median_cv, accuracy = 0.1))

# Plot the CV ranges
fig_cvs <- ggplot(df.plate_means, aes(plate_id, cv_intensity)) +
  scale_fill_npg() +
  scale_color_npg() +
  scale_y_continuous(label = percent) +
  coord_cartesian(ylim = c(0,1.5)) +
  geom_boxplot(aes(color = sample_type), show.legend = FALSE, width = 0.7, outlier.shape = 2, outlier.alpha = 0.5, outlier.size = 0.5) +
  geom_line(data = filter(df.plate_means, lipid_group == "Spike"), aes(group = lipid_id), color = "lightgray", size = 0.25) +
  geom_point(data = filter(df.plate_means, lipid_group == "Spike"), aes(fill = sample_type), shape = 21, alpha = 0.3, show.legend = FALSE) +
  geom_label(data = df.plate_means_medians, aes(y = median_cv, label = median_label), size = 2.5) +
  facet_grid(sample_type ~.) +
  labs(x = "Plate ID",
       y = "CV",
       title = "Distribution of Metabolite Precision",
       subtitle = "CV of ln intensities by plate are plotted (cv = sqrt(exp(sd^2)-1))",
       caption = "Spike metabolites (QRESS) are highlighted as circles and connected by lines; Plot limited to 150% CV for clarity, medians unaffected")

# Save the plot
my_ggsave("fig_cvs")

# Plot the CV ranges
fig_cvs <- ggplot(df.plate_norm_means, aes(plate_id, cv_intensity)) +
  scale_fill_npg() +
  scale_color_npg() +
  scale_y_continuous(label = percent) +
  coord_cartesian(ylim = c(0,1.5)) +
  geom_boxplot(aes(color = sample_type), show.legend = FALSE, width = 0.7, outlier.shape = 2, outlier.alpha = 0.5, outlier.size = 0.5) +
  geom_line(data = filter(df.plate_norm_means, lipid_group == "Spike"), aes(group = lipid_id), color = "lightgray", size = 0.25) +
  geom_point(data = filter(df.plate_norm_means, lipid_group == "Spike"), aes(fill = sample_type), shape = 21, alpha = 0.3, show.legend = FALSE) +
  geom_label(data = df.plate_norm_means_medians, aes(y = median_cv, label = median_label), size = 2.5) +
  facet_grid(sample_type ~.) +
  labs(x = "Plate ID",
       y = "CV",
       title = "Normalized Distribution of Metabolite Precision",
       subtitle = "CV of ln intensities by plate are plotted (cv = sqrt(exp(sd^2)-1))",
       caption = "Spike metabolites (QRESS) are highlighted as circles and connected by lines; Plot limited to 150% CV for clarity, medians unaffected")

# Save the plot
my_ggsave("fig_norm_cvs")

#Plot density plot 
#fig_density <- ggplot(df.data, aes(plate_id, l_intensity, fill=sample_type, color = plate_id)) +
#  geom_violin(alpha=0.5, outlier.colour="transparent", trim = FALSE) + geom_boxplot(width=0.1) +
#  facet_grid(sample_type ~.)
#my_ggsave("fig_violin")