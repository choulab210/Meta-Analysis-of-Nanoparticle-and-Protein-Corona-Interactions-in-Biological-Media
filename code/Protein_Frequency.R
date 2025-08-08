############################################################
# 1. Load required libraries --------------------------------
############################################################
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("scales", quietly = TRUE)) install.packages("scales")

library(readxl)
library(dplyr)
library(ggplot2)
library(scales)

############################################################
# 2. Read dataset -------------------------------------------
############################################################
file_path <- "Dataset-for-protein-freq-analysis-v2.xlsx"  # Adjust path if needed
df <- read_excel(file_path, sheet = 1)                    # Read first sheet

############################################################
# 3. Identify protein columns --------------------------------
#    Assumes first two columns are IDs (e.g., NP and study)
############################################################
protein_cols <- colnames(df)[-(1:2)]  # Select all columns except first two

############################################################
# 4. Compute detection rate ---------------------------------
############################################################
det_rate <- df %>%
  mutate(across(all_of(protein_cols), ~ as.integer(.x > 0))) %>%  # 1 = detected, 0 = not detected
  summarise(across(all_of(protein_cols), mean, na.rm = TRUE)) %>%
  tidyr::pivot_longer(cols = everything(),
                      names_to = "protein",
                      values_to = "detection_rate") %>%
  arrange(desc(detection_rate))

############################################################
# 5. Display top detected proteins --------------------------
############################################################
print(head(det_rate, 20))  # Show top 20 proteins by detection rate

############################################################
# 6. (Optional) Export detection rates ----------------------
############################################################
write.csv(det_rate, "/content/protein_detection_rates.csv", row.names = FALSE)
# In Colab, download with:
# from google.colab import files
# files.download('/content/protein_detection_rates.csv')

############################################################
# 7. Visualize top 20 proteins ------------------------------
############################################################
top20 <- det_rate %>%
  slice_max(detection_rate, n = 20)

ggplot(top20, aes(x = reorder(protein, detection_rate), y = detection_rate)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Top 20 Most Frequently Detected Proteins",
    x = "Protein",
    y = "Detection Proportion"
  ) +
  theme_minimal(base_size = 12)