# R script to normalize phosphoproteome based on whole proteome #
# Author: Indrani Bera dated January 9/1/2025 #
# Adjust steps as per the experimental needs such as imputation/filtering rows among others #


# Load required libraries
library(tidyverse)
library(preprocessCore)
library(data.table)

# Step 1: Read proteomics and phosphoproteomics files
proteome <- read_tsv("report.pg_matrix(10).tsv")
phosphoproteome <- read_tsv("report.phosphosites_90.tsv")
write.csv(proteome, "proteome_initial.csv", row.names = FALSE)
write.csv(phosphoproteome, "phosphoproteome_initial.csv", row.names = FALSE)

# Step 2: Check the dimensions of the dataframes
cat("Dimensions of proteome file: ", dim(proteome), "\n")
cat("Dimensions of phosphoproteome file: ", dim(phosphoproteome), "\n")

# Step 3: Select relevant columns from phosphoproteome file
phosphoproteome_selected <- phosphoproteome %>%
  select(Protein, Protein.Names, Residue, Site, Sequence, where(is.numeric))
write.csv(phosphoproteome_selected, "phosphoproteome_selected.csv", row.names = FALSE)

# Step 4: Assign categorical and numerical variables for phosphoproteome file
phosphoproteome_selected <- phosphoproteome_selected %>%
  mutate(across(c(Protein, Protein.Names, Residue, Site, Sequence), as.factor))

proteome_selected <- proteome %>%
  select(Protein.Group, Protein.Names, Genes, where(is.numeric))
write.csv(proteome_selected, "proteome_selected.csv", row.names = FALSE)
proteome_selected <- proteome_selected %>%
  mutate(across(c(Protein.Group, Protein.Names, Genes), as.factor))

# Step 5: Check percentage of missing values in numerical columns
# Discuss if missing values higher than 50%

# Subset only numerical columns
numerical_cols_phospho <- phosphoproteome_selected[sapply(phosphoproteome_selected, is.numeric)]
numerical_cols_proteo <- proteome_selected[sapply(proteome_selected, is.numeric)]

# Calculate the total percentage of zeros in numerical columns
total_zero_percentage_phos <- sum(numerical_cols_phospho == 0) / 
  (nrow(numerical_cols_phospho) * ncol(numerical_cols_phospho)) * 100
# Calculate the total percentage of zeros in numerical columns
total_na_percentage_proteo <- sum(is.na(numerical_cols_proteo)) / 
  (nrow(numerical_cols_proteo) * ncol(numerical_cols_proteo)) * 100

# Print the result
cat("Total percentage of missing values (represented as 0):", total_zero_percentage_phos, "%\n")
cat("Total percentage of missing values (represented as 0):", total_na_percentage_proteo, "%\n")


if (any(total_zero_percentage_phos > 50) | any(total_na_percentage_proteo > 50)) {
  stop("More than 50% missing values detected. Please discuss before proceeding.")
}

# Step 6: Check numerical column counts and calculate mean of technical replicates
if (ncol(phosphoproteome_selected[, sapply(phosphoproteome_selected, is.numeric)]) !=
    ncol(proteome_selected[, sapply(proteome_selected, is.numeric)])) {
  stop("The number of numerical columns in proteome and phosphoproteome files are not the same.")
}

process_technical_replicates <- function(data) {
  # Select only numerical columns
  numeric_data <- data[, sapply(data, is.numeric)]
  
  # Check if the number of columns is even
  if (ncol(numeric_data) %% 2 != 0) {
    stop("Number of numerical columns is not even. Cannot process technical replicates.")
  }
  
  # Calculate biological means by averaging every two consecutive columns (technical replicates)
  num_pairs <- ncol(numeric_data) / 2
  biological_means <- data.frame(matrix(ncol = num_pairs, nrow = nrow(numeric_data)))
  
  for (i in 1:num_pairs) {
    biological_means[, i] <- rowMeans(numeric_data[, (2 * i - 1):(2 * i)], na.rm = TRUE)
  }
  
  # Assign new column names for clarity (e.g., Bio_1, Bio_2, etc.)
  colnames(biological_means) <- paste0("Bio_", 1:num_pairs)
  
  # Return the resulting dataframe
  return(biological_means)
}

# Process the phosphoproteome and proteome data
phosphoproteome_bio <- process_technical_replicates(phosphoproteome_selected)
proteome_bio <- process_technical_replicates(proteome_selected)

write.csv(phosphoproteome_bio, "phosphoproteome_bio.csv", row.names = FALSE)
write.csv(proteome_bio, "proteome_bio.csv", row.names = FALSE)

# Step 9: Log2 transform numerical columns and keep only the transformed columns
phosphoproteome_log2 <- phosphoproteome_bio %>%
  mutate(across(where(is.numeric), log2)) %>%
  rename_with(~ paste0("log2_", .), everything())  # Add 'log2_' prefix to column names

proteome_log2 <- proteome_bio %>%
  mutate(across(where(is.numeric), log2)) %>%
  rename_with(~ paste0("log2_", .), everything())  # Add 'log2_' prefix to column names

# Save the log2-transformed dataframes to CSV
write.csv(phosphoproteome_log2, "phosphoproteome_log2.csv", row.names = FALSE)
write.csv(proteome_log2, "proteome_log2.csv", row.names = FALSE)

# Step 10: Impute missing values (NA, NaN, -Inf) with 0.05
phosphoproteome_imputed <- phosphoproteome_log2 %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.) | is.nan(.) | is.infinite(.), 0.05)))

proteome_imputed <- proteome_log2 %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.) | is.nan(.) | is.infinite(.), 0.05)))

# Save the imputed dataframes
write.csv(phosphoproteome_imputed, "phosphoproteome_imputed.csv", row.names = FALSE)
write.csv(proteome_imputed, "proteome_imputed.csv", row.names = FALSE)

# Step 11: Quantile normalization
#phosphoproteome_quantile <- as.data.frame(normalize.quantiles(as.matrix(phosphoproteome_imputed[, sapply(phosphoproteome_imputed, is.numeric)])))
#proteome_quantile <- as.data.frame(normalize.quantiles(as.matrix(proteome_imputed[, sapply(proteome_imputed, is.numeric)])))

#write.csv(phosphoproteome_quantile, "phosphoproteome_quantile.csv", row.names = FALSE)
#write.csv(proteome_quantile, "proteome_quantile.csv", row.names = FALSE)

# Step 11: Attach categorical columns back to the imputed dataframes
phosphoproteome_final <- bind_cols(
  phosphoproteome_selected %>% select(where(is.character) | where(is.factor)),
  phosphoproteome_imputed
)

proteome_final <- bind_cols(
  proteome_selected %>% select(where(is.character) | where(is.factor)),
  proteome_imputed
)

# Save the final dataframes with both categorical and numerical columns
write.csv(phosphoproteome_final, "phosphoproteome_final.csv", row.names = FALSE)
write.csv(proteome_final, "proteome_final.csv", row.names = FALSE)

phosphoproteome_final$Protein.Names <- trimws(tolower(phosphoproteome_final$Protein.Names))
proteome_final$Protein.Names <- trimws(tolower(proteome_final$Protein.Names))


# Filter phospho_data to retain only rows where Protein.Ids match Protein.ID in proteo_data
filtered_phospho_data <- phosphoproteome_final %>%
  filter(Protein.Names %in% proteome_final$Protein.Names)

filtered_proteome_data <- filtered_phospho_data %>%
  select(Protein.Names) %>%  # Keep only the key column for matching
  left_join(proteome_final, by = "Protein.Names")

# Ensure the filtered_proteome_data has the same row order as filtered_phospho_data
filtered_proteome_data <- filtered_proteome_data %>%
  slice(match(filtered_phospho_data$Protein.Names, Protein.Names))
########################################################################
# Select only numeric columns from both datasets
phospho_numeric <- filtered_phospho_data %>%
  select(where(is.numeric))

proteo_numeric <- filtered_proteome_data %>%
  select(where(is.numeric))

# Check if the row counts match
if (nrow(phospho_numeric) != nrow(proteo_numeric)) {
  stop("Row counts between phospho_numeric and proteo_numeric do not match!")
}

# Initialize normalized phospho data as a numeric matrix
normalized_phospho_data <- matrix(NA, nrow = nrow(phospho_numeric), ncol = ncol(phospho_numeric))
colnames(normalized_phospho_data) <- colnames(phospho_numeric)

# Perform linear modeling row by row
for (i in 1:nrow(phospho_numeric)) {
  # Extract current row of numeric phosphoproteomics and proteomics data
  phospho_row <- as.numeric(phospho_numeric[i, ])
  proteo_row <- as.numeric(proteo_numeric[i, ])
  
  # Combine the rows into a dataframe for regression
  data <- data.frame(phospho = phospho_row, proteo = proteo_row)
  
  # Fit a linear model
  model <- lm(phospho ~ proteo, data = data)
  normalized_phospho_data[i, ] <- residuals(model)
}

# Convert the matrix back to a dataframe
normalized_phospho_data <- as.data.frame(normalized_phospho_data)
# Rename columns starting with "log2_Bio_"
colnames(normalized_phospho_data) <- gsub("^log2_Bio_", "Residual_Bio_", colnames(normalized_phospho_data))

# Identify categorical columns in filtered_phospho_data
categorical_columns <- filtered_phospho_data[, sapply(filtered_phospho_data, is.character) | 
                                               sapply(filtered_phospho_data, is.factor)]

combined_phospho_data <- cbind(categorical_columns, normalized_phospho_data)
combined_phospho_data <- as.data.frame(combined_phospho_data)
numeric <- combined_phospho_data[, sapply(combined_phospho_data, is.numeric)]

# Create a logical vector for filtering rows
filter_condition <- apply(numeric, 1, function(row) any(row >= 2 | row <= -2))

# Filter rows based on the condition
filtered_combined_phospho_data <- combined_phospho_data[filter_condition, ]
write.csv(filtered_combined_phospho_data, "filtered_combined_phospho_data.csv", row.names = FALSE)


filtered_phospho_data$Protein.Names <- rownames(filtered_phospho_data)
filtered_combined_phospho_data$Protein.Names <- rownames(filtered_combined_phospho_data)

# Merge the dataframes, keeping Protein.Names only once
filtered_phospho_data_norm <- merge(
  filtered_phospho_data,
  filtered_combined_phospho_data,
  by = "Protein.Names",
  suffixes = c("", "_combined")  # Prevent .x and .y suffixes
)

filtered_phospho_data_norm <- filtered_phospho_data_norm[, !grepl("_combined$", names(filtered_phospho_data_norm))]
# Save the filtered dataframe
write.csv(filtered_phospho_data_norm, "filtered_phospho_data_norm.csv", row.names = FALSE)


##########################################################################

#### perform differential expression based on group information ####
# Load required libraries

remove.packages("EnhancedVolcano")
install.packages("BiocManager")
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
library(limma)
library(ggplot2)
library(ggrepel)
library(preprocessCore)
# Read the filtered_phospho_data_norm dataframe into phospho_data
phospho_data <- filtered_phospho_data_norm
# Merge the columns Protein, Residue, and Site into a new column Protein_Psite
phospho_data$Protein_Psite <- paste(phospho_data$Protein, phospho_data$Residue, phospho_data$Site, sep = "_")

# Select the first 8 numerical columns beginning with "log2"
numeric_data <- phospho_data[, grep("^log2", names(phospho_data))[1:8]]
rownames(numeric_data) <- phospho_data$Protein_Psite  # Replace 'Protein' with the correct column name

# Select non-numeric columns (categorical columns)
categorical_columns <- phospho_data[, !sapply(phospho_data, is.numeric)]

# Quantile normalize the numeric data
numeric_data_qn <- normalize.quantiles(as.matrix(numeric_data))
# Replace row names in the normalized data with Protein_Psite
rownames(numeric_data_qn) <- rownames(numeric_data)
# Convert the normalized matrix back to a dataframe
numeric_data_qn <- as.data.frame(numeric_data_qn)
# Restore the original column names
colnames(numeric_data_qn) <- colnames(numeric_data)


# Create a group factor (Light vs Dark)
group <- factor(c(rep("Light", 4), rep("Dark", 4)))


# Ensure numeric data is properly prepared
# 'numeric_data_qn' should already be quantile normalized

# Create a design matrix for the group factor
design <- model.matrix(~0 + group)  # "~0" ensures no intercept is included
colnames(design) <- levels(group)   # Rename columns to "Light" and "Dark"

# Fit the linear model
fit <- lmFit(numeric_data_qn, design)

# Create contrast matrix for "Light vs Dark"
contrast_matrix <- makeContrasts(Light_vs_Dark = Light - Dark, levels = design)

# Apply the contrast matrix to the fitted model
fit2 <- contrasts.fit(fit, contrast_matrix)

# Compute moderated t-statistics, F-statistics, and log-odds
fit2 <- eBayes(fit2)

# Extract the top differentially expressed proteins
results_proteome <- topTable(fit2, coef = "Light_vs_Dark", number = Inf)

#write.csv(top_table, "differentially_expressed_proteins.csv", row.names = TRUE)

# Define significance and condition
results_proteome$Significance <- "NO"
# if log2Foldchange > 1.5 and pvalue < 0.05, set as "UP" 
results_proteome$Significance[results_proteome$logFC > 1.5 & results_proteome$P.Value < 0.05] <- "Up-regulated"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
results_proteome$Significance[results_proteome$logFC < -1.5 & results_proteome$P.Value < 0.05] <- "Down-regulated"

results_proteome$delabel <- NA
results_proteome$delabel[results_proteome$Significance != "NO"] <- rownames(results_proteome)[results_proteome$Significance != "NO"]


#results_proteome$Significance <- ifelse(results_proteome$P.Value < 0.05 & abs(results_proteome$logFC) > 1.5, "Significant", "Not Significant")
#results_proteome$Condition <- ifelse(results_proteome$logFC > 0, "Light", "Dark")
#results_proteome$color <- ifelse(results_proteome$Significance == "Not Significant", "black",
 #                                ifelse(results_proteome$Condition == "Light", "green", "red"))

# Create the volcano plot
volcano_plot <- ggplot(results_proteome, aes(x = logFC, y = -log10(P.Value), col=Significance, label=delabel)) +
  geom_point(size = 2.6) +  # Adjust size and transparency
  scale_color_manual(values=c("red", "black", "green")) +
  labs(title = "Volcano Plot Light vs Dark") + 
  xlab(bquote(Log[2]*FC))+ 
  ylab(bquote(-Log[10]*P.value)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_text_repel(size = 3.8) +
  theme_classic(base_size = 8) + 
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),  # Larger title, centered
    axis.title = element_text(size = 18, face = "bold"),  # Larger axis labels
    axis.text = element_text(size = 14, face = "bold")   # Larger tick labels
  )
 # Increase plot title size

# Save the volcano plot to a high-resolution JPEG file
jpeg(file = "volcano_plot_light_vs_dark2.jpeg", width = 9, height = 10, units = "in", res = 600)
print(volcano_plot)
dev.off()

# Save results with differential expression and volcano plot labels
write.csv(results_proteome, file = "differential_expression_results_light_vs_dark.csv", row.names = FALSE)

################################################

# Load the pheatmap package
library(pheatmap)

# Save as JPEG
jpeg("heatmap_numeric_data.jpeg", width = 9, height = 10, units = "in", res = 900)  # Adjust resolution and size
p<-pheatmap(
  numeric_data, 
  color = colorRampPalette(c("blue", "white", "red"))(50), # Blue-to-red gradient
  scale = "row",       # Normalize rows to emphasize patterns
  cluster_rows = TRUE, # Cluster rows
  cluster_cols = TRUE, # Cluster columns
  show_rownames = TRUE,  # Show row names
  show_colnames = TRUE,  # Show column names
  fontsize_row = 1,     # Adjust font size for row names
  fontsize_col = 10,    # Adjust font size for column names
  #main = "Heatmap" # Add a title
)
dev.off()  # Close the graphics device
#########################################################
#########################################################


# Load necessary libraries
install.packages("FactoMineR")
library("FactoMineR")
# Install factoextra package if not installed
#install.packages("factoextra")

# Load the factoextra package
library(factoextra)

groups <- factor(c(rep("Light", 4), rep("Dark", 4)))
groups
# Assuming 'numeric_data' is your data with the first 8 numerical columns
numeric_data_8 <- numeric_data[, 1:8]
#numeric_data_8 <- t(numeric_data_8)
data_normalized <- scale(numeric_data_8)
head(data_normalized)

data.pca <- prcomp(data_normalized)
summary(data.pca)
data.pca$loadings[, 1:2]
fviz_eig(data.pca, addlabels = TRUE)
fviz_pca_var(data.pca, col.var = "black")
fviz_pca_var(data.pca, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)
fviz_cos2(data.pca, choice = "var", axes = 1:2)

library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)
str(data.pca$Sample)


# Create PCA plot
# Plot only the scores (samples)
g <- ggbiplot(data.pca,
              obs.scale = 1,
              var.scale = 1,
              groups = groups,
              ellipse = TRUE,
              circle = TRUE,
              ellipse.prob = 0.68,
              var.axes = FALSE)  # Disable variable arrows

# Customize the plot
g <- g + scale_color_discrete(name = "Group") +
  theme(legend.direction = 'horizontal',
        legend.position = 'top')
print(g)

