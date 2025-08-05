library(readxl)
library(dplyr)
library(circacompare)

#cell type to analyze
cell_type <- "SCH Gaba"

#define paths for the excel files from metacycle. in raw folder. rename to the desired cell type
path_circa4 <- "./Results/2025-07-21_circa4_celltype/Raw/SCH_Gaba_cyc_analysis.xlsx"
path_sd1 <- "./Results/2025-07-21_SD1_celltype/Raw/SCH_Gaba_cyc_analysis.xlsx"

#find shared significantly cycling genes from circa4 and SD1
#filter by BH.Q < 0.05
sig_genes_circa4 <- read_excel(path_circa4, sheet = "sig_cyl_gene") %>%
  filter(meta2d_pvalue < 0.05) %>%
  pull(CycID)

sig_genes_sd1 <- read_excel(path_sd1, sheet = "sig_cyl_gene") %>%
  filter(meta2d_pvalue < 0.05) %>%
  pull(CycID)

shared_significant_genes <- intersect(sig_genes_circa4, sig_genes_sd1)

message("Found ", length(shared_significant_genes), " shared significant (BH-corrected) genes for ", cell_type, ".")
print(shared_significant_genes)

#get average measurement data from both
avg_data_circa4 <- read_excel(path_circa4, sheet = "averagesheet")
avg_data_sd1 <- read_excel(path_sd1, sheet = "averagesheet")

#filter both tables to keep only the rows for our shared genes
measurements_circa4 <- avg_data_circa4 %>% filter(gene %in% shared_significant_genes)
measurements_sd1 <- avg_data_sd1 %>% filter(gene %in% shared_significant_genes)

#get ready for circacompare
timepoints <- c(1, 5, 9, 13, 17, 21)

#empty list to store results
all_results <- list()

#circacompare results path
output_path <- "./Results/CircaCompare_Results/"
plot_output_path <- file.path(output_path, paste0(gsub(" ", "_", cell_type), "_plots"))
dir.create(plot_output_path, showWarnings = FALSE, recursive = TRUE)

#loop through each shared significant gene
for (gene_name in shared_significant_genes) {
  
  #row of measurements for this gene from both conditions
  gene_measures_circa4 <- measurements_circa4 %>% filter(gene == gene_name) %>% select(matches("ZT")) %>% as.numeric()
  gene_measures_sd1 <- measurements_sd1 %>% filter(gene == gene_name) %>% select(matches("ZT")) %>% as.numeric()
  
  #dataframe for this iteration
  x <- data.frame(
    time = rep(timepoints, 2),
    group = c(rep("circa4", length(timepoints)), rep("SD1", length(timepoints))),
    measure = c(gene_measures_circa4, gene_measures_sd1),
    alpha_threshold = 0.99 #to bypass circacompare cycling check, using all genes that pass metacycle check
  )
  
  #run circacompare and store the result in our list
  tryCatch({
    cc_result <- circacompare(x = x, col_time = "time", col_group = "group", col_outcome = "measure")
    all_results[[gene_name]] <- cc_result
    
    plot_filename <- file.path(plot_output_path, paste0(gene_name, "_circacompare_plot.png"))
    png(plot_filename, width = 8, height = 6, units = "in", res = 300)
    print(cc_result$plot)
    dev.off()
    
  }, error = function(e) {
    message("  Skipping gene '", gene_name, "' due to an error during circacompare: ", e$message)
  })
}

#combine the summary statistics from all genes into one table
summary_list <- list() 

for (gene in names(all_results)) {
  res <- all_results[[gene]]
  
  if (!is.null(res$summary) && nrow(res$summary) > 0) {
    
    wide_summary <- res$summary %>%
      tidyr::pivot_wider(names_from = parameter, values_from = value)
    
    summary_list[[gene]] <- wide_summary %>%
      mutate(gene = gene, .before = 1)
    
  } else {
    message("  Skipping gene '", gene, "' from summary file because it has no summary table.")
  }
}

summary_df <- bind_rows(summary_list)

#Define the output path and save the file
dir.create(output_path, showWarnings = FALSE)
output_filename <- file.path(output_path, paste0(gsub(" ", "_", cell_type), "_circacompare_results.csv"))
write.csv(summary_df, output_filename, row.names = FALSE)
rds_filename <- file.path(output_path, paste0(gsub(" ", "_", cell_type), "_circacompare_full_object.rds"))
saveRDS(all_results, file = rds_filename)

#will have summary csv, plot file, and R object saved.