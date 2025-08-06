library(readxl)
library(dplyr)
library(circacompare)

#cell type to analyze
path_to_R = '/media/volume/volume_spatial/hugo/R/'

list_celltype = c('SMC', 'Sncg_Gaba', 'STR_Gaba', 'PAL_STR_Gaba_Chol',
        'L2_3_PIR_ENTl_Glut', 'MEA_Glut', 'OB_STR_CTX_IMN', 'L5_NP_CTX_Glut',
        'COAa_PAA_MEA_Glut', 'Vip_Gaba', 'Lamp5_Gaba', 'DG_Glut', 'ABC',
       'STR_PAL_Gaba', 'PVH_Glut', 'L6b_CTX_Glut', 'AHN_Glut', 'SCH_Gaba',
       'Sst_Gaba', 'VLMC', 'L5_ET_CTX_Glut', 'Pericyte', 'LHA_Glut',
       'Ependymal', 'Pvalb_Gaba', 'Microglia', 'OPC', 'L2_3_IT_PIR_ENTl_Glut',
       'STR_D1_Gaba', 'L6_CT_CTX_Glut', 'Endothelial', 'Oligodendrocyte',
       'L4_5_IT_CTX_Glut', 'STR_D2_Gaba', 'L6_IT_CTX_Glut', 'Astro_TE')

circacompare_loop = function(path_to_R = '/media/volume/volume_spatial/hugo/R/',
                            cell_type = cell_type){

#define paths for the excel files from metacycle. in raw folder. rename to the desired cell type
path_circa4 <- paste0(path_to_R,"Results/2025-07-21_circa4_celltype/Raw/",cell_type,"_cyc_analysis.xlsx")
path_sd1 <- paste0(path_to_R,"Results/2025-07-21_SD1_celltype/Raw/",cell_type,"_cyc_analysis.xlsx")

#find shared significantly cycling genes from circa4 and SD1
#filter by pval< 0.05
sig_genes_circa4 <- read_excel(path_circa4, sheet = "sig_cyl_gene") %>%
  filter(meta2d_pvalue < 0.05) %>%
  pull(CycID)

sig_genes_sd1 <- read_excel(path_sd1, sheet = "sig_cyl_gene") %>%
  filter(meta2d_pvalue < 0.05) %>%
  pull(CycID)

shared_significant_genes <- intersect(sig_genes_circa4, sig_genes_sd1)

message("Found ", length(shared_significant_genes), " shared significant genes for ", cell_type, ".")
print(shared_significant_genes)

#get average measurement data from both
avg_data_circa4 <- read_excel(path_circa4, sheet = "averagesheet")
avg_data_sd1 <- read_excel(path_sd1, sheet = "averagesheet")

#filter both tables to keep only the rows for our shared genes
measurements_circa4 <- avg_data_circa4 %>% filter(gene %in% shared_significant_genes)
measurements_sd1 <- avg_data_sd1 %>% filter(gene %in% shared_significant_genes)

#get ready for circacompare
timepoints <- c(1,1,1,5,5,5,9,9,9,13,13,13,17,17,17,21,21,21)

#empty list to store results
all_results <- list()

#circacompare results path
output_path <- paste0(path_to_R,"./Results/CircaCompare_Results/")
plot_output_path <- file.path(output_path, paste0(gsub(" ", "_", cell_type), "_plots")) 
dir.create(plot_output_path, showWarnings = FALSE, recursive = TRUE)

Sys.time()
#loop through each shared significant gene
for (gene_name in shared_significant_genes) {
  
  #row of measurements for this gene from both conditions
  gene_measures_circa4 <- measurements_circa4 %>% filter(gene == gene_name) %>% select(matches("ZT")) %>% as.numeric()
  gene_measures_sd1 <- measurements_sd1 %>% filter(gene == gene_name) %>% select(matches("ZT")) %>% as.numeric()
  
  #dataframe for this iteration
  x <- data.frame(
    time = rep(timepoints, 2),
    group = c(rep("circa4", length(timepoints)), rep("SD1", length(timepoints))),
    measure = c(gene_measures_circa4, gene_measures_sd1)
  )

  #run circacompare and store the result in our list
  tryCatch({
    cc_result <- circacompare(x = x, col_time = "time", col_group = "group", col_outcome = "measure", alpha_threshold = 0.99)
    all_results[[gene_name]] <- cc_result

    plot_filename <- file.path(plot_output_path, paste0(gene_name, "_circacompare_plot.png"))
    png(plot_filename, width = 8, height = 6, units = "in", res = 300)
    print(cc_result$plot)
    dev.off()
    
  }, error = function(e) {
    message("  Skipping gene '", gene_name, "' due to an error during circacompare: ", e$message)
  })
}
Sys.time()

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
}
#will have summary csv, plot file, and R object saved.