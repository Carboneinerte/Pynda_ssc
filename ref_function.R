### Preprocessing_data
library(dplyr)
library(arrow)
library(MetaCycle)

PreProcessingData = function(run_name, circascore='all'){
  print("Start importing data")
  path_to_data = paste0('/media/volume/volume_spatial/hugo/R/data/',run_name,'_norm_combined.parquet')
  data=read_parquet(path_to_data)
  
  data = CircaFilter(data, run_name, circascore)
  
  return(data)
}

### Circascore filter
CircaFilter = function(data, run_name,circascore){
  
  if (circascore == 'zero'){
    data = dplyr::filter(data, data$circascore == 0)
  } else if (circascore == 'non-zero'){
    data = dplyr::filter(data, data$circascore != 0)
  } else if (circascore == 'all') {print('no circascore')}
  print('Data Preprocessing done')
  
  return (data)
}

# data = CircaFilter(run, circascore) ### Usage


### Valid cells
GetValidCells <- function(run_name){
  if (run_name == 'circa4'){
    valid_cells <- c("ABC","AHN Glut","Astro TE","CLA EPd CTX Glut","COAa PAA MEA Glut","Choroid","DG Glut",                     "Endothelial","Ependymal",
                     "L2 3 IT PIR ENTl Glut","L2 3 PIR ENTl Glut","L4 5 IT CTX Glut","L5 ET CTX Glut","L5 NP CTX Glut",
                     "L6 CT CTX Glut","L6 IT CTX Glut","L6b CTX Glut","LHA Glut","Lamp5 Gaba","MEA Glut","Microglia","OB STR CTX IMN","OPC",
                     "Oligodendrocyte","PAL STR Gaba Chol","PVH Glut","PVT Glut","Pericyte","Pvalb Gaba","SCH Gaba",                     "SMC",                     "SPA Glut",
                     "STR D1 Gaba","STR D2 Gaba","STR Gaba","STR PAL Gaba","Sncg Gaba",
                     "Sst Gaba","VLMC","Vip Gaba"
    )}
  else if (run_name == 'SD1'){
    valid_cells <- c('ABC','Astro TE','Choroid','CLA EPd CTX Glut','COAa PAA MEA Glut','DG Glut','Endothelial',
                     'Ependymal','HY Glut','L2 3 IT PIR ENTl Glut','L2 3 PIR ENTl Glut','L4 5 IT CTX Glut',
                     'L5 ET CTX Glut','L5 NP CTX Glut','L6 CT CTX Glut','L6 IT CTX Glut','L6b CTX Glut',
                     'Lamp5 Gaba','Microglia','OB STR CTX IMN','Oligodendrocyte','OPC','PAL STR Gaba Chol',
                     'Pericyte','Pvalb Gaba','SCH Gaba','SMC','Sncg Gaba','Sst Gaba','STR D1 Gaba','STR D2 Gaba',
                     'STR Gaba','STR PAL Gaba','Vip Gaba','VLMC'
    )}
  print('Valid celltypes names: Loaded')
  return (valid_cells)
}

# celltype = GetValidCells(run) ### Usage

### Valid Regions
GetValidRegions = function(run_name){
  if (run_name == 'circa4'){
    valid_regions = c("CTX","Ependymal","HIPP","HY","PVT","SCH","STR","VLMC","WM")
  }
  else if(run_name == 'SD1'){
    valid_regions = c("CTX","Ependymal","HIPP","HY","SCH","STR","VLMC","WM")
  }
  print('Valid region names: Loaded')
  return (valid_regions)
}

### Clock genes
clockgenelist = c("Arntl","Clock", "Cry1","Cry2", "Npas2","Nr1d1", "Per1",  "Per2", "Per3", "Rora", "Rorb","Rorc")


MetaCycleAnalysis <- function(data, condition, run_name, path_to_save, date, 
                              use_gene_lists = FALSE, gene_list_path = NULL) {
  
  celltype_col <- "cell_type_final"
  region_col   <- "region_automap_name"

  if (use_gene_lists && is.null(gene_list_path)) {
    stop("If 'use_gene_lists' is TRUE, you must provide a path in 'gene_list_path'.")
  }
  
  if (length(condition) == 2 && all(c("celltype", "region") %in% condition)) {
    message("===== Starting COMBINED analysis for Cell Type + Region =====")
    if (!celltype_col %in% names(data) || !region_col %in% names(data)) {
      stop(paste("Error: One or both required columns not found in data:", celltype_col, ",", region_col))
    }
    analysis_data <- data %>%
      filter(.data[[celltype_col]] %in% GetValidCells(run_name), .data[[region_col]] %in% GetValidRegions(run_name)) %>%
      mutate(combined_group = paste(.data[[celltype_col]], .data[[region_col]], sep = "_in_"))
    items_to_process <- unique(analysis_data$combined_group)
    column_to_filter <- "combined_group"
    analysis_by <- "combined_group"
    data_to_use <- analysis_data
    message(paste("Found", length(items_to_process), "unique cell type/region combinations to analyze."))
  } else if (length(condition) == 1 && condition %in% c("celltype", "region")) {
    analysis_by <- condition[1]
    message(paste0("\n\n===== Starting SINGLE analysis by: ", toupper(analysis_by), " ====="))
    if (analysis_by == "celltype") {
      items_to_process <- GetValidCells(run_name)
      column_to_filter <- celltype_col
    } else {
      items_to_process <- GetValidRegions(run_name)
      column_to_filter <- region_col
    }
    data_to_use <- data
  } else {
    stop("Invalid 'condition' variable. Please use c('celltype'), c('region'), or c('celltype', 'region').")
  }
  
  siglist <- list()
  siglistBH <- list()
  
  for (item in items_to_process) {
    message("Processing group: ", item)
    
    data1 <- data_to_use %>% filter(.data[[column_to_filter]] == item)
    
    genes_to_analyze <- NULL
    if (use_gene_lists) {
      gene_file <- "" 
      
      if(analysis_by == "combined_group") {
        celltype_for_list <- sub("_in_.*", "", item)
        region_for_list <- sub(".*_in_", "", item)
        gene_file <- file.path(gene_list_path, region_for_list, paste0(celltype_for_list, ".csv"))
      } else {
        gene_file <- file.path(gene_list_path, paste0(item, ".csv"))
      }
      
      if (file.exists(gene_file)) {
        message("  Found gene list. Reading and cleaning...")
        genes_raw <- read_csv(gene_file)
        
        genes_to_analyze = genes_raw$`0`
        # genes_to_analyze <- gsub('\\"', '', genes_raw)
        # genes_to_analyze <- gsub('^\\d+\\s+', '', genes_to_analyze)
        # genes_to_analyze <- trimws(genes_to_analyze)
        # genes_to_analyze <- genes_to_analyze[genes_to_analyze != "x" & nchar(genes_to_analyze) > 0]
        
        matching_genes <- intersect(genes_to_analyze, names(data1))
        
        if(length(matching_genes) == 0){
          message("  WARNING: No genes from the list matched the data columns for '", item, "'. Skipping.")
          next
        }
        
        message("  Found ", length(matching_genes), " matching genes to analyze.")
        meta_cols <- names(data1)[!grepl("^[A-Z]", names(data1))] 
        cols_to_keep <- c(meta_cols, matching_genes)
        data1 <- data1[, cols_to_keep]
        
      } else {
        message("  WARNING: Gene list file not found for '", item, "'. Path checked: ", gene_file)
        next 
      }
    }
    
    if (nrow(data1) == 0) { message("  Skipping '", item, "' as it contains no data."); next }
    
    datasamplezt <- split.data.frame(data1, data1$sample)
    dataaveragelist <- list(); datasdlist <- list()
    
    for (ZT in names(datasamplezt)) {
      df <- datasamplezt[[ZT]]
      group_labels <- rep(1:3, length.out = nrow(df)); set.seed(123); df$group <- sample(group_labels)
      cols_for_summary <- if(!is.null(genes_to_analyze)) intersect(genes_to_analyze, names(df)) else names(df)[1:5006]
      if(length(cols_for_summary) > 0) {
        dfaverage <- df %>% group_by(group) %>% summarize(across(all_of(cols_for_summary), mean, .names = "{.col}"))
        dfsd <- df %>% group_by(group) %>% summarize(across(all_of(cols_for_summary), sd, .names = "{.col}"))
        dataaveragelist[[ZT]] <- dfaverage; datasdlist[[ZT]] <- dfsd
      }
    }
    
    dftoexport <- do.call(rbind, dataaveragelist)
    tdftoexport <- as.data.frame(t(dftoexport))
    if("group" %in% rownames(tdftoexport)) { tdftoexport <- tdftoexport[-which(rownames(tdftoexport) == "group"), ] }
    
    if (!is.data.frame(tdftoexport) && !is.matrix(tdftoexport)) {
      message("  Skipping '", item, "' because its data is 1-dimensional (likely from a single cell/sample)."); next
    }
    
    tdftoexport <- tdftoexport[rowSums(tdftoexport == 0, na.rm=TRUE) <= ncol(tdftoexport) / 3, ]
    
    safe_item_name <- gsub("[ /]", "_", item)
    outfile_prefix <- paste0(path_to_save, "/Raw/", safe_item_name)
    write.csv(tdftoexport, paste0(outfile_prefix, "_data.csv"))
    
    tryCatch({
      if (nrow(tdftoexport) < 1 || ncol(tdftoexport) < 2) {
        message("  Skipping MetaCycle for '", item, "' due to insufficient data after filtering.")
      } else {
        timepointstolookat <- as.numeric(gsub("ZT", "", gsub(".*(ZT\\d+).*", "\\1", rownames(dftoexport))))
        if(any(is.na(timepointstolookat))){ message("  Skipping '", item, "' because it contains non-ZT samples."); next }
        method_to_use <- if (length(timepointstolookat) >= 18) c("LS", "JTK") else "LS"
        d <- meta2d(infile = paste0(outfile_prefix, "_data.csv"), outdir = path_to_save,
                    filestyle = "csv", timepoints = timepointstolookat, minper = 20, 
                    maxper = 28, cycMethod = method_to_use, outputFile = FALSE)
        
        if("meta2d_pvalue" %in% names(d$meta)) {
          d$meta$meta2d_pvalue <- as.numeric(d$meta$meta2d_pvalue)
          sigcyclegene <- d$meta %>% filter(meta2d_pvalue < 0.05)
          d$meta$meta2d_BH.Q <- as.numeric(d$meta$meta2d_BH.Q)
          sigcyclegeneBH <- d$meta %>% filter(meta2d_BH.Q < 0.05)
        } else {
          d$meta$LS_pvalue <- as.numeric(d$meta$LS_pvalue)
          sigcyclegene <- d$meta %>% filter(LS_pvalue < 0.05)
        }

        dftosd <- as.data.frame(t(do.call(rbind, datasdlist)))[-1, ]; dftosd$gene <- rownames(dftosd) 
        tdftoexport$gene <- rownames(tdftoexport); d[["averagesheet"]] <- tdftoexport; d[["sdsheet"]] <- dftosd
        d[["sig_cyl_gene"]] <- sigcyclegene
        writexl::write_xlsx(Filter(Negate(is.null), d), paste0(outfile_prefix, "_cyc_analysis.xlsx"))
        if(nrow(sigcyclegene) > 0) { siglist[[item]] <- sigcyclegene }
        if(nrow(sigcyclegeneBH) > 0) { siglistBH[[item]] <- sigcyclegeneBH }
      }
    }, error = function(e) { message("  Skipping '", item, "' due to an error: ", e$message) })
  }
  
  if(length(siglist) > 0) {
    message("\nAnalysis complete. Writing summary files...")
    writexl::write_xlsx(siglist, paste0(path_to_save, "/Summary/", date, "_", run_name, "_cyc_siggene_analysis.xlsx"))
    writexl::write_xlsx(siglistBH, paste0(path_to_save, "/Summary/", date, "_", run_name, "_cyc_siggeneBH_analysis.xlsx"))
    
    summary_df <- data.frame(group = names(siglist), cycling_gene_count = sapply(siglist, nrow))
    names(summary_df)[1] <- analysis_by
    write.csv(summary_df, paste0(path_to_save, "/Summary/", date, "_", run_name, "_cycling_gene_per_group.csv"), row.names = FALSE)
    
    summary_df_BH <- data.frame(group = names(siglistBH), cycling_gene_count = sapply(siglist, nrow))
    names(summary_df_BH)[1] <- analysis_by
    write.csv(summary_df_BH, paste0(path_to_save, "/Summary/", date, "_", run_name, "_cycling_gene_per_group.csv"), row.names = FALSE)
    
    all_cycling_genes <- do.call(rbind, lapply(names(siglist), function(name) data.frame(gene=siglist[[name]]$CycID, item_name=name)))
    if (!is.null(all_cycling_genes) && nrow(all_cycling_genes) > 0) {
      gene_item_counts <- all_cycling_genes %>% group_by(gene) %>% summarize(group_count = n(), groups = paste(item_name, collapse = " | ")) %>% arrange(desc(group_count))
      names(gene_item_counts) <- c("gene", paste0(analysis_by, "_count"), paste0(analysis_by, "s"))
      write.csv(gene_item_counts, paste0(path_to_save, "/Summary/", date, "_", run_name, "_group_per_cycling_gene.csv"), row.names = FALSE)
    }
    message("...Summary files written successfully.\n")
  } else { message("No significant cycling genes found.") }
  
  message("===== All analyses complete. =====")
}
