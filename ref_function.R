### Preprocessing_data
library(dplyr)
library(arrow)
library(MetaCycle)
library(readr)
library(knitr)

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
    valid_cells <- c('ABC','AHN Glut','Astro TE','CEA Gaba','CLA EPd CTX Glut','COAa PAA MEA Glut','Choroid',
'Endothelial','Ependymal','L2 3 IT PIR ENTl Glut','L4 5 IT CTX Glut','L5 ET CTX Glut',
'L5 NP CTX Glut','L6 CT CTX Glut','L6 IT CTX Glut','L6b CTX Glut','LHA Glut','Lamp5 Gaba',
'MEA Glut','Microglia','OB STR CTX IMN','OPC','Oligodendrocyte','PAL STR Gaba Chol','PVH Glut','PVT Glut','Pericyte','Pvalb Gaba','SCH Gaba','SMC','STR D1 Gaba','STR D2 Gaba','STR Gaba','STR PAL Gaba','Sncg Gaba','Sst Gaba','VLMC','Vip Gaba'
)}
  else if (run_name == 'SD1'){
    valid_cells <- c("ABC","AHN Glut","Astro TE","CEA Gaba","CLA EPd CTX Glut","COAa PAA MEA Glut",
"Choroid","Endothelial","Ependymal","L2 3 IT PIR ENTl Glut","L4 5 IT CTX Glut",
"L5 ET CTX Glut","L5 NP CTX Glut","L6 CT CTX Glut","L6 IT CTX Glut","L6b CTX Glut",
"LHA Glut","Lamp5 Gaba","MEA Glut","MPO Glut","Microglia","OB STR CTX IMN","OPC",
"Oligodendrocyte","PAL STR Gaba Chol","PVH Glut","Pericyte","Pvalb Gaba","SCH Gaba",
"SMC","SPA Glut","STR D1 Gaba","STR D2 Gaba","STR Gaba","STR PAL Gaba","Sncg Gaba",
"Sst Gaba","VLMC","Vip Gaba")}
  print('Valid celltypes names: Loaded')
  return (valid_cells)
}

# celltype = GetValidCells(run) ### Usage

### Valid Regions
GetValidRegions = function(run_name){
  if (run_name == 'circa4'){
    valid_regions = c("AHN", "AMY", "CTX", "LHA", "PVH", "PVT", "SCH", "STR", "WM")
  }
  else if(run_name == 'SD1'){
    valid_regions = c("AHN", "AMY", "CTX", "LHA", "PVH", "SCH", "STR")
  }
  print('Valid region names: Loaded')
  return (valid_regions)
}

### Valid Cell Classes
GetValidCellClasses <- function(run_name){
  if (run_name == 'circa4'){
    valid_classes <- c("Glial", "Neuronal", "Vascular","Ependymal") 
  } else if (run_name == 'SD1'){
    valid_classes <- c("Glial", "Neuronal", "Vascular","Ependymal")
  }
  print('Valid cell class names: Loaded')
  return (valid_classes)
}

### Valid Neurotransmitters
GetValidNeurotransmitters <- function(run_name){
  if (run_name == 'circa4'){
    valid_neurotransmitters <- c("Glutamate", "Gaba", "Acetylcholine") 
  } else if (run_name == 'SD1'){
    valid_neurotransmitters <- c("Glutamate", "Gaba", "Acetylcholine")
  }
  print('Valid neurotransmitter names: Loaded')
  return (valid_neurotransmitters)
}

### Clock genes
#clockgenelist = c("Arntl","Clock", "Cry1","Cry2", "Npas2","Nr1d1", "Per1",  "Per2", "Per3", "Rora", "Rorb","Rorc")


MetaCycleAnalysis <- function(data, condition, run_name, path_to_save, date, 
                              use_gene_lists = FALSE, gene_list_path = NULL) {
  
  col_map <- list(
    celltype = "cell_type_final",
    region = "region_automap_name",
    class = "cell_class",
    neurotransmitter = "neurotransmitter"
  )

  if (length(condition) == 2 && "region" %in% condition) {
    #will identify the condition other than region below
    primary_condition <- condition[condition != "region"]
    
    message(paste("===== Starting COMBINED analysis for", primary_condition, "+ Region ====="))
    
    primary_col <- col_map[[primary_condition]]
    region_col  <- col_map[["region"]]
    
    #check that all cols are correctly named
    if (is.null(primary_col) || !primary_col %in% names(data) || !region_col %in% names(data)) {
      stop("One or both required columns not found in data for the combined analysis.")
    }
    
    #to select right function for getvalid___
    get_valid_items_func <- switch(primary_condition,
                                   "celltype" = GetValidCells,
                                   "class" = GetValidCellClasses,
                                   "neurotransmitter" = GetValidNeurotransmitters,
                                   stop("Invalid primary condition for combined analysis.")
    )
    valid_primary_items <- get_valid_items_func(run_name)
    
    analysis_data <- data %>%
      filter(.data[[primary_col]] %in% valid_primary_items, .data[[region_col]] %in% GetValidRegions(run_name)) %>%
      mutate(combined_group = paste(.data[[primary_col]], .data[[region_col]], sep = "_in_"))
    
    items_to_process <- unique(analysis_data$combined_group)
    column_to_filter <- "combined_group"
    analysis_by <- "combined_group"
    data_to_use <- analysis_data
    
  } else if (length(condition) == 1 && condition %in% names(col_map)) {
    analysis_by <- condition[1]
    message(paste0("\n\n===== Starting SINGLE analysis by: ", toupper(analysis_by), " ====="))
    
    get_valid_items_func <- switch(analysis_by,
                                   "celltype" = GetValidCells,
                                   "region" = GetValidRegions,
                                   "class" = GetValidCellClasses,
                                   "neurotransmitter" = GetValidNeurotransmitters,
                                   stop("Invalid analysis type specified.")
    )
    items_to_process <- get_valid_items_func(run_name)
    column_to_filter <- col_map[[analysis_by]]
    data_to_use <- data
    
  } else {
    stop("Invalid 'condition' variable. Please use a single valid condition or a combination of one condition + 'region'.")
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
        genes_raw <- read_csv(gene_file, col_names = FALSE, col_types = c('c','c'))
        genes_raw = na.omit(genes_raw)
        genes_to_analyze = genes_raw$X2
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
        meta_cols <- names(data1)[!names(data1) %in% matching_genes]
        cols_to_keep <- unique(c(meta_cols, matching_genes))
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
                    maxper = 28, cycMethod = method_to_use, outputFile = FALSE,
                    parallelize = FALSE, nCores = 40)
        
        if("meta2d_pvalue" %in% names(d$meta)) {
          d$meta$meta2d_pvalue <- as.numeric(d$meta$meta2d_pvalue)
          sigcyclegene <- d$meta %>% filter(meta2d_pvalue < 0.05)
          d$meta$meta2d_BH.Q <- as.numeric(d$meta$meta2d_BH.Q)
          sigcyclegeneBH <- d$meta %>% filter(meta2d_BH.Q < 0.05)
        } else {
          d$meta$LS_pvalue <- as.numeric(d$meta$LS_pvalue)
          sigcyclegene <- d$meta %>% filter(LS_pvalue < 0.05)
          d$meta$LS_BH.Q <- as.numeric(d$meta$LS_BH.Q)
          sigcyclegeneBH <- d$meta %>% filter(LS_BH.Q < 0.05)
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
    
    summary_df_BH <- data.frame(group = names(siglistBH), cycling_gene_count = sapply(siglistBH, nrow))
    names(summary_df_BH)[1] <- analysis_by
    write.csv(summary_df_BH, paste0(path_to_save, "/Summary/", date, "_", run_name, "_cycling_gene_per_group_BH.csv"), row.names = FALSE)
    
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
