# run_ = 'circa4'
# circascore = 'all' 

### Preprocessing_data
library(dplyr)
library(arrow)
library(MetaCycle)

PreProcessingData = function(run_name, circascore='all'){
  print("Start importing data")
  path_to_data = paste0('./data/',run_name,'_norm_combined.parquet')
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
    valid_cells <- c("ABC","Astro TE","Choroid","CLA EPd CTX Glut","COAa PAA MEA Glut","DG Glut","Endothelial",
                     "Ependymal","HY Glut","L2 3 IT PIR ENTl Glut","L2 3 PIR ENTl Glut","L4 5 IT CTX Glut",
                     "L5 ET CTX Glut","L5 NP CTX Glut","L6 CT CTX Glut","L6 IT CTX Glut","L6b CTX Glut",
                     "Lamp5 Gaba","Microglia","OB STR CTX IMN","Oligodendrocyte","OPC","PAL STR Gaba Chol",
                     "Pericyte","Pvalb Gaba","PVT PT Glut","SCH Gaba","SMC","Sncg Gaba","Sst Gaba","STR D1 Gaba",
                     "STR D2 Gaba","STR Gaba","STR PAL Gaba","Vip Gaba","VLMC"
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
clockgenelist = c("Arntl","Clock", "Cry1","Cry2", "Npas2","Nr1d1", "Per1",  "Per2", "Per3", "Rora", "Rorb",
                  "Rorc", "Dbp", "Nfil3", "Hlf", "Ciart")


MetaCycleAnalysis <- function(data, condition, run_name, path_to_save, date) {
  
  celltype_col <- "cell_type_final"
  region_col   <- "region_automap_name"

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
    } else { # region
      items_to_process <- GetValidRegions(run_name)
      column_to_filter <- region_col
    }
    data_to_use <- data
    
  } else {
    stop("Invalid 'condition' variable. Please use c('celltype'), c('region'), or c('celltype', 'region').")
  }
  
  siglist <- list()
  
  for (item in items_to_process) {
    message("Processing group: ", item)
    
    data1 <- data_to_use %>% filter(.data[[column_to_filter]] == item)
    
    if (nrow(data1) == 0) {
      message("  Skipping '", item, "' as it contains no data.")
      next
    }
    
    datasamplezt <- split.data.frame(data1, data1$sample)
    dataaveragelist <- list()
    datasdlist <- list()
    
    for (ZT in names(datasamplezt)) {
      df <- datasamplezt[[ZT]]
      group_labels <- rep(1:3, length.out = nrow(df))
      set.seed(123)
      df$group <- sample(group_labels)
      dfaverage <- df %>% group_by(group) %>% summarize(across(1:5006, mean, .names = "{.col}"))
      dfsd <- df %>% group_by(group) %>% summarize(across(1:5006, sd, .names = "{.col}"))
      dataaveragelist[[ZT]] <- dfaverage
      datasdlist[[ZT]] <- dfsd
    }
    
    dftoexport <- do.call(rbind, dataaveragelist)
    tdftoexport <- as.data.frame(t(dftoexport))[-1, ]
    tdftoexport <- tdftoexport[rowSums(tdftoexport == 0) <= ncol(tdftoexport) / 3, ]
    
    safe_item_name <- gsub("[ /]", "_", item)
    outfile_prefix <- paste0(path_to_save, "/Raw/", safe_item_name)
    write.csv(tdftoexport, paste0(outfile_prefix, "_data.csv"))
    
    tryCatch({
      if (nrow(tdftoexport) < 1 || ncol(tdftoexport) < 2) {
        message("  Skipping MetaCycle for '", item, "' due to insufficient data.")
      } else {
        timepointstolookat <- as.numeric(gsub("ZT", "", gsub(".*(ZT\\d{2})\\..*", "\\1", rownames(dftoexport))))
        method_to_use <- if (length(timepointstolookat) >= 18) c("LS", "JTK") else "LS"
        
        d <- meta2d(
          infile = paste0(outfile_prefix, "_data.csv"),
          outdir = path_to_save,
          filestyle = "csv", timepoints = timepointstolookat,
          minper = 20, maxper = 28, cycMethod = method_to_use,
          outputFile = FALSE
        )
        sigcyclegene <- if ("meta2d_pvalue" %in% names(d$meta)) d$meta %>% filter(meta2d_pvalue < 0.05) else d$meta %>% filter(LS_pvalue < 0.05)
        
        dftosd <- as.data.frame(t(do.call(rbind, datasdlist)))[-1, ]
        dftosd$gene <- rownames(dftosd) 
        tdftoexport$gene <- rownames(tdftoexport)
        d[["averagesheet"]] <- tdftoexport
        d[["sdsheet"]] <- dftosd
        d[["sig_cyl_gene"]] <- sigcyclegene
        
        writexl::write_xlsx(Filter(Negate(is.null), d), paste0(outfile_prefix, "_cyc_analysis.xlsx"))
        if(nrow(sigcyclegene) > 0) {
          siglist[[item]] <- sigcyclegene
        }
      }
    }, error = function(e) {
      message("  Skipping '", item, "' due to an error: ", e$message)
    })
  }
  
  if(length(siglist) > 0) {
    message("\nAnalysis complete. Writing summary files...")
    
    writexl::write_xlsx(siglist, paste0(path_to_save, "/Summary/", date, "_", run_name, "_cyc_siggene_by_", analysis_by, ".xlsx"))
    
    summary_df <- data.frame(group = names(siglist), cycling_gene_count = sapply(siglist, nrow))
    names(summary_df)[1] <- analysis_by
    write.csv(summary_df, paste0(path_to_save, "/Summary/", date, "_", run_name, "_cycling_gene_count_by_", analysis_by, ".csv"), row.names = FALSE)
    
    all_cycling_genes <- do.call(rbind, lapply(names(siglist), function(name) data.frame(gene=siglist[[name]]$CycID, item_name=name)))
    if (!is.null(all_cycling_genes) && nrow(all_cycling_genes) > 0) {
      gene_item_counts <- all_cycling_genes %>% 
        group_by(gene) %>% 
        summarize(group_count = n(), groups = paste(item_name, collapse = " | ")) %>% 
        arrange(desc(group_count))
      names(gene_item_counts) <- c("gene", paste0(analysis_by, "_count"), paste0(analysis_by, "s"))
      write.csv(gene_item_counts, paste0(path_to_save, "/Summary/", date, "_", run_name, "_groups_per_cycling_gene.csv"), row.names = FALSE)
    }
    message("...Summary files written successfully.\n")
  } else {
    message("No significant cycling genes found.")
  }
  
  message("===== All analyses complete. =====")
}
