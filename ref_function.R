# run_ = 'circa4'
# circascore = 'all' 

### Preprocessing_data
library(dplyr)
library(arrow)
library(MetaCycle)


PreProcessingData = function(run_name, path_to_data, circascore='all'){
  print("Start importing data")
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
  return(valid_regions)
}

### Clock genes
clockgenelist = c("Arntl","Clock", "Cry1","Cry2", "Npas2","Nr1d1", "Per1",  "Per2", "Per3", "Rora", "Rorb",
                "Rorc", "Dbp", "Nfil3", "Hlf", "Ciart")


MetaCycleCelltype = function(data, condition){
  if ("celltype" %in% condition){
    celltypes = GetValidCells(run_name)
  }
  for (cell in celltypes) {
    message("Processing: ", cell)
    data1 <- data[data$cell_type_final == cell, ]
    datasamplezt <- split.data.frame(data1, data1$ZT)
    
    dataaveragelist <- list()
    datasdlist <- list()
    
    for (ZT in names(datasamplezt)) {
      message("  Processing sample: ", ZT)
      df <- datasamplezt[[ZT]]
      group_labels <- rep(1:3, length.out = nrow(df))
      set.seed(123)
      df$group <- sample(group_labels)
      dfaverage <- df %>% group_by(group) %>% summarize(across(1:5006, mean))
      dfsd <- df %>% group_by(group) %>% summarize(across(1:5006, sd))
      dataaveragelist[[ZT]] <- dfaverage
      datasdlist[[ZT]] <- dfsd
    }
    
    dftoexport <- do.call(rbind, dataaveragelist)
    tdftoexport <- as.data.frame(t(dftoexport))
    tdftoexport <- tdftoexport[-1, ]
    #filter: drop genes with too many zeros
    tdftoexport <- tdftoexport[rowSums(tdftoexport == 0) <= ncol(tdftoexport) / 3, ]
    
    dftosd <- do.call(rbind, datasdlist)
    dftosd <- as.data.frame(t(dftosd))
    dftosd <- dftosd[-1, ]
    
    outfile_prefix <- paste0(path_to_save,"/Raw/", cell)
    write.csv(tdftoexport, paste0(outfile_prefix, "_data.csv"))
    write.csv(dftosd, paste0(outfile_prefix, "_sd.csv"))
    
    timepointstolookat <- rownames(dftoexport)
    timepointstolookat <- gsub(".*(ZT\\d{2})\\..*", "\\1", timepointstolookat)
    timepointstolookat <- gsub("circa2-ZT", "", timepointstolookat)
    timepointstolookat <- gsub("ZT", "", timepointstolookat)
    timepointstolookat <- as.numeric(timepointstolookat)
    
    tryCatch({
      message("  Timepoints length: ", length(timepointstolookat))
      message("  Data columns: ", ncol(tdftoexport) - 1)  # Subtract 1 for the 'gene' column
      if (length(timepointstolookat) == 18) {
        d <- meta2d(
          infile = paste0(outfile_prefix, "_data.csv"),
          outdir = paste0(path_to_save),
          filestyle = "csv", timepoints = timepointstolookat,
          minper = 20, maxper = 28, cycMethod = c("LS", "JTK"),
          analysisStrategy = "auto", outputFile = FALSE,
          outIntegration = "both", adjustPhase = "predictedPer",
          combinePvalue = "fisher", weightedPerPha = FALSE, ARSmle = "auto",
          ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE,
          outSymbol = "", parallelize = TRUE, nCores = 64, inDF = NULL
        )
        sigcyclegene <- d$meta %>% filter(meta2d_pvalue < 0.05)
      } else {
        d <- meta2d(
          infile = paste0(outfile_prefix, "_data.csv"),
          outdir = paste0(path_to_save),
          filestyle = "csv", timepoints = timepointstolookat,
          minper = 20, maxper = 28, cycMethod = c("LS"),
          analysisStrategy = "auto", outputFile = FALSE,
          outIntegration = "both", adjustPhase = "predictedPer",
          combinePvalue = "fisher", weightedPerPha = FALSE, ARSmle = "auto",
          ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE,
          outSymbol = "", parallelize = FALSE, nCores = 16, inDF = NULL
        )
        sigcyclegene <- d$meta %>% filter(LS_pvalue < 0.05)
      }
      
      tdftoexport$gene <- rownames(tdftoexport)
      dftosd$gene <- rownames(dftosd)
      d[["averagesheet"]] <- tdftoexport
      d[["sdsheet"]] <- dftosd
      d[["sig_cyl_gene"]] <- sigcyclegene
      my_list_clean <- Filter(Negate(is.null), d)
      
      writexl::write_xlsx(my_list_clean, paste0(outfile_prefix, "_cyc_analysis.xlsx"))
      siglist[[cell]] <- sigcyclegene
      
    }, error = function(e) {
      message("Skipping cell type ", cell, " due to error: ", e$message)
    })
  }
  
  #contains all significantly cycling genes for each cell type
  writexl::write_xlsx(siglist, paste0(path_to_save,"/Summary/",date,"_",run_name,"_cyc_siggene_analysis.xlsx"))
  
  #create summary of cycling genes per cell type
  summary_df <- data.frame(
    cell_type = names(siglist),
    cycling_gene_count = sapply(siglist, nrow)
  )
  
  #write summary to CSV
  write.csv(summary_df, paste0(path_to_save,"/Summary/",date,"_",run_name,"_cycling_gene_per_celltype.csv"), row.names = FALSE)
  
  #create a list of all significantly cycling genes with their cell types
  all_cycling_genes <- list()
  
  for (cell in names(siglist)) {
    if (nrow(siglist[[cell]]) > 0) {
      # Extract gene names from each cell type's significant gene list
      genes <- siglist[[cell]]$CycID
      
      # Add each gene with its cell type to our list
      for (gene in genes) {
        if (gene %in% names(all_cycling_genes)) {
          all_cycling_genes[[gene]] <- c(all_cycling_genes[[gene]], cell)
        } else {
          all_cycling_genes[[gene]] <- cell
        }
      }
    }
  }
  
  #create a data frame showing each gene and the count of cell types it appears in
  gene_celltype_counts <- data.frame(
    gene = names(all_cycling_genes),
    cell_type_count = sapply(all_cycling_genes, length),
    cell_types = sapply(all_cycling_genes, function(x) paste(x, collapse = ", "))
  )
  
  #sort by number of cell types (descending)
  gene_celltype_counts <- gene_celltype_counts[order(-gene_celltype_counts$cell_type_count), ]
  
  #write to CSV
  write.csv(gene_celltype_counts, paste0(path_to_save,"/Summary/",date,"_",run_name,"_","celltype_per_cycling_gene.csv"), row.names = FALSE)
}

MetaCycleRegion = function(data, condition){
  if ("region" %in% condition){
    regions = GetValidRegions(run_name)
  }
  
  for (cell in regions) {
    message("Processing: ", cell)
    data1 <- data[data$region_automap_name == cell, ]
    datasamplezt <- split.data.frame(data1, data1$sample)
    
    dataaveragelist <- list()
    datasdlist <- list()
    
    for (ZT in names(datasamplezt)) {
      message("  Processing sample: ", ZT)
      df <- datasamplezt[[ZT]]
      group_labels <- rep(1:3, length.out = nrow(df))
      set.seed(123)
      df$group <- sample(group_labels)
      dfaverage <- df %>% group_by(group) %>% summarize(across(1:5006, mean))
      dfsd <- df %>% group_by(group) %>% summarize(across(1:5006, sd))
      dataaveragelist[[ZT]] <- dfaverage
      datasdlist[[ZT]] <- dfsd
    }
    
    dftoexport <- do.call(rbind, dataaveragelist)
    tdftoexport <- as.data.frame(t(dftoexport))
    tdftoexport <- tdftoexport[-1, ]
    #filter: drop genes with too many zeros
    tdftoexport <- tdftoexport[rowSums(tdftoexport == 0) <= ncol(tdftoexport) / 3, ]
    
    dftosd <- do.call(rbind, datasdlist)
    dftosd <- as.data.frame(t(dftosd))
    dftosd <- dftosd[-1, ]
    
    outfile_prefix <- paste0(path_to_save,"/Raw/", cell)
    write.csv(tdftoexport, paste0(outfile_prefix, "_data.csv"))
    write.csv(dftosd, paste0(outfile_prefix, "_sd.csv"))
    
    timepointstolookat <- rownames(dftoexport)
    timepointstolookat <- gsub(".*(ZT\\d{2})\\..*", "\\1", timepointstolookat)
    timepointstolookat <- gsub("circa2-ZT", "", timepointstolookat)
    timepointstolookat <- gsub("ZT", "", timepointstolookat)
    timepointstolookat <- as.numeric(timepointstolookat)
    
    tryCatch({
      message("  Timepoints length: ", length(timepointstolookat))
      message("  Data columns: ", ncol(tdftoexport) - 1)  # Subtract 1 for the 'gene' column
      if (length(timepointstolookat) == 18) {
        d <- meta2d(
          infile = paste0(outfile_prefix, "_data.csv"),
          outdir = paste0(path_to_save),
          filestyle = "csv", timepoints = timepointstolookat,
          minper = 20, maxper = 28, cycMethod = c("LS", "JTK"),
          analysisStrategy = "auto", outputFile = FALSE,
          outIntegration = "both", adjustPhase = "predictedPer",
          combinePvalue = "fisher", weightedPerPha = FALSE, ARSmle = "auto",
          ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE,
          outSymbol = "", parallelize = FALSE, nCores = 16, inDF = NULL
        )
        sigcyclegene <- d$meta %>% filter(meta2d_pvalue < 0.05)
      } else {
        d <- meta2d(
          infile = paste0(outfile_prefix, "_data.csv"),
          outdir = paste0(path_to_save),
          filestyle = "csv", timepoints = timepointstolookat,
          minper = 20, maxper = 28, cycMethod = c("LS"),
          analysisStrategy = "auto", outputFile = FALSE,
          outIntegration = "both", adjustPhase = "predictedPer",
          combinePvalue = "fisher", weightedPerPha = FALSE, ARSmle = "auto",
          ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE,
          outSymbol = "", parallelize = FALSE, nCores = 16, inDF = NULL
        )
        sigcyclegene <- d$meta %>% filter(LS_pvalue < 0.05)
      }
      
      tdftoexport$gene <- rownames(tdftoexport)
      dftosd$gene <- rownames(dftosd)
      d[["averagesheet"]] <- tdftoexport
      d[["sdsheet"]] <- dftosd
      d[["sig_cyl_gene"]] <- sigcyclegene
      my_list_clean <- Filter(Negate(is.null), d)
      
      writexl::write_xlsx(my_list_clean, paste0(outfile_prefix, "_cyc_analysis.xlsx"))
      siglist[[cell]] <- sigcyclegene
      
    }, error = function(e) {
      message("Skipping cell type ", cell, " due to error: ", e$message)
    })
  }
  
  #contains all significantly cycling genes for each cell type
  writexl::write_xlsx(siglist, paste0(path_to_save,"/Summary/",date,"_",run_name,"_cyc_siggene_analysis.xlsx"))
  
  #create summary of cycling genes per cell type
  summary_df <- data.frame(
    cell_type = names(siglist),
    cycling_gene_count = sapply(siglist, nrow)
  )
  
  #write summary to CSV
  write.csv(summary_df, paste0(path_to_save,"/Summary/",date,"_",run_name,"_cycling_gene_per_celltype.csv"), row.names = FALSE)
  
  #create a list of all significantly cycling genes with their cell types
  all_cycling_genes <- list()
  
  for (cell in names(siglist)) {
    if (nrow(siglist[[cell]]) > 0) {
      # Extract gene names from each cell type's significant gene list
      genes <- siglist[[cell]]$CycID
      
      # Add each gene with its cell type to our list
      for (gene in genes) {
        if (gene %in% names(all_cycling_genes)) {
          all_cycling_genes[[gene]] <- c(all_cycling_genes[[gene]], cell)
        } else {
          all_cycling_genes[[gene]] <- cell
        }
      }
    }
  }
  
  #create a data frame showing each gene and the count of cell types it appears in
  gene_celltype_counts <- data.frame(
    gene = names(all_cycling_genes),
    cell_type_count = sapply(all_cycling_genes, length),
    cell_types = sapply(all_cycling_genes, function(x) paste(x, collapse = ", "))
  )
  
  #sort by number of cell types (descending)
  gene_celltype_counts <- gene_celltype_counts[order(-gene_celltype_counts$cell_type_count), ]
  
  #write to CSV
  write.csv(gene_celltype_counts, paste0(path_to_save,"/Summary/",date,"_",run_name,"_","celltype_per_cycling_gene.csv"), row.names = FALSE)
}