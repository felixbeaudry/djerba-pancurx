#! /usr/bin/env Rscript

library(tidyr)
library(dplyr)
library(optparse)
library(immunedeconv)
library(data.table)
library(ggplot2)

#### TAKE IN FILES AND SET CONFIG ####

  option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, help="RSEM .genes.results input file", metavar="PATH"),
    make_option(c("-g", "--genes"), type="character", default=NULL, help="Ensembl gene name conversion file", metavar="PATH"),
    make_option(c("-c", "--comparison"), type="character", default=NULL, help="comparison matrix file", metavar="PATH"),
    make_option(c("-l", "--loadings"), type="character", default=NULL, help="comparison matrix file", metavar="PATH"),
    make_option(c("-p", "--cibersort"), type="character", default=NULL, help="comparison matrix file", metavar="PATH"),
    make_option(c("-o", "--output"), type="character", default=NULL, help="cibersort scores and percentiles", metavar="PATH")
  )
  
  opt_parser <- OptionParser(option_list=option_list, add_help_option=TRUE)
  opt <- parse_args(opt_parser)
  
  sample_tpm_file_path <- opt$input
  output_file_path  = opt$output
  
  # CIBERSORT.R file path
  set_cibersort_binary(paste0(opt$cibersort, '/CIBERSORT.R'))
  
  # Cibersort loadings matrix file, eg. LM22.txt 
  set_cibersort_mat(opt$loadings) 
  
  #load the control matrix # Path to matrix of comparator TPMs used to rank the results 
  comparison_matrix = read.csv(opt$comparison, sep="\t", header=TRUE)

#### PROCESS AND CLEAN DATA ####
  
  if(file.exists(sample_tpm_file_path)){

    TPMs <- fread(sample_tpm_file_path)
     
  } else {
    cat("rna not found \n")
    
  }
  
  #load the sample to be processed and join with the control matrix
  TPMs_full =  as.data.frame(inner_join(TPMs[,-3], comparison_matrix))
  
  rownames(TPMs_full) <- TPMs_full$gene_list
  TPMs_full <- TPMs_full[,-1] 

  #some genes missing from analysis, may impact results
  TPMs_full <- TPMs_full[!is.na(TPMs_full$this_sample),]
  
#### RUN ANALYSIS IF POSSIBLE ####
  
  if (nrow(TPMs_full[ TPMs_full$this_sample > 0,]) == 0) {
    
    # No data available, we need to produce an empty file
    cat("No Data in TPM file, writing empty file\n")
    results_dataframe = data.frame(matrix(c(0,0,0), nrow = 1))
    colnames(res) = c("cell_type","Score","ptile")
    
  } else {
    
    # process with immunedeconv
    results_dataframe = immunedeconv::deconvolute(TPMs_full, "cibersort_abs")
  
  }
  
  this_result <- results_dataframe %>% dplyr::select(cell_type, this_sample)
  cohort_file <- results_dataframe %>% dplyr::select(-this_sample)
  
  this_result$this_percentile <- NA
  for(this_cell_type in this_result$cell_type){
    ## assumes genes missing from the cohort will be missing in the sample
    if(length(this_result$this_sample[this_result$cell_type == this_cell_type]) > 0){
      percentile_function <- ecdf(as.numeric(unlist(cohort_file[cohort_file$cell_type == this_cell_type,])[-1]))
      this_result$this_percentile[this_result$cell_type == this_cell_type] <- percentile_function(this_result$this_sample[this_result$cell_type == this_cell_type])
    }
  }
  
#### OUTPUT ####
  
  write.table(
    this_result,
    file = output_file_path,
    row.names = FALSE, quote = FALSE, sep = "\t"
  )
  
  # couldn't get immunedeconv to turn off writing file to working directory
  file.remove(paste0('./CIBERSORT-Results.txt'))
  file.remove(paste0(opt$cibersort, '/CIBERSORT-Results.txt'))
  

