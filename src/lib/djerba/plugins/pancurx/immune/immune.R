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
  names(comparison_matrix)[1] <- c('HGNC')
  
  #gene_set should be 'hgnc_complete_set.txt' from https://www.genenames.org/download/archive/
  gene_set_data <-  read.csv(opt$genes, sep="\t", header=TRUE, check.names=FALSE)[]
  
#### PROCESS AND CLEAN DATA ####
  
  if(file.exists(sample_tpm_file_path)){
  
    gene_list <- gene_set_data[gene_set_data$locus_group == "protein-coding gene",c("symbol","ensembl_gene_id" )]
    
    #input_data <- read.csv(opt$input, sep="\t", header=TRUE, check.names=FALSE)
    tpm_file <- fread(sample_tpm_file_path)[,c("Gene Name", "Coverage")]
  
    tpm_file <- left_join(gene_list, tpm_file, by=c("symbol" = "Gene Name"), multiple = "all")
    names(tpm_file)[1] <- c('HGNC')
    
    TPMs <- tpm_file %>%
      group_by(HGNC) %>%
      summarise(this_sample = round(mean(Coverage), 2))
    
    TPMs$this_sample[is.na(TPMs$this_sample)] <- 0
    
  } else {
    cat("rna not found \n")
    
  }
  
  #load the sample to be processed and join with the control matrix
  TPMs_full =  as.data.frame(inner_join(TPMs, comparison_matrix))
  
  rownames(TPMs_full) <- TPMs_full$HGNC
  TPMs_full <- TPMs_full[,-1] #remove "HGNC" column

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
  

#### OUTPUT ####
  
  write.table(
    results_dataframe,
    file = output_file_path,
    row.names = FALSE, quote = FALSE, sep = "\t"
  )
  
  # couldn't get immunedeconv to turn off writing file to working directory
  file.remove(paste0('./CIBERSORT-Results.txt'))
  file.remove(paste0(opt$cibersort, '/CIBERSORT-Results.txt'))
  
#### PLOT ####  
  

  cells_by_sample <-   transpose(results_dataframe[-1])
  names(cells_by_sample) <- gsub(pattern = " ", replacement = "_", x = results_dataframe$cell_type)
  names(cells_by_sample) <- gsub(pattern = "\\+", replacement = "pos", x = names(cells_by_sample))
  
  this_equals_zero = 1e-05
  cells_by_sample[cells_by_sample == 0] <- this_equals_zero
  
  svg(paste0(output_file_path,".svg"), width=3.65 * 2, height=1.56 *2)

    #ADD DISTRIBUTION TO PLOT
    par(mar=c(0,1,0,1)+0.1, fig=c(0,1,0.78,1))
    boxplot(log(cells_by_sample$T_cell_CD8pos[-1]), ylim=c(log(this_equals_zero), log(0.5)), horizontal=T,   ylab="Indel Count", xaxt="n",xaxs="i", frame=F, width=0.5, pch=20)
    points(log(cells_by_sample$T_cell_CD8pos[1]), 1, cex=2, pch=4)
    mtext(expression(paste( plain("CD8+ T"))),side=2, cex=0.75)
    
    par(mar=c(0,1,0,1)+0.1, fig=c(0,1,0.56,0.78),new=T)
    boxplot(log(cells_by_sample$T_cell_CD4pos_naive[-1]), ylim=c(log(this_equals_zero), log(0.5)), horizontal=T,   ylab="Indel Count", xaxt="n",xaxs="i", frame=F, width=0.5, pch=20)
    points(log(cells_by_sample$T_cell_CD4pos_naive[1]), 1, cex=2, pch=4)
    mtext(expression(paste( plain("CD4+ T"))),side=2, cex=0.75)
    
    par(mar=c(0,1,0,1)+0.1, fig=c(0,1,0.34,0.56),new=T)
    boxplot(log(cells_by_sample$NK_cell_activated[-1]), ylim=c(log(this_equals_zero), log(0.5)), horizontal=T,   ylab="Indel Count", xaxt="n",xaxs="i", frame=F, width=0.5, pch=20)
    points(log(cells_by_sample$NK_cell_activated[1]), 1, cex=2, pch=4)
    mtext(expression(paste( plain("NK T"))),side=2, cex=0.75)
    
    par(mar=c(2,1,0,1)+0.1,fig=c(0,1,0,0.34),new=T)
    boxplot(log(cells_by_sample$Macrophage_M1[-1]), ylim=c(log(this_equals_zero), log(0.5)), horizontal=T,   xaxs="i",  frame=F,width=0.5, pch=20)
    points(log(cells_by_sample$Macrophage_M1[1]), 1, cex=2, pch=4)
    mtext(expression(paste( plain("M1 Mac."))),side=2, cex=0.75)
    

  dev.off()
  
  

  cat("success!\n")
