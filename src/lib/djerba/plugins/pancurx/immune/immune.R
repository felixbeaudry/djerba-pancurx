#! /usr/bin/env Rscript

library(tidyr)
library(dplyr)
library(optparse)
library(immunedeconv)
library(data.table)

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

# CIBERSORT.R file path
set_cibersort_binary(paste0(opt$cibersort, '/CIBERSORT.R'))
# Cibersort loadings matrix file, eg. LM22.txt 
set_cibersort_mat(opt$loadings) 

#load the control matrix # Path to matrix of comparator TPMs used to rank the results 
comparison_matrix = read.csv(opt$comparison, sep=" ", header=TRUE)

#gene_set should be 'hgnc_complete_set.txt' from https://www.genenames.org/download/archive/
gene_set_data <-  read.csv(opt$genes, sep="\t", header=TRUE, check.names=FALSE)[]
gene_list <- gene_set_data[gene_set_data$locus_group == "protein-coding gene",c("symbol","ensembl_gene_id" )]

sample_tpm_file_path <- opt$input
output_file_path  = opt$output

if(file.exists(sample_tpm_file_path)){

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

if (nrow(TPMs_full[ TPMs_full$this_sample > 0,]) == 0) {
  
  # No data available, we need to produce an empty file
  cat("No Data in TPM file, writing empty file\n")
  res = data.frame(matrix(c(0,0,0), nrow = 1))
  colnames(res) = c("cell_type","Score","ptile")
  
} else {
  
  # process with immunedeconv
  res = immunedeconv::deconvolute(TPMs_full, "cibersort_abs")

}

write.table(
  res,
  file = output_file_path,
  row.names = FALSE, quote = FALSE, sep = "\t"
)

# couldn't get immunedeconv to turn off writing file to working directory
file.remove(paste0('./CIBERSORT-Results.txt'))
file.remove(paste0(opt$cibersort, '/CIBERSORT-Results.txt'))
cat("success!\n")
