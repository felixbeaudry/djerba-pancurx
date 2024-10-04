#! /usr/bin/env Rscript

library(dplyr)
library(data.table)
library(optparse)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="input tpm file", metavar="PATH"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output processed tpm file", metavar="PATH"),
  make_option(c("-c", "--cohort"), type="character", default=NULL, help="cohort comparison file", metavar="PATH")
)

opt_parser <- OptionParser(option_list=option_list, add_help_option=TRUE)
opt <- parse_args(opt_parser)

file_path <- opt$input
cohort_path <- opt$cohort
output_path <- opt$output

#donor_id='BTC0079'
#sample_id='BTC_0079_Lv_P_526'
#file_path= paste0("/Volumes/pcsi/pipeline/hg38/", donor_id,"/", sample_id,"/rna/star/2.7.4a/stringtie/2.0.6/",sample_id, "_stringtie_abundance.txt")
#cohort_path=paste0("/Volumes/pcsi/users/fbeaudry/reports/djerba-pancurx/src/lib/djerba/data/pancurx/LBR.immune_rna.txt")

if(file.exists(file_path) & file.exists(cohort_path)){
  
  cohort_file <- fread(cohort_path)
  gene_file <- fread(file_path)[,c("Gene Name", "TPM")]
  names(gene_file)[1] <- c('gene_name')
  
  gene_group <- gene_file %>%
    group_by(gene_name) %>%
    summarise(TPM = sum(TPM))
  
  gene_filtered <- gene_group %>% filter( gene_name %in% cohort_file$gene_list)
  gene_filtered$this_sample <- gene_filtered$TPM / sum(gene_filtered$TPM) * 1000000

  gene_filtered$this_percentile <- NA
  for(this_gene in cohort_file$gene_list){
    ## assumes genes missing from the cohort will be missing in the sample
    if(length(gene_filtered$this_sample[gene_filtered$gene_name == this_gene]) > 0){
      percentile_function <- ecdf(as.numeric(unlist(cohort_file[cohort_file$gene_list == this_gene,])[-1]))
      gene_filtered$this_percentile[gene_filtered$gene_name == this_gene] <- percentile_function(gene_filtered$this_sample[gene_filtered$gene_name == this_gene])
    }
  }
  
  joined_file <- left_join(cohort_file, gene_filtered[,c('gene_name','this_sample','this_percentile')], by=c("gene_list"="gene_name"))
  
  write.table( joined_file[,c('gene_list','this_sample','this_percentile')],  file = output_path, row.names = F, quote = FALSE, sep = "\t", col.names = T)
  
} else { cat("rna not found for ",file_path, "\n") 
  write.table( c(),  file = output_path, row.names = T, quote = FALSE, sep = "\t", col.names = F)}
    
  

