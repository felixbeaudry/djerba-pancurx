
#sample_tpm_file_path <- "/Volumes/pcsi/pipeline/hg38/BTC0023/BTC_0023_Lv_M_526/rna/star/2.7.4a/stringtie/2.0.6/BTC_0023_Lv_M_526_stringtie_abundance.txt"
sample_tpm_file_path = '/Volumes/pcsi/users/fbeaudry/reports/BTC0053/tpm.txt'

#genes_path = '/Volumes/pcsi/users/fbeaudry/djerba-pancurx/src/lib/djerba/data/pancurx/hgnc_complete_set.txt'
#gene_set_data <-  read.csv(genes_path, sep="\t", header=TRUE, check.names=FALSE) 

MATRIX_file = '/Volumes/pcsi/references/cohort_lists/LBR.immune_rna.txt'

comparison_matrix = read.csv(MATRIX_file, sep="\t", header=TRUE)


CIBERSORT_PATH = '/Volumes/PCSI/users/fbeaudry/reports/djerba_dev/lib/python3.10/site-packages/djerba/plugins/pancurx/immune/CIBERSORT.R'
CIBERSORT_MAT_PATH = '/Volumes/pcsi/references/cohort_lists/LM22.txt'

# CIBERSORT.R file path
set_cibersort_binary(CIBERSORT_PATH) 
# Cibersort loadings matrix file LM22.txt 
set_cibersort_mat(CIBERSORT_MAT_PATH) 
