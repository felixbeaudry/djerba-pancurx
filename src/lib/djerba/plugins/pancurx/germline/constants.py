GERM_VARIANT_COUNT = "germ_variant_count" 
GERM_NONSILENT_COUNT = "germ_nonsilent_count"
GERM_NONSIL_SUBSET_COUNT = 'germ_nonsil_genes'
GERM_NONSIL_SUBSET_RARE_COUNT = 'germ_nonsil_genes_rare'
GERM_PATHOGENIC_COUNT = 'germ_pathogenic'

EXCLUDED_GERMLINE_PATHOGENIC_VARIANTS = ['rs11571833']

GERMLINE_GENE_ORDER = ["POLE", "POLD1", "EPCAM", "MLH1", "MSH2", "MSH6", "PMS1", "PMS2", "BRCA1", "BRCA2", "PALB2", "ATM", "APC", "MUTYH", "CDKN2A", "STK11", "TP53", "SMAD4", "RAD51C", "CHEK2" ,"BRIP1", "RNF43"]

NONSILENT_CHANGES =  [
	"frameshift" ,
	"nonframeshift" ,
	"nonsynonymous" ,
	"stoploss" ,
	"stopgain" ,
	"splicing" ,
	"strong amplification" ,
	"homozygous deletion",
	"deletion breakpoint" ,
	"duplication breakpoint" ,
	"inversion breakpoint" ,
	"translocation breakpoint" 
]