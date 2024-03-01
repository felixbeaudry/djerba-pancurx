# FILE PATH VARIABLES
DONOR = 'donor'
TUMOUR_SAMPLE_ID = 'tumour'
NORMAL_SAMPLE_ID = 'normal'
SEQTYPE = 'seqtype'
ALIGNER = 'aligner'
ALIGNER_VERSION = 'aligner_version'

# FILE PATH DEFAULTS
DEFAULT_ROOTPATH = '/.mounts/labs/PCSI/pipeline/hg38'
DEFAULT_SEQTYPE = 'wgs'
DEFAULT_ALIGNER = 'bwa'
DEFAULT_ALIGNER_VERSION = '0.7.17'
DEFAULT_NEOPATHS = "netMHC/pan-2.8/polysolver/1.0/"
#TODO: replace wildcard!
DEFAULT_XENOME_PATHS = "xenome/1.0.1-r/*.log" 

# DATA FILES
DEFAULT_COSMIC_CENSUS_FILE_PATH = "/.mounts/labs/PCSI/users/rdenroche/phoenix/phoenixParser/cosmic_census_20161115.tsv"
DEFAULT_DSBR_LIST_FILE_PATH = "/.mounts/labs/PCSI/users/rdenroche/lists/dsbr_list.tsv"
DEFAULT_MMR_LIST_FILE_PATH = "/.mounts/labs/PCSI/users/rdenroche/lists/mmr_list.tsv"
DEFAULT_HUGO_SYN_FILE_PATH = "/.mounts/labs/PCSI/raw_data/hugo/HUGO_synonyms_171213.txt"
DEFAULT_REFSEQ_FILE_PATH = "/.mounts/labs/PCSI/references/GRCh38/rna/refSeq/hg38_refSeq_genes.txt"

#TODO: does not exist
DEFAULT_RNA_SIG_FILE_PATH = "/.mounts/labs/PCSI/users/azhang/RNA/RNA-subtype.txt"

COSMIC_MUTATION_TYPES = {
			"A" : "strong amplification",
			"D" : "homozygous deletion",
			"F" : "frameshift",
			"Mis" : "nonsynonymous,nonframeshift",
			"N" : "stopgain",
			"S" : "splicing",
			"T" : "translocation",
        }