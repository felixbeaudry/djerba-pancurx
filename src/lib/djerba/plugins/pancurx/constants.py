#defaults
CURRENT_REPORT_VERSION = '1.0'
DIPLOID_CUTOFF = 2.3
DEFAULT_DELETION_CUTOFF = 0.5
DEFAULT_AMP_MULTIPLIER_CUTOFF = 4
DEFAULT_DEL_LOH_CUTOFF = 1.5
DEFAULT_GAIN_ADDEND_CUTOFF = 1

# FILE PATH DEFAULTS
DEFAULT_CORE_LIMS_URL = "http://pinery.gsi.oicr.on.ca/samples?type=Identity"
DEFAULT_ROOTPATH = '/.mounts/labs/PCSI/pipeline/hg38'
#TEMP: plot path is different
DEFAULT_PLOTPATH = '/.mounts/labs/PCSI/users/fbeaudry/reports/btc_plots'
DEFAULT_DATA_LOCATION = '/.mounts/labs/PCSI/users/fbeaudry/reports/djerba-run-data/'
DEFAULT_ARCHIVEPATH = 'https://www.hpc.oicr.on.ca/archive/projects/PCSI/pipeline/hg38'
DEFAULT_SEQTYPE = 'wgs'
DEFAULT_ALIGNER = 'bwa'
DEFAULT_ALIGNER_VERSION = '0.7.17'
DEFAULT_CELLULOID_VERSION = "v0.11.7"
DEFAULT_BLURB_FILE = 'blurb_template.txt'
DEFAULT_MANE_FILE = 'MANE.GRCh38.v1.3.summary.txt'

#DEFAULT_CIBERSORT_PATH='CIBERSORT.R'
DEFAULT_CIBERSORT_LOADINGS_PATH='pancurx/LM22.txt'
DEFAULT_CIBERSORT_GENES_PATH='pancurx/hgnc_complete_set.txt'
DEFAULT_CIBERSORT_COMPARISON_PATH='/.mounts/labs/PCSI/users/fbeaudry/reports/djerba-run-data/LBR.immune_rna.txt'


DEFAULT_RNA_SEQTYPE = 'rna'
DEFAULT_RNA_ALIGNER = 'star'
DEFAULT_RNA_ALIGNER_VERSION = '2.7.4a'

# DATA FILES
DEFAULT_HUGO_SYN_FILE_PATH = "/.mounts/labs/PCSI/raw_data/hugo/HUGO_synonyms_171213.txt"
DEFAULT_COSMIC_CENSUS_FILE_PATH = "/.mounts/labs/PCSI/users/rdenroche/phoenix/phoenixParser/cosmic_census_20161115.tsv" #TODO: check use?
DEFAULT_REFSEQ_FILE_PATH = "/.mounts/labs/PCSI/references/GRCh38/rna/refSeq/hg38_refSeq_genes.txt"
DEFAULT_REFERENCE_COHORT = '/.mounts/labs/PCSI/users/fbeaudry/PanCuRx_Analysis_Pipeline/scripts/plot/data/pdac_20170926_matrix.csv'

DEFAULT_GENE_FILE = 'genes.txt'
DEFAULT_GERMLINE_GENE_FILE = 'germline.genes.txt'
DEFAULT_SIGNATURES_FILE = 'signatures.txt'

# report attributes
REPORT_VERSION = "report_version"
SLIDE_DATE = 'slide_date'
REPORT_DATE = "report_date"

# case attributes
TUMOUR_ID = "tumour"
TUMOUR_SAMPLE_ID = "tumour"
DONOR = "donor"
DONOR_ID = "donor"
NORMAL_ID = "normal"
NORMAL_SAMPLE_ID = "normal"
EXTERNAL_IDS = "external_ids"
TUMOUR_TISSUE = "tumour_tissue"
NORMAL_TISSUE = "normal_tissue"
SEQTYPE = 'seqtype'
ALIGNER = 'aligner'
ALIGNER_VERSION = 'aligner_version'
SAMPLE_TYPE = 'sample_type'
LOCATION_SUBTYPE = 'location_subtype'

# sample attributes
SUMMARY_FILE = 'summary_file'
CELLULARITY = "cellularity"
TUMOUR_COVERAGE = "tumour_coverage"
NORMAL_COVERAGE = "normal_coverage"
TMB = 'tmb'
SNV_LOAD = 'snv_load'
INDEL_LOAD = 'indel_load'
SV_LOAD = 'sv_load'
PLOIDY = 'ploidy'
MOFFIT = 'moffit'
HRD = 'hrd'
MENGHI = 'menghi_tdp'
MMRD = 'mmrd'
GERM_VARIANT_COUNT = "germ_variant_count" 
GERM_NONSILENT_COUNT = "germ_nonsilent_count"
GERM_NONSIL_SUBSET_COUNT = 'germ_nonsil_genes'
GERM_NONSIL_SUBSET_RARE_COUNT = 'germ_nonsil_genes_rare'
GERM_PATHOGENIC_COUNT = 'germ_pathogenic'
ALEXANDROV_CLASS = "alexandrov_class"
COLLISSON_CLASS = "collisson_class"
WADDELL_CLASS = "waddell_class"
MOFFITT_CLASS = "moffitt_class"
INFERRED_SEX = "inferred_sex"
HLA_TYPES = "hla_types"

# strings
NONE_SPECIFIED = "NONE_SPECIFIED"
VALUE = "value"
ABOVE_CUTOFF = "above_cutoff"

#paths

ONCOSLIDE_DRIVER_PLOT  = "oncoslide_driver_plot"
ONCOSLIDE_SNV_PLOT  = "oncoslide_snv_plot"
ONCOSLIDE_INDEL_PLOT  = "oncoslide_indel_plot"
ONCOSLIDE_SV_PLOT  = "oncoslide_sv_plot"
ONCOSLIDE_CNV_PLOT  = "oncoslide_cnv_plot"
PARAM_PATH = "param_path"
CELLULOID_PARAM_PATH = "param_path"
TDP_PATH = "tdp_path"
COSMIC_SIGNNLS_PATH = "cosmic_signnls_path"
SUMMARY_FILE_PATH = "summary_file_path"
CELLULOID_PLOT = "celluloid_plot"
CELLULOID_DIR = 'celluloid_dir'
TUMOUR_COVERAGE_PATH = "tumour_coverage_path"
NORMAL_COVERAGE_PATH = "normal_coverage_path"
TUMOUR_BAM_PATH = "tumour_coverage_path"
NORMAL_BAM_PATH = "normal_coverage_path"
SEX_PATH = "sex_path"
GERMLINE_PATH = "germline_path"
GERMLINE_ANNOVAR_PATH = "germline_annovar_path"
DSBR_SCORE_BAR = "dsbr_score_image_path"
MMR_SCORE_BAR = "mmr_score_image_path"
SAMPLE_VARIANTS_FILE = "sample_variants"
WHOLE_GENOME_PLOT = "whole_genome_plot_path"
INDEL_VAF_PLOT = "indel_vaf_plot"
SNV_VAF_PLOT = "snv_vaf_plot"
INDEL_BIN_PLOT = "indel_bin_plot"
SV_BIN_PLOT = 'sv_bin_plot'
SNV_CONTEXT_PLOT = "snv_context_plot"
COSMIC_STACK_PLOT = 'cosmic_stack_plot'
TMB_STACK_PLOT = 'tmb_stack_plot'
SV_STACK_PLOT = 'sv_stack_plot'
HISTBOX_INDEL = "histbox_indel"
HISTBOX_SV = "histbox_sv"
HISTBOX_SNV = "histbox_snv"
WHOLE_CHROMOSOME_PLOTS = "whole_chromosome_plots"
DSBR_RESULTS = "DSBR_results"
MMR_RESULTS = "MMR_results"
TDP_RESULTS = "TDP_results"
REPORTING_NAME = "reporting_name"
MAVIS_FUSIONS_PATH = "mavis_fusions"
TPM_PATH = "tpm_path"
STAR_QC_PATH = "star_qc"

CHROMOSOME_1_PLOT = 'chromosome_1_plot'
CHROMOSOME_2_PLOT = 'chromosome_2_plot' 
CHROMOSOME_3_PLOT = 'chromosome_3_plot'
CHROMOSOME_4_PLOT = 'chromosome_4_plot'
CHROMOSOME_5_PLOT = 'chromosome_5_plot'
CHROMOSOME_6_PLOT = 'chromosome_6_plot'
CHROMOSOME_7_8_PLOT = 'chromosome_7_8_plot'
CHROMOSOME_9_10_PLOT = 'chromosome_9_10_plot'
CHROMOSOME_11_12_PLOT = 'chromosome_11_12_plot'
CHROMOSOME_13_14_PLOT = 'chromosome_13_14_plot'
CHROMOSOME_15_17_PLOT = 'chromosome_15_17_plot'
CHROMOSOME_18_22_PLOT = 'chromosome_18_22_plot'
CHROMOSOME_X_Y_PLOT = 'chromosome_x_y_plot'

CHROMOSOME_PLOTS = {"chr1" : CHROMOSOME_1_PLOT, 
                    "chr2" : CHROMOSOME_2_PLOT, 
                    "chr3" : CHROMOSOME_3_PLOT, 
                    "chr4" : CHROMOSOME_4_PLOT, 
                    "chr5" : CHROMOSOME_5_PLOT, 
                    "chr6" : CHROMOSOME_6_PLOT, 
                    "chr7-8" : CHROMOSOME_7_8_PLOT, 
                    "chr9-10" : CHROMOSOME_9_10_PLOT, 
                    "chr11-12" : CHROMOSOME_11_12_PLOT,
                    "chr13-14" : CHROMOSOME_13_14_PLOT, 
                    "chr15-17" : CHROMOSOME_15_17_PLOT, 
                    "chr18-22" : CHROMOSOME_18_22_PLOT, 
                    "chrX-Y" : CHROMOSOME_X_Y_PLOT}

EXCLUDED_GERMLINE_PATHOGENIC_VARIANTS = ['rs11571833']

DSBR_GENES = ["BRCA1", "BRCA2", "PALB2", "RAD51", "RAD51B", "RAD51C", "RAD51D", "XRCC2", "XRCC3", "DMC1", "BRIP1"]
MMR_GENES = ["MLH1", "MSH2", "MSH6", "PMS2", "EPCAM"]

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


NONSILENT_CHANGES_RANK =  {
    "nonsynonymous"  : 0,
    "stopgain" : 1,
	"frameshift" : 2,
    "splicing" : 3,
	"nonframeshift" : 4 ,
	"stoploss" : 5,
}


COSMIC_MUTATION_TYPES = {
	"A" : "strong amplification",
	"D" : "homozygous deletion",
	"F" : "frameshift",
	"Mis" : "nonsynonymous,nonframeshift",
	"N" : "stopgain",
	"S" : "splicing",
	"T" : "translocation",
}

#MOVED TO signatures.txt Aug 21st, 2024 by Fe
# COSMIC_SIGNATURE_SET = {
# 	"Signature.1": "SBS1",
# 	"Signature.2": "SBS2",
# 	"Signature.3": "SBS3",
# 	"Signature.5": "SBS5",
# 	"Signature.6": "SBS6",
# 	"Signature.8": "SBS8",
# 	"Signature.13": "SBS13",
# 	"Signature.17a": "SBS17a",
#     "Signature.17b": "SBS17b",
# 	"Signature.18": "SBS18",
# 	"Signature.20": "SBS20",
#     "Signature.24": "SBS24",
# 	"Signature.26": "SBS26",
#     "Signature.30": "SBS30"
# }

DSBR_SNV_LOAD = "dsbr_snv_load"
DSBR_CT_RATIO = "dsbr_ct_ratio"
DSBR_DEL4_LOAD = "dsbr_del4_load" 
DSBR_DEL4_RATIO = "dsbr_del4_ratio" 
DSBR_DELSV_LOAD =	"dsbr_delsv_load" 
DSBR_DELSV_RATIO = 	"dsbr_delsv_ratio" 
DSBR_DUP_LOAD = "dsbr_dup_load" 
DSBR_SV_LOAD = 	"dsbr_sv_load" 
DSBR_FIRST_GENES =	"dsbr_first_genes"
DSBR_SECOND_GENES = "dsbr_second_genes"
MMR_SNV_LOAD = "mmr_snv_load"
MMR_INDEL_LOAD = "mmr_indel_load"
MMR_FIRST_GENES = "mmr_first_genes"
MMR_SECOND_GENES = "mmr_second_genes"

DSBR_DEFAULT_HALLMARK_CUTOFFS = {
	DSBR_SNV_LOAD : 12000,
	DSBR_CT_RATIO : 0.3,
	DSBR_DEL4_LOAD: 200,
	DSBR_DEL4_RATIO: 0.4,
	DSBR_DELSV_LOAD: 50,
	DSBR_DELSV_RATIO: 0.45,
	DSBR_DUP_LOAD: 75,
	DSBR_SV_LOAD: 200
}

MMR_DEFAULT_HALLMARK_CUTOFFS = {
	MMR_SNV_LOAD : 35000,
	MMR_INDEL_LOAD : 3000
}

DSBR_HTML_HEADERS = {
	DSBR_SNV_LOAD : "SNV Load >",
	DSBR_CT_RATIO: "SNV C>T Ratio <",
	DSBR_DEL4_LOAD: "4bp+ Deletion Load >",
	DSBR_DEL4_RATIO: "4bp+ Deletion Ratio >",
	DSBR_DELSV_LOAD: "100-10kbp Deletion Load >",
	DSBR_DELSV_RATIO: "100-10kbp Deletion Ratio >",
	DSBR_DUP_LOAD: "10k-1mbp Duplication Load >",
	DSBR_SV_LOAD: "Structural Variant Load >",
	DSBR_FIRST_GENES: "First Gene Hit(s)",
	DSBR_SECOND_GENES: "Second Gene Hit(s)"
}


MMR_HTML_HEADERS = {
	MMR_SNV_LOAD : "SNV Load >",
	MMR_INDEL_LOAD : "Indel Load >",
	MMR_FIRST_GENES : "First Gene Hit(s)",
	MMR_SECOND_GENES : "Second Gene Hit(s)"
}

TMB_RANGE_CUTOFF = {
	'high': 10,
	'elevated': 5,
	'typical': 1,
	'low': 0
}

SNV_RANGE_CUTOFF = {
	'high': 20000,
	'elevated': 10000,
	'typical': 2500,
	'low': 0
}

INDEL_RANGE_CUTOFF = {
	'high': 2000,
	'elevated': 1000,
	'typical': 250,
	'low': 0
}

SV_RANGE_CUTOFF = {
	'high': 280,
	'elevated': 140,
	'typical': 35,
	'low': 0
}



#updated 2024.04.03
LIMS_TISSUE_CODES = {
	'Ab':'Abdomen',
	'Ad':'Adipose',
	'Ae':'Adnexa',
	'Ag':'Adrenal',
	'An':'Anus',
	'Ao':'Anorectal',
	'Ap':'Appendix',
	'As':'Ascites',
	'At':'Astrocytoma',
	'Av':'Ampulla',
	'Ax':'Axillary',
	'Ba':'Back',
	'Bd':'Bile',
	'Bi':'Biliary',
	'Bl':'Bladder',
	'Bm':'Bone',
	'Bn':'Brain',
	'Bo':'Bone',
	'Br':'Breast',
	'Bu':'Buccal',
	'Bw':'Bowel',
	'Ca':'Carotid',
	'Cb':'Cord',
	'Cc':'Cecum',
	'Ce':'Cervix',
	'Cf':'Cell-Free',
	'Ch':'Chest',
	'Cj':'Conjunctiva',
	'Ck':'Cheek',
	'Cn':'Central',
	'Co':'Colon',
	'Cr':'Colorectal',
	'Cs':'Cul-de-sac',
	'Ct':'Circulating',
	'Di':'Diaphragm',
	'Du':'Duodenum',
	'Ea':'Ear',
	'En':'Endometrial',
	'Ep':'Epidural',
	'Es':'Esophagus',
	'Ey':'Eye',
	'Fa':'Fallopian',
	'Fb':'Fibroid',
	'Fs':'Foreskin',
	'Ft':'Foot',
	'Ga':'Gastric',
	'Gb':'Gallbladder',
	'Ge':'Gastroesophageal',
	'Gi':'Gastrointestinal',
	'Gj':'Gastrojejunal',
	'Gl':'Glomus',
	'Gn':'Gingiva',
	'Gt':'Genital',
	'Hp':'Hypopharynx',
	'Hr':'Heart',
	'Ic':'ileocecum',
	'Il':'Ileum',
	'Ki':'Kidney',
	'La':'Lacrimal',
	'Lb':'Limb',
	'Le':'Leukocyte',
	'Lg':'Leg',
	'Li':'Large',
	'Ln':'Lymph',
	'Lp':'Lymphoblast',
	'Lt':'Placenta',
	'Lu':'Lung',
	'Lv':'Liver',
	'Lx':'Larynx',
	'Ly':'Lymphocyte',
	'Md':'Mediastinum',
	'Me':'Mesenchyme',
	'Mn':'Mandible',
	'Mo':'Mouth',
	'Ms':'Mesentary',
	'Mu':'Muscle',
	'Mx':'Maxilla',
	'Nk':'Neck',
	'nn':'Unknown',
	'No':'Nose',
	'Np':'Nasopharynx',
	'Oc':'Oral',
	'Om':'Omentum',
	'Or':'Orbit',
	'Ov':'Ovary',
	'Pa':'Pancreas',
	'Pb':'Peripheral',
	'Pc':'Pancreatobiliary',
	'Pd':'Parathyroid',
	'Pe':'Pelvic',
	'Pg':'Parotid',
	'Ph':'Paratracheal',
	'Pi':'Penis',
	'Pl':'Plasma',
	'Pm':'Peritoneum',
	'Pn':'Peripheral',
	'Po':'Peri-aorta',
	'Pr':'Prostate',
	'Pt':'Palate',
	'Pu':'Pleura',
	'Py':'periampullary',
	'Ra':'Right',
	'Rc':'Rectosigmoid',
	'Re':'Rectum',
	'Ri':'Rib',
	'Rp':'Retroperitoneum',
	'Sa':'Saliva',
	'Sb':'Small',
	'Sc':'Scalp',
	'Se':'Serum',
	'Sg':'Salivary',
	'Si':'Small',
	'Sk':'Skin',
	'Sm':'Skeletal',
	'Sn':'Spine',
	'So':'Soft',
	'Sp':'Spleen',
	'Sr':'Serosa',
	'Ss':'Sinus',
	'St':'Stomach',
	'Su':'Sternum',
	'Ta':'Tail',
	'Te':'Testes',
	'Tg':'Thymic',
	'Th':'Thymus',
	'Tn':'Tonsil',
	'To':'Throat',
	'Tr':'Trachea',
	'Tu':'Tongue',
	'Ty':'Thyroid',
	'Uc':'Urachus',
	'Ue':'Ureter',
	'Um':'Umbilical',
	'Up':'Urine',
	'Ur':'Urethra',
	'Us':'Urine',
	'Ut':'Uterus',
	'Uw':'Urine',
	'Vg':'Vagina',
	'Vu':'Vulva',
	'Wm':'Worm'
}
