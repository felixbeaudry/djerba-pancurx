PLOIDY = "ploidy"
ALEXANDROV_CLASS = "alexandrov_class"
COLLISSON_CLASS = "collisson_class"
WADDELL_CLASS = "waddell_class"
MOFFITT_CLASS = "moffitt_class"
INFERRED_SEX = "inferred_sex"
HLA_TYPES = "hla_types"
DIPLOID_CUTOFF = 2.3
DSBR_SCORE_BAR = "dsbr_score_image_path"
MMR_SCORE_BAR = "mmr_score_image_path"
DSBR_RESULTS = "DSBR_results"
MMR_RESULTS = "MMR_results"
REPORTING_NAME = "reporting_name"
VALUE = "value"
ABOVE_CUTOFF = "above_cutoff"

COSMIC_SIGNATURE_SET = {
	"Signature.1": "csnnls_sig1",
	"Signature.2": "csnnls_sig2",
	"Signature.3": "csnnls_sig3",
	"Signature.5": "csnnls_sig5",
	"Signature.6": "csnnls_sig6",
	"Signature.8": "csnnls_sig8",
	"Signature.13": "csnnls_sig13",
	"Signature.17": "csnnls_sig17",
	"Signature.18": "csnnls_sig18",
	"Signature.20": "csnnls_sig20",
	"Signature.26": "csnnls_sig26",
	"Residuals": "csnnls_residuals",
	"N.Mutations": "csnnls_n.mutations",
}

DSBR_GENES = ["BRCA1", "BRCA2", "PALB2", "RAD51", "RAD51B", "RAD51C", "RAD51D", "XRCC2", "XRCC3", "DMC1", "BRIP1"]

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

MMR_GENES = ["MLH1", "MSH2", "MSH6", "PMS2", "EPCAM"]



MMR_SNV_LOAD = "mmr_snv_load"
MMR_INDEL_LOAD = "mmr_indel_load"
MMR_FIRST_GENES = "mmr_first_genes"
MMR_SECOND_GENES = "mmr_second_genes"

MMR_DEFAULT_HALLMARK_CUTOFFS = {
	MMR_SNV_LOAD : 35000,
	MMR_INDEL_LOAD : 3000
}

MMR_HTML_HEADERS = {
	MMR_SNV_LOAD : "SNV Load >",
	MMR_INDEL_LOAD : "Indel Load >",
	MMR_FIRST_GENES : "First Gene Hit(s)",
	MMR_SECOND_GENES : "Second Gene Hit(s)"
}
