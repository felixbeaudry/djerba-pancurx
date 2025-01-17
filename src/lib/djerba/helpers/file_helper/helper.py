"""
Helper for writing a subset of files to the shared workspace
"""

import os
import csv
import gzip
import logging
import requests
import json
import re

import djerba.core.constants as core_constants
import djerba.util.ini_fields as ini 
from djerba.helpers.base import helper_base
import djerba.plugins.pancurx.constants as phe
from djerba.util.logger import logger
from djerba.util.subprocess_runner import subprocess_runner
import djerba.plugins.pancurx.tools as tools

class main(helper_base):

    PRIORITY = 20


    def configure(self, config):
        """
        Writes a subset of provenance, and informative JSON files, to the workspace
        """
        config = self.apply_defaults(config)
        wrapper = self.get_config_wrapper(config)

        wrapper = tools.fill_file_if_null(self, wrapper, phe.TUMOUR_ID, phe.TUMOUR_ID, core_constants.DEFAULT_SAMPLE_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.DONOR_ID, phe.DONOR_ID, core_constants.DEFAULT_SAMPLE_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.NORMAL_ID, phe.NORMAL_ID, core_constants.DEFAULT_SAMPLE_INFO)

        donor = wrapper.get_my_string(phe.DONOR)
        sample = wrapper.get_my_string(phe.TUMOUR_SAMPLE_ID)
        normal = wrapper.get_my_string(phe.NORMAL_SAMPLE_ID)
        
        if wrapper.my_param_is_null('rootpath'):
            wrapper.set_my_param('rootpath', phe.DEFAULT_ROOTPATH)

        if wrapper.my_param_is_null('envpath'):
            wrapper.set_my_param('envpath', phe.DEFAULT_ENV_PATH)

        if wrapper.my_param_is_null('binpath'):
            wrapper.set_my_param('binpath', phe.DEFAULT_DJERBA_BIN_PATH)

        path_finder = assemble_file_paths(wrapper.get_my_string('rootpath'), donor, sample, normal)

        if wrapper.my_param_is_null(phe.COSMIC_SIGNNLS_PATH):
            wrapper.set_my_param(phe.COSMIC_SIGNNLS_PATH, path_finder.make_cosmic_signals_path())

        if wrapper.my_param_is_null(phe.SEX_PATH):
            wrapper.set_my_param(phe.SEX_PATH, path_finder.make_sex_chromosome_path())

        if wrapper.my_param_is_null(phe.GERMLINE_ANNOVAR_PATH):
            wrapper.set_my_param(phe.GERMLINE_ANNOVAR_PATH, path_finder.make_germline_annovar())

        if wrapper.my_param_is_null(phe.SUMMARY_FILE_PATH):
            wrapper.set_my_param(phe.SUMMARY_FILE_PATH, path_finder.make_sample_summary_path())

        if wrapper.my_param_is_null(phe.SAMPLE_VARIANTS_FILE):
            wrapper.set_my_param(phe.SAMPLE_VARIANTS_FILE, path_finder.make_sample_variants_path())

        if wrapper.my_param_is_null(phe.CELLULOID_PLOT):
            wrapper.set_my_param(phe.CELLULOID_PLOT, path_finder.make_celluloid_plot_path())

        if wrapper.my_param_is_null(phe.TDP_PATH):
            wrapper.set_my_param(phe.TDP_PATH, path_finder.make_tdp_path())

        if wrapper.my_param_is_null(phe.SV_FILE):
            wrapper.set_my_param(phe.SV_FILE, path_finder.make_structural_path())

        DATA_LOCATION = phe.DEFAULT_DATA_LOCATION
        wrapper = tools.fill_file_if_null(self, wrapper, 'template_type', 'template_type', core_constants.DEFAULT_SAMPLE_INFO)

        if wrapper.my_param_is_null('comparison_cohort_file'):
            file_name = '.'.join( (wrapper.get_my_string('template_type'), 'summary.csv'))
            file_path = os.path.join(DATA_LOCATION, file_name)
            wrapper.set_my_param('comparison_cohort_file', file_path)

        if wrapper.my_param_is_null(phe.BLURB_FILE):
            file_path = os.path.join(DATA_LOCATION, phe.DEFAULT_BLURB_FILE)
            wrapper.set_my_param(phe.BLURB_FILE, file_path)

        if wrapper.my_param_is_null('genes_of_interest_file'):
            file_name = '.'.join( (wrapper.get_my_string('template_type'), phe.DEFAULT_GENE_FILE))
            file_path = os.path.join(DATA_LOCATION, file_name)
            wrapper.set_my_param('genes_of_interest_file', file_path)

        if wrapper.my_param_is_null('signatures_of_interest_file'):
            file_name = '.'.join( (wrapper.get_my_string('template_type'), phe.DEFAULT_SIGNATURES_FILE))
            file_path = os.path.join(DATA_LOCATION, file_name)
            wrapper.set_my_param('signatures_of_interest_file', file_path)

        if wrapper.my_param_is_null('sites_of_interest_file'):
            file_name = '.'.join( (wrapper.get_my_string('template_type'), phe.DEFAULT_SITES_FILE))
            file_path = os.path.join(DATA_LOCATION, file_name)
            wrapper.set_my_param('sites_of_interest_file', file_path)

        if wrapper.my_param_is_null('germline_genes_of_interest_file'):
            file_name = '.'.join( (wrapper.get_my_string('template_type'), phe.DEFAULT_GERMLINE_GENE_FILE))
            file_path = os.path.join(DATA_LOCATION, file_name)
            wrapper.set_my_param('germline_genes_of_interest_file', file_path)

        if wrapper.my_param_is_null('immune_genes_of_interest_file'):
            file_name = '.'.join( (wrapper.get_my_string('template_type'), 'immune.genes.txt'))
            file_path = os.path.join(DATA_LOCATION, file_name)
            wrapper.set_my_param('immune_genes_of_interest_file', file_path)

        if wrapper.my_param_is_null('immune_cells_of_interest_file'):
            file_name = '.'.join( (wrapper.get_my_string('template_type'), 'immune.cells.txt'))
            file_path = os.path.join(DATA_LOCATION, file_name)
            wrapper.set_my_param('immune_cells_of_interest_file', file_path)

        if wrapper.my_param_is_null(phe.COVERAGE_PATHS):
            wrapper.set_my_param(phe.COVERAGE_PATHS, path_finder.make_coverage_paths())

        if wrapper.my_param_is_null("bam_paths"):
            wrapper.set_my_param("bam_paths", path_finder.make_bam_paths())

        if wrapper.my_param_is_null(phe.CELLULOID_PARAM_PATH):
            wrapper.set_my_param(phe.CELLULOID_PARAM_PATH, path_finder.make_celluloid_param_path())

        if wrapper.my_param_is_null("seg_path"):
            wrapper.set_my_param("seg_path", path_finder.make_celluloid_seg_path())

        if wrapper.my_param_is_null(phe.CELLULOID_DIR):
            wrapper.set_my_param(phe.CELLULOID_DIR, path_finder.make_celluloid_dir_path())

        if wrapper.my_param_is_null(phe.GERMLINE_PATH):
            wrapper.set_my_param(phe.GERMLINE_PATH, path_finder.make_germline_path())

        if wrapper.my_param_is_null(phe.SNV_PATH):
            wrapper.set_my_param(phe.SNV_PATH, path_finder.make_snv_path())

        if wrapper.my_param_is_null("indel_path"):
            wrapper.set_my_param("indel_path", path_finder.make_indel_path())

        if wrapper.my_param_is_null("snv_annovar"):
            wrapper.set_my_param("snv_annovar", path_finder.make_snv_annovar())

        if wrapper.my_param_is_null("indel_annovar"):
            wrapper.set_my_param("indel_annovar", path_finder.make_indel_annovar())

        if wrapper.my_param_is_null(phe.STAR_QC_PATH):
            wrapper.set_my_param(phe.STAR_QC_PATH, path_finder.make_star_qc())

        if wrapper.my_param_is_null(phe.MAVIS_FUSIONS_PATH):
            wrapper.set_my_param(phe.MAVIS_FUSIONS_PATH, path_finder.make_mavis_fusion())

        if wrapper.my_param_is_null(phe.TPM_PATH):
            wrapper.set_my_param(phe.TPM_PATH, path_finder.make_stringtie())

        if wrapper.my_param_is_null(phe.WHOLE_GENOME_PLOT):
            wrapper.set_my_param(phe.WHOLE_GENOME_PLOT, path_finder.make_whole_genome_plot_path())

        if wrapper.my_param_is_null(phe.WHOLE_CHROMOSOME_PLOTS):
            wrapper.set_my_param(phe.WHOLE_CHROMOSOME_PLOTS, path_finder.make_whole_chromosome_plots_path())

        if wrapper.my_param_is_null('svg_plots'):
            wrapper.set_my_param('svg_plots', path_finder.make_svg_plots())

        all_paths = {

            phe.COSMIC_SIGNNLS_PATH : wrapper.get_my_string(phe.COSMIC_SIGNNLS_PATH),
            phe.SEX_PATH: wrapper.get_my_string(phe.SEX_PATH),
            phe.GERMLINE_ANNOVAR_PATH: wrapper.get_my_string(phe.GERMLINE_ANNOVAR_PATH),
            phe.SUMMARY_FILE_PATH: wrapper.get_my_string(phe.SUMMARY_FILE_PATH),
            phe.SAMPLE_VARIANTS_FILE: wrapper.get_my_string(phe.SAMPLE_VARIANTS_FILE),

            phe.CELLULOID_PLOT: wrapper.get_my_string(phe.CELLULOID_PLOT),
            phe.TDP_PATH : wrapper.get_my_string(phe.TDP_PATH),
            phe.SV_FILE : wrapper.get_my_string(phe.SV_FILE),

            'comparison_cohort_file': wrapper.get_my_string('comparison_cohort_file'),
            'genes_of_interest_file': wrapper.get_my_string('genes_of_interest_file'),
            'germline_genes_of_interest_file': wrapper.get_my_string('germline_genes_of_interest_file'),
            'immune_genes_of_interest_file': wrapper.get_my_string('immune_genes_of_interest_file'),
            'signatures_of_interest_file': wrapper.get_my_string('signatures_of_interest_file'),
            'sites_of_interest_file': wrapper.get_my_string('sites_of_interest_file'),
            'immune_cells_of_interest_file': wrapper.get_my_string('immune_cells_of_interest_file'),
            phe.BLURB_FILE: wrapper.get_my_string(phe.BLURB_FILE),

            phe.COVERAGE_PATHS: wrapper.get_my_string(phe.COVERAGE_PATHS) ,
            "bam_paths": wrapper.get_my_string("bam_paths")  ,
            
            phe.CELLULOID_PARAM_PATH : wrapper.get_my_string(phe.CELLULOID_PARAM_PATH) ,
            "seg_path" : wrapper.get_my_string("seg_path") ,
            phe.CELLULOID_DIR: wrapper.get_my_string(phe.CELLULOID_DIR) ,

            phe.GERMLINE_PATH : wrapper.get_my_string(phe.GERMLINE_PATH) ,
            phe.SNV_PATH : wrapper.get_my_string(phe.SNV_PATH) ,
            "indel_path" : wrapper.get_my_string("indel_path") ,
            "snv_annovar": wrapper.get_my_string("snv_annovar") ,
            "indel_annovar": wrapper.get_my_string("indel_annovar") ,

            phe.STAR_QC_PATH: wrapper.get_my_string(phe.STAR_QC_PATH),
            phe.MAVIS_FUSIONS_PATH: wrapper.get_my_string(phe.MAVIS_FUSIONS_PATH) ,
            phe.TPM_PATH: wrapper.get_my_string(phe.TPM_PATH) ,

            # plots
            phe.WHOLE_GENOME_PLOT: wrapper.get_my_string(phe.WHOLE_GENOME_PLOT) ,
            phe.WHOLE_CHROMOSOME_PLOTS: wrapper.get_my_string(phe.WHOLE_CHROMOSOME_PLOTS) ,

            'svg_plots': wrapper.get_my_string('svg_plots') ,
            'envpath': wrapper.get_my_string('envpath') ,
            'binpath': wrapper.get_my_string('binpath') ,

        }
        self.write_path_info(all_paths)
        return wrapper.get_config()

    def extract(self, config):
        wrapper = self.get_config_wrapper(config)

    def specify_params(self):
        self.logger.debug("Specifying params for provenance helper")
        self.set_priority_defaults(self.PRIORITY)
        discovered = [
            phe.TUMOUR_ID,
            phe.DONOR_ID,
            phe.NORMAL_ID,
            phe.COSMIC_SIGNNLS_PATH,
            phe.SEX_PATH  ,
            phe.GERMLINE_ANNOVAR_PATH,
            phe.SUMMARY_FILE_PATH,
            phe.SAMPLE_VARIANTS_FILE,
            phe.BLURB_FILE,
            'comparison_cohort_file',
            'immune_genes_of_interest_file',
            'germline_genes_of_interest_file',
            'genes_of_interest_file',
            'signatures_of_interest_file',
            'sites_of_interest_file',
            'template_type',
            phe.CELLULOID_PLOT,
            phe.TDP_PATH,
            phe.SV_FILE,
            phe.COVERAGE_PATHS,
            "bam_paths",
            phe.CELLULOID_PARAM_PATH,
            "seg_path",
            phe.CELLULOID_DIR,
            phe.GERMLINE_PATH,
            phe.SNV_PATH,
            "indel_path" ,
            "snv_annovar",
            "indel_annovar",
            phe.STAR_QC_PATH,
            phe.MAVIS_FUSIONS_PATH,
            phe.TPM_PATH,
            phe.WHOLE_GENOME_PLOT,
            phe.WHOLE_CHROMOSOME_PLOTS,
            'svg_plots',
            'rootpath',
            'binpath',
            'envpath',
            'immune_cells_of_interest_file'
        ]
        for key in discovered:
            self.add_ini_discovered(key)

    def write_path_info(self, path_info):
        self.workspace.write_json(core_constants.DEFAULT_PATH_INFO, path_info)
        self.logger.debug("Wrote path info to workspace: {0}".format(core_constants.DEFAULT_PATH_INFO))

class assemble_file_paths(logger):
    
    
    def __init__(self, rootpath, donor, sample, normal_sample, log_level=logging.WARNING, log_path=None):
        self.log_level = log_level
        self.log_path = log_path
        self.logger = self.get_logger(log_level, __name__, log_path)
        self.donor = donor
        self.sample = sample
        self.normal_sample = normal_sample
        seq_path = os.path.join(phe.DEFAULT_SEQTYPE, phe.DEFAULT_ALIGNER, phe.DEFAULT_ALIGNER_VERSION)
        rna_seq_path = os.path.join(phe.DEFAULT_RNA_SEQTYPE, phe.DEFAULT_RNA_ALIGNER, phe.DEFAULT_RNA_ALIGNER_VERSION)

        self.root_extended = os.path.join(rootpath, donor, sample, seq_path) 
        self.rna_root_extended = os.path.join(rootpath, donor, sample, rna_seq_path)
        self.normal_root_extended = os.path.join(rootpath, donor, normal_sample, seq_path)

        self.archive_extended = os.path.join(phe.DEFAULT_ARCHIVEPATH, donor, sample, seq_path)
        self.normal_archive_extended = os.path.join(phe.DEFAULT_ARCHIVEPATH, donor, normal_sample, seq_path)

    def make_stringtie(self):
        file_name = "".join((self.sample,"_stringtie_abundance.txt"))
        file_path = os.path.join(self.rna_root_extended, "stringtie/2.0.6/", file_name)
        return(file_path)

    def make_mavis_fusion(self):
        file_name = "".join(("mavis_summary_all_",self.sample,".tab"))
        file_path = os.path.join(self.rna_root_extended, "starfusion/1.9.0/mavis/summary", file_name)
        return(file_path)

    def make_star_qc(self):
        file_name = "".join(("Log.final.out"))
        file_path = os.path.join(self.rna_root_extended, "collapsed", file_name)
        return(file_path)
       
    def make_bam_paths(self):
        file_name = "".join((self.sample, '.bam'))
        normal_file_name = "".join((self.normal_sample, '.bam'))
        file_paths = {
            "tumour": os.path.join(self.archive_extended, "collapsed", file_name),
            "normal": os.path.join(self.normal_archive_extended, "collapsed", normal_file_name)
        }
        return(file_paths)

    def make_celluloid_param_path(self):
        file_name = "".join(("parameters_",self.sample,".txt"))
        file_path = os.path.join(self.root_extended, "celluloidXY", phe.DEFAULT_CELLULOID_VERSION, "solution", file_name)
        return(file_path)

    def make_celluloid_seg_path(self):
        file_name = "".join(("segments_",self.sample,".txt.sorted.bed.gz"))
        file_path = os.path.join(self.root_extended, "celluloidXY", phe.DEFAULT_CELLULOID_VERSION, "solution", file_name)
        return(file_path)

    def make_celluloid_dir_path(self):
        file_path = os.path.join(self.archive_extended, "celluloidXY", phe.DEFAULT_CELLULOID_VERSION)
        return(file_path)

    def make_cosmic_signals_path(self):
        file_name = "".join((self.sample,"_SBS_signatures.txt"))
        file_path = os.path.join(self.root_extended, "cosmicSigNNLS", file_name)
        return(file_path)

    def make_coverage_paths(self):
        FILE_EXTENSION = "_coverage_collapsed.metrics"
        file_name = "".join((self.sample, FILE_EXTENSION))
        normal_file_name = "".join((self.normal_sample, FILE_EXTENSION))
        file_paths = {
            "tumour": os.path.join(self.root_extended, "coverage", file_name),
            "normal": os.path.join(self.normal_root_extended, "coverage", normal_file_name)
        }
        return(file_paths)

    def make_germline_annovar_old(self):
        #deprecated by snakemake
        file_name = "".join((self.sample,"_germline_final.hg38_multianno.txt.gz"))
        file_path = os.path.join(self.root_extended, "final_germline_variants/annovar/20170716", file_name)
        snv_annovar = {
            "data": file_path,
            "index": "".join((file_path,".tbi"))
        }
        return(file_path)

    def make_germline_annovar(self):
        file_name = "".join((self.sample,"_germline_final.hg38_multianno.txt.gz"))
        file_path = os.path.join(self.root_extended, "final_germline_variants", file_name)
        snv_annovar = {
            "data": file_path,
            "index": "".join((file_path,".tbi"))
        }
        return(file_path)

    def make_germline_path(self):
        file_name = "".join((self.sample,"_germline_final.vcf.gz"))
        file_path = os.path.join(self.root_extended, "final_germline_variants", file_name)
        return(file_path)

    def make_germline_path(self):
        file_name = "".join((self.sample,"_germline_final.vcf.gz"))
        file_path = os.path.join(self.root_extended, "final_germline_variants", file_name)
        return(file_path)

    def make_hrdetect_path(self):
        file_name = "".join((self.sample,".hrdetect_score.csv"))
        file_path = os.path.join(self.root_extended, "hrdetect", file_name)
        return(file_path)

    def make_indel_annovar(self):
        file_name = "".join((self.sample,"_merged_final_indels.hg38_multianno.txt.gz"))
        file_path = os.path.join(self.root_extended, "final_indel/annovar/20170716", file_name)
        snv_annovar = {
            "data": file_path,
            "index": "".join((file_path,".tbi"))
        }
        return(snv_annovar)

    def make_indel_path(self):
        file_name = "".join((self.sample,"_merged_final_indels.vcf.gz"))
        file_path = os.path.join(self.root_extended, "final_indel", file_name)
        return(file_path)

    def make_sex_chromosome_path_old(self):
        # univa pipeline location, deprecated
        file_name = "".join((self.donor,".final_gender.txt"))
        file_path = os.path.join(self.root_extended, "gendertype/final", file_name)
        return(file_path)

    def make_sex_chromosome_path(self):
        file_name = "".join((self.sample,".xyr_gender.txt"))
        file_path = os.path.join(self.root_extended, "gendertype", file_name)
        return(file_path)

    def make_snv_path(self):
        file_name = "".join((self.sample,"_merged_final.vcf.gz"))
        file_path = os.path.join(self.root_extended, "final_strelka2_mutect2", file_name)
        return(file_path)

    def make_structural_path(self):
        file_name = "".join((self.sample,"_finalannotatedSV.txt"))
        file_path = os.path.join(self.root_extended, "mergeSV", file_name)
        return(file_path)

    def make_tdp_path(self):
        file_name = "".join((self.sample,".tdp_score.csv"))
        file_path = os.path.join(self.root_extended, "tdp_score", file_name)
        return(file_path)

    def make_snv_annovar(self):
        file_name = "".join((self.sample,"_merged_final.hg38_multianno.txt.gz"))
        file_path = os.path.join(self.root_extended, "final_strelka2_mutect2/annovar/20170716", file_name)
        snv_annovar = {
            "data": file_path,
            "index": "".join((file_path,".tbi"))
        }
        return(snv_annovar)

    def make_sample_summary_path(self):
        file_name = "".join((self.sample,".summary.csv"))
        file_path = os.path.join(self.root_extended, "results", file_name)
        return(file_path)

    def make_sample_variants_path(self):
        file_name = "".join((self.sample,".variants.csv"))
        file_path = os.path.join(self.root_extended, "results", file_name)
        return(file_path)

    ## COPY OF CELLULOID PLOT

    def make_celluloid_plot_path(self):
        file_name = "".join((self.sample,"-celluloid_contour.png"))
        file_path = os.path.join(self.root_extended, "results/plots", file_name)
        return(file_path)

    ## PLOTTING PATHS

    def make_svg_plots(self):
        plot_dict = {}
        plot_names = ['drivers',
                      'snv.bar_count',
                      'snv.histbox_count',
                      'snv.vaf',
                      'indel.bar_count',
                      'indel.histbox_count',
                      'indel.vaf',
                      'indel.stack_count',
                      'sv.bar_count',
                      'sv.histbox_count',
                      'sv.bins',
                      'sv.stack_count',
                      'SBS_context_bar',
                      'COSMIC_stack',
                      'stack-tmb',
                      'oncoplot.SVs',
                      ]
       
        core_file_path = os.path.join(self.root_extended, "results/djerba_plots")

        for this_plot_name in plot_names:

            file_name = ".".join((self.sample,this_plot_name,"svg"))
            file_path = os.path.join(core_file_path, file_name)
            plot_dict[this_plot_name] = file_path

        return(plot_dict)

    def make_whole_genome_plot_path(self):
        file_name = "".join((self.sample,".whole_genome.png"))
        file_path = os.path.join(self.root_extended, "results/djerba_plots")
        file_path = os.path.join(file_path, file_name)
        return(file_path)

    def make_whole_chromosome_plots_path(self):
        whole_chromosome_plots_paths = {}
        core_file_path = os.path.join(self.root_extended, "results/djerba_plots")
        for this_chromosome in phe.CHROMOSOME_PLOTS.keys():
            file_name = "".join((self.sample,".whole_",this_chromosome,".png"))
            file_path = os.path.join(core_file_path, file_name)
            whole_chromosome_plots_paths[phe.CHROMOSOME_PLOTS[this_chromosome]] = file_path
        return(whole_chromosome_plots_paths)

