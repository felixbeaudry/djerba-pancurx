"""
Helper for writing a subset of cardea to the shared workspace
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
import djerba.helpers.phoenix_helper.constants as phe
from djerba.helpers.phoenix_helper.phoenix_tools import phoenix_processor as phoenix_tools
from djerba.util.logger import logger
from djerba.util.subprocess_runner import subprocess_runner

class main(helper_base):

    PRIORITY = 20
    CELLULOID_PARAM_PATH = "param_path"
    SEX_PATH = "sex_path"
    TDP_PATH = "tdp_path"
    COSMIC_SIGNNLS_PATH = "cosmic_signnls_path"
    GERMLINE_PATH = "germline_path"
    GERMLINE_ANNOVAR_PATH = "germline_annovar_path"
    ONCOSLIDE_PLOT = "oncoslide_plot"
    SUMMARY_FILE_PATH = "summary_file_path"
    DSBR_SCORE_BAR = "dsbr_score_image_path"
    MMR_SCORE_BAR = "mmr_score_image_path"
    SAMPLE_VARIANTS_FILE = "sample_variants"

    def configure(self, config):
        """
        Writes a subset of provenance, and informative JSON files, to the workspace
        """
        config = self.apply_defaults(config)
        wrapper = self.get_config_wrapper(config)
        donor = wrapper.get_my_string(phe.DONOR)
        sample = wrapper.get_my_string(phe.TUMOUR_SAMPLE_ID)
        normal = wrapper.get_my_string(phe.NORMAL_SAMPLE_ID)
        seqtype = wrapper.get_my_string(phe.SEQTYPE)
        aligner = wrapper.get_my_string(phe.ALIGNER)
        aligner_version = wrapper.get_my_string(phe.ALIGNER_VERSION)
        SEQ_PATH = os.path.join(seqtype, aligner, aligner_version)
        path_finder = assemble_file_paths(SEQ_PATH, donor, sample, normal)
        all_paths = {

            self.COSMIC_SIGNNLS_PATH : path_finder.make_cosmic_signals_path(),
            "coveragePaths": path_finder.make_coverage_paths(),
            "indelPath" : path_finder.make_indel_path(),
            self.CELLULOID_PARAM_PATH : path_finder.make_celluloid_param_path(),
            "snvPath" : path_finder.make_snv_path(),
            self.GERMLINE_PATH : path_finder.make_germline_path(),
            "svPath" : path_finder.make_structural_path(),
            "segPath" : path_finder.make_celluloid_seg_path(),
            self.SEX_PATH: path_finder.make_sex_chromosome_path(),
            self.TDP_PATH : path_finder.make_tdp_path(),
            "snv_annovar": path_finder.make_snv_annovar(),
            "indel_annovar": path_finder.make_indel_annovar(),
            self.GERMLINE_ANNOVAR_PATH: path_finder.make_germline_annovar(),

            # eventually the summary and figure scripts will be ported into djerba
            # and djerba won't need to look for these files

            self.SUMMARY_FILE_PATH: path_finder.make_sample_summary_path(),
            self.SAMPLE_VARIANTS_FILE: path_finder.make_sample_variants_path(),

            self.ONCOSLIDE_PLOT: path_finder.make_oncoslice_path(),
            self.DSBR_SCORE_BAR: path_finder.make_dsbr_score_path(),
            self.MMR_SCORE_BAR: path_finder.make_mmr_score_path()

        }
        self.write_path_info(all_paths)
        return wrapper.get_config()

    def extract(self, config):
        wrapper = self.get_config_wrapper(config)
        sample_data = {
            phe.DONOR : wrapper.get_my_string(phe.DONOR),
            phe.TUMOUR_SAMPLE_ID : wrapper.get_my_string(phe.TUMOUR_SAMPLE_ID),
            phe.NORMAL_SAMPLE_ID : wrapper.get_my_string(phe.NORMAL_SAMPLE_ID),
            'seq_type' : wrapper.get_my_string(phe.SEQTYPE),
        }
        # self.parse_germline(sample_data)

    def specify_params(self):
        self.logger.debug("Specifying params for provenance helper")
        self.set_priority_defaults(self.PRIORITY)
        self.set_ini_default(phe.SEQTYPE, phe.DEFAULT_SEQTYPE)
        self.set_ini_default(phe.ALIGNER, phe.DEFAULT_ALIGNER)
        self.set_ini_default(phe.ALIGNER_VERSION, phe.DEFAULT_ALIGNER_VERSION)
        self.add_ini_required(phe.DONOR)
        self.add_ini_required(phe.TUMOUR_SAMPLE_ID)
        self.add_ini_required(phe.NORMAL_SAMPLE_ID)

    def write_path_info(self, path_info):
        self.workspace.write_json(core_constants.DEFAULT_PATH_INFO, path_info)
        self.logger.debug("Wrote path info to workspace: {0}".format(core_constants.DEFAULT_PATH_INFO))

class assemble_file_paths(logger):

    CELLULOID_VERSION = "v0.11.7"

    def __init__(self, SEQ_PATH, donor, sample, normal_sample, log_level=logging.WARNING, log_path=None):
        self.log_level = log_level
        self.log_path = log_path
        self.logger = self.get_logger(log_level, __name__, log_path)
        self.donor = donor
        self.sample = sample
        self.normal_sample = normal_sample
        self.root_extended = os.path.join(phe.DEFAULT_ROOTPATH, donor, sample, SEQ_PATH)
        self.normal_root_extended = os.path.join(phe.DEFAULT_ROOTPATH, donor, normal_sample, SEQ_PATH)

    def make_celluloid_param_path(self):
        file_name = "".join(("parameters_",self.sample,".txt"))
        file_path = os.path.join(self.root_extended, "celluloidXY", self.CELLULOID_VERSION, "solution", file_name)
        return(file_path)

    def make_celluloid_seg_path(self):
        file_name = "".join(("segments_",self.sample,".txt.sorted.bed.gz"))
        file_path = os.path.join(self.root_extended, "celluloidXY", self.CELLULOID_VERSION, "solution", file_name)
        return(file_path)

    def make_cosmic_signals_path(self):
        file_name = "".join((self.sample,"_signatures.txt"))
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

    def make_dsbr_score_path(self):
        file_name = "".join((self.sample,"-DSBR_score_bar.png"))
        file_path = os.path.join(self.root_extended, "results/plots", file_name)
        return(file_path)
    
    def make_mmr_score_path(self):
        file_name = "".join((self.sample,"-MMR_score_bar.png"))
        file_path = os.path.join(self.root_extended, "results/plots", file_name)
        return(file_path)

    def make_germline_annovar(self):
        file_name = "".join((self.sample,"_germline_final.hg38_multianno.txt.gz"))
        file_path = os.path.join(self.root_extended, "final_germline_variants/annovar/20170716", file_name)
        snv_annovar = {
            "data": file_path,
            "index": "".join((file_path,".tbi"))
        }
        return(snv_annovar)

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

    def make_sex_chromosome_path(self):
        file_name = "".join((self.donor,".final_gender.txt"))
        file_path = os.path.join(self.root_extended, "gendertype/final", file_name)
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

    def make_oncoslice_path(self):
        file_name = "".join((self.sample,"-onco_slice-1000_no_name.png"))
        file_path = os.path.join(self.root_extended, "results/plots", file_name)
        return(file_path)
 