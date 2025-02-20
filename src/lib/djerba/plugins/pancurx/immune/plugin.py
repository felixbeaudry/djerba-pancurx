"""Plugin to generate the "Gremline" report section"""

import logging
from time import strptime
import csv
import os
import re
import collections
import subprocess
import djerba.core.constants as core_constants
from djerba.plugins.base import plugin_base
from djerba.util.render_mako import mako_renderer
import djerba.plugins.pancurx.constants as phe
import djerba.plugins.pancurx.tools as tools
from djerba.util.subprocess_runner import subprocess_runner
from djerba.util.environment import directory_finder

class main(plugin_base):

    PRIORITY = 650
    PLUGIN_VERSION = '1.0.0'
    TEMPLATE_NAME = 'immune_template.html'


    def configure(self, config):
        config = self.apply_defaults(config)
        wrapper = self.get_config_wrapper(config)
        wrapper = tools.fill_file_if_null(self, wrapper, 'immune_genes_of_interest_file', 'immune_genes_of_interest_file', core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, 'immune_cells_of_interest_file', 'immune_cells_of_interest_file', core_constants.DEFAULT_PATH_INFO)

        wrapper = tools.fill_file_if_null(self, wrapper, phe.SAMPLE_VARIANTS_FILE, phe.SAMPLE_VARIANTS_FILE, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.PARAM_PATH, phe.PARAM_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.PARAM_PATH, phe.PARAM_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.TPM_PATH, phe.TPM_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, 'binpath', 'binpath', core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, 'envpath', 'envpath', core_constants.DEFAULT_PATH_INFO)

        return wrapper.get_config()

    def extract(self, config):
        wrapper = self.get_config_wrapper(config)
        attributes = wrapper.get_my_attributes()
        self.check_attributes_known(attributes)
        envpath = wrapper.get_my_string('envpath')
        binpath = wrapper.get_my_string('binpath')


        data = self.get_starting_plugin_data(wrapper, self.PLUGIN_VERSION)
        genes_of_interest = tools.get_genes_of_interest(self, wrapper.get_my_string('immune_genes_of_interest_file'))

        somatic_variants = tools.parse_somatic_variants(self, wrapper.get_my_string(phe.SAMPLE_VARIANTS_FILE), genes_of_interest)

        mane_transcript_path = os.path.join(phe.DEFAULT_DATA_LOCATION, phe.DEFAULT_MANE_FILE)
        mane_transcripts = tools.parse_mane_transcript(self, mane_transcript_path)

        data[core_constants.RESULTS]['reportable_variants'] = tools.get_subset_of_somatic_variants(self, somatic_variants, genes_of_interest, mane_transcripts)
        data[core_constants.RESULTS]['gene_expression'] = tools.get_gene_expression(self, genes_of_interest, wrapper.get_my_string(phe.TPM_PATH),  phe.DEFAULT_CIBERSORT_COMPARISON_PATH, envpath, binpath)
        ploidy = tools.parse_celluloid_params(self, wrapper.get_my_string(phe.PARAM_PATH), "ploidy_numeric")
        data[core_constants.RESULTS][phe.PLOIDY] = ploidy

        cnvs_and_abs = self.workspace.read_json('cnvs_and_abs.json')
        data[core_constants.RESULTS]['immune_markers'] = tools.filter_immune(self,  genes_of_interest, cnvs_and_abs, data[core_constants.RESULTS]['gene_expression'])

        work_dir = self.workspace.get_work_dir()
        self.run_immunedeconv( work_dir, binpath, envpath ) 

        ##parse_cosmic_signatures
        data[core_constants.RESULTS]['immune_cells'] = tools.parse_immune_cells(self, wrapper.get_my_string('immune_cells_of_interest_file'))

        return data

    def render(self, data):
        renderer = mako_renderer(self.get_module_dir())
        return renderer.render_name(self.TEMPLATE_NAME, data)

    def specify_params(self):
        discovered = [
            phe.PARAM_PATH,
            phe.SAMPLE_VARIANTS_FILE,
            'immune_genes_of_interest_file',
            phe.TPM_PATH,
            'binpath',
            'envpath',
            'immune_cells_of_interest_file'
        ]
        for key in discovered:
            self.add_ini_discovered(key)
        self.set_ini_default(core_constants.ATTRIBUTES, 'research')
        self.set_priority_defaults(self.PRIORITY)

    def run_immunedeconv(self,  work_dir, binpath, envpath):
        if self.workspace.has_file('immunedeconv.txt'):
            msg = "CIBERSORT file found"
            self.logger.info(msg)
        else:
            
            r_script = os.path.join(binpath, "pancurx/immune/immune.R")
            r_command = os.path.join(envpath, "lib/R/bin/Rscript")

            cmd = [
                r_command, r_script ,
                '--cibersort', os.path.join(binpath, "pancurx/immune/"),
                '--loadings',  phe.DEFAULT_CIBERSORT_LOADINGS_PATH,
                '--genes',  phe.DEFAULT_CIBERSORT_GENES_PATH,
                '--comparison',  phe.DEFAULT_CIBERSORT_COMPARISON_PATH,
                '--input', os.path.join(work_dir, 'tpm.txt'),
                '--output', os.path.join(work_dir, 'immunedeconv.txt'),
            ]

            runner = subprocess_runner()
            runner.run(cmd, "immunedeconv R script")
