"""Plugin to generate the "Gremline" report section"""

import logging
from time import strptime
import csv
import os
import re
import collections
import djerba.core.constants as core_constants
from djerba.plugins.base import plugin_base
from djerba.util.render_mako import mako_renderer
import djerba.plugins.pancurx.constants as phe
import djerba.plugins.pancurx.tools as tools
from djerba.util.environment import directory_finder
class main(plugin_base):

    PRIORITY = 125
    PLUGIN_VERSION = '1.0.0'
    TEMPLATE_NAME = 'fusions_template.html'


    def configure(self, config):
        config = self.apply_defaults(config)
        wrapper = self.get_config_wrapper(config)
        wrapper = tools.fill_file_if_null(self, wrapper, 'genes_of_interest_file', 'genes_of_interest_file', core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.MAVIS_FUSIONS_PATH, phe.MAVIS_FUSIONS_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.STAR_QC_PATH, phe.STAR_QC_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.TPM_PATH, phe.TPM_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.PARAM_PATH, phe.PARAM_PATH, core_constants.DEFAULT_PATH_INFO)

        return wrapper.get_config()

    def extract(self, config):
        wrapper = self.get_config_wrapper(config)
        attributes = wrapper.get_my_attributes()
        self.check_attributes_known(attributes)
        data = self.get_starting_plugin_data(wrapper, self.PLUGIN_VERSION)
        all_fusions = tools.parse_fusions(self, wrapper.get_my_string(phe.MAVIS_FUSIONS_PATH))
        genes_of_interest = tools.get_genes_of_interest(self, wrapper.get_my_string('genes_of_interest_file'))


        cnvs_and_abs = self.workspace.read_json('cnvs_and_abs.json')
        data[core_constants.RESULTS]['gene_expression'] = tools.get_gene_expression(self, genes_of_interest, wrapper.get_my_string(phe.TPM_PATH),  phe.DEFAULT_CIBERSORT_COMPARISON_PATH)

        fusions, fusion_count = tools.filter_fusions(self, all_fusions, genes_of_interest, cnvs_and_abs, data[core_constants.RESULTS]['gene_expression'])
        data[core_constants.RESULTS]['fusions'] = fusions
        data[core_constants.RESULTS]['fusion_count'] = fusion_count
        data[core_constants.RESULTS]['star_qc'] = tools.parse_star_qc(self, wrapper.get_my_string(phe.STAR_QC_PATH))

        ploidy = tools.parse_celluloid_params(self, wrapper.get_my_string(phe.PARAM_PATH), "ploidy_numeric")
        data[core_constants.RESULTS][phe.PLOIDY] = ploidy
        return data

    def render(self, data):
        renderer = mako_renderer(self.get_module_dir())
        return renderer.render_name(self.TEMPLATE_NAME, data)

    def specify_params(self):
        discovered = [
            phe.MAVIS_FUSIONS_PATH,
            phe.STAR_QC_PATH,
            phe.TPM_PATH,
            phe.PARAM_PATH,
            'genes_of_interest_file',
        ]
        for key in discovered:
            self.add_ini_discovered(key)
        self.set_ini_default(core_constants.ATTRIBUTES, 'research')
        self.set_priority_defaults(self.PRIORITY)

