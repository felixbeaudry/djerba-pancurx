"""
Helper for writing a subset of cardea to the shared workspace
"""
import os
import djerba.core.constants as core_constants
import djerba.util.ini_fields as ini 
from djerba.helpers.base import helper_base
import djerba.plugins.pancurx.constants as phe
from djerba.util.logger import logger
from djerba.util.subprocess_runner import subprocess_runner

class main(helper_base):

    PRIORITY = 10

    def configure(self, config):
        """
        Writes a subset of provenance, and informative JSON files, to the workspace
        """
        config = self.apply_defaults(config)
        wrapper = self.get_config_wrapper(config)
        if wrapper.my_param_is_null('template_type'):
            wrapper.set_my_param('template_type', 'PCX')

        DATA_LOCATION = '/.mounts/labs/PCSI/users/fbeaudry/djerba-pancurx/src/lib/djerba/data/pancurx/'
        if wrapper.my_param_is_null('comparison_cohort_file'):
            file_name = '.'.join( (wrapper.get_my_string('template_type'), 'summary.csv'))
            file_path = os.path.join(DATA_LOCATION, file_name)
            wrapper.set_my_param('comparison_cohort_file', file_path)

        if wrapper.my_param_is_null('genes_of_interest_file'):
            file_name = '.'.join( (wrapper.get_my_string('template_type'), phe.DEFAULT_GENE_FILE))
            file_path = os.path.join(DATA_LOCATION, file_name)
            wrapper.set_my_param('genes_of_interest_file', file_path)

        if wrapper.my_param_is_null('germline_genes_of_interest_file'):
            file_name = '.'.join( (wrapper.get_my_string('template_type'), phe.DEFAULT_GERMLINE_GENE_FILE))
            file_path = os.path.join(DATA_LOCATION, file_name)
            wrapper.set_my_param('germline_genes_of_interest_file', file_path)

        if wrapper.my_param_is_null('immune_genes_of_interest_file'):
            file_name = '.'.join( (wrapper.get_my_string('template_type'), 'immune.genes.txt'))
            file_path = os.path.join(DATA_LOCATION, file_name)
            wrapper.set_my_param('immune_genes_of_interest_file', file_path)

        all_params = {
            phe.DONOR :  wrapper.get_my_string(phe.DONOR),
            phe.TUMOUR_SAMPLE_ID : wrapper.get_my_string(phe.TUMOUR_SAMPLE_ID),
            phe.NORMAL_SAMPLE_ID : wrapper.get_my_string(phe.NORMAL_SAMPLE_ID),
            'comparison_cohort_file': wrapper.get_my_string('comparison_cohort_file'),
            'genes_of_interest_file': wrapper.get_my_string('genes_of_interest_file'),
            'germline_genes_of_interest_file': wrapper.get_my_string('germline_genes_of_interest_file'),
            'immune_genes_of_interest_file': wrapper.get_my_string('immune_genes_of_interest_file'),
            'template_type': wrapper.get_my_string('template_type'),
        }
        self.write_sample_info(all_params)
        return wrapper.get_config()

    def extract(self, config):
        wrapper = self.get_config_wrapper(config)

    def specify_params(self):
        self.logger.debug("Specifying params for sample helper")
        self.set_priority_defaults(self.PRIORITY)
        self.add_ini_required(phe.DONOR)
        self.add_ini_required(phe.TUMOUR_SAMPLE_ID)
        self.add_ini_required(phe.NORMAL_SAMPLE_ID)
        discovered = [
            'comparison_cohort_file',
            'immune_genes_of_interest_file',
            'germline_genes_of_interest_file',
            'genes_of_interest_file',
            'template_type',
        ]
        for key in discovered:
            self.add_ini_discovered(key)

    def write_sample_info(self, path_info):
        self.workspace.write_json(core_constants.DEFAULT_SAMPLE_INFO, path_info)
        self.logger.debug("Wrote path info to workspace: {0}".format(core_constants.DEFAULT_SAMPLE_INFO))
