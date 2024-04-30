"""Plugin to generate the "Patient & Physician Info" report section"""

import logging
from time import strptime
import csv
import os
import re
import djerba.core.constants as core_constants
from djerba.plugins.base import plugin_base
from djerba.util.render_mako import mako_renderer
from djerba.util.image_to_base64 import converter
import djerba.plugins.pancurx.constants as phe
import djerba.plugins.pancurx.tools as tools

class main(plugin_base):

    PRIORITY = 200
    PLUGIN_VERSION = '1.0.0'
    TEMPLATE_NAME = 'classification_template.html'

    def configure(self, config):
        config = self.apply_defaults(config)
        wrapper = self.get_config_wrapper(config)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.PARAM_PATH, phe.PARAM_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.SEX_PATH, phe.SEX_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.TDP_PATH, phe.TDP_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.SUMMARY_FILE_PATH, phe.SUMMARY_FILE_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, 'template_type', 'template_type', core_constants.DEFAULT_SAMPLE_INFO)

        class_params = [
            phe.ALEXANDROV_CLASS,
            phe.COLLISSON_CLASS,
            phe.WADDELL_CLASS,
            phe.MOFFITT_CLASS,
            phe.HLA_TYPES,
        ]
        for key in class_params:
            if wrapper.my_param_is_null(key):
                wrapper.set_my_param(key, "NA")
        return wrapper.get_config()

    def extract(self, config):
        wrapper = self.get_config_wrapper(config)
        attributes = wrapper.get_my_attributes()
        self.check_attributes_known(attributes)
        data = self.get_starting_plugin_data(wrapper, self.PLUGIN_VERSION)
        results_keys = [
            
            phe.ALEXANDROV_CLASS,
            phe.COLLISSON_CLASS,
            phe.MOFFITT_CLASS,
            phe.HLA_TYPES
        ]
        data[core_constants.RESULTS] = {k: wrapper.get_my_string(k) for k in results_keys}
        data[core_constants.RESULTS][phe.INFERRED_SEX] = self.parse_sex(wrapper.get_my_string(phe.SEX_PATH), wrapper.get_my_string('template_type'))
        ploidy = tools.parse_celluloid_params(self, wrapper.get_my_string(phe.PARAM_PATH), phe.PLOIDY)
        data[core_constants.RESULTS][phe.PLOIDY] = ploidy
        data[core_constants.RESULTS]['tandem_duplicator_phenotype_score'] = tools.parse_TDP(self, wrapper.get_my_string(phe.TDP_PATH))
        summary_results = tools.parse_summary_file(self, wrapper.get_my_string(phe.SUMMARY_FILE_PATH))
        dsbr_results, dsbr_tally = tools.parse_multifactor_marker(self, summary_results, phe.DSBR_HTML_HEADERS, phe.DSBR_DEFAULT_HALLMARK_CUTOFFS)
        data[core_constants.RESULTS][phe.DSBR_RESULTS] = dsbr_results
        data[core_constants.RESULTS]['dsbr_tally'] = dsbr_tally
        mmr_results, mmr_tally = tools.parse_multifactor_marker(self, summary_results, phe.MMR_HTML_HEADERS, phe.MMR_DEFAULT_HALLMARK_CUTOFFS)
        data[core_constants.RESULTS][phe.MMR_RESULTS] = mmr_results
        data[core_constants.RESULTS]['mmr_tally'] = mmr_tally
        data[core_constants.RESULTS]['template_type'] = '_'.join((wrapper.get_my_string('template_type'), self.TEMPLATE_NAME))
        #TODO: replace waddell pull by calculation from SV counts
        data[core_constants.RESULTS][phe.WADDELL_CLASS] = summary_results['waddell']
        return data

    def render(self, data):
        renderer = mako_renderer(self.get_module_dir())
        template_name = data[core_constants.RESULTS]['template_type'] 
        return renderer.render_name(template_name, data)

    def specify_params(self):
        discovered = [
            phe.PARAM_PATH,
            phe.SEX_PATH,
            phe.TDP_PATH,
            phe.ALEXANDROV_CLASS,
            phe.COLLISSON_CLASS,
            phe.WADDELL_CLASS,
            phe.MOFFITT_CLASS,
            phe.HLA_TYPES,
            phe.SUMMARY_FILE_PATH,
            'template_type'
        ]
        for key in discovered:
            self.add_ini_discovered(key)
        self.set_ini_default(core_constants.ATTRIBUTES, 'research')
        self.set_priority_defaults(self.PRIORITY)

    def parse_sex(self, sex_file_path, template="PCX"):
        sex_file_path = tools.check_path_exists(self, sex_file_path)
        with open(sex_file_path) as sex_file:
            for row in sex_file:
                short_sex = row.strip()
                if template == "PCX":
                    if short_sex == 'M':
                        inferredSex = "Male"
                    elif short_sex == 'F':
                        inferredSex = "Female"
                    else:
                        msg = "Cannot infer sex {0} from file: {1}".format(short_sex, sex_file_path)
                        self.logger.error(msg)
                        raise RuntimeError(msg)
                else:
                    if short_sex == 'M':
                        inferredSex = "XY"
                    elif short_sex == 'F':
                        inferredSex = "XX"
                    else:
                        msg = "Cannot infer sex {0} from file: {1}".format(short_sex, sex_file_path)
                        self.logger.error(msg)
                        raise RuntimeError(msg)

        return(inferredSex)
            
	