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
import djerba.plugins.pancurx.tools as tools

class main(helper_base):

    PRIORITY = 10

    def configure(self, config):
        """
        Writes a subset of provenance, and informative JSON files, to the workspace
        """
        config = self.apply_defaults(config)
        wrapper = self.get_config_wrapper(config)

        if wrapper.my_param_is_null(phe.EXTERNAL_IDS):
            if self.workspace.has_file(core_constants.DEFAULT_SAMPLE_INFO):
                wrapper = tools.fill_file_if_null(self, wrapper, phe.EXTERNAL_IDS, phe.EXTERNAL_IDS, core_constants.DEFAULT_SAMPLE_INFO)
            else:
                external_ids = tools.parse_lims(self, wrapper.get_my_string(phe.DONOR))
                wrapper.set_my_param(phe.EXTERNAL_IDS, external_ids)

        if self.workspace.has_file(core_constants.DEFAULT_SAMPLE_INFO):
            self.logger.debug("sample info {0} exists already".format(core_constants.DEFAULT_SAMPLE_INFO))
        else:
            all_params = {
                phe.DONOR :  wrapper.get_my_string(phe.DONOR),
                phe.TUMOUR_SAMPLE_ID : wrapper.get_my_string(phe.TUMOUR_SAMPLE_ID),
                phe.NORMAL_SAMPLE_ID : wrapper.get_my_string(phe.NORMAL_SAMPLE_ID),
                phe.EXTERNAL_IDS : wrapper.get_my_string(phe.EXTERNAL_IDS),

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
        self.add_ini_required('template_type')
        discovered = [
            phe.EXTERNAL_IDS,
 
        ]
        for key in discovered:
            self.add_ini_discovered(key)

    def write_sample_info(self,sample_info):
        self.workspace.write_json(core_constants.DEFAULT_SAMPLE_INFO, sample_info)
        self.logger.debug("Wrote sample info to workspace: {0}".format(core_constants.DEFAULT_SAMPLE_INFO))
