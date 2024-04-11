"""Plugin to generate slide"""

import logging
from time import strftime
import re
import json
import os
import csv
import subprocess
import djerba.core.constants as core_constants
from djerba.plugins.base import plugin_base
from djerba.util.render_mako import mako_renderer
import djerba.plugins.pancurx.constants as phe
import djerba.plugins.pancurx.tools as tools

class main(plugin_base):

    PRIORITY = 1000
    PLUGIN_VERSION = '1.0.0'
    TEMPLATE_NAME = 'slide_template.html'

    def configure(self, config):
        config = self.apply_defaults(config)
        wrapper = self.get_config_wrapper(config)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.CELLULOID_PLOT, phe.CELLULOID_PLOT, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.PARAM_PATH, phe.PARAM_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.TDP_PATH, phe.TDP_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.COSMIC_SIGNNLS_PATH, phe.COSMIC_SIGNNLS_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.SUMMARY_FILE_PATH, phe.SUMMARY_FILE_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.INDEL_BIN_PLOT, phe.INDEL_BIN_PLOT, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.SV_BIN_PLOT, phe.SV_BIN_PLOT, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.SNV_CONTEXT_PLOT, phe.SNV_CONTEXT_PLOT, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.COSMIC_STACK_PLOT, phe.COSMIC_STACK_PLOT, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper,phe.SV_STACK_PLOT, phe.SV_STACK_PLOT, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper,phe.TMB_STACK_PLOT, phe.TMB_STACK_PLOT, core_constants.DEFAULT_PATH_INFO)

        if wrapper.my_param_is_null(phe.SLIDE_DATE):
            wrapper.set_my_param(phe.SLIDE_DATE, phe.NONE_SPECIFIED)

        if wrapper.my_param_is_null(phe.EXTERNAL_IDS):
            wrapper.set_my_param(phe.EXTERNAL_IDS, "external_ids")
        if wrapper.my_param_is_null(phe.TUMOUR_TISSUE):
            tissue = tools.get_tissue_from_sample_id(wrapper.get_my_string(phe.TUMOUR_ID))
            wrapper.set_my_param(phe.TUMOUR_TISSUE, tissue)
        if wrapper.my_param_is_null(phe.NORMAL_TISSUE):
            tissue = tools.get_tissue_from_sample_id(wrapper.get_my_string(phe.NORMAL_ID))
            wrapper.set_my_param(phe.NORMAL_TISSUE, tissue)

        if wrapper.my_param_is_null(phe.TUMOUR_COVERAGE_PATH):
          if self.workspace.has_file(core_constants.DEFAULT_PATH_INFO):
              path_info = self.workspace.read_json(core_constants.DEFAULT_PATH_INFO)
              workflow_paths = path_info.get('coveragePaths')
              workflow_path = workflow_paths['tumour']
              if workflow_path == None:
                  msg = 'Cannot find {0}'.format(phe.TUMOUR_COVERAGE_PATH)
                  self.logger.error(msg)
                  raise RuntimeError(msg)
              wrapper.set_my_param(phe.TUMOUR_COVERAGE_PATH, workflow_path)
        if wrapper.my_param_is_null(phe.NORMAL_COVERAGE_PATH):
          if self.workspace.has_file(core_constants.DEFAULT_PATH_INFO):
              path_info = self.workspace.read_json(core_constants.DEFAULT_PATH_INFO)
              workflow_paths = path_info.get('coveragePaths')
              workflow_path = workflow_paths['normal']
              if workflow_path == None:
                  msg = 'Cannot find {0}'.format(phe.NORMAL_COVERAGE_PATH)
                  self.logger.error(msg)
                  raise RuntimeError(msg)
              wrapper.set_my_param(phe.NORMAL_COVERAGE_PATH, workflow_path)

        return wrapper.get_config()

    def extract(self, config):
        wrapper = self.get_config_wrapper(config)
        attributes = wrapper.get_my_attributes()
        self.check_attributes_known(attributes)
        data = self.get_starting_plugin_data(wrapper, self.PLUGIN_VERSION)
        results_keys = [
            phe.TUMOUR_ID ,
            phe.EXTERNAL_IDS,
            phe.TUMOUR_TISSUE ,
            phe.NORMAL_ID ,
            phe.NORMAL_TISSUE ,
            phe.TMB,
            phe.SNV_LOAD ,
            phe.INDEL_LOAD ,
            phe.SV_LOAD ,
            phe.MOFFIT ,
            phe.HRD ,
            phe.MENGHI ,
            phe.MMRD 
        ]
        data[core_constants.RESULTS] = {k: wrapper.get_my_string(k) for k in results_keys}
        data[core_constants.RESULTS][phe.CELLULOID_PLOT] = tools.convert_plot(self, wrapper.get_my_string(phe.CELLULOID_PLOT), phe.CELLULOID_PLOT)
        data[core_constants.RESULTS][phe.TUMOUR_COVERAGE] = tools.parse_coverage(self, wrapper.get_my_string(phe.TUMOUR_COVERAGE_PATH))
        data[core_constants.RESULTS][phe.NORMAL_COVERAGE] = tools.parse_coverage(self, wrapper.get_my_string(phe.NORMAL_COVERAGE_PATH))
        data[core_constants.RESULTS][phe.INDEL_BIN_PLOT] = tools.convert_svg_plot(self, wrapper.get_my_string(phe.INDEL_BIN_PLOT), phe.INDEL_BIN_PLOT)
        data[core_constants.RESULTS][phe.SV_BIN_PLOT] = tools.convert_svg_plot(self, wrapper.get_my_string(phe.SV_BIN_PLOT), phe.SV_BIN_PLOT)
        data[core_constants.RESULTS][phe.COSMIC_STACK_PLOT] = tools.convert_svg_plot(self, wrapper.get_my_string(phe.COSMIC_STACK_PLOT), phe.COSMIC_STACK_PLOT)
        data[core_constants.RESULTS][phe.SNV_CONTEXT_PLOT] = tools.convert_svg_plot(self, wrapper.get_my_string(phe.SNV_CONTEXT_PLOT), phe.SNV_CONTEXT_PLOT)
        data[core_constants.RESULTS][phe.SV_STACK_PLOT] = tools.convert_svg_plot(self, wrapper.get_my_string(phe.SV_STACK_PLOT), phe.SV_STACK_PLOT)
        data[core_constants.RESULTS][phe.TMB_STACK_PLOT] = tools.convert_svg_plot(self, wrapper.get_my_string(phe.TMB_STACK_PLOT), phe.TMB_STACK_PLOT)


        if wrapper.get_my_string(phe.SLIDE_DATE) == phe.NONE_SPECIFIED:
            data[core_constants.RESULTS][phe.SLIDE_DATE] = strftime("%a %b %d %H:%M:%S %Y")
        else:
            data[core_constants.RESULTS][phe.SLIDE_DATE] = wrapper.get_my_string(phe.SLIDE_DATE)
        ploidy = tools.parse_celluloid_params(self, wrapper.get_my_string(phe.PARAM_PATH), phe.PLOIDY)
        cellularity = tools.parse_celluloid_params(self, wrapper.get_my_string(phe.PARAM_PATH), phe.CELLULARITY)
        data[core_constants.RESULTS][phe.PLOIDY] = ploidy
        data[core_constants.RESULTS][phe.CELLULARITY] = cellularity
        return(data)

    def render(self, data):
        renderer = mako_renderer(self.get_module_dir())
        return renderer.render_name(self.TEMPLATE_NAME, data)

    def specify_params(self):
        discovered = [
            phe.PARAM_PATH,
            phe.COSMIC_SIGNNLS_PATH,
            phe.TDP_PATH,
            phe.SUMMARY_FILE_PATH,
            phe.TUMOUR_COVERAGE_PATH,
            phe.NORMAL_COVERAGE_PATH,
            phe.CELLULOID_PLOT,
            phe.SLIDE_DATE ,
            phe.TUMOUR_ID ,
            phe.EXTERNAL_IDS,
            phe.TUMOUR_TISSUE ,
            phe.NORMAL_ID ,
            phe.NORMAL_TISSUE ,
            phe.TMB,
            phe.SNV_LOAD ,
            phe.INDEL_LOAD ,
            phe.SV_LOAD ,
            phe.MOFFIT ,
            phe.HRD ,
            phe.MENGHI ,
            phe.MMRD ,
            phe.INDEL_BIN_PLOT,
            phe.SV_BIN_PLOT,
            phe.SNV_CONTEXT_PLOT,
            phe.COSMIC_STACK_PLOT,
            phe.SV_STACK_PLOT,
            phe.TMB_STACK_PLOT
        ]
        for key in discovered:
            self.add_ini_discovered(key)
        self.set_ini_default(core_constants.ATTRIBUTES, 'slide')
        self.set_priority_defaults(self.PRIORITY)
