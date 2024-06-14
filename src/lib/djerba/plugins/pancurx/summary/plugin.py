"""Plugin to generate the "Patient & Physician Info" report section"""

import logging
from time import strftime
import re
import json
import os
import csv
import subprocess
import requests
import djerba.core.constants as core_constants
from djerba.plugins.base import plugin_base
from djerba.util.render_mako import mako_renderer
import djerba.plugins.pancurx.constants as phe
import djerba.plugins.pancurx.tools as tools
from djerba.util.image_to_base64 import converter

class main(plugin_base):

    PRIORITY = 100
    PLUGIN_VERSION = '1.0.0'
    TEMPLATE_NAME = 'summary_template.html'
    
    def configure(self, config):
        config = self.apply_defaults(config)
        wrapper = self.get_config_wrapper(config)

        wrapper = tools.fill_file_if_null(self, wrapper, phe.TUMOUR_ID, phe.TUMOUR_ID, core_constants.DEFAULT_SAMPLE_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.DONOR_ID, phe.DONOR_ID, core_constants.DEFAULT_SAMPLE_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.NORMAL_ID, phe.NORMAL_ID, core_constants.DEFAULT_SAMPLE_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, 'template_type', 'template_type', core_constants.DEFAULT_SAMPLE_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.EXTERNAL_IDS, phe.EXTERNAL_IDS, core_constants.DEFAULT_SAMPLE_INFO)

        if wrapper.my_param_is_null(phe.REPORT_DATE):
            wrapper.set_my_param(phe.REPORT_DATE, phe.NONE_SPECIFIED)

        if wrapper.my_param_is_null(phe.SAMPLE_TYPE):
            wrapper.set_my_param(phe.SAMPLE_TYPE, phe.NONE_SPECIFIED)
        if wrapper.my_param_is_null(phe.LOCATION_SUBTYPE):
            wrapper.set_my_param(phe.LOCATION_SUBTYPE, phe.NONE_SPECIFIED)
        
        if wrapper.my_param_is_null(phe.TUMOUR_TISSUE):
            tissue = tools.get_tissue_from_sample_id(wrapper.get_my_string(phe.TUMOUR_ID))
            wrapper.set_my_param(phe.TUMOUR_TISSUE, tissue)
        if wrapper.my_param_is_null(phe.NORMAL_TISSUE):
            tissue = tools.get_tissue_from_sample_id(wrapper.get_my_string(phe.NORMAL_ID))
            wrapper.set_my_param(phe.NORMAL_TISSUE, tissue)

        wrapper = tools.fill_file_if_null(self, wrapper, phe.PARAM_PATH, phe.PARAM_PATH, core_constants.DEFAULT_PATH_INFO)
        
        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'tumour', phe.TUMOUR_COVERAGE_PATH, core_constants.DEFAULT_PATH_INFO, 'coveragePaths')
        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'normal', phe.NORMAL_COVERAGE_PATH, core_constants.DEFAULT_PATH_INFO, 'coveragePaths')

        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'drivers', phe.ONCOSLIDE_DRIVER_PLOT, core_constants.DEFAULT_PATH_INFO, 'svg_plots')
        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'snv.bar_count', phe.ONCOSLIDE_SNV_PLOT, core_constants.DEFAULT_PATH_INFO, 'svg_plots')
        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'indel.bar_count', phe.ONCOSLIDE_INDEL_PLOT, core_constants.DEFAULT_PATH_INFO, 'svg_plots')
        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'sv.bar_count', phe.ONCOSLIDE_SV_PLOT, core_constants.DEFAULT_PATH_INFO, 'svg_plots')
        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'oncoplot.SVs', phe.ONCOSLIDE_CNV_PLOT, core_constants.DEFAULT_PATH_INFO, 'svg_plots')

        return wrapper.get_config()

    def extract(self, config):
        wrapper = self.get_config_wrapper(config)
        attributes = wrapper.get_my_attributes()
        self.check_attributes_known(attributes)
        data = self.get_starting_plugin_data(wrapper, self.PLUGIN_VERSION)
        results_keys = [
            phe.REPORT_VERSION,
            phe.DONOR_ID,
            phe.TUMOUR_ID,
            phe.NORMAL_ID,
            phe.EXTERNAL_IDS,
            phe.TUMOUR_TISSUE,
            phe.NORMAL_TISSUE,
            phe.SAMPLE_TYPE,
            phe.LOCATION_SUBTYPE,
        ]
        data[core_constants.RESULTS] = {k: wrapper.get_my_string(k) for k in results_keys}
        if wrapper.get_my_string(phe.REPORT_DATE) == phe.NONE_SPECIFIED:
            data[core_constants.RESULTS][phe.REPORT_DATE] = strftime("%a %b %d %H:%M:%S %Y")
        else:
            data[core_constants.RESULTS][phe.REPORT_DATE] = wrapper.get_my_string(phe.REPORT_DATE)
        cellularity = tools.parse_celluloid_params(self, wrapper.get_my_string(phe.PARAM_PATH), phe.CELLULARITY)
        data[core_constants.RESULTS][phe.CELLULARITY] = cellularity
        data[core_constants.RESULTS][phe.TUMOUR_COVERAGE] = tools.parse_coverage(self, wrapper.get_my_string(phe.TUMOUR_COVERAGE_PATH))
        data[core_constants.RESULTS][phe.NORMAL_COVERAGE] = tools.parse_coverage(self, wrapper.get_my_string(phe.NORMAL_COVERAGE_PATH))

        data[core_constants.RESULTS][phe.ONCOSLIDE_DRIVER_PLOT] = tools.convert_svg_plot(self, wrapper.get_my_string(phe.ONCOSLIDE_DRIVER_PLOT),   phe.ONCOSLIDE_DRIVER_PLOT)
        data[core_constants.RESULTS][phe.ONCOSLIDE_SNV_PLOT] = tools.convert_svg_plot(self, wrapper.get_my_string(phe.ONCOSLIDE_SNV_PLOT),   phe.ONCOSLIDE_SNV_PLOT)
        data[core_constants.RESULTS][phe.ONCOSLIDE_INDEL_PLOT] = tools.convert_svg_plot(self, wrapper.get_my_string(phe.ONCOSLIDE_INDEL_PLOT),   phe.ONCOSLIDE_INDEL_PLOT)
        data[core_constants.RESULTS][phe.ONCOSLIDE_SV_PLOT] = tools.convert_svg_plot(self, wrapper.get_my_string(phe.ONCOSLIDE_SV_PLOT),   phe.ONCOSLIDE_SV_PLOT)
        data[core_constants.RESULTS][phe.ONCOSLIDE_CNV_PLOT] = tools.convert_svg_plot(self, wrapper.get_my_string(phe.ONCOSLIDE_CNV_PLOT),   phe.ONCOSLIDE_CNV_PLOT)
        data[core_constants.RESULTS]['template_type'] = '_'.join((wrapper.get_my_string('template_type'), self.TEMPLATE_NAME))

        return(data)

    def render(self, data):
        renderer = mako_renderer(self.get_module_dir())
        template_name = data[core_constants.RESULTS]['template_type'] 
        return renderer.render_name(template_name, data)

    def specify_params(self):
        discovered = [
            phe.TUMOUR_ID,
            phe.DONOR_ID,
            phe.NORMAL_ID,

            phe.ONCOSLIDE_DRIVER_PLOT,
            phe.ONCOSLIDE_SNV_PLOT,
            phe.ONCOSLIDE_INDEL_PLOT,
            phe.ONCOSLIDE_SV_PLOT,
            phe.ONCOSLIDE_CNV_PLOT,
            phe.TUMOUR_COVERAGE_PATH,
            phe.NORMAL_COVERAGE_PATH,
            phe.PARAM_PATH,
            phe.REPORT_DATE,
            phe.EXTERNAL_IDS,
            phe.TUMOUR_TISSUE,
            phe.NORMAL_TISSUE,
            'template_type',
            phe.SAMPLE_TYPE,
            phe.LOCATION_SUBTYPE 

            
        ]
        for key in discovered:
            self.add_ini_discovered(key)
        self.set_ini_default(core_constants.ATTRIBUTES, 'research')
        self.set_ini_default(phe.REPORT_VERSION, phe.CURRENT_REPORT_VERSION)
        self.set_priority_defaults(self.PRIORITY)
        self.set_ini_default('render_priority', 30)
