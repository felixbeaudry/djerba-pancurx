"""Plugin to generate the "Gremline" report section"""

import logging
from time import strptime, strftime
import csv
import os
import re
import djerba.core.constants as core_constants
from djerba.plugins.base import plugin_base
from djerba.util.render_mako import mako_renderer
from djerba.util.image_to_base64 import converter
import djerba.plugins.pancurx.pdo.pdo_constants as phe
import djerba.plugins.pancurx.tools as tools
from djerba.util.environment import directory_finder

class main(plugin_base):

    PRIORITY = 100
    PLUGIN_VERSION = '1.0.0'
    TEMPLATE_NAME = 'pdo_template.html'


    def configure(self, config):
        config = self.apply_defaults(config)
        wrapper = self.get_config_wrapper(config)

        if wrapper.my_param_is_null(phe.REPORT_DATE):
            wrapper.set_my_param(phe.REPORT_DATE, phe.NONE_SPECIFIED)

        if wrapper.my_param_is_null(phe.PDO_TYPE):
            pdo_type = check_pdo_name_for_met(wrapper.get_my_string(phe.PDO_ID))
            wrapper.set_my_param(phe.PDO_TYPE, pdo_type)

        #psp_id = wrapper.get_my_string(phe.PSP_ID)
        pdo_id = wrapper.get_my_string(phe.PDO_ID)

        if wrapper.my_param_is_null(phe.VIABILITY_SCORE_PLOT):
            plot_file_name = '.'.join(( pdo_id,'qc.svg'))
            plot_file_path = os.path.join(self.workspace.print_location(),  plot_file_name)
            wrapper.set_my_param(phe.VIABILITY_SCORE_PLOT, plot_file_path)

        if wrapper.my_param_is_null(phe.SENSITIVITIES_TABLE):
            file_name = '.'.join(( pdo_id,'sensivities.txt'))
            file_path = os.path.join(self.workspace.print_location(),  file_name)
            wrapper.set_my_param(phe.SENSITIVITIES_TABLE, file_path)

        if wrapper.my_param_is_null(phe.VIABILITY_SCORE):
            file_name = '.'.join(( pdo_id,'zprime.txt'))
            file_path = os.path.join(self.workspace.print_location(),  file_name)
            viability_score = tools.parse_zprime(self, file_path)
            wrapper.set_my_param(phe.VIABILITY_SCORE, viability_score)

        if wrapper.my_param_is_null(phe.CURVE_DIR):
            plot_dir_path = os.path.join(self.workspace.print_location(),  'curves')
            wrapper.set_my_param(phe.CURVE_DIR, plot_dir_path)

        return wrapper.get_config()

    def extract(self, config):
        wrapper = self.get_config_wrapper(config)
        attributes = wrapper.get_my_attributes()
        self.check_attributes_known(attributes)
        data = self.get_starting_plugin_data(wrapper, self.PLUGIN_VERSION)
        
        results_keys = [
            phe.PDO_ID,
            phe.PSP_ID,
            phe.PANX_ID,
            phe.ADOPT_ID,
            phe.PDO_TYPE
            ]
        
        data[core_constants.RESULTS] = {k: wrapper.get_my_string(k) for k in results_keys}        
        if wrapper.get_my_string(phe.REPORT_DATE) == phe.NONE_SPECIFIED:
            data[core_constants.RESULTS][phe.REPORT_DATE] = strftime("%a %b %d %H:%M:%S %Y")
        else:
            data[core_constants.RESULTS][phe.REPORT_DATE] = wrapper.get_my_string(phe.REPORT_DATE)
        
        data[core_constants.RESULTS][phe.VIABILITY_SCORE] = float(wrapper.get_my_string(phe.VIABILITY_SCORE)) 
        data[core_constants.RESULTS][phe.VIABILITY_SCORE_PLOT] = tools.convert_svg_plot(self, wrapper.get_my_string(phe.VIABILITY_SCORE_PLOT),   phe.VIABILITY_SCORE_PLOT)
        data[core_constants.RESULTS][phe.SENSITIVITIES_TABLE] = tools.parse_sensitivities(self, wrapper.get_my_string(phe.SENSITIVITIES_TABLE))
                 

        data[core_constants.RESULTS]['template_type'] = self.TEMPLATE_NAME
        
        drug_curve_plots = []

        # List all files in the directory and pass them to the function
        drug_count = 0
        #plot_files = os.listdir(wrapper.get_my_string(phe.CURVE_DIR))
        plot_files = [f for f in os.listdir(wrapper.get_my_string(phe.CURVE_DIR)) if os.path.isfile(os.path.join(wrapper.get_my_string(phe.CURVE_DIR), f))]
        plot_files.sort()
        for filename in plot_files:
            file_path = os.path.join(wrapper.get_my_string(phe.CURVE_DIR), filename)
            if os.path.isfile(file_path):
                plot_count = '_'.join(('drug_plot_', str(drug_count)))
                result = tools.convert_svg_plot(self, file_path, plot_count)
                drug_curve_plots.append(result)
                drug_count = drug_count + 1

        data[core_constants.RESULTS]['drug_curve_plots'] = drug_curve_plots

        return data

    def render(self, data):
        renderer = mako_renderer(self.get_module_dir())
        template_name = data[core_constants.RESULTS]['template_type'] 
        return renderer.render_name(template_name, data)

    def specify_params(self):
        discovered = [
            phe.PDO_ID,
            phe.PSP_ID,
            phe.PANX_ID,
            phe.ADOPT_ID,
            phe.PDO_TYPE,
            phe.REPORT_DATE,
            phe.VIABILITY_SCORE,
            phe.VIABILITY_SCORE_PLOT,
            phe.SENSITIVITIES_TABLE,
            phe.CURVE_DIR
        ]
        for key in discovered:
            self.add_ini_discovered(key)
        self.set_ini_default(core_constants.ATTRIBUTES, 'research')
        self.set_priority_defaults(self.PRIORITY)

def check_pdo_name_for_met(input_string):
    substrings = input_string.split('_')
    if len(substrings) >= 3 and substrings[2] == 'M':
        return 'Metastasis'
    elif len(substrings) >= 3 and substrings[2] == 'P':
        return 'Primary'
    return None

