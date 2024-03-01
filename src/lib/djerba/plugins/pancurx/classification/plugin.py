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
import djerba.plugins.pancurx.classification.constants as phe

class main(plugin_base):

    PRIORITY = 200
    PLUGIN_VERSION = '1.0.0'
    TEMPLATE_NAME = 'classification_template.html'
    PARAM_PATH = "param_path"
    SEX_PATH = "sex_path"
    TDP_PATH = "tdp_path"
    COSMIC_SIGNNLS_PATH = "cosmic_signnls_path"
    SUMMARY_FILE_PATH = "summary_file_path"


    def configure(self, config):
        config = self.apply_defaults(config)
        wrapper = self.get_config_wrapper(config)
        wrapper = self.fill_file_if_null(wrapper, self.PARAM_PATH, self.PARAM_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = self.fill_file_if_null(wrapper, self.SEX_PATH, self.SEX_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = self.fill_file_if_null(wrapper, self.TDP_PATH, self.TDP_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = self.fill_file_if_null(wrapper, self.COSMIC_SIGNNLS_PATH, self.COSMIC_SIGNNLS_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = self.fill_file_if_null(wrapper, self.SUMMARY_FILE_PATH, self.SUMMARY_FILE_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = self.fill_file_if_null(wrapper, phe.DSBR_SCORE_BAR, phe.DSBR_SCORE_BAR, core_constants.DEFAULT_PATH_INFO)
        wrapper = self.fill_file_if_null(wrapper, phe.MMR_SCORE_BAR, phe.MMR_SCORE_BAR, core_constants.DEFAULT_PATH_INFO)
        return wrapper.get_config()

    def extract(self, config):
        wrapper = self.get_config_wrapper(config)
        attributes = wrapper.get_my_attributes()
        self.check_attributes_known(attributes)
        data = self.get_starting_plugin_data(wrapper, self.PLUGIN_VERSION)
        results_keys = [
            phe.WADDELL_CLASS, #TODO: from SV counts and positions
            phe.ALEXANDROV_CLASS,
            phe.COLLISSON_CLASS,
            phe.MOFFITT_CLASS,
            phe.HLA_TYPES
        ]
        data[core_constants.RESULTS] = {k: wrapper.get_my_string(k) for k in results_keys}
        data[core_constants.RESULTS][phe.INFERRED_SEX] = self.parse_sex(wrapper.get_my_string(self.SEX_PATH))
        ploidy = self.parse_celluloid_params(wrapper.get_my_string(self.PARAM_PATH), phe.PLOIDY)
        data[core_constants.RESULTS][phe.PLOIDY] = ploidy
        data[core_constants.RESULTS]['cosmic_signatures'] = self.parse_cosmic_signatures(wrapper.get_my_string(self.COSMIC_SIGNNLS_PATH))
        data[core_constants.RESULTS]['tandem_duplicator_phenotype_score'] = self.parse_TDP(wrapper.get_my_string(self.TDP_PATH))
        summary_results = self.parse_summary_file(wrapper.get_my_string(self.SUMMARY_FILE_PATH))
        data[core_constants.RESULTS][phe.DSBR_RESULTS] = self.parse_multifactor_marker(summary_results, phe.DSBR_HTML_HEADERS, phe.DSBR_DEFAULT_HALLMARK_CUTOFFS)
        data[core_constants.RESULTS][phe.MMR_RESULTS] = self.parse_multifactor_marker(summary_results, phe.MMR_HTML_HEADERS, phe.MMR_DEFAULT_HALLMARK_CUTOFFS)
        data[core_constants.RESULTS][phe.DSBR_SCORE_BAR] = self.convert_plot(wrapper.get_my_string(phe.DSBR_SCORE_BAR), phe.DSBR_SCORE_BAR)
        data[core_constants.RESULTS][phe.MMR_SCORE_BAR] = self.convert_plot(wrapper.get_my_string(phe.MMR_SCORE_BAR), phe.MMR_SCORE_BAR)
        return data

    def render(self, data):
        renderer = mako_renderer(self.get_module_dir())
        return renderer.render_name(self.TEMPLATE_NAME, data)

    def specify_params(self):
        discovered = [
            self.PARAM_PATH,
            self.SEX_PATH,
            self.COSMIC_SIGNNLS_PATH,
            self.TDP_PATH,
            phe.ALEXANDROV_CLASS,
            phe.COLLISSON_CLASS,
            phe.WADDELL_CLASS,
            phe.MOFFITT_CLASS,
            phe.HLA_TYPES,
            self.SUMMARY_FILE_PATH,
            phe.MMR_SCORE_BAR,
            phe.DSBR_SCORE_BAR
        ]
        for key in discovered:
            self.add_ini_discovered(key)
        self.set_ini_default(core_constants.ATTRIBUTES, 'clinical')
        self.set_priority_defaults(self.PRIORITY)

    def convert_plot(self, plot_path, plot_name):
        """Read VAF plot from file and return as a base64 string"""
        image_converter = converter(self.log_level, self.log_path)
        converted_plot = image_converter.convert_png(plot_path, plot_name)
        return converted_plot

    def parse_sex(self, sex_file_path):
        sex_file_path = self.check_path_exists(sex_file_path)
        with open(sex_file_path) as sex_file:
            for row in sex_file:
                short_sex = row.strip()
                if short_sex == 'M':
                    inferredSex = "Male"
                elif short_sex == 'F':
                    inferredSex = "Female"
                else:
                    msg = "Cannot infer sex {0} from file: {1}".format(short_sex, sex_file_path)
                    self.logger.error(msg)
                    raise RuntimeError(msg)
        return(inferredSex)
    
    def parse_celluloid_params(self, celluloid_params_file_path, ploidy_or_cellularity):
        celluloid_params_file_path = self.check_path_exists(celluloid_params_file_path)
        with open(celluloid_params_file_path, 'r') as celluloid_params_file:
            for row in csv.DictReader(celluloid_params_file, delimiter=" "):
                cellularity = "{:.1f}".format(float(row["T1"])*100)
                ploidy = "{:.3f}".format(float(row["Ploidy"]))
        if ploidy_or_cellularity == "ploidy":
            ploidy_string = ''
            if float(ploidy) > phe.DIPLOID_CUTOFF:
                ploidy_string = "polyploid ({0})".format(ploidy)
            elif float(ploidy) <= phe.DIPLOID_CUTOFF:
                ploidy_string = "diploid ({0})".format(ploidy)
            else:
                msg = "Ploidy {0} not a number".format(ploidy)
                self.logger.error(msg)
                raise RuntimeError(msg)
            return(ploidy_string)
        else:
            return(cellularity)
            
    def parse_cosmic_signatures(self, cosmic_signatures_file_path):
        signature_results = {}
        row = {}
        cosmic_signatures_file_path = self.check_path_exists(cosmic_signatures_file_path)
        with open(cosmic_signatures_file_path, 'r') as cosmic_signatures_file:
            line = cosmic_signatures_file.readline()
            header = line.split()
            line = cosmic_signatures_file.readline()
            for header_position, signature_value in enumerate(line.split()):
                row[header[header_position]] = signature_value
        for signature_name in phe.COSMIC_SIGNATURE_SET:
            if signature_name in row:
                signature_results[phe.COSMIC_SIGNATURE_SET[signature_name]] = row[signature_name]
        return(signature_results)
	
    def parse_TDP(self, TDP_file_path):
        row = {}
        TDP_file_path = self.check_path_exists(TDP_file_path)
        with open(TDP_file_path, 'r') as TDP_file:
            line = TDP_file.readline().strip()
            header = line.split(',')
            line = TDP_file.readline().strip()
            row = dict(zip(header, line.split(',')))
        return(row["score"])

    def fill_file_if_null(self, wrapper, workflow_name, ini_param, path_info):
      if wrapper.my_param_is_null(ini_param):
          if self.workspace.has_file(path_info):
              path_info = self.workspace.read_json(path_info)
              workflow_path = path_info.get(workflow_name)
              if workflow_path == None:
                  msg = 'Cannot find {0}'.format(ini_param)
                  self.logger.error(msg)
                  raise RuntimeError(msg)
              wrapper.set_my_param(ini_param, workflow_path)
      return(wrapper)
    
    def check_path_exists(self, path):
        if not os.path.exists(path):
            msg = "Cannot find file: {0}".format(path)
            self.logger.error(msg)
            raise RuntimeError(msg)
        else:
            return(path)

    def parse_multifactor_marker(self, summary_results, html_headers, marker_cutoffs):
        result = []
        for this_key in summary_results:
            this_reporting_name = html_headers.get(this_key)
            this_value = summary_results.get(this_key)
            if this_key in marker_cutoffs:
                this_value = float(this_value)
                this_cutoff = marker_cutoffs.get(this_key)
                if this_value > float(this_cutoff):
                    above_cutoff = True
                else:
                    above_cutoff = False
                this_value = round(this_value, 2)
                this_reporting_name = " ".join((this_reporting_name, str(this_cutoff), ":"))
                result_tmp = {
                    phe.REPORTING_NAME: this_reporting_name,
                    phe.VALUE : this_value,
                    phe.ABOVE_CUTOFF: above_cutoff
                }
                result.append(result_tmp)
            elif this_key in html_headers:
                if this_value == "":
                    above_cutoff = False
                else:
                    above_cutoff = True
                    this_value = re.sub(r'\|', ', ', this_value)
                result_tmp = {
                    phe.REPORTING_NAME: this_reporting_name,
                    phe.VALUE : this_value,
                    phe.ABOVE_CUTOFF: above_cutoff
                }
                result.append(result_tmp)
        return(result)

    def parse_summary_file(self, summary_file_path):
        row = {}
        summary_file_path = self.check_path_exists(summary_file_path)
        with open(summary_file_path, 'r') as summary_file:
            line = summary_file.readline().strip()
            header = line.split(',')
            line = summary_file.readline().strip()
            row = dict(zip(header, line.split(',')))
        return(row)