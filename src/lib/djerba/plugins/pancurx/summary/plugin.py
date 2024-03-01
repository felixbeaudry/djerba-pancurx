"""Plugin to generate the "Patient & Physician Info" report section"""

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
import djerba.plugins.pancurx.summary.constants as phe
from djerba.util.image_to_base64 import converter

class main(plugin_base):

    PRIORITY = 100
    PLUGIN_VERSION = '1.0.0'
    TEMPLATE_NAME = 'summary_template.html'
    PARAM_PATH = "param_path"
    TUMOUR_COVERAGE_PATH = "tumour_coverage_path"
    NORMAL_COVERAGE_PATH = "normal_coverage_path"

    def configure(self, config):
        config = self.apply_defaults(config)
        wrapper = self.get_config_wrapper(config)
        if wrapper.my_param_is_null(phe.REPORT_DATE):
            wrapper.set_my_param(phe.REPORT_DATE, phe.NONE_SPECIFIED)
        if wrapper.my_param_is_null(phe.EXTERNAL_IDS):
            external_ids = self.parse_LIMS(wrapper.get_my_string(phe.DONOR_ID), wrapper.get_my_string(phe.TUMOUR_ID))
            wrapper.set_my_param(phe.EXTERNAL_IDS, external_ids)
        if wrapper.my_param_is_null(phe.TUMOUR_TISSUE):
            tissue = self.get_tissue_from_sample_id(wrapper.get_my_string(phe.TUMOUR_ID))
            wrapper.set_my_param(phe.TUMOUR_TISSUE, tissue)
        if wrapper.my_param_is_null(phe.NORMAL_TISSUE):
            tissue = self.get_tissue_from_sample_id(wrapper.get_my_string(phe.NORMAL_ID))
            wrapper.set_my_param(phe.NORMAL_TISSUE, tissue)
        wrapper = self.fill_file_if_null(wrapper, phe.ONCOSLIDE_PLOT, phe.ONCOSLIDE_PLOT, core_constants.DEFAULT_PATH_INFO)
        wrapper = self.fill_file_if_null(wrapper, self.PARAM_PATH, self.PARAM_PATH, core_constants.DEFAULT_PATH_INFO)

        if wrapper.my_param_is_null(self.TUMOUR_COVERAGE_PATH):
          if self.workspace.has_file(core_constants.DEFAULT_PATH_INFO):
              path_info = self.workspace.read_json(core_constants.DEFAULT_PATH_INFO)
              workflow_paths = path_info.get('coveragePaths')
              workflow_path = workflow_paths['tumour']
              if workflow_path == None:
                  msg = 'Cannot find {0}'.format(self.TUMOUR_COVERAGE_PATH)
                  self.logger.error(msg)
                  raise RuntimeError(msg)
              wrapper.set_my_param(self.TUMOUR_COVERAGE_PATH, workflow_path)
        if wrapper.my_param_is_null(self.NORMAL_COVERAGE_PATH):
          if self.workspace.has_file(core_constants.DEFAULT_PATH_INFO):
              path_info = self.workspace.read_json(core_constants.DEFAULT_PATH_INFO)
              workflow_paths = path_info.get('coveragePaths')
              workflow_path = workflow_paths['normal']
              if workflow_path == None:
                  msg = 'Cannot find {0}'.format(self.NORMAL_COVERAGE_PATH)
                  self.logger.error(msg)
                  raise RuntimeError(msg)
              wrapper.set_my_param(self.NORMAL_COVERAGE_PATH, workflow_path)

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
        ]
        data[core_constants.RESULTS] = {k: wrapper.get_my_string(k) for k in results_keys}
        if wrapper.get_my_string(phe.REPORT_DATE) == phe.NONE_SPECIFIED:
            data[core_constants.RESULTS][phe.REPORT_DATE] = strftime("%a %b %d %H:%M:%S %Y")
        else:
            data[core_constants.RESULTS][phe.REPORT_DATE] = wrapper.get_my_string(phe.REPORT_DATE)
        cellularity = self.parse_celluloid_params(wrapper.get_my_string(self.PARAM_PATH), phe.CELLULARITY)
        data[core_constants.RESULTS][phe.CELLULARITY] = cellularity
        data[core_constants.RESULTS][phe.TUMOUR_COVERAGE] = self.parse_coverage(wrapper.get_my_string(self.TUMOUR_COVERAGE_PATH))
        data[core_constants.RESULTS][phe.NORMAL_COVERAGE] = self.parse_coverage(wrapper.get_my_string(self.NORMAL_COVERAGE_PATH))
        data[core_constants.RESULTS][phe.ONCOSLIDE_PLOT] = self.convert_plot(wrapper.get_my_string(phe.ONCOSLIDE_PLOT), 
                                                                                   phe.ONCOSLIDE_PLOT)
        return(data)

    def render(self, data):
        renderer = mako_renderer(self.get_module_dir())
        return renderer.render_name(self.TEMPLATE_NAME, data)

    def specify_params(self):
        discovered = [
            self.TUMOUR_COVERAGE_PATH,
            self.NORMAL_COVERAGE_PATH,
            self.PARAM_PATH,
            phe.TUMOUR_ID,
            phe.DONOR_ID,
            phe.NORMAL_ID,
            phe.REPORT_DATE,
            phe.EXTERNAL_IDS,
            phe.TUMOUR_TISSUE,
            phe.NORMAL_TISSUE,
            phe.ONCOSLIDE_PLOT
        ]
        for key in discovered:
            self.add_ini_discovered(key)
        self.set_ini_default(core_constants.ATTRIBUTES, 'clinical')
        self.set_ini_default(phe.REPORT_VERSION, phe.CURRENT_REPORT_VERSION)
        self.set_priority_defaults(self.PRIORITY)
        self.set_ini_default('render_priority', 30)

    def add_underscore_to_donor(self, donor):
        if "EPPIC" in donor:
            donor = donor.replace("EPPIC", "EPPIC_")
        else:
            donor = re.sub(r'^([A-Z0-9]+)(....)', r'\1_\2', donor)
        return(donor)

    def find_external_id_in_json_dict(self, lims_dict, this_donor):
        external_id = "NA"
        for donor in range(len(lims_dict)):
            if lims_dict[donor]['name'] == this_donor:
                for attribute in range(len(lims_dict[donor]['attributes'])):
                    if lims_dict[donor]['attributes'][attribute]['name'] == "External Name":
                        external_id = lims_dict[donor]['attributes'][attribute]['value']
                        external_id = external_id.replace(',', ' ')
        return(external_id)

    def parse_LIMS_dump(self, donor, lims_dump_path):
        with open(lims_dump_path, 'r') as file:
            lims_file = ''.join(file.readlines())
        lims_dict = json.loads(lims_file)
        external_id = self.find_external_id_in_json_dict(lims_dict, donor)
        return(external_id)

    def parse_LIMS_URL(self, donor, lims_url, project="BTC"):
        CURL_COMMAND = "curl -X GET"
        lims_curl_command = " ".join((CURL_COMMAND, lims_url))
        if project in donor:
            lims_curl = "".join((lims_curl_command, '&project=', project))
            json_file = subprocess.getoutput(lims_curl)
        else:
            json_file = subprocess.getoutput(lims_curl_command)
        #print(json_file)
        lims_dict = json.loads(json_file)
        external_id = self.find_external_id_in_json_dict(lims_dict, donor)
        return(external_id)

    def parse_LIMS(self, donor, tumour):
        donor = self.add_underscore_to_donor(donor)
        external_id = "NA"
        external_id = self.parse_LIMS_dump(donor, phe.DEFAULT_LIMS_DUMP)
        if tumour in external_id:
            external_id = external_id[tumour]     
        if external_id == "NA":
            #external_id = self.parse_LIMS_URL(donor, phe.DEFAULT_CORE_LIMS_URL)
            external_id = "Pinery W.I.P."
        return(external_id)

    def get_tissue_from_sample_id(self, sample_name):
        tissue = "NA"
        if re.match(r'^[A-Z0-9]+_...._(..)_', sample_name):
            match_group = re.match(r'^[A-Z0-9]+_...._(..)_', sample_name).group(1)
            tissue = phe.LIMS_TISSUE_CODES.get(match_group, match_group)
        return(tissue)
    
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


    def parse_coverage(self, coverage_file_path):
        coverage_file_path = self.check_path_exists(coverage_file_path)
        with open(coverage_file_path, 'r') as file:
            ## Number of columns is variable among rows
            lines = file.read().splitlines()
            for column_number in range(len(lines)):
                column = lines[column_number]
                if "GENOME_TERRITORY" in column:
                    genome, mean_coverage, sd, median_coverage, *extra = lines[column_number + 1].split('\t')
                    median_coverage = float(median_coverage)
        return(median_coverage)


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
            msg = "Cannot find hugo file: {0}".format(path)
            self.logger.error(msg)
            raise RuntimeError(msg)
        else:
            return(path)

    def convert_plot(self, plot_path, plot_name):
        """Read VAF plot from file and return as a base64 string"""
        image_converter = converter(self.log_level, self.log_path)
        converted_plot = image_converter.convert_png(plot_path, plot_name)
        return converted_plot
