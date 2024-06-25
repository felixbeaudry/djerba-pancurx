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
from djerba.util.image_to_base64 import converter
import djerba.plugins.pancurx.constants as phe
import djerba.plugins.pancurx.tools as tools

class main(plugin_base):

    PRIORITY = 700
    PLUGIN_VERSION = '1.0.0'
    TEMPLATE_NAME = 'all_genes_template.html'


    def configure(self, config):
        config = self.apply_defaults(config)
        wrapper = self.get_config_wrapper(config)

        wrapper = tools.fill_file_if_null(self, wrapper, phe.SAMPLE_VARIANTS_FILE, phe.SAMPLE_VARIANTS_FILE, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.PARAM_PATH, phe.PARAM_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.SEX_PATH, phe.SEX_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, 'template_type', 'template_type', core_constants.DEFAULT_SAMPLE_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.CELLULOID_PLOT, phe.CELLULOID_PLOT, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'indel.vaf', phe.INDEL_VAF_PLOT, core_constants.DEFAULT_PATH_INFO, 'svg_plots')
        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'snv.vaf', phe.SNV_VAF_PLOT, core_constants.DEFAULT_PATH_INFO, 'svg_plots')
        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'indel.stack_count', phe.INDEL_BIN_PLOT, core_constants.DEFAULT_PATH_INFO, 'svg_plots')
        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'sv.bins', phe.SV_BIN_PLOT, core_constants.DEFAULT_PATH_INFO, 'svg_plots')
        wrapper = tools.fill_file_if_null(self, wrapper, phe.WHOLE_GENOME_PLOT, phe.WHOLE_GENOME_PLOT, core_constants.DEFAULT_PATH_INFO)

        
        return wrapper.get_config()

    def extract(self, config):
        wrapper = self.get_config_wrapper(config)
        attributes = wrapper.get_my_attributes()
        self.check_attributes_known(attributes)
        data = self.get_starting_plugin_data(wrapper, self.PLUGIN_VERSION)
        inferred_sex_chromosomes = tools.parse_sex(self, wrapper.get_my_string(phe.SEX_PATH), 'LBR')

        somatic_variants = tools.parse_somatic_variants(self, wrapper.get_my_string(phe.SAMPLE_VARIANTS_FILE))
        data[core_constants.RESULTS]['all_variants'] = tools.get_all_somatic_variants(self, somatic_variants, inferred_sex_chromosomes)
        ploidy = tools.parse_celluloid_params(self, wrapper.get_my_string(phe.PARAM_PATH), "ploidy_numeric")
        data[core_constants.RESULTS][phe.PLOIDY] = ploidy
        data[core_constants.RESULTS][phe.CELLULOID_PLOT] = tools.convert_plot(self, wrapper.get_my_string(phe.CELLULOID_PLOT), phe.CELLULOID_PLOT)
        data[core_constants.RESULTS][phe.INDEL_VAF_PLOT] = tools.convert_svg_plot(self, wrapper.get_my_string(phe.INDEL_VAF_PLOT), phe.INDEL_VAF_PLOT)
        data[core_constants.RESULTS][phe.SNV_VAF_PLOT] = tools.convert_svg_plot(self, wrapper.get_my_string(phe.SNV_VAF_PLOT), phe.SNV_VAF_PLOT)
        data[core_constants.RESULTS][phe.INDEL_BIN_PLOT] = tools.convert_svg_plot(self, wrapper.get_my_string(phe.INDEL_BIN_PLOT), phe.INDEL_BIN_PLOT)
        data[core_constants.RESULTS][phe.SV_BIN_PLOT] = tools.convert_svg_plot(self, wrapper.get_my_string(phe.SV_BIN_PLOT), phe.SV_BIN_PLOT)
        data[core_constants.RESULTS][phe.WHOLE_GENOME_PLOT] = tools.convert_plot(self, wrapper.get_my_string(phe.WHOLE_GENOME_PLOT), phe.WHOLE_GENOME_PLOT)
        data[core_constants.RESULTS]['template_type'] = '_'.join((wrapper.get_my_string('template_type'), self.TEMPLATE_NAME))
        
        ploidy_long = tools.parse_celluloid_params(self, wrapper.get_my_string(phe.PARAM_PATH), phe.PLOIDY)
        data[core_constants.RESULTS]['ploidy_long'] = ploidy_long
        data[core_constants.RESULTS][phe.INFERRED_SEX] = tools.parse_sex(self, wrapper.get_my_string(phe.SEX_PATH), wrapper.get_my_string('template_type'))

        return data

    def render(self, data):
        renderer = mako_renderer(self.get_module_dir())
        template_name = data[core_constants.RESULTS]['template_type'] 
        return renderer.render_name(template_name, data)
    
    def specify_params(self):
        discovered = [
            phe.WHOLE_GENOME_PLOT,
            phe.CELLULOID_PLOT,
            phe.INDEL_VAF_PLOT,
            phe.INDEL_BIN_PLOT,
            phe.SV_BIN_PLOT,
            phe.SNV_VAF_PLOT,
            phe.SEX_PATH,
            phe.PARAM_PATH,
            phe.SAMPLE_VARIANTS_FILE,
            'template_type'

        ]
        for key in discovered:
            self.add_ini_discovered(key)
        self.set_ini_default(core_constants.ATTRIBUTES, 'research')
        self.set_priority_defaults(self.PRIORITY)

