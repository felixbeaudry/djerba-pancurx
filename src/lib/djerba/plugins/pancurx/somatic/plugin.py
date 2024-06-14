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

    PRIORITY = 300
    PLUGIN_VERSION = '1.0.0'
    TEMPLATE_NAME = 'somatic_template.html'


    def configure(self, config):
        config = self.apply_defaults(config)
        wrapper = self.get_config_wrapper(config)

        wrapper = tools.fill_file_if_null(self, wrapper, 'genes_of_interest_file', 'genes_of_interest_file', core_constants.DEFAULT_SAMPLE_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, 'template_type', 'template_type', core_constants.DEFAULT_SAMPLE_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, 'comparison_cohort_file', 'comparison_cohort_file', core_constants.DEFAULT_SAMPLE_INFO)

        wrapper = tools.fill_file_if_null(self, wrapper, phe.SAMPLE_VARIANTS_FILE, phe.SAMPLE_VARIANTS_FILE, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.SUMMARY_FILE_PATH, phe.SUMMARY_FILE_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.PARAM_PATH, phe.PARAM_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.CELLULOID_DIR, phe.CELLULOID_DIR, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.COSMIC_SIGNNLS_PATH, phe.COSMIC_SIGNNLS_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.SEX_PATH, phe.SEX_PATH, core_constants.DEFAULT_PATH_INFO)

        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'SBS_context_bar', phe.SNV_CONTEXT_PLOT, core_constants.DEFAULT_PATH_INFO, 'svg_plots')
        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'indel.histbox_count', phe.HISTBOX_INDEL, core_constants.DEFAULT_PATH_INFO, 'svg_plots')
        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'sv.histbox_count', phe.HISTBOX_SV, core_constants.DEFAULT_PATH_INFO, 'svg_plots')
        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'snv.histbox_count', phe.HISTBOX_SNV, core_constants.DEFAULT_PATH_INFO, 'svg_plots')
        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'oncoplot.SVs', phe.ONCOSLIDE_CNV_PLOT, core_constants.DEFAULT_PATH_INFO, 'svg_plots')
        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'indel.vaf', phe.INDEL_VAF_PLOT, core_constants.DEFAULT_PATH_INFO, 'svg_plots')
        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'snv.vaf', phe.SNV_VAF_PLOT, core_constants.DEFAULT_PATH_INFO, 'svg_plots')
        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'indel.stack_count', phe.INDEL_BIN_PLOT, core_constants.DEFAULT_PATH_INFO, 'svg_plots')
        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'sv.bins', phe.SV_BIN_PLOT, core_constants.DEFAULT_PATH_INFO, 'svg_plots')

        wrapper = tools.fill_file_if_null(self, wrapper, phe.WHOLE_GENOME_PLOT, phe.WHOLE_GENOME_PLOT, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.CELLULOID_PLOT, phe.CELLULOID_PLOT, core_constants.DEFAULT_PATH_INFO)

        return wrapper.get_config()

    def extract(self, config):
        wrapper = self.get_config_wrapper(config)
        attributes = wrapper.get_my_attributes()
        self.check_attributes_known(attributes)
        data = self.get_starting_plugin_data(wrapper, self.PLUGIN_VERSION)
        
        genes_of_interest = tools.get_genes_of_interest(self, wrapper.get_my_string('genes_of_interest_file'))
        somatic_variants = tools.parse_somatic_variants(self, wrapper.get_my_string(phe.SAMPLE_VARIANTS_FILE), genes_of_interest)
        data[core_constants.RESULTS]['reportable_variants'] = tools.get_subset_of_somatic_variants(self, somatic_variants, genes_of_interest)
        summary_results = tools.parse_summary_file(self, wrapper.get_my_string(phe.SUMMARY_FILE_PATH))
        data[core_constants.RESULTS]['loads'] = tools.get_loads_from_summary(summary_results)
        data[core_constants.RESULTS]['sigs'] = tools.parse_cosmic_signatures(self, wrapper.get_my_string(phe.COSMIC_SIGNNLS_PATH))
        data[core_constants.RESULTS][phe.CELLULOID_DIR] =  wrapper.get_my_string(phe.CELLULOID_DIR)
        ploidy = tools.parse_celluloid_params(self, wrapper.get_my_string(phe.PARAM_PATH), "ploidy_numeric")
        data[core_constants.RESULTS][phe.PLOIDY] = ploidy
        ploidy_long = tools.parse_celluloid_params(self, wrapper.get_my_string(phe.PARAM_PATH), phe.PLOIDY)
        data[core_constants.RESULTS]['ploidy_long'] = ploidy_long
        data[core_constants.RESULTS][phe.INFERRED_SEX] = tools.parse_sex(self, wrapper.get_my_string(phe.SEX_PATH), wrapper.get_my_string('template_type'))


        data[core_constants.RESULTS][phe.SNV_CONTEXT_PLOT] = tools.convert_svg_plot(self, wrapper.get_my_string(phe.SNV_CONTEXT_PLOT), phe.SNV_CONTEXT_PLOT)
        data[core_constants.RESULTS][phe.HISTBOX_INDEL] = tools.convert_svg_plot(self, wrapper.get_my_string(phe.HISTBOX_INDEL), phe.HISTBOX_INDEL)
        data[core_constants.RESULTS][phe.HISTBOX_SV] = tools.convert_svg_plot(self, wrapper.get_my_string(phe.HISTBOX_SV), phe.HISTBOX_SV)
        data[core_constants.RESULTS][phe.HISTBOX_SNV] = tools.convert_svg_plot(self, wrapper.get_my_string(phe.HISTBOX_SNV), phe.HISTBOX_SNV)
        data[core_constants.RESULTS][phe.WHOLE_GENOME_PLOT] = tools.convert_plot(self, wrapper.get_my_string(phe.WHOLE_GENOME_PLOT), phe.WHOLE_GENOME_PLOT)
        data[core_constants.RESULTS][phe.CELLULOID_PLOT] = tools.convert_plot(self, wrapper.get_my_string(phe.CELLULOID_PLOT), phe.CELLULOID_PLOT)
        data[core_constants.RESULTS][phe.INDEL_VAF_PLOT] = tools.convert_svg_plot(self, wrapper.get_my_string(phe.INDEL_VAF_PLOT), phe.INDEL_VAF_PLOT)
        data[core_constants.RESULTS][phe.SNV_VAF_PLOT] = tools.convert_svg_plot(self, wrapper.get_my_string(phe.SNV_VAF_PLOT), phe.SNV_VAF_PLOT)
        data[core_constants.RESULTS][phe.INDEL_BIN_PLOT] = tools.convert_svg_plot(self, wrapper.get_my_string(phe.INDEL_BIN_PLOT), phe.INDEL_BIN_PLOT)
        data[core_constants.RESULTS][phe.SV_BIN_PLOT] = tools.convert_svg_plot(self, wrapper.get_my_string(phe.SV_BIN_PLOT), phe.SV_BIN_PLOT)

        data[core_constants.RESULTS]['loads']['snv_percentile'] = self.get_percentile(data[core_constants.RESULTS]['loads']['snv_count'],  wrapper.get_my_string('comparison_cohort_file'), 'snv_count')
        data[core_constants.RESULTS]['loads']['sv_percentile'] = self.get_percentile(data[core_constants.RESULTS]['loads']['sv_count'],  wrapper.get_my_string('comparison_cohort_file'), 'sv_count')
        data[core_constants.RESULTS]['loads']['indel_percentile'] = self.get_percentile(data[core_constants.RESULTS]['loads']['indel_count'], wrapper.get_my_string('comparison_cohort_file'), 'indel_count')
        data[core_constants.RESULTS]['template_type'] = '_'.join((wrapper.get_my_string('template_type'), self.TEMPLATE_NAME))
        data[core_constants.RESULTS][phe.ONCOSLIDE_CNV_PLOT] = tools.convert_svg_plot(self, wrapper.get_my_string(phe.ONCOSLIDE_CNV_PLOT),   phe.ONCOSLIDE_CNV_PLOT)

        return data

    def render(self, data):
        renderer = mako_renderer(self.get_module_dir())
        template_name = data[core_constants.RESULTS]['template_type'] 
        return renderer.render_name(template_name, data)

    def specify_params(self):
        discovered = [
            phe.SEX_PATH,
            phe.COSMIC_SIGNNLS_PATH,
            phe.PARAM_PATH,
            phe.SAMPLE_VARIANTS_FILE,
            phe.WHOLE_GENOME_PLOT,
            phe.CELLULOID_PLOT,
            phe.INDEL_VAF_PLOT,
            phe.INDEL_BIN_PLOT,
            phe.SV_BIN_PLOT,
            phe.SNV_VAF_PLOT,
            phe.ONCOSLIDE_CNV_PLOT,
            phe.SNV_CONTEXT_PLOT,
            phe.HISTBOX_INDEL,
            phe.HISTBOX_SV,
            phe.HISTBOX_SNV,
            phe.SUMMARY_FILE_PATH,
            phe.CELLULOID_DIR,
            'genes_of_interest_file',
            'template_type',
            'comparison_cohort_file'
        ]
        for key in discovered:
            self.add_ini_discovered(key)
        self.set_ini_default(core_constants.ATTRIBUTES, 'research')
        self.set_priority_defaults(self.PRIORITY)



    def get_percentile(self, sample_variant_count, cohort_path, variant_column_name):
        these_counts = []
        cohort_path = tools.check_path_exists(self, cohort_path)
        with open(cohort_path, 'r') as cohort_file:
            for row in csv.DictReader(cohort_file, delimiter=","):
                these_counts.append(int(row[variant_column_name]))
        cohort_count = {}
        for this_count in these_counts:
            if this_count in cohort_count.keys():
                cohort_count[this_count] = cohort_count[this_count] + 1
            else:
                cohort_count[this_count] = 1
        sorted_cohort_count = collections.OrderedDict(sorted(cohort_count.items()))
        
        lowerCount = 0
        equalCount = 0
        for this_count in sorted_cohort_count.keys():
            if int(sample_variant_count) > this_count:
                lowerCount = lowerCount + sorted_cohort_count[this_count]
            elif int(sample_variant_count) == this_count:
                equalCount = equalCount + sorted_cohort_count[this_count]
        rank = round(((lowerCount + (0.5 * equalCount)) / (len(these_counts) + 1)*100))
        return(rank )

