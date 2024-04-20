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
        wrapper = tools.fill_file_if_null(self, wrapper, phe.TUMOUR_ID, phe.TUMOUR_ID, core_constants.DEFAULT_SAMPLE_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.DONOR_ID, phe.DONOR_ID, core_constants.DEFAULT_SAMPLE_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.NORMAL_ID, phe.NORMAL_ID, core_constants.DEFAULT_SAMPLE_INFO)

        wrapper = tools.fill_file_if_null(self, wrapper, phe.SAMPLE_VARIANTS_FILE, phe.SAMPLE_VARIANTS_FILE, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.PARAM_PATH, phe.PARAM_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.TDP_PATH, phe.TDP_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.COSMIC_SIGNNLS_PATH, phe.COSMIC_SIGNNLS_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.SUMMARY_FILE_PATH, phe.SUMMARY_FILE_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.CELLULOID_PLOT, phe.CELLULOID_PLOT, core_constants.DEFAULT_PATH_INFO)

        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'indel.stack_count', phe.INDEL_BIN_PLOT, core_constants.DEFAULT_PATH_INFO, 'svg_plots')
        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'sv.bins', phe.SV_BIN_PLOT, core_constants.DEFAULT_PATH_INFO, 'svg_plots')
        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'SBS_context_bar', phe.SNV_CONTEXT_PLOT, core_constants.DEFAULT_PATH_INFO, 'svg_plots')
        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'COSMIC_stack', phe.COSMIC_STACK_PLOT, core_constants.DEFAULT_PATH_INFO, 'svg_plots')
        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'sv.stack_count', phe.SV_STACK_PLOT, core_constants.DEFAULT_PATH_INFO, 'svg_plots')
        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'stack-tmb', phe.TMB_STACK_PLOT, core_constants.DEFAULT_PATH_INFO, 'svg_plots')
        wrapper = tools.fill_file_if_null(self, wrapper, 'genes_of_interest_file', 'genes_of_interest_file', core_constants.DEFAULT_SAMPLE_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, 'germline_genes_of_interest_file', 'germline_genes_of_interest_file', core_constants.DEFAULT_SAMPLE_INFO)

        if wrapper.my_param_is_null(phe.SLIDE_DATE):
            wrapper.set_my_param(phe.SLIDE_DATE, phe.NONE_SPECIFIED)
        if wrapper.my_param_is_null(phe.MOFFIT):
            wrapper.set_my_param(phe.MOFFIT, 'NA')

        if wrapper.my_param_is_null(phe.EXTERNAL_IDS):
            external_ids = tools.parse_lims(self, wrapper.get_my_string(phe.DONOR_ID))
            wrapper.set_my_param(phe.EXTERNAL_IDS, external_ids)

        if wrapper.my_param_is_null(phe.TUMOUR_TISSUE):
            tissue = tools.get_tissue_from_sample_id(wrapper.get_my_string(phe.TUMOUR_ID))
            wrapper.set_my_param(phe.TUMOUR_TISSUE, tissue)
        if wrapper.my_param_is_null(phe.NORMAL_TISSUE):
            tissue = tools.get_tissue_from_sample_id(wrapper.get_my_string(phe.NORMAL_ID))
            wrapper.set_my_param(phe.NORMAL_TISSUE, tissue)

        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'tumour', phe.TUMOUR_COVERAGE_PATH, core_constants.DEFAULT_PATH_INFO, 'coveragePaths')
        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'normal', phe.NORMAL_COVERAGE_PATH, core_constants.DEFAULT_PATH_INFO, 'coveragePaths')
        return wrapper.get_config()

    def extract(self, config):
        wrapper = self.get_config_wrapper(config)
        attributes = wrapper.get_my_attributes()
        self.check_attributes_known(attributes)
        data = self.get_starting_plugin_data(wrapper, self.PLUGIN_VERSION)
        results_keys = [
            phe.DONOR_ID,
            phe.TUMOUR_ID ,
            phe.EXTERNAL_IDS,
            phe.TUMOUR_TISSUE ,
            phe.NORMAL_ID ,
            phe.NORMAL_TISSUE ,
            phe.MOFFIT ,
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
        summary_results = tools.parse_summary_file(self, wrapper.get_my_string(phe.SUMMARY_FILE_PATH))
        data[core_constants.RESULTS]['loads'] = tools.get_loads_from_summary(summary_results)
        data[core_constants.RESULTS]['sigs'] = tools.get_sigs_from_summary(summary_results)

        if wrapper.get_my_string(phe.SLIDE_DATE) == phe.NONE_SPECIFIED:
            data[core_constants.RESULTS][phe.SLIDE_DATE] = strftime("%a %b %d %H:%M:%S %Y")
        else:
            data[core_constants.RESULTS][phe.SLIDE_DATE] = wrapper.get_my_string(phe.SLIDE_DATE)
        ploidy = tools.parse_celluloid_params(self, wrapper.get_my_string(phe.PARAM_PATH), phe.PLOIDY)
        cellularity = tools.parse_celluloid_params(self, wrapper.get_my_string(phe.PARAM_PATH), phe.CELLULARITY)
        data[core_constants.RESULTS][phe.PLOIDY] = ploidy
        data[core_constants.RESULTS][phe.CELLULARITY] = float(cellularity)
        ploidy = tools.parse_celluloid_params(self, wrapper.get_my_string(phe.PARAM_PATH), "ploidy_numeric")
        data[core_constants.RESULTS]['ploidy_numeric'] = ploidy

        data[core_constants.RESULTS]['tandem_duplicator_phenotype_score'] = tools.parse_TDP(self, wrapper.get_my_string(phe.TDP_PATH))
        summary_results = tools.parse_summary_file(self, wrapper.get_my_string(phe.SUMMARY_FILE_PATH))
        dsbr_results, dsbr_tally = tools.parse_multifactor_marker(self, summary_results, phe.DSBR_HTML_HEADERS, phe.DSBR_DEFAULT_HALLMARK_CUTOFFS)
        data[core_constants.RESULTS][phe.DSBR_RESULTS] = dsbr_results
        data[core_constants.RESULTS]['dsbr_tally'] = dsbr_tally
        mmr_results, mmr_tally = tools.parse_multifactor_marker(self, summary_results, phe.MMR_HTML_HEADERS, phe.MMR_DEFAULT_HALLMARK_CUTOFFS)
        data[core_constants.RESULTS][phe.MMR_RESULTS] = mmr_results
        data[core_constants.RESULTS]['mmr_tally'] = mmr_tally

        genes_of_interest = tools.get_genes_of_interest(self, wrapper.get_my_string('genes_of_interest_file'))
        somatic_variants = tools.parse_somatic_variants(self, wrapper.get_my_string(phe.SAMPLE_VARIANTS_FILE), genes_of_interest)
        data[core_constants.RESULTS]['reportable_variants'] = tools.get_subset_of_somatic_variants(self, somatic_variants, genes_of_interest)

        genes_of_interest = tools.get_genes_of_interest(self, wrapper.get_my_string('germline_genes_of_interest_file'))
        germline_variants = tools.parse_germline_variants(self, wrapper.get_my_string(phe.SAMPLE_VARIANTS_FILE))
        germ_nonsil_genes, germ_nonsil_genes_rare, germ_pathogenic, reportable_germline_variants = tools.get_subset_of_germline_variants(germline_variants, genes_of_interest)

        data[core_constants.RESULTS]['reportable_germline_variants'] = reportable_germline_variants
        data[core_constants.RESULTS][phe.GERM_NONSIL_SUBSET_RARE_COUNT] = germ_nonsil_genes_rare


        return(data)

    def render(self, data):
        renderer = mako_renderer(self.get_module_dir())
        return renderer.render_name(self.TEMPLATE_NAME, data)

    def specify_params(self):
        discovered = [
            phe.DONOR_ID,
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
            phe.MOFFIT ,
            phe.INDEL_BIN_PLOT,
            phe.SV_BIN_PLOT,
            phe.SNV_CONTEXT_PLOT,
            phe.COSMIC_STACK_PLOT,
            phe.SV_STACK_PLOT,
            phe.TMB_STACK_PLOT,
            phe.SAMPLE_VARIANTS_FILE,
            'genes_of_interest_file',
            'germline_genes_of_interest_file'
        ]
        for key in discovered:
            self.add_ini_discovered(key)
        self.set_ini_default(core_constants.ATTRIBUTES, 'slide')
        self.set_priority_defaults(self.PRIORITY)
