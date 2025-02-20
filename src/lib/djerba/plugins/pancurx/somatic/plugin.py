"""Plugin to generate the "Gremline" report section"""

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

    PRIORITY = 100
    PLUGIN_VERSION = '1.0.0'
    TEMPLATE_NAME = 'somatic_template.html'


    def configure(self, config):
        config = self.apply_defaults(config)
        wrapper = self.get_config_wrapper(config)

        wrapper = tools.fill_file_if_null(self, wrapper, 'template_type', 'template_type', core_constants.DEFAULT_SAMPLE_INFO)

        wrapper = tools.fill_file_if_null(self, wrapper, 'genes_of_interest_file', 'genes_of_interest_file', core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, 'comparison_cohort_file', 'comparison_cohort_file', core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, 'signatures_of_interest_file', 'signatures_of_interest_file', core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, 'germline_genes_of_interest_file', 'germline_genes_of_interest_file', core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, 'sites_of_interest_file', 'sites_of_interest_file', core_constants.DEFAULT_PATH_INFO)

        wrapper = tools.fill_file_if_null(self, wrapper, phe.SAMPLE_VARIANTS_FILE, phe.SAMPLE_VARIANTS_FILE, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.SUMMARY_FILE_PATH, phe.SUMMARY_FILE_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.PARAM_PATH, phe.PARAM_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.CELLULOID_DIR, phe.CELLULOID_DIR, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.COSMIC_SIGNNLS_PATH, phe.COSMIC_SIGNNLS_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.SEX_PATH, phe.SEX_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.SNV_PATH, phe.SNV_PATH, core_constants.DEFAULT_PATH_INFO)

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
        germline_genes_of_interest = tools.get_genes_of_interest(self, wrapper.get_my_string('germline_genes_of_interest_file'))
        all_genes_of_interest = genes_of_interest | germline_genes_of_interest

        extra_sites_call = tools.parse_extra_sites(self, wrapper.get_my_string(phe.SNV_PATH), wrapper.get_my_string('sites_of_interest_file'))

        somatic_variants = tools.parse_somatic_variants(self, wrapper.get_my_string(phe.SAMPLE_VARIANTS_FILE), all_genes_of_interest, extra_sites_call)

        mane_transcript_path = os.path.join(phe.DEFAULT_DATA_LOCATION, phe.DEFAULT_MANE_FILE)
        mane_transcripts = tools.parse_mane_transcript(self, mane_transcript_path)
        

        reportable_variants, tier_mode = tools.get_subset_of_somatic_variants(self, somatic_variants, genes_of_interest, mane_transcripts)
        data[core_constants.RESULTS]['reportable_variants'] = reportable_variants
        data[core_constants.RESULTS]['tier_mode'] = tier_mode
        summary_results = tools.parse_summary_file(self, wrapper.get_my_string(phe.SUMMARY_FILE_PATH))
        data[core_constants.RESULTS]['loads'] = tools.get_loads_from_summary(summary_results)
        data[core_constants.RESULTS]['sigs'] = tools.parse_cosmic_signatures(self, wrapper.get_my_string(phe.COSMIC_SIGNNLS_PATH), wrapper.get_my_string('signatures_of_interest_file'))
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

        data[core_constants.RESULTS]['loads']['snv_percentile'] = tools.get_percentile(self, data[core_constants.RESULTS]['loads']['snv_count'],  wrapper.get_my_string('comparison_cohort_file'), 'snv_count')
        data[core_constants.RESULTS]['loads']['sv_percentile'] = tools.get_percentile(self, data[core_constants.RESULTS]['loads']['sv_count'],  wrapper.get_my_string('comparison_cohort_file'), 'sv_count')
        data[core_constants.RESULTS]['loads']['indel_percentile'] = tools.get_percentile(self, data[core_constants.RESULTS]['loads']['indel_count'], wrapper.get_my_string('comparison_cohort_file'), 'indel_count')
        data[core_constants.RESULTS]['template_type'] = '_'.join((wrapper.get_my_string('template_type'), self.TEMPLATE_NAME))
        data[core_constants.RESULTS][phe.ONCOSLIDE_CNV_PLOT] = tools.convert_svg_plot(self, wrapper.get_my_string(phe.ONCOSLIDE_CNV_PLOT),   phe.ONCOSLIDE_CNV_PLOT)

        cnvs_and_abs = tools.subset_and_deduplicate(somatic_variants)
        self.workspace.write_json('cnvs_and_abs.json', cnvs_and_abs)

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
            phe.SNV_PATH,
            'genes_of_interest_file',
            'template_type',
            'comparison_cohort_file',
            'signatures_of_interest_file',
            'germline_genes_of_interest_file',
            'sites_of_interest_file'
        ]
        for key in discovered:
            self.add_ini_discovered(key)
        self.set_ini_default(core_constants.ATTRIBUTES, 'research')
        self.set_priority_defaults(self.PRIORITY)


