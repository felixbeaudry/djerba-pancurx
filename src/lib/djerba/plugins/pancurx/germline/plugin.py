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

    PRIORITY = 300
    PLUGIN_VERSION = '1.0.0'
    TEMPLATE_NAME = 'germline_template.html'


    def configure(self, config):
        config = self.apply_defaults(config)
        wrapper = self.get_config_wrapper(config)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.SUMMARY_FILE_PATH, phe.SUMMARY_FILE_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.SAMPLE_VARIANTS_FILE, phe.SAMPLE_VARIANTS_FILE, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, phe.PARAM_PATH, phe.PARAM_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, 'template_type', 'template_type', core_constants.DEFAULT_SAMPLE_INFO)
        wrapper = tools.fill_file_if_null(self, wrapper, 'germline_genes_of_interest_file', 'germline_genes_of_interest_file', core_constants.DEFAULT_SAMPLE_INFO)

        return wrapper.get_config()

    def extract(self, config):
        wrapper = self.get_config_wrapper(config)
        attributes = wrapper.get_my_attributes()
        self.check_attributes_known(attributes)
        data = self.get_starting_plugin_data(wrapper, self.PLUGIN_VERSION)
        summary_results = tools.parse_summary_file(self, wrapper.get_my_string(phe.SUMMARY_FILE_PATH))
        germ_variant_count, germ_nonsilent_count = self.get_germline_variant_counts(summary_results)

        genes_of_interest = tools.get_genes_of_interest(self, wrapper.get_my_string('germline_genes_of_interest_file'))
        germline_variants = self.parse_germline_variants(wrapper.get_my_string(phe.SAMPLE_VARIANTS_FILE))
        germ_nonsil_genes, germ_nonsil_genes_rare, germ_pathogenic, reportable_germline_variants = self.get_subset_of_germline_variants(germline_variants, genes_of_interest)
        ploidy = tools.parse_celluloid_params(self, wrapper.get_my_string(phe.PARAM_PATH), "ploidy_numeric")

        results = {
            phe.PLOIDY: ploidy,
            phe.GERM_VARIANT_COUNT : germ_variant_count,
            phe.GERM_NONSILENT_COUNT: germ_nonsilent_count,
            phe.GERM_NONSIL_SUBSET_COUNT: germ_nonsil_genes,
            phe.GERM_NONSIL_SUBSET_RARE_COUNT: germ_nonsil_genes_rare,
            phe.GERM_PATHOGENIC_COUNT: germ_pathogenic,
            'reportable_germline_variants': reportable_germline_variants
        }
        data[core_constants.RESULTS] = results
        data[core_constants.RESULTS]['template_type'] = '_'.join((wrapper.get_my_string('template_type'), self.TEMPLATE_NAME))

        return data

    def render(self, data):
        renderer = mako_renderer(self.get_module_dir())
        template_name = data[core_constants.RESULTS]['template_type'] 
        return renderer.render_name(template_name, data)

    def specify_params(self):
        discovered = [
            phe.PARAM_PATH,
            phe.SAMPLE_VARIANTS_FILE,
            phe.SUMMARY_FILE_PATH,
            'template_type',
            'germline_genes_of_interest_file'
        ]
        for key in discovered:
            self.add_ini_discovered(key)
        self.set_ini_default(core_constants.ATTRIBUTES, 'research')
        self.set_priority_defaults(self.PRIORITY)

    def parse_germline_variants(self, sample_variants_file_path):
        sample_variants = {}
        sample_variants_file_path = tools.check_path_exists(self, sample_variants_file_path)
        with open(sample_variants_file_path, 'r') as sample_variants_file:
            for row in csv.DictReader(sample_variants_file, delimiter=","):
                ## only saving non-silent germline variants to save on memory
                if 'germline' in row['mutation_class'] and row['mutation_type'] in phe.NONSILENT_CHANGES:
                    gene = row['gene']
                    variant_key = f"{row['mutation_type']},{row['position']},{row['base_change']}"

                    variant_characteristics = {
                        'gene': gene,
                        'mutation_type': row['mutation_type'],
                        'mutation_class': row['mutation_class'],
                        'rarity': row['rarity'],
                        'clinvar': row['clinvar'],
                        'dbsnp': row['dbsnp'],
                        'copy_number': row['copy_number'],
                        'ab_counts': row['ab_counts'],
                        'cosmic_census_flag': row['cosmic_census_flag'],
                        'nuc_context': row['nuc_context'],
                        'aa_context': row['aa_context'],
                    }
                    if gene in sample_variants:
                        sample_variants[gene][variant_key] = variant_characteristics
                    else:
                        sample_variants[gene] = {}
                        sample_variants[gene][variant_key] = variant_characteristics
        return(sample_variants) 

    def get_subset_of_germline_variants(self, sample_variants, gene_order):
        germline_nonsilent_gene_count = 0
        germ_nonsil_genes_rare = 0
        germ_pathogenic = 0
        reportable_germline_variants = []
        for gene in gene_order:
            if gene in sample_variants:
                for variant in sample_variants[gene]:
                    germline_nonsilent_gene_count = germline_nonsilent_gene_count + 1
                    if sample_variants[gene][variant]['rarity'] != "common":
                        reportable_germline_variants.append(sample_variants[gene][variant])
                        germ_nonsil_genes_rare += 1
                        mutation_type = sample_variants[gene][variant]['mutation_type']
                        clinvar = sample_variants[gene][variant]['clinvar'] 
                        if (mutation_type in ["frameshift", "stopgain"] or clinvar.startswith('CLINSIG=pathogenic')) and \
                            sample_variants[gene][variant]['dbsnp'] != "NA" and \
                            sample_variants[gene][variant]['dbsnp'] not in phe.EXCLUDED_GERMLINE_PATHOGENIC_VARIANTS :
                            germ_pathogenic += 1
        return(germline_nonsilent_gene_count, germ_nonsil_genes_rare, germ_pathogenic, reportable_germline_variants)

    def get_germline_variant_counts(self, summary_results):
        germ_variant_count = int(summary_results.get("germline_snv_count")) + int(summary_results.get("germline_indel_count"))
        germ_nonsilent_count = int(summary_results.get("germline_missense_count")) + int(summary_results.get("germline_nonsense_count"))
        return(germ_variant_count, germ_nonsilent_count)

