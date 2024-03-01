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
import djerba.plugins.pancurx.germline.constants as phe

class main(plugin_base):

    PRIORITY = 300
    PLUGIN_VERSION = '1.0.0'
    TEMPLATE_NAME = 'germline_template.html'
    SAMPLE_VARIANTS_FILE = "sample_variants"
    SUMMARY_FILE_PATH = "summary_file_path"

    def configure(self, config):
        config = self.apply_defaults(config)
        wrapper = self.get_config_wrapper(config)
        wrapper = self.fill_file_if_null(wrapper, self.SUMMARY_FILE_PATH, self.SUMMARY_FILE_PATH, core_constants.DEFAULT_PATH_INFO)
        wrapper = self.fill_file_if_null(wrapper, self.SAMPLE_VARIANTS_FILE, self.SAMPLE_VARIANTS_FILE, core_constants.DEFAULT_PATH_INFO)
        return wrapper.get_config()

    def extract(self, config):
        wrapper = self.get_config_wrapper(config)
        attributes = wrapper.get_my_attributes()
        self.check_attributes_known(attributes)
        data = self.get_starting_plugin_data(wrapper, self.PLUGIN_VERSION)
        summary_results = self.parse_summary_file(wrapper.get_my_string(self.SUMMARY_FILE_PATH))
        germ_variant_count, germ_nonsilent_count = self.get_germline_variant_counts(summary_results)
        germline_variants = self.parse_germline_variants(wrapper.get_my_string(self.SAMPLE_VARIANTS_FILE))
        germ_nonsil_genes, germ_nonsil_genes_rare, germ_pathogenic, reportable_germline_variants = self.get_subset_of_germline_variants(germline_variants)
        results = {
            phe.GERM_VARIANT_COUNT : germ_variant_count,
            phe.GERM_NONSILENT_COUNT: germ_nonsilent_count,
            phe.GERM_NONSIL_SUBSET_COUNT: germ_nonsil_genes,
            phe.GERM_NONSIL_SUBSET_RARE_COUNT: germ_nonsil_genes_rare,
            phe.GERM_PATHOGENIC_COUNT: germ_pathogenic,
            'reportable_germline_variants': reportable_germline_variants
        }
        data[core_constants.RESULTS] = results
        return data

    def render(self, data):
        renderer = mako_renderer(self.get_module_dir())
        return renderer.render_name(self.TEMPLATE_NAME, data)

    def specify_params(self):
        discovered = [
            self.SAMPLE_VARIANTS_FILE,
            self.SUMMARY_FILE_PATH
        ]
        for key in discovered:
            self.add_ini_discovered(key)
        self.set_ini_default(core_constants.ATTRIBUTES, 'clinical')
        self.set_priority_defaults(self.PRIORITY)

    def parse_germline_variants(self, sample_variants_file_path):
        sample_variants = {}
        sample_variants_file_path = self.check_path_exists(sample_variants_file_path)
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

    def get_subset_of_germline_variants(self, sample_variants):
        germline_nonsilent_gene_count = 0
        germ_nonsil_genes_rare = 0
        germ_pathogenic = 0
        reportable_germline_variants = []
        for gene in phe.GERMLINE_GENE_ORDER:
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

    def parse_summary_file(self, summary_file_path):
        row = {}
        summary_file_path = self.check_path_exists(summary_file_path)
        with open(summary_file_path, 'r') as summary_file:
            line = summary_file.readline().strip()
            header = line.split(',')
            line = summary_file.readline().strip()
            row = dict(zip(header, line.split(',')))
        return(row)
    
    def check_path_exists(self, path):
        if not os.path.exists(path):
            msg = "Cannot find file: {0}".format(path)
            self.logger.error(msg)
            raise RuntimeError(msg)
        else:
            return(path)

