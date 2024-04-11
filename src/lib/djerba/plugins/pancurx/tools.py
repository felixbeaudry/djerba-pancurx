"""pancurx supporting functions"""
import csv
from decimal import Decimal
import logging
import re
import os
from djerba.util.image_to_base64 import converter
import djerba.plugins.pancurx.constants as phe


def add_underscore_to_donor(donor):
    if "EPPIC" in donor:
        donor = donor.replace("EPPIC", "EPPIC_")
    else:
        donor = re.sub(r'^([A-Z0-9]+)(....)', r'\1_\2', donor)
    return(donor)

def check_path_exists(self, path):
    if not os.path.exists(path):
        msg = "Cannot find file: {0}".format(path)
        self.logger.error(msg)
        raise RuntimeError(msg)
    else:
        return(path)
    
def convert_plot(self, plot_path, plot_name):
    """Read VAF plot from file and return as a base64 string"""
    image_converter = converter(self.log_level, self.log_path)
    converted_plot = image_converter.convert_png(plot_path, plot_name)
    return converted_plot

def convert_svg_plot(self, plot_path, plot_name):
    """Read VAF plot from file and return as a base64 string"""
    image_converter = converter(self.log_level, self.log_path)
    converted_plot = image_converter.convert_svg(plot_path, plot_name)
    return converted_plot

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

def fill_categorized_file_if_null(self, wrapper, file_type_name, ini_param, path_info, category):
    if wrapper.my_param_is_null(ini_param):
        if self.workspace.has_file(path_info):
            path_info = self.workspace.read_json(path_info)
            workflow_paths = path_info.get(category)
            workflow_path = workflow_paths[file_type_name]
            if workflow_path == None:
                msg = 'Cannot find {0}'.format(ini_param)
                self.logger.error(msg)
                raise RuntimeError(msg)
            wrapper.set_my_param(ini_param, workflow_path)
    return(wrapper)

def get_subset_of_somatic_variants(self, sample_variants, subset_order):
    reportable_variants = []
    for gene in subset_order :
        if gene in sample_variants:
            for variant in sample_variants[gene]:
                reportable_variants.append(sample_variants[gene][variant])

    return(reportable_variants)

def get_all_somatic_variants(self, sample_variants):
    all_variants = []
    sample_variants_sorted = dict(sorted(sample_variants.items()))
    for gene in sample_variants_sorted:
        for variant in sample_variants_sorted[gene]:
            if sample_variants_sorted[gene][variant]['mutation_type'] != 'NA':
                all_variants.append(sample_variants_sorted[gene][variant])
    return(all_variants)

def get_tissue_from_sample_id(sample_name):
    tissue = "NA"
    if re.match(r'^[A-Z0-9]+_...._(..)_', sample_name):
        match_group = re.match(r'^[A-Z0-9]+_...._(..)_', sample_name).group(1)
        tissue = phe.LIMS_TISSUE_CODES.get(match_group, match_group)
    return(tissue)

def parse_celluloid_params(self, celluloid_params_file_path, ploidy_or_cellularity):
    celluloid_params_file_path = check_path_exists(self, celluloid_params_file_path)
    with open(celluloid_params_file_path, 'r') as celluloid_params_file:
        for row in csv.DictReader(celluloid_params_file, delimiter=" "):
            cellularity = "{:.1f}".format(float(row["T1"])*100)
            ploidy = "{:.3f}".format(float(row["Ploidy"]))
    ploidy_string = ''
    if float(ploidy) > phe.DIPLOID_CUTOFF:
        ploidy_string = "polyploid ({0})".format(ploidy)
    elif float(ploidy) <= phe.DIPLOID_CUTOFF:
        ploidy_string = "diploid ({0})".format(ploidy)
    else:
        msg = "Ploidy {0} not a number".format(ploidy)
        self.logger.error(msg)
        raise RuntimeError(msg)
    if ploidy_or_cellularity == "ploidy":
        return(ploidy_string)
    elif ploidy_or_cellularity == "ploidy_numeric":
        return(float(ploidy))
    else:
        return(cellularity)

def parse_cosmic_signatures(self, cosmic_signatures_file_path):
    signature_results = {}
    row = {}
    cosmic_signatures_file_path = check_path_exists(self, cosmic_signatures_file_path)
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

def parse_coverage(self, coverage_file_path):
    coverage_file_path = check_path_exists(self, coverage_file_path)
    with open(coverage_file_path, 'r') as file:
        ## Number of columns is variable among rows
        lines = file.read().splitlines()
        for column_number in range(len(lines)):
            column = lines[column_number]
            if "GENOME_TERRITORY" in column:
                genome, mean_coverage, sd, median_coverage, *extra = lines[column_number + 1].split('\t')
                median_coverage = float(median_coverage)
    return(median_coverage)

def get_genes_of_interest(self, genes_of_interest_file_path):
    genes_of_interest = []
    genes_of_interest_file_path = check_path_exists(self, genes_of_interest_file_path)
    with open(genes_of_interest_file_path, 'r') as genes_of_interest_file:
        for row in csv.reader(genes_of_interest_file, delimiter="\t"):
            genes_of_interest.append(row[0])
    return(genes_of_interest)

def parse_somatic_variants(self, sample_variants_file_path, gene_list = []):
    sample_variants = {}
    sample_variants_file_path = check_path_exists(self, sample_variants_file_path)
    with open(sample_variants_file_path, 'r') as sample_variants_file:
        for row in csv.DictReader(sample_variants_file, delimiter=","):
            ## only saving non-silent germline variants to save on memory
            if 'somatic' in row['mutation_class']:
               
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
                    'position': row['position']
                }
                # add variant to gene or create gene
                if gene in sample_variants:
                    sample_variants[gene][variant_key] = variant_characteristics
                else:
                    sample_variants[gene] = {}
                    sample_variants[gene][variant_key] = variant_characteristics
            elif row['gene'] in gene_list and ('NA' in row['mutation_class'] or 'germline' in row['mutation_class']):
                gene = row['gene']
                variant_key = "No Variant"

                variant_characteristics = {
                    'gene': gene,
                    'mutation_type': 'NA',
                    'mutation_class': 'NA',
                    'rarity': 'NA',
                    'clinvar': 'NA',
                    'dbsnp': 'NA',
                    'copy_number': row['copy_number'],
                    'ab_counts': row['ab_counts'],
                    'cosmic_census_flag': 'NA',
                    'nuc_context': 'NA',
                    'aa_context': 'NA',
                    'position': 'NA'
                }

                sample_variants[gene] = {}
                sample_variants[gene][variant_key] = variant_characteristics

    return(sample_variants) 

def parse_summary_file(self, summary_file_path):
    row = {}
    summary_file_path = check_path_exists(self, summary_file_path)
    with open(summary_file_path, 'r') as summary_file:
        line = summary_file.readline().strip()
        header = line.split(',')
        line = summary_file.readline().strip()
        row = dict(zip(header, line.split(',')))
    return(row)

    

    