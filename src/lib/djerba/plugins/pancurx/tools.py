"""pancurx supporting functions"""
import csv
from decimal import Decimal
import logging
import re
import os
import requests
import json
from djerba.util.image_to_base64 import converter
import djerba.plugins.pancurx.constants as phe
import shutil

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

def copy_if_not_exists(source_path, destination_path):
    if not os.path.exists(destination_path):
        shutil.copyfile(source_path, destination_path)

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

def find_external_id_in_json_dict(lims_dict, this_donor):
    external_id = "NA"
    for donor in range(len(lims_dict)):
        if lims_dict[donor]['name'] == this_donor:
            for attribute in range(len(lims_dict[donor]['attributes'])):
                if lims_dict[donor]['attributes'][attribute]['name'] == "External Name":
                    external_id = lims_dict[donor]['attributes'][attribute]['value']
                    external_id = external_id.replace(',', ' ')
    return(external_id)

def get_all_somatic_variants(self, sample_variants):
    all_variants = []
    sample_variants_sorted = dict(sorted(sample_variants.items()))
    for gene in sample_variants_sorted:
        for variant in sample_variants_sorted[gene]:
            if sample_variants_sorted[gene][variant]['mutation_type'] != 'NA':
                all_variants.append(sample_variants_sorted[gene][variant])
    return(all_variants)

def get_genes_of_interest(self, genes_of_interest_file_path):
    genes_of_interest = []
    genes_of_interest_file_path = check_path_exists(self, genes_of_interest_file_path)
    with open(genes_of_interest_file_path, 'r') as genes_of_interest_file:
        for row in csv.reader(genes_of_interest_file, delimiter="\t"):
            genes_of_interest.append(row[0])
    return(genes_of_interest)

def get_germline_variant_counts( summary_results):
    germ_variant_count = int(summary_results.get("germline_snv_count")) + int(summary_results.get("germline_indel_count"))
    germ_nonsilent_count = int(summary_results.get("germline_missense_count")) + int(summary_results.get("germline_nonsense_count"))
    return(germ_variant_count, germ_nonsilent_count)

def get_loads_from_summary(summary_file):
    loads = {
        'snv_count' : int(summary_file.get("snv_count")),
        'indel_count' : int(summary_file.get("indel_count")),
        'sv_count' : int(summary_file.get("sv_count")),
        'TMB' : round((int(summary_file.get("snv_count")) + int(summary_file.get("indel_count"))) / 3000, 2),
        'nonsyn_count' : summary_file.get('nonsyn_count'),
        'del_frameshift_count' : summary_file.get('del_frameshift_count'),
        'sv_del_bp_gene_count': summary_file.get('sv_del_bp_gene_count'),
        'stopgain_count': summary_file.get('stopgain_count'),
        'del_nonframeshift_count': summary_file.get('del_nonframeshift_count'),
        'sv_dup_bp_gene_count': summary_file.get('sv_dup_bp_gene_count'),
        'stoploss_count': summary_file.get('stoploss_count'),
        'ins_frameshift_count': summary_file.get('ins_frameshift_count'),
        'sv_inv_bp_gene_count': summary_file.get('sv_inv_bp_gene_count'),
        'splice_count': summary_file.get('splice_count'),
        'ins_nonframeshift_count': summary_file.get('ins_nonframeshift_count'),
        'sv_tra_bp_gene_count': summary_file.get('sv_tra_bp_gene_count'),
        'total_del_snv' :int(summary_file.get('splice_count')) + 
                            int(summary_file.get('stoploss_count')) +
                                int(summary_file.get('stopgain_count')) + 
                                int(summary_file.get('nonsyn_count')),
        'total_del_indel' : int(summary_file.get('ins_nonframeshift_count') )+ 
                            int(summary_file.get('ins_frameshift_count') )+ 
                            int( summary_file.get('del_nonframeshift_count')) + 
                            int(summary_file.get('del_frameshift_count')),
        'total_del_sv': int(summary_file.get('sv_dup_bp_gene_count') )+ 
                        int(summary_file.get('sv_del_bp_gene_count')) + 
                        int(summary_file.get('sv_inv_bp_gene_count')) + 
                        int(summary_file.get('sv_tra_bp_gene_count'))
    }
    loads['TMB_quantified'] = get_load_quantified(loads['TMB'], phe.TMB_RANGE_CUTOFF)
    loads['snv_quantified'] = get_load_quantified(loads['snv_count'], phe.SNV_RANGE_CUTOFF)
    loads['indel_quantified'] = get_load_quantified(loads['indel_count'], phe.INDEL_RANGE_CUTOFF)
    loads['sv_quantified'] = get_load_quantified(loads['sv_count'], phe.SV_RANGE_CUTOFF)

    return(loads)

def get_load_quantified(sample_value, quantification_dictionary):
    sample_quantification = ''
    for this_quantification, this_cutoff in quantification_dictionary.items():
        if this_cutoff < sample_value:
            sample_quantification = this_quantification
            break
    quantification_string = "".join((str(sample_value), " (",sample_quantification,")"))
    return(quantification_string)


def get_subset_of_germline_variants(sample_variants, gene_order):
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
                    if (clinvar.startswith('CLINSIG=pathogenic')) and \
                        sample_variants[gene][variant]['dbsnp'] != "NA" and \
                        sample_variants[gene][variant]['dbsnp'] not in phe.EXCLUDED_GERMLINE_PATHOGENIC_VARIANTS :
                        print("----germline pathogenic----", gene, variant)
                        germ_pathogenic += 1
    return(germline_nonsilent_gene_count, germ_nonsil_genes_rare, germ_pathogenic, reportable_germline_variants)

def get_subset_of_somatic_variants(self, sample_variants, subset_order):
    reportable_variants = []
    for gene in subset_order :
        if gene in sample_variants:
            for variant in sample_variants[gene]:
                reportable_variants.append(sample_variants[gene][variant])

    return(reportable_variants)

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
            ploidy = "{:.2f}".format(float(row["Ploidy"]))
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
            signature_results[phe.COSMIC_SIGNATURE_SET[signature_name]] = float(row[signature_name])
    sigs_total = sum(signature_results.values(), )
    sigs = {key: round( value / sigs_total , 2) for key, value in signature_results.items()}

    return(sigs)

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

def parse_fusions(self, fusions_file_path, min_split_reads = 8):
    fusion_count = 0
    sample_fusions = []
    fusions_file_path = check_path_exists(self, fusions_file_path)
    with open(fusions_file_path, 'r') as fusions_file:
        for row in csv.DictReader(fusions_file, delimiter="\t"):
            if int(row['break1_split_reads']) >= min_split_reads and \
                int(row['break2_split_reads']) >= min_split_reads and \
                row['gene1_aliases'] != 'None' and \
                row['gene2_aliases'] != 'None' :  
                fusion_characteristics = {
                    'gene1_aliases': row['gene1_aliases'],
                    'break1_chromosome': row['break1_chromosome'],
                    'gene2_aliases': row['gene2_aliases'],
                    'break2_chromosome': row['break2_chromosome'],
                    'fusion_product': "::".join((row['gene1_aliases'], row['gene2_aliases'])),
                    'gene_product_type': row['gene_product_type'],
                    'event_type': row['event_type'],
                    'fusion_splicing_pattern': row['fusion_splicing_pattern'],
                    'break1_split_reads': row['break1_split_reads'],
                    'break2_split_reads': row['break2_split_reads'],
                    'linking_split_reads': row['linking_split_reads'],                    
                }
                fusion_count = fusion_count + 1
                sample_fusions.append(fusion_characteristics)
    return(sample_fusions, fusion_count)


def parse_germline_variants(self, sample_variants_file_path):
    sample_variants = {}
    sample_variants_file_path = check_path_exists(self, sample_variants_file_path)
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

def parse_lims(self, donor):
    donor = add_underscore_to_donor(donor)
    r = requests.get(phe.DEFAULT_CORE_LIMS_URL, allow_redirects=True)
    if r.status_code == 404:
        msg = "Trouble pulling LIMS"
        raise MissingLIMSError(msg)
    else:
        lims_json_dict = json.loads(r.text)
        external_ids = find_external_id_in_json_dict(lims_json_dict , donor)
    return(external_ids)

def parse_multifactor_marker(self, summary_results, html_headers, marker_cutoffs):
    result = []
    hallmark_tally = 0
    for this_key in summary_results:
        this_reporting_name = html_headers.get(this_key)
        this_value = summary_results.get(this_key)
        if this_key in marker_cutoffs:
            this_value = float(this_value)
            this_cutoff = marker_cutoffs.get(this_key)
            if this_reporting_name == 'SNV C>T Ratio <':
                if this_value < float(this_cutoff):
                    above_cutoff = True
                    hallmark_tally = hallmark_tally + 1
                else:
                    above_cutoff = False
            elif this_value > float(this_cutoff):
                above_cutoff = True
                hallmark_tally = hallmark_tally + 1
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
                hallmark_tally = hallmark_tally + 1
                this_value = re.sub(r'\|', ', ', this_value)
            result_tmp = {
                phe.REPORTING_NAME: this_reporting_name,
                phe.VALUE : this_value,
                phe.ABOVE_CUTOFF: above_cutoff
            }
            result.append(result_tmp)
    return(result, hallmark_tally)

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

def parse_TDP(self, TDP_file_path):
    row = {}
    TDP_file_path = check_path_exists(self, TDP_file_path)
    with open(TDP_file_path, 'r') as TDP_file:
        line = TDP_file.readline().strip()
        header = line.split(',')
        line = TDP_file.readline().strip()
        row = dict(zip(header, line.split(',')))
    return(round(float(row["score"]),2))


def try_two_null_files(self, wrapper, workflow_name, ini_param, path_info, first_file):
    if wrapper.my_param_is_null(ini_param):
        this_file = os.path.join(self.workspace.print_location(), first_file )
        if self.workspace.has_file(this_file):
            wrapper.set_my_param(ini_param, this_file)
        elif self.workspace.has_file(path_info):
            path_info = self.workspace.read_json(path_info)
            workflow_path = path_info.get(workflow_name)
            if workflow_path == None:
                msg = 'Cannot find {0}'.format(ini_param)
                self.logger.error(msg)
                raise RuntimeError(msg)
            wrapper.set_my_param(ini_param, workflow_path)
    return(wrapper)



class MissingLIMSError(Exception):
    pass
