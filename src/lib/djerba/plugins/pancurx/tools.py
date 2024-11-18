"""pancurx supporting functions"""
import csv
from decimal import Decimal
import logging
import re
import os
import gzip
import requests
import json
from djerba.util.image_to_base64 import converter
import djerba.plugins.pancurx.constants as phe
from djerba.util.environment import directory_finder
from djerba.util.subprocess_runner import subprocess_runner
import shutil
from scipy import stats
import collections
import time

def add_underscore_to_donor(donor):
    if "EPPIC" in donor:
        donor = donor.replace("EPPIC", "EPPIC_")
    else:
        donor = re.sub(r'^([A-Z0-9]+)(....)', r'\1_\2', donor)
    return(donor)

def calculate_TDP_status(sdf):
    ###CONDITIONS(As discussed with Faiyaz at DART in May 2024)
    # Menghi Score > 0 and duplication ratio > 0.205

    tdp_tally = 0
    if sdf["tdp_score"] != 'NA':
        tdp_score = round(float(sdf["tdp_score"]), 2)
    else: 
        tdp_score= 0

    if tdp_score > 0:
        tdp_score_call = True
        tdp_tally = tdp_tally + 1
    else:
        tdp_score_call = False
    tdp_score_dict = {
        "reporting_name": "TDP score > 0 :", 
        "value": tdp_score, 
        "above_cutoff": tdp_score_call
    }
    
    dup_count = sdf["sv_dup_count"]
    sv_count = sdf["sv_count"]
    dup_ratio = round(int(dup_count)/int(sv_count), 3)
    
    if dup_ratio > 0.205:
        dup_ratio_call = True
        tdp_tally = tdp_tally + 1
    else:
        dup_ratio_call = False

    dup_ratio_dict = {
        "reporting_name": "Dup. ratio > 0.205 :", 
        "value": dup_ratio, 
        "above_cutoff": dup_ratio_call
    }
    
    tdp_results = [tdp_score_dict, dup_ratio_dict]
    return(tdp_results, tdp_tally)

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

def custom_sort_key(this_dictionary):
    this_rank = None
    if this_dictionary['mutation_type'] in phe.NONSILENT_CHANGES_RANK:
        this_rank = phe.NONSILENT_CHANGES_RANK[this_dictionary['mutation_type']]
    else:
        this_rank = len(phe.NONSILENT_CHANGES_RANK)
    return(this_rank)

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
    return(external_id)

def get_all_somatic_variants(self, sample_variants, inferred_sex):
    all_variants = []
    sample_variants_sorted = dict(sorted(sample_variants.items()))
    for gene in sample_variants_sorted:
        for variant in sample_variants_sorted[gene]:
            if sample_variants_sorted[gene][variant]['mutation_type'] != 'NA':
                this_cytoband = sample_variants_sorted[gene][variant]['gene_chr']
                if this_cytoband != '':
                    this_chromosome = list(this_cytoband)[0]
                    if inferred_sex == "XY":
                        all_variants.append(sample_variants_sorted[gene][variant])
                    elif this_chromosome != 'Y':
                        all_variants.append(sample_variants_sorted[gene][variant])
                    else:
                        msg = "Skipping Y chromosome in XX"
                        self.logger.info(msg)
                else:
                    msg = "Skipping unplaced cytoband"
                    self.logger.info(msg)
    sorted_list_of_dicts = sorted(all_variants, key=custom_sort_key)
    return(sorted_list_of_dicts)

def get_genes_of_interest(self, genes_of_interest_file_path):
    genes_of_interest = {}
    genes_of_interest_file_path = check_path_exists(self, genes_of_interest_file_path)
    with open(genes_of_interest_file_path, 'r') as genes_of_interest_file:
        for row in csv.reader(genes_of_interest_file, delimiter="\t"):
            genes_of_interest[row[0]] = row[1]
    return(genes_of_interest)

def get_gene_expression(self, genes_of_interest, input_tpm_path, comparison_cohort_path):
    expression_dict = {}
    input_tpm_path = check_path_exists(self, input_tpm_path)
    comparison_cohort_path = check_path_exists(self, comparison_cohort_path)

    work_dir = self.workspace.get_work_dir()

    if self.workspace.has_file('tpm.txt'):
        msg = "TPM file found"
        self.logger.info(msg)
    else:
        r_command = os.path.join(phe.DEFAULT_ENV_PATH, "lib/R/bin/Rscript")
        r_script = os.path.join(phe.DEFAULT_DJERBA_BIN_PATH, "pancurx/fusions/get_tpm.R")

        cmd = [
            r_command, r_script ,
            '--cohort', comparison_cohort_path,
            '--input', input_tpm_path,
            '--output', os.path.join(work_dir, 'tpm.txt'),
        ]

        runner = subprocess_runner()
        runner.run(cmd, "tpm R script")

    input_tpm_path = os.path.join(work_dir, 'tpm.txt')
    with open(input_tpm_path, 'r') as input_tpm_file:
        for row in csv.DictReader(input_tpm_file, delimiter="\t", fieldnames=['gene','TPM','percentile'] ):
            if row['gene'] in genes_of_interest:
                data = {
                    'gene' : row['gene'],
                    'input_tpm': row['TPM'],
                    'percentile': row['percentile'],
                }
                expression_dict[row['gene']] = data
    with open(comparison_cohort_path, 'r') as comparison_cohort_file:
        for row in csv.DictReader(comparison_cohort_file, delimiter="\t"):
            if row['gene_list'] in expression_dict:
                dict_items = list(row.values())
                cohort_expression_list = dict_items[1:]
                this_tpm = expression_dict[row['gene_list']]['input_tpm']
                expression_dict[row['gene_list']]['cohort_perc'] = round(float(expression_dict[row['gene_list']]['percentile'] ) * 100, 1)
    return(expression_dict)

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


def get_percentile(self, sample_variant_count, cohort_path, variant_column_name):
    these_counts = []
    cohort_path = check_path_exists(self, cohort_path)
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



def get_subset_of_germline_variants(sample_variants, gene_order, chosen_transcripts, cnvs_and_abs):
    germline_nonsilent_gene_count = 0
    germ_nonsil_genes_rare = 0
    germ_pathogenic = 0
    reportable_germline_variants = []
    for gene in gene_order.keys():
        if gene in sample_variants:
            for variant in sample_variants[gene]:
                germline_nonsilent_gene_count = germline_nonsilent_gene_count + 1
                if sample_variants[gene][variant]['rarity'] != "common" and \
                    (sample_variants[gene][variant]['clinvar'] not in [ "Benign", "Benign/Likely_benign", "Likely_benign" ] or gene == 'DPYD'):
                    this_chosen_transcript = chosen_transcripts[gene]
                    if this_chosen_transcript in sample_variants[gene][variant]['aa_change']:
                        chosen_frame = sample_variants[gene][variant]['aa_change'][this_chosen_transcript]
                        sample_variants[gene][variant]['nuc_context'] = chosen_frame['dna_seq']
                        sample_variants[gene][variant]['aa_context'] = chosen_frame['prot_seq']
                    else:
                        for this_transcript in sample_variants[gene][variant]['aa_change']:
                            sample_variants[gene][variant]['nuc_context'] = sample_variants[gene][variant]['aa_change'][this_transcript]['dna_seq']
                            sample_variants[gene][variant]['aa_context'] = sample_variants[gene][variant]['aa_change'][this_transcript]['prot_seq']
                            continue
                        
                    sample_variants[gene][variant]['copy_number'] =  cnvs_and_abs[gene]['copy_number']
                    sample_variants[gene][variant]['ab_counts'] =  cnvs_and_abs[gene]['ab_counts']

                    sample_variants[gene][variant]['cosmic_census_flag'] =  "NA"

                    reportable_germline_variants.append(sample_variants[gene][variant])
                    germ_nonsil_genes_rare += 1
                    clinvar = sample_variants[gene][variant]['clinvar'] 
                    if (clinvar.startswith('Pathogenic')) and \
                        sample_variants[gene][variant]['dbsnp'] != "NA" and \
                        sample_variants[gene][variant]['dbsnp'] not in phe.EXCLUDED_GERMLINE_PATHOGENIC_VARIANTS :
                        print("----germline pathogenic----", gene, variant)
                        germ_pathogenic += 1
    return(germline_nonsilent_gene_count, germ_nonsil_genes_rare, germ_pathogenic, reportable_germline_variants)

def get_subset_of_somatic_variants(self, sample_variants, subset_order, mane_transcripts):
    reportable_variants = []
    core_variant_count = 0
    for gene in subset_order.keys():
        if gene in sample_variants:
            for variant in sample_variants[gene]:
                if variant == "No Variant" and subset_order[gene] in ['driver','discovery']:
                    pass
                else:
                    if variant != "No Variant" and subset_order[gene] in ['driver','action']:
                        core_variant_count += 1
                    if sample_variants[gene][variant]['mutation_class'] != 'EXTRASITE':
                        sample_variants = get_mane_context(sample_variants, gene, variant, mane_transcripts, 'nuc_context')
                        sample_variants = get_mane_context(sample_variants, gene, variant, mane_transcripts, 'aa_context')
                    sample_variants[gene][variant]['tier'] = subset_order[gene]
                    reportable_variants.append(sample_variants[gene][variant])
    if core_variant_count > 6:
        tier_mode = 'driver'
    else:
        tier_mode = 'discovery'
    return(reportable_variants, tier_mode)

def get_mane_context(sample_variants, gene, variant, mane_transcripts, context_key):
    if sample_variants[gene][variant][context_key] not in ['NA',''] and \
        sample_variants[gene][variant] != "No Variant"   :
        
        if gene in mane_transcripts and mane_transcripts[gene] in sample_variants[gene][variant][context_key]  :
            sample_variants[gene][variant][context_key] =   sample_variants[gene][variant][context_key][mane_transcripts[gene]]
        else:
            
            for this_context in sample_variants[gene][variant][context_key]:
                sample_variants[gene][variant][context_key] = sample_variants[gene][variant][context_key][this_context]
                break
    return(sample_variants)


def get_mane_somatic_variants(self, sample_variants,  mane_transcripts):

    for gene in sample_variants:

        for variant in sample_variants[gene]:
            sample_variants = get_mane_context(sample_variants, gene, variant, mane_transcripts, 'nuc_context')
            sample_variants = get_mane_context(sample_variants, gene, variant, mane_transcripts, 'aa_context')

    return(sample_variants)

def get_tissue_from_sample_id(sample_name):
    tissue = "NA"
    if re.match(r'^[A-Z0-9]+_...._(..)_', sample_name):
        match_group = re.match(r'^[A-Z0-9]+_...._(..)_', sample_name).group(1)
        tissue = phe.LIMS_TISSUE_CODES.get(match_group, match_group)
    return(tissue)

    
def parse_annovar_germline_variants(self, sample_variants_file_path):
    RARE_FREQUENCY = 0.01

    sample_variants = {}
    sample_variants_file_path = check_path_exists(self, sample_variants_file_path)
    with gzip.open(sample_variants_file_path, 'rt') as sample_variants_file:
        for row in csv.DictReader(sample_variants_file, delimiter="\t"):
            if row['Func.refGene'] in ['exonic',  'splicing', 'exonic;splicing'] \
                and row['ExonicFunc.refGene'] != 'synonymous SNV':

                    #order of estimators is legacy, could consider mean
                    for frequency_estimator in [ "AF", "ExAC_ALL", "1000g2015aug_all", ]:
                        
                        if frequency_estimator in row:
                            if row[frequency_estimator] in ["NA","."]:
                                pop_frequency = 0
                            else:
                                pop_frequency = float(row[frequency_estimator])
                            if pop_frequency == 0:
                                pop_frequency_string = "Novel"
                            elif pop_frequency > RARE_FREQUENCY:
                                pop_frequency_string = "Common"
                            elif pop_frequency < RARE_FREQUENCY:
                                pop_frequency_string = "Rare"
                            continue
                    if  pop_frequency_string in ["Rare", "Novel"] :      
                        snp_id = '.'.join((row['Gene.refGene'], row['Start']))

                        if row['Func.refGene'] in ['exonic', 'exonic;splicing']:
                            aa_changes = row[ 'AAChange.refGene'].split(",")
                            
                            aa_dict = {}

                            for this_change in aa_changes:
                                this_change_split = this_change.split(':')
                                if len(this_change_split) == 5 and this_change_split[0] == row['Gene.refGene'] :
                                    aa_dict[this_change_split[1]] = {
                                        'exon' : this_change_split[2],
                                        'dna_seq' : this_change_split[3],
                                        'prot_seq' :  this_change_split[4],
                                    }
                            exonic_function = row['ExonicFunc.refGene'].split(" ")
                            mutation_type = exonic_function[0]
                        else:
                            mutation_type = 'splice'
                            aa_changes = row[ 'GeneDetail.refGene'].split(":")
                            
                            aa_dict = {}

                            for this_change in aa_changes:
                                this_change_split = this_change.split(':')
                                if len(this_change_split) == 3 :
                                    aa_dict[this_change_split[1]] = {
                                        'exon' : this_change_split[2],
                                        'dna_seq' : this_change_split[3],
                                        'prot_seq' : 'NA',
                                    }
                        if(row['Start'] == row['End']):
                            mutation_class = ' '.join(('germline', 'snp'))
                        else:
                            mutation_class = ' '.join(('germline', 'indel'))

                        if row[ 'avsnp150']  == 'NA' and  pop_frequency_string == "Novel":
                            dbsnp = pop_frequency_string
                        else:
                            dbsnp = row[ 'avsnp150']

                        results = {
                                'chr' : row['#Chr'],
                                'gene' : row['Gene.refGene'],
                                'rarity' : pop_frequency_string,
                                'start' : row['Start'],
                                'end' : row['End'],
                                'mutation_class' : mutation_class,
                                'dbsnp' : dbsnp,
                                'clinvar' : row['CLNSIG'],
                                'mutation_type' : mutation_type,
                                'aa_change' : aa_dict,
                                
                            }
                        if row['Gene.refGene'] in sample_variants:
                            sample_variants[row['Gene.refGene']][snp_id] = results
                        else:
                            sample_variants[row['Gene.refGene']] = {}
                            sample_variants[row['Gene.refGene']][snp_id] = results



    return(sample_variants)                                  

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

def get_signatures_of_interest(self, cosmic_signature_set_path):
    signatures_of_interest = {}
    cosmic_signature_set_path = check_path_exists(self, cosmic_signature_set_path)
    with open(cosmic_signature_set_path, 'r') as cosmic_signature_set:
        for row in csv.reader(cosmic_signature_set, delimiter="\t"):
            signatures_of_interest[row[0]] = {
                    'sbs' : row[1],
                    'aetiology': row[2]
                }
    return(signatures_of_interest)

def parse_cosmic_signatures(self, cosmic_signatures_file_path, cosmic_signature_set_path):

    row = {}
    cosmic_signatures_file_path = check_path_exists(self, cosmic_signatures_file_path)
    with open(cosmic_signatures_file_path, 'r') as cosmic_signatures_file:
        line = cosmic_signatures_file.readline()
        header = line.split()
        line = cosmic_signatures_file.readline()
        for header_position, signature_value in enumerate(line.split()):
            if header[header_position] not in ['sampleName', 'Residuals', 'AbsResiduals', 'RSS', 'N.Mutations']:
                row[header[header_position]] = float(signature_value)
    sigs_total = sum(row.values())
    signatures_of_interest = get_signatures_of_interest(self, cosmic_signature_set_path)
    for signature_name in signatures_of_interest:
        if signature_name in row:
            signatures_of_interest[signature_name]['count'] = row[signature_name]
            signatures_of_interest[signature_name]['proportion'] = round( row[signature_name] / sigs_total , 2)
    return(signatures_of_interest)

def parse_extra_sites(self, snv_file_path, site_set_path):
    extra_sites_set = []
    with open(site_set_path, 'r') as site_set:
        for row in csv.DictReader(site_set, delimiter="\t"):
            this_site = {
                'chrom'	: row['chrom'],
                'pos'	:row['pos'],
                'type'	:row['type'],
                'variant':row['variant'],
                'element': row['element'],
            }
            extra_sites_set.append(this_site)


    snv_file_path = check_path_exists(self, snv_file_path)
    extra_sites_call = []
    with gzip.open(snv_file_path, 'rt') as snv_file:
        for line in snv_file:
            if not line.startswith('#'):
                columns = line.strip().split('\t')
                for this_site in extra_sites_set:
                    if columns[0] == this_site['chrom'] and columns[1] == this_site['pos']:
                        #take max AF as more likely from tumour
                        first_genotype = float(columns[10].split(':')[2])
                        second_genotype = float(columns[10].split(':')[2])
                        this_site['vaf'] = max((first_genotype, second_genotype))
                        extra_sites_call.append(this_site)
    return(extra_sites_call)

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
                mean_coverage = round(float(mean_coverage), 1)
    return(median_coverage, mean_coverage)

def parse_star_qc(self, star_qc_file_path):
    qc_metrics = {}
    total_unmapped = 0
    fieldnames = ['metric', 'value']
    star_qc_file_path = check_path_exists(self, star_qc_file_path)
    with open(star_qc_file_path, 'r') as star_qc_file:
        for row in csv.DictReader(star_qc_file, delimiter="|", fieldnames=fieldnames):
            if row['value'] != None:
                qc_metrics[row['metric'].strip()] = row['value'].strip()
            if row['metric'].strip() in ['Number of reads unmapped: too many mismatches', 'Number of reads unmapped: too short', 'Number of reads unmapped: other']:
                total_unmapped = total_unmapped + int(row['value'])
    qc_metrics['% of reads unmapped'] = round(( total_unmapped / int(qc_metrics['Number of input reads']) ) * 100, 1)
    if qc_metrics['% of reads unmapped'] <= 50 and int(qc_metrics['Uniquely mapped reads number']) >= 3000000:
        qc_metrics['qc_final'] = 'good'
    else:
        qc_metrics['qc_final'] = 'poor'
    return(qc_metrics)

def filter_fusions(self, all_fusions, subset_order, cnvs_and_abs, gene_expression, min_split_reads = 7):
    reportable_variants = []
    fusion_count = 0

    for this_gene in subset_order.keys():
        fusion_not_found = True
        if this_gene in gene_expression and this_gene in cnvs_and_abs:
            this_variant = {
                'gene' : this_gene,
                'gene_chr' : cnvs_and_abs[this_gene]['gene_chr'],
                'input_tpm' : round(float(gene_expression[this_gene]['input_tpm']), 1),
                'cohort_perc': gene_expression[this_gene]['cohort_perc'], 
                'copy_number' : cnvs_and_abs[this_gene]['copy_number'],
            }
            for this_fusion in all_fusions:
                if this_gene == this_fusion['gene1_aliases'] or this_gene == this_fusion['gene2_aliases']:
                    if int(this_fusion['break1_split_reads']) >= min_split_reads and int(this_fusion['break2_split_reads']) >= min_split_reads:
                        fusion_count = fusion_count + 1
                        this_variant['variant'] =  ''.join((this_fusion['fusion_product'], ' (', this_fusion['gene_product_type'],' ',this_fusion['event_type'],  ')'))
                        reportable_variants.append(this_variant)
                        fusion_not_found = False
            if fusion_not_found and subset_order[this_gene] in ['action','driver']:
                if float(gene_expression[this_gene]['cohort_perc']) < 5:
                    this_variant['variant'] =  'downregulated'
                    reportable_variants.append(this_variant)
                elif float(gene_expression[this_gene]['cohort_perc']) > 95:
                    this_variant['variant'] =  'upregulated'
                    reportable_variants.append(this_variant)
    return(reportable_variants, fusion_count)


def parse_fusions(self, fusions_file_path):
    sample_fusions = []
    fusions_file_path = check_path_exists(self, fusions_file_path)
    with open(fusions_file_path, 'r') as fusions_file:
        for row in csv.DictReader(fusions_file, delimiter="\t"):

            if row['break1_homologous_seq'] != 'None' and \
                row['transcript1'] != 'None'  and \
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
                sample_fusions.append(fusion_characteristics)
    return(sample_fusions)


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
                    'clinvar': row['clinvar'],
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
        try:
            lims_json_dict = json.loads(r.text)
        except json.decoder.JSONDecodeError:
            msg = 'LIMS file is currently writing, waiting for 45 seconds'
            self.logger.info(msg)
            time.sleep(45)

        external_ids = find_external_id_in_json_dict(lims_json_dict , donor)
    return(external_ids)

def parse_mane_transcript(self, mane_transcript_path):
    genes_and_transcripts = {}
    mane_transcript_path = check_path_exists(self, mane_transcript_path)
    with open(mane_transcript_path, 'rt') as mane_transcript_file:
        for row in csv.DictReader(mane_transcript_file, delimiter="\t"):
            if row['MANE_status'] == 'MANE Select':
                refseq_nuc = row['RefSeq_nuc'].split('.')
                genes_and_transcripts[row['symbol']] = refseq_nuc[0]
    return(genes_and_transcripts)

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

def parse_sex(self, sex_file_path, template="PCX"):
    sex_file_path = check_path_exists(self, sex_file_path)
    with open(sex_file_path) as sex_file:
        for row in sex_file:
            short_sex = row.strip()
            if template == "PCX":
                if short_sex == 'M':
                    inferredSex = "Male"
                elif short_sex == 'F':
                    inferredSex = "Female"
                else:
                    msg = "Cannot infer sex {0} from file: {1}".format(short_sex, sex_file_path)
                    self.logger.error(msg)
                    raise RuntimeError(msg)
            else:
                if short_sex == 'M':
                    inferredSex = "XY"
                elif short_sex == 'F':
                    inferredSex = "XX"
                else:
                    msg = "Cannot infer sex {0} from file: {1}".format(short_sex, sex_file_path)
                    self.logger.error(msg)
                    raise RuntimeError(msg)

    return(inferredSex)



def parse_somatic_variants(self, sample_variants_file_path, gene_list = {}, extra_sites = []):
    sample_variants = {}
    sample_variants_file_path = check_path_exists(self, sample_variants_file_path)
    with open(sample_variants_file_path, 'r') as sample_variants_file:
        for row in csv.DictReader(sample_variants_file, delimiter=","):
            gene = row['gene']
            if 'somatic' in row['mutation_class']:
               
                variant_key = f"{row['mutation_type']},{row['position']},{row['base_change']}"
                
                variant_characteristics = {
                    'gene': gene,
                    'gene_chr': row['gene_chr'],
                    'mutation_type': row['mutation_type'],
                    'mutation_class': row['mutation_class'],
                    'rarity': row['rarity'],
                    'clinvar': row['clinvar'],
                    'dbsnp': row['dbsnp'],
                    'copy_number': row['copy_number'],
                    'ab_counts': row['ab_counts'],
                    'cosmic_census_flag': row['cosmic_census_flag'],
                    'nuc_context': split_sequence_context_by_transcript(row['nuc_context']),
                    'aa_context': split_sequence_context_by_transcript(row['aa_context']),
                    'position': row['position'],
                    'tumour_freq': row['tumour_freq']
                }
                # add variant to gene or create gene
                if gene in sample_variants:
                    sample_variants[gene][variant_key] = variant_characteristics
                else:
                    sample_variants[gene] = {}
                    sample_variants[gene][variant_key] = variant_characteristics
            for this_site in extra_sites:
                if row['gene'] == this_site['element']:

                    variant_key = f"{this_site['chrom']},{this_site['pos']}"
                    
                    variant_characteristics = {
                        'gene': this_site['element'],
                        'gene_chr': row['gene_chr'],
                        'position': this_site['pos'],
                        'nuc_context': this_site['variant'],
                        'mutation_type': this_site['type'],
                        'mutation_class': 'EXTRASITE',
                        'rarity': 'NA',
                        'clinvar': 'NA',
                        'dbsnp': 'NA',
                        'copy_number': row['copy_number'],
                        'ab_counts': row['ab_counts'],
                        'cosmic_census_flag': 'NA',
                        'aa_context': 'NA',
                        'tumour_freq': this_site['vaf']
                    }

                    # add variant to gene or create gene
                    if gene in sample_variants:
                        sample_variants[gene][variant_key] = variant_characteristics
                    else:
                        sample_variants[gene] = {}
                        sample_variants[gene][variant_key] = variant_characteristics


    #run through again to get CNs for genes without variants
    with open(sample_variants_file_path, 'r') as sample_variants_file:
        for row in csv.DictReader(sample_variants_file, delimiter=","):
            gene = row['gene']
            if (row['gene'] in gene_list.keys() ) \
                and ('NA' in row['mutation_class'] or 'germline' in row['mutation_class']) \
                and gene not in sample_variants \
                and 'somatic' not in row['mutation_class']:
                    gene = row['gene']
                    variant_key = "No Variant"

                    variant_characteristics = {
                        'gene': gene,
                        'gene_chr': row['gene_chr'],
                        'copy_number': row['copy_number'],
                        'ab_counts': row['ab_counts'],
                        'mutation_type': 'NA',
                        'mutation_class': 'NA',
                        'rarity': 'NA',
                        'clinvar': 'NA',
                        'dbsnp': 'NA',
                        'cosmic_census_flag': 'NA',
                        'nuc_context': 'NA',
                        'aa_context': 'NA',
                        'position': 'NA',
                        'tumour_freq':  'NA'
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

def split_sequence_context_by_transcript(context_string):
    full_nuc_context = {}
    if context_string not in ["NA", ""]:
        split_nuc_context = context_string.split("|")
        for this_nuc_context in split_nuc_context:
            this_nuc_context_split = this_nuc_context.split(":")
            full_nuc_context[this_nuc_context_split[0]] = this_nuc_context_split[1]
        return(full_nuc_context)
    else:
        return("NA")


def subset_and_deduplicate(data ):
    subset_data = {}
    for this_gene in data:
        for this_variant in data[this_gene]:
            if this_gene not in subset_data: 
                subset_data[this_gene] = {
                    'copy_number' : data[this_gene][this_variant]['copy_number'],
                    'ab_counts' : data[this_gene][this_variant]['ab_counts'],
                    'gene_chr' : data[this_gene][this_variant]['gene_chr'],
                }
    return(subset_data)

def try_two_null_files(self, wrapper, workflow_name, ini_param, path_info, first_file):
    if wrapper.my_param_is_null(ini_param):
        
        if self.workspace.has_file(first_file):
            msg = '{0} found in workspace {1} '.format(first_file, self.workspace.print_location() )
            self.logger.info(msg)
            this_file = os.path.join(self.workspace.print_location(), first_file )
            wrapper.set_my_param(ini_param, this_file)
        elif self.workspace.has_file(path_info):
            msg = '{0} not found in workspace {1}, checking {2} '.format(first_file, self.workspace.print_location(), path_info )
            self.logger.info(msg)


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
