"""phoenix supporting functions"""
import csv
from decimal import Decimal
import logging
import re
import os

import djerba.helpers.phoenix_helper.constants as phe
from djerba.util.logger import logger

class get_germline_variants(logger):

	NORMAL = "normal"
	TUMOUR = "tumour"
	GERMLINE_INDEL = "germline indel"
	GERMLINE_SNP = "germline snp"
	TRANSITIONS = ["AG", "GA", "CT", "TC",]
	RARE_MUTATION_CUTOFF = 0.01
	RARE = "rare"
	COMMON = "common"
	NOVEL = "novel"
	VARIANTS = 'variants'
	ANNOVAR_HEADER = ["Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "GeneDetail.refGene", 
				"ExonicFunc.refGene", "AAChange.refGene", "avsnp150", "1000g2015aug_all", "1000g2015aug_afr", 
				"1000g2015aug_eas", "1000g2015aug_eur", "1000g2015aug_sas", "ExAC_ALL", "ExAC_AFR", "ExAC_AMR", 
				"ExAC_EAS", "ExAC_FIN", "ExAC_NFE", "ExAC_OTH", "ExAC_SAS", "AF", "AF_popmax", "AF_male", "AF_female", 
				"AF_raw", "AF_afr", "AF_sas", "AF_amr", "AF_eas", "AF_nfe", "AF_fin", "AF_asj", "AF_oth", 
				"non_topmed_AF_popmax", "non_neuro_AF_popmax", "non_cancer_AF_popmax", "controls_AF_popmax", 
				"CLNALLELEID", "CLNDN", "CLNDISDB", "CLNREVSTAT", "CLNSIG", "Otherinfo"]
	WANTED_ANNOVAR_ANNOTATIONS = ['1000g2015aug_all', 'ExAC_ALL', 'AF']
	PATHOGENIC_VARIANTS = ["rs80357389"]

	def __init__(self, germline_file_path, germline_annovar_file_path, normal_id, tumour_id , gene_dict,  log_level=logging.WARNING, log_path=None):
		self.log_level = log_level
		self.log_path = log_path
		self.logger = self.get_logger(log_level, __name__, log_path)
		self.normal_id = normal_id
		self.tumour_id = tumour_id
		# self.synonym_dict = synonym_dict
		self.gene_dict = gene_dict
		self.counts_dict = {
			"snvCount" : 0,
			"insCount" : 0,
			"delCount" : 0,
			"commonCount" : 0,
			"rareCount" : 0,
			"novelCount" : 0,
			"titvRatio" : "NA",
			"tiCount" : 0,
			"missenseCount" : 0,
			"nonsenseCount" : 0,
			"noncodingCount" : 0,
			"spliceCount" : 0,
		}
		self.first_gt = ''
		germline_file_path = self.check_path_exists(germline_file_path)
		self.parse_germline_variants(germline_file_path)
		self.germline_annovar_file_path = germline_annovar_file_path

	def add_annovar_consequences_to_counts(self, consequence):
		if consequence == "splicing":
			self.counts_dict["splice_count"] += 1
		elif re.match(r'^(nonframeshift|nonsynonymous|stoploss)', consequence):
			self.counts_dict["missense_count"] += 1
		elif re.match(r'^(stopgain|frameshift)', consequence):
			self.counts_dict["nonsense_count"]  += 1

	def add_semicolon_if_none(field):
		if not field.endswith(";"):
			field += ";"
		return(field)
	
	def add_to_population_counts(self, population_frequency, wanted_annovar_annotations_found):
		rarity = ''
		if population_frequency > self.RARE_MUTATION_CUTOFF:
			self.counts_dict["commonCount"] += 1
			rarity = self.COMMON
		elif population_frequency > 0:
			self.counts_dict["rareCount"] += 1
			rarity = self.RARE
		elif wanted_annovar_annotations_found:
			self.counts_dict["novelCount"] += 1
			rarity = self.NOVEL
		return(rarity)

	def check_for_clinical_variants(self, annovar_annotation):
		clinvar_headers = sorted(key for key in annovar_annotation.keys() if key.startswith("CLNSIG"))
		clinical_variants = ""
		for clinvar_header in clinvar_headers:
			clin_info = annovar_annotation[clinvar_header]
			clinical_variants += f"{clinvar_header}={clin_info};"
		return(clinical_variants)

	def check_for_cosmic_consensus(self, gene, consequence):
		COSMIC_CENSUS_TYPES = "cosmic_census_types"
		cosmic_flag = "NA"
		if gene in self.gene_dict and COSMIC_CENSUS_TYPES in self.gene_dict[gene] and consequence in self.gene_dict[gene][COSMIC_CENSUS_TYPES]:
			cosmic_flag = "cosmic_mutation"
		return(cosmic_flag)

	def check_for_cosmic_in_info(info):
		cosmic_match = re.search(r'COSMIC=(.*?);', info)
		if cosmic_match:
			cosmic_found = cosmic_match.group(1)
		return(cosmic_found)

	def check_for_enhancer_or_promoter(self, info):
		mutation_type = ""
		if re.search(r'GENE=(.*?)\(Promoter\)', info):
			gene = re.search(r'GENE=(.*?)\(Promoter\)', info).group(1)
			self.counts_dict["noncoding_count"] += 1
			mutation_type = "altered promoter"

		if re.search(r'NCENC=(.*)', info):
			for g in re.findall(r'(.*?)\((.*)\|.*\)', re.search(r'NCENC=(.*)', info).group(1)):
				if g[0] == "Enhancer":
					gene = g[1]
					self.counts_dict["noncoding_count"] += 1
					mutation_type = "altered enhancer"
		return(gene, mutation_type)
	
	def check_path_exists(self, path):
		if not os.path.exists(path):
			msg = "Cannot find hugo file: {0}".format(path)
			self.logger.error(msg)
			raise RuntimeError(msg)
		else:
			return(path)

	def check_which_sample_is_first(self, row):
		if self.normal_id in row and self.tumour_id in row:
			first_gt = self.NORMAL
		else:
			first_gt = self.TUMOUR
		return(first_gt)

	def do_annovar_lookup(self, chromosome, pos, ref_base, alt_base):
		DEL_SYMBOL = '-' #(?)
		tabix = self.check_path_exists(self.germline_annovar_file_path)
		annovar_annotation = []
		#unclear where the numbers 5 and 10 come from
		upstream_interval = pos - 5
		downstream_interval = pos + 10
		tabix_iterval = tabix.fetch(chromosome, upstream_interval, downstream_interval)
		if DEL_SYMBOL in tabix_iterval:
			line = tabix_iterval.next()
			split_row = line.strip().split('\t')
			if len(alt_base) > len(ref_base) and ref_base != DEL_SYMBOL :
				ref_base = DEL_SYMBOL
				alt_base = alt_base[1:]
			if len(ref_base) > len(alt_base) and alt_base != DEL_SYMBOL:
				alt_base = DEL_SYMBOL
				ref_base = ref_base[1:]
				pos = pos + 1
			if len(ref_base) > len(alt_base) and ref_base[0] == alt_base[0]:
				alt_base = DEL_SYMBOL
				ref_base = ref_base[1:]
				pos = pos + 1
			if split_row[0] == chr and int(split_row[1]) == pos and \
				(split_row[3] == ref_base or split_row[52] == alt_base) and \
					(split_row[4] == alt_base or split_row[53] == alt_base):
				for header_index in range(5, len(self.ANNOVAR_HEADER)):
					if split_row[header_index] == ".":
						split_row[header_index] = "NA"
					split_row[header_index] = split_row[header_index].replace(",", ".")
					annovar_annotation[self.ANNOVAR_HEADER[header_index]] = split_row[header_index]
		else:
			print(f"did not find variant {chromosome}, {pos}, {ref_base}, {alt_base} \
		 			when searching from {upstream_interval} to {downstream_interval}")
			raise ValueError
		return(annovar_annotation)

	def find_population_frequency(self, annotation):
		population_frequency = -1
		wanted_annovar_annotations_found = False
		for annotation_column in self.WANTED_ANNOVAR_ANNOTATIONS:
			if annotation_column in annotation:
				wanted_annovar_annotations_found = True
				if annotation[annotation_column] != "NA":
					if float(annotation[annotation_column]) > population_frequency:
						population_frequency = float(annotation[annotation_column])				
		return(population_frequency, wanted_annovar_annotations_found)

	def get_depth_and_frequency(self, genotype):
		match = re.match(r'^.*:.*,(.*):(.*):.*:.*', genotype)
		if match:
			depth = int(match.group(2))
			freq = "."
			if depth > 0:
				freq = float(match.group(1)) / depth
		return(depth, freq)

	def get_mutation_type(self, ref, alt):
		if len(ref) == len(alt):
			mut_class = self.GERMLINE_SNP
			self.counts_dict["snv_count"] += 1
			if f"{ref}{alt}" in self.TRANSITIONS:
				ti_count += 1
		elif len(ref) > len(alt):
			mut_class = self.GERMLINE_INDEL
			self.counts_dict["del_count"] += 1
		else:
			mut_class = self.GERMLINE_INDEL
			self.counts_dict["ins_count"] += 1
		return(mut_class)

	def get_variant_results(self, annovar_annotation, id, chr, pos, ref, alt, info, normal_genotype, tumour_genotype, best_gene_name, mutation_type, annovar_dict):
		normal_depth, normal_freq = self.get_depth_and_frequency(normal_genotype)
		tumour_depth, tumour_freq = self.get_depth_and_frequency(tumour_genotype)
		clinical_variants = self.check_for_clinical_variants(annovar_annotation)

		population_frequency, wanted_annovar_annotations_found = self.find_population_frequency(annovar_annotation)
		base_change = f"{ref}>{alt}"
		chr_pos = f"{chr}:{pos}"
		var_name = f"SGV {mutation_type} {chr_pos} {base_change}"
		results = { var_name : { "dbsnp": id,
					"position": chr_pos,
					"base_change": base_change,
					"mutation_type": mutation_type,
					"mutation_class": self.get_mutation_type(ref, alt),
					"nuc_context": annovar_dict["nuc"],
					"aa_context": annovar_dict["aa"],
					"tumour_freq": "{:.3f}".format(tumour_freq),
					"tumour_depth": tumour_depth,
					"normal_freq": "{:.3f}".format(normal_freq),
					"normal_depth":  normal_depth,
					"cosmic": self.check_for_cosmic_in_info(info),
					"rarity": self.add_to_population_counts(population_frequency, wanted_annovar_annotations_found),
					"1000G_all": annovar_annotation["1000g2015aug_all"],
					"ExAC_all": annovar_annotation["ExAC_ALL"],
					"gnomad_all": annovar_annotation["AF"],
					"cadd_phred": annovar_annotation["CADD_phred"],
					clinical_variants: annovar_annotation["CLINSIG"],
					"cosmic_census_flag": self.check_for_cosmic_consensus(best_gene_name, mutation_type)
				} }
		return(results)
	
	def parse_germline_variants(self, germline_file_path, germline_annovar_file_path):
		with open(germline_file_path, 'r') as germline_file:
			for row in csv.DictReader(germline_file, delimiter="\t"):
				if row.startswith("#CHROM"):
					first_gt = self.check_which_sample_is_first(row)
				elif row.startswith("#"):
					pass
				else:
					if first_gt == self.NORMAL:
						chr, pos, id, ref, alts, qual, filter, info, format, normal_genotype, tumour_genotype = row
					else:
						chr, pos, id, ref, alts, qual, filter, info, format, tumour_genotype, normal_genotype = row
					
					info = self.add_semicolon_if_none(info) ## needed??
					annovar_dict = self.parse_annovar(info)

					for alt in alts.split(','):
						for t in sorted(annovar_dict.keys()):
							for gene in sorted(annovar_dict[t].keys()):
								consequence = annovar_dict[t][gene]["consequence"] #mut_type = consequence
								self.add_annovar_consequences_to_counts(consequence)
								if consequence == "" | alt == '*':
									pass
								else:
									#TODO: make synonym dict functional
									best_gene_name = gene # self.find_best_gene_symbol(gene, chr, self.gene_dict, self.synonym_dict, "ANNOVAR (SGV)", warnings=True)
									annovar_annotation = self.do_annovar_lookup( chr, pos, ref, alt, germline_annovar_file_path, row)
									if id == self.PATHOGENIC_VARIANTS:
										annovar_annotation["CLINSIG"] = "CLINSIG=pathogenic"       
									self.gene_dict[gene][self.VARIANTS] = self.get_variant_results(annovar_annotation, id, chr, pos, ref, alt, info, 
																	  normal_genotype, tumour_genotype, best_gene_name, consequence , annovar_dict[t][gene])

						# additional variants if also in another gene's enhancer (?)
						gene, mutation_type = self.check_for_enhancer_or_promoter(info)
						if mutation_type != "" and gene != "chmm/segway" and gene != "drm" and gene in self.gene_dict:
							annovar_annotation = self.do_annovar_lookup( chr, pos, ref, alt, germline_annovar_file_path, row)
							self.gene_dict[gene][self.VARIANTS] = self.get_variant_results(annovar_annotation, id, chr, pos, ref, alt, info, 
																	  normal_genotype, tumour_genotype, gene, mutation_type , annovar_dict[t][gene])

		self.counts_dict["germline_indel_count"] = self.counts_dict["delCount"] + self.counts_dict["insCount"]
		if self.counts_dict["snvCount"] - self.counts_dict["tiCount"] > 0:
			self.counts_dict["germline_titv_ratio"] = self.counts_dict["tiCount"] / (self.counts_dict["snvCount"] - self.counts_dict["tiCount"] )

def parse_annovar(info, base_change=None):
	
    info_hash = {}
    gene = "."
    consequence = "."
    g, nuc_context, aa_context =  None, None, None
    types, genes = [], []
    anno_hash = {}

    if re.search(r'intergenic', info, re.IGNORECASE) or re.search(r'ANNOVAR=(upstream|downstream)', info):
        return info_hash
    elif match := re.search(r'ANNOVAR=([a-zA-Z0-9]+\|?[a-z]*\|?)\|([a-zA-Z0-9-|()\+\->_:\.]*)', info):
        types = match.group(1).split('|')
        annotation_part2 = match.group(2)

        if all(this_type in ['UTR', 'intergenic', 'intronic'] for this_type in types):
            return info_hash
        elif len(types) == 1:
            gene = annotation_part2
            variants = annotation_part2.split('|')
            if all(this_variant in ['UTR', 'intergenic', 'intronic'] for this_variant in variants):
                return info_hash
            elif types[0] == 'exonic':
                genes = gene.split('|')
            elif types[0] == 'splicing':
                genes = [g + ')' if '(' in g else g for g in gene.split(')')]
            else:
                genes = [gene]
        elif len(types) >= 2:
            if match2 := re.search(r'(([a-zA-Z0-9\-]+\|)+)([a-zA-Z0-9>\+\-_\.:|()]+)$', annotation_part2):
                genes = match2.group(1).split('|') + [g + ')' if '(' in g else g for g in match2.group(3).split(')')]
            else:
                print(f"parse error: {types}, {annotation_part2}, {info}")
                raise ValueError("Parse Error")

        for i in range(len(types)):
            t = types[i]
            if t == "exonic":
                for g in [g for g in genes[i].split('|') if g not in ['splicing'] and '(' not in g]:
                    anno_hash.setdefault(t, {}).setdefault(g, 0)
                    anno_hash[t][g] += 1
            elif t == "splicing":
                if genes[i] is None:
                    print(f"{types}, {genes}, {info}")
                    raise ValueError("Error")
                if match3 := re.match(r'^(.*?)\((.*?)\)', genes[i]):
                    gene = match3.group(1)
                    context = match3.group(2)
                    anno_hash[t][gene] = context
                    genes[i] = re.sub(r'^.*?\(.*?\)', '', genes[i]).lstrip(',')
                elif base_change is not None:
                    anno_hash[t][genes[i]] = base_change

        for t in sorted(anno_hash.keys()):
            if t == "splicing":
                for g in sorted(anno_hash[t].keys()):
                    nuc_context = "|".join(f"{m.group(1)}:{m.group(2)}" for m in re.finditer(r'(.*?):.*?:(.*?)$', anno_hash[t][g]))
                    nuc_context = nuc_context.lstrip('|')
                    aa_context = "NA"

                    info_hash.setdefault(t, {}).setdefault(g, {})
                    info_hash[t][g]['consequence'] = "splicing"
                    info_hash[t][g]['nuc'] = nuc_context
                    info_hash[t][g]['aa'] = aa_context
            elif t == "exonic":
                for g in sorted(anno_hash[t].keys()):
                    nuc_context, aa_context, consequence = "", "", "."

                    if match4 := re.search(r'ANNOVAR_EXONIC=([a-z]+-?([a-zA-Z]+)?)\|(.*:NM_[0-9]+:.*);ANNOVAR', info):
                        consequence = f"{match4.group(1)},{match4.group(3)}"
                        consequence = re.sub(r'\|$', '', consequence)

                        for isoform in consequence.split('|'):
                            if m := re.match(r'^([a-z]+-?[a-zA-Z]+,)?' + re.escape(g) + r':(NM_[0-9]+):.*:(.*?):(.*?)\|?$', isoform):
                                nuc_context += f"|{m.group(2)}:{m.group(3)}"
                                aa_context += f"|{m.group(2)}:{m.group(4)}"
                            elif m := re.match(r'^([a-z]+-?[a-zA-Z]+,)?' + re.escape(g) + r':(NM_[0-9]+):wholegene$', isoform):
                                nuc_context += f"|{m.group(2)}:wholegene"
                                aa_context += f"|{m.group(2)}:wholegene"
                            else:
                                other_genes = [gene for gene in genes if gene != g]
                                if len(other_genes) > 0:
                                    for gene in other_genes:
                                        if m := re.match(r'^([a-z]+-?[a-zA-Z]+,)?' + re.escape(gene) + r':(NM_[0-9]+):.*:(.*?):(.*?)(\|)?$', isoform):
                                            nuc_context += f"|{m.group(2)}:{m.group(3)}"
                                            aa_context += f"|{m.group(2)}:{m.group(4)}"
                                        elif m := re.match(r'^([a-z]+-?[a-zA-Z]+,)?' + re.escape(gene) + r':(NM_[0-9]+):wholegene$', isoform):
                                            nuc_context += f"|{m.group(2)}:wholegene"
                                            aa_context += f"|{m.group(2)}:wholegene"

                                    if not any(char.isalnum() for char in nuc_context):
                                        print(f"isoform mismatch1: {isoform}\t{consequence}\t{info}")
                                        print(f"{g}: {other_genes}, types: {types}")
                                        for gene in other_genes:
                                            if m := re.match(r'^([a-z]+-?[a-zA-Z]+,)?' + re.escape(gene) + r':(NM_[0-9]+):.*:(.*?):(.*?)(\|)?$', isoform):
                                                nuc_context += f"|{m.group(2)}:{m.group(3)}"
                                                aa_context += f"|{m.group(2)}:{m.group(4)}"
                                            else:
                                                print(f"{gene} not the right one")
                                        raise ValueError("Error")
                                else:
                                    print(f"isoform mismatch2: {isoform}, {consequence}, {info}")
                                    print(f"{g}: {genes}, types: {types}")
                                    raise ValueError("Error")
                    else:
                        info_hash = {}
        return info_hash
    else:
        if re.search(r'ANNOVAR=', info) and 'ncRNA' not in info and 'intronic' not in info:
            print(f"gene mismatch: {info}")
        info_hash = {}
        return info_hash