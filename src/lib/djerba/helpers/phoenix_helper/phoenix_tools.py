"""phoenix supporting functions"""
import csv
from decimal import Decimal
import logging
import re
import os

import djerba.helpers.phoenix_helper.constants as phe
from djerba.util.logger import logger


def find_best_gene_symbol(gene, chromosome, gene_dict, synonym_dict, origin, warnings=True):
	if gene in gene_dict:
		return gene
	elif gene in synonym_dict:
		for synonym_chromosome in synonym_dict[gene]:
			if synonym_chromosome == chromosome:
				if ';' not in synonym_dict[gene][synonym_chromosome]:
					return synonym_dict[gene][synonym_chromosome]
				else:
					if warnings:
						msg = "{0}: multiple same chromosome synonyms for {1} on {2}: \
							{synonym_dict[gene][synonym_chromosome]}".format(origin, gene, chromosome)
						print(msg)
					return synonym_dict[gene][synonym_chromosome]
		else:
			if warnings:
				msg = "{0}: no synonyms for {1} on {2} using original name".format(origin, gene, chromosome)
				print(msg)
			return gene
	else:
		if warnings:
			msg = "{0}: no match or synonym entry for {1}, using original name".format(origin, gene)
			print(msg)
		return gene
		
class phoenix_processor(logger):

	def __init__(self, sample_data, log_level=logging.WARNING, log_path=None):
		self.log_level = log_level
		self.log_path = log_path
		self.logger = self.get_logger(log_level, __name__, log_path)
		self.gene_dict = {} #formerly `vars`
		self.synonym_dict = {} #formerly `synonym`
		self.sample_data = sample_data #formerly `data`
		self.parse_hugo()
		self.gene_positions = self.look_up_gene_positions()
		self.parse_cosmic_census()

	def add_X_to_Y_synonyms(self, row, cytoband):
		this_gene = row['Approved Symbol']
		for chromosome in row['Chromosome'].split(' and '):
			if this_match := re.match(r'^(.*?)[pq]', chromosome):
				cytoband = f'chr{this_match.group(1)}'
			else:
				cytoband = f'chr{chromosome}'
			if len(self.gene_dict[this_gene]['synonyms'].split(';')) > 0:
				for interchromosome_synonym in self.gene_dict[this_gene]['synonyms'].split(';'):
					if interchromosome_synonym in self.synonym_dict:
						if cytoband in self.synonym_dict[interchromosome_synonym]:
							self.synonym_dict[interchromosome_synonym][cytoband] += f';{this_gene}'
						else:
							self.synonym_dict[interchromosome_synonym][cytoband] = this_gene
		return(cytoband)

	def look_up_gene_positions(self):
		refseq_file_path = phe.DEFAULT_REFSEQ_FILE_PATH
		if not os.path.exists(refseq_file_path):
			msg = "Cannot find hugo file: {0}".format(refseq_file_path)
			self.logger.error(msg)
			raise RuntimeError(msg)
		gene_positions = {}
		with open(refseq_file_path, 'r') as refseq_file:
			for row in refseq_file:
				chromosome, start, end, gene = row.strip().split('\t')
				if gene in self.gene_dict:
					gene = self.find_best_gene_symbol(gene, chromosome,  self.gene_dict, self.synonym_dict, refseq_file_path, warnings=False)
					gene_positions[gene] = f"{chromosome}:{start}-{end}"
					# self.gene_dict[gene]['gene'] = gene # not sure why you'd want to overwrite the gene name??
					self.gene_dict[gene]['gene_position'] = f"{chromosome}:{start}-{end}"
		return(gene_positions)
			
	def merge_synonyms(self, row):
		if len(row['Previous Symbols']) > 0 & len(row['Synonyms']) > 0:
			synonyms = f"{row['Previous Symbols']}, {row['Synonyms']}"
			synonyms.replace(',', ';').strip(';')
		elif len(row['Previous Symbols']) > 0 & len(row['Synonyms']) == 0:
			synonyms = row['Previous Symbols'].strip(';')
		elif len(row['Previous Symbols']) == 0 & len(row['Synonyms']) > 0:
			synonyms = row['Synonyms'].strip(';')
		else:
			synonyms = ""
		return(synonyms)

	def parse_hugo(self):
		hugo_file_path = phe.DEFAULT_HUGO_SYN_FILE_PATH
		if not os.path.exists(hugo_file_path):
			msg = "Cannot find hugo file: {0}".format(hugo_file_path)
			self.logger.error(msg)
			raise RuntimeError(msg)
		cytoband = ""
		with open(hugo_file_path) as hugo_file:
			for row in csv.DictReader(hugo_file, delimiter="\t"):
				if row['Status'] == "Symbol Withdrawn":
					pass
				else:
					gene_name = row['Approved Symbol']
					gene_info = {
						'gene_name': gene_name,
						'full_name': row['Approved Name'].replace(",", ";"),
						'gene_chr': row['Chromosome'].replace(",", ";"),
						'entrez_id': row['Entrez Gene ID'],
						'ensembl_id': row['Ensembl Gene ID'],
						'synonyms': self.merge_synonyms(row),
						'synonym_cytoband': cytoband
					}
					if row['RefSeq IDs']:
						gene_info['refseq_id'] = row['RefSeq IDs'].replace(",", ";")
					if row['Pubmed IDs']:
						gene_info['pubmed_ids'] = row['Pubmed IDs'].replace(",", ";")
					if row['Gene Family Name']:
						gene_info['gene_fam_id'] = row['Gene Family ID']
						gene_info['gene_fam_name'] = row['Gene Family Name'].replace(",", ";")
					self.gene_dict[gene_name] = gene_info
					cytoband = self.add_X_to_Y_synonyms(row, cytoband) 
			
		for gene_symbol, gene_info in self.gene_dict.items():
			if gene_symbol in self.synonym_dict:
				self.synonym_dict[gene_symbol][gene_info['synonym_cytoband']] = gene_symbol

	def parse_cosmic_census(self):
		cosmic_census_file_path = phe.DEFAULT_COSMIC_CENSUS_FILE_PATH
		if not os.path.exists(cosmic_census_file_path):
			msg = "Cannot find hugo file: {0}".format(cosmic_census_file_path)
			self.logger.error(msg)
			raise RuntimeError(msg)

		with open(cosmic_census_file_path) as cosmic_census_file:
			for row in csv.DictReader(cosmic_census_file, delimiter="\t"):
				gene_found_in_cosmic = False
				synonyms = [row['Gene Symbol']]
				#TODO: add synonyms
				# synonyms = {row['Gene Symbol']: 1}

				# if 'Synonyms' in row:
				#	for name in row['Synonyms'].split(','):
				#		synonyms[name] = 1

				chromosome = row['Genome Location'].split(':')[0]
				chromosome = "chr" + chromosome
		
				for gene_name in synonyms:
					best_gene_name = self.find_best_gene_symbol(gene_name, chromosome, self.gene_dict, self.synonym_dict, cosmic_census_file_path, warnings=False)
					if best_gene_name in self.gene_dict:
						gene_found_in_cosmic = True
						self.gene_dict = row

				# 		# add cosmic data to gene
				# 		$vars->{$name}{cosmic_census_data} = "";
				# 		for my $h (@header)
				# 		{
				# 			if (exists $row{$h})
				# 			{
				# 				$vars->{$name}{cosmic_census_data} .= "$h:$row{$h}|";
				# 			}
				# 		}
				# 		$vars->{$name}{cosmic_census_data} =~ s/\|$//;
				# 		$vars->{$name}{cosmic_census_data} =~ s/,/;/g;

				# 		for my $t (split(/, /, $row{'Mutation Types'}))
				# 		{
				# 			if (exists $mutTypes{$t})
				# 			{
				# 				for my $type (split(/,/, $mutTypes{$t}))
				# 				{
				# 					$vars->{$name}{cosmic_census_types}{$type}++;
				# 				}
				# 			}
				# 		}

				# 		for my $g (split(/, /, $row{'Translocation Partner'}))
				# 		{
				# 			$vars->{$name}{cosmic_census_fusions}{$g}++;
				# 		}


				if gene_found_in_cosmic is False:
				# 	if exists $row{Synonyms}:
				# 		warn "Couldn't find HUGO gene for Cosmic census gene with the following names: $row{'Gene Symbol'},$row{Synonyms}\n";
				# 	else
				#	msg = "Cannot find hugo gene {0} in file: {1}".format(cosmic_census_file_path, row['Gene Symbol'])
				#	self.logger.error(msg)
				# warn "Couldn't find HUGO gene for Cosmic census gene with the following names: $row{'Gene Symbol'}\n";
					pass


