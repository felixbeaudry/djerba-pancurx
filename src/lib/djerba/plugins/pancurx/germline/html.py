"""Methods to generate SNV/indel html"""

import djerba.core.constants as core_constants
import djerba.plugins.pancurx.germline.constants as phe
from djerba.util.html import html_builder as hb



def make_table_rows(mutation_info):
    row_fields = mutation_info['reportable_germline_variants']
    rows = []
    for row in row_fields:
        nuc_context = row['nuc_context'].split(":")
        aa_context = row['aa_context'].split(":")
        mutation_string = ''.join((row['mutation_type']," (",nuc_context[1],';',aa_context[1],')'))
        if row['cosmic_census_flag'] == 'cosmic_mutation':
            cosmic_census_flag = "COSMIC census hit"
        additional_string = '; '.join((row['dbsnp'],row['rarity'],cosmic_census_flag))
        cells = [
            hb.td(row['gene'], italic=True),
            hb.td(mutation_string),
            hb.td(row['copy_number']),
            hb.td(row['ab_counts']),
            hb.td(additional_string)
        ]
        rows.append(hb.tr(cells))
    return rows

