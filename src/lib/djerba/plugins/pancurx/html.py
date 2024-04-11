"""Methods to generate SNV/indel html"""

import djerba.core.constants as core_constants
import djerba.plugins.pancurx.constants as phe
from djerba.util.html import html_builder as hb

def make_hallmark_tally_table(hallmark_tally, hallmark_length, color):
    rows = []
    hallmark_index = 0
    this_width = 100 / hallmark_length
    while hallmark_index <= hallmark_length:
        if hallmark_index <= (hallmark_tally -1):
            this_row = '<td style="background-color: {0}; width:{1}%">|</td>'.format(color, this_width)
        else:
            this_row = '<td style="width:{0}%">|</td>'.format(this_width)
        rows.append(this_row)
        hallmark_index = hallmark_index + 1
    return(rows)


def make_class_table_rows(multifactor_marker_results, type):
    row_fields = multifactor_marker_results 
    rows = []
    row_count = 1
    row_hold = ""
    for row in row_fields:
        cells = [
            td_class(row[phe.ABOVE_CUTOFF], type),
            hb.td(''.join(("<strong>",row[phe.REPORTING_NAME],"<strong>"))),
            hb.td(row[phe.VALUE]),
        ]
        cells = ''.join(cells)
        if row_count % 2 == 0:
            row_join = ''.join((row_hold, cells))
            rows.append(hb.tr(row_join))
        else: 
            row_hold = cells
        row_count += 1
    return rows

@staticmethod
def td_class(level, type):
    # make a table cell with an OncoKB level symbol
    # permitted levels must have a format defined in style.css
    if level == True:
        shape = '&check;'
    elif level == False:
        shape = '&#x2715;'
    div = '<div class="circle {0}-{1}">{2}</div>'.format(type, level, shape)
    return hb.td(div)

def make_somatic_table_rows(row_fields, ploidy):
    rows = []
    for row in row_fields:
        mutation = process_gene_context(row)
        additional_string = process_additional(row)
        cn_class, cn = process_cn(row['copy_number'], ploidy)
        ab_class, ab = process_ab(row['ab_counts'])
        cells = [
            hb.td(row['gene'], italic=True),
            hb.td(mutation),
            '<td {0}>{1}</td>'.format(cn_class, cn),
            '<td {0}>{1}</td>'.format(ab_class, ab),
            hb.td(additional_string)
        ]
        rows.append(hb.tr(cells))
    return rows

def process_gene_context(row):

    if row['mutation_type'] == 'NA':
        mutation_string = ''
    else:
        mutation_string = row['mutation_type']
        if row['mutation_class'] in ("somatic snv", "somatic indel", "germline snp", "germline indel") and row['nuc_context'] != '':
        
            nuc_context = row['nuc_context'].split("|")
            nuc_context = nuc_context[0].split(":")
            mutation_string = ''.join((row['mutation_type']," (",nuc_context[1],')'))
            
            if row['aa_context'] != 'NA':
                aa_context = row['aa_context'].split("|")
                aa_context = aa_context[0].split(":")
                mutation_string = ''.join((row['mutation_type']," (",nuc_context[1],';',aa_context[1],')'))
                if len(mutation_string) > 40:
                    mutation_string = ''.join((row['mutation_type']," (",nuc_context[1],')'))

        elif row['mutation_class'] in ("somatic sv", "somatic cnv"):
            if row['position'] != 'NA':
                mutation_string = ''.join((row['mutation_type']," (",row['position'],')'))


    return(mutation_string)

def process_additional(row):
    additional_string = []
    if row['dbsnp'] != 'NA':
        additional_string.append(row['dbsnp']) 
    if row['rarity'] != 'NA':
        additional_string.append(row['rarity']) 
    if row['cosmic_census_flag'] == 'cosmic_mutation':
        additional_string.append("COSMIC census hit") 
    return('; '.join(additional_string))

def process_ab(ab_full):
    ab_split = ab_full.split('|')
    ab = 'NA'
    for this_ab in ab_split:
        if this_ab == 'NA':
            next
        else:
            ab = this_ab
            break
    a_b_list = ab.split('.')
    ab_class = ''
    if a_b_list[0] == "0":
        ab_class = 'class="cn_loh"'
    if len(ab_split) > 1:
        ab = "".join((ab,'*'))
    return(ab_class, ab)

def process_cn(cn_full, ploidy):
    cn_split = cn_full.split('|')
    cn = 'NA'
    for this_cn in cn_split:
        if this_cn == 'NA':
            next
        else:
            cn = this_cn
            break
    try:
        float(cn)
        cn_class = get_cn_cell_class(cn, ploidy)
    except ValueError:
        cn_class = ''
    if len(cn_split) > 1:
        cn = "".join((cn,'*'))
    return(cn_class, cn)

def get_cn_cell_class(cn, ploidy):
    ploidy = float(ploidy)
    cn_class = ''
    if float(cn) < phe.DEFAULT_DELETION_CUTOFF:
        cn_class = 'class="cn_loss"'
    elif float(cn) < phe.DEFAULT_DEL_LOH_CUTOFF:
        cn_class = 'class="cn_loh"'
    elif float(cn) >= (ploidy * phe.DEFAULT_AMP_MULTIPLIER_CUTOFF):
        cn_class = 'class="cn_very_gain"'
    elif float(cn) > (ploidy + phe.DEFAULT_GAIN_ADDEND_CUTOFF):
        cn_class = 'class="cn_gain"'
    return(cn_class)
