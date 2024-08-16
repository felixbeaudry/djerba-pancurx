"""Methods to generate SNV/indel html"""

import djerba.core.constants as core_constants
import djerba.plugins.pancurx.constants as phe
from djerba.util.html import html_builder as hb

def k_comma_format(value):
    value_formatted = f'{value:,}'
    return(value_formatted)

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


def make_fusion_table_rows(row_fields):
    rows = []
    for row in row_fields:

        cells = [
            hb.td(row['fusion_product'], italic=True),
            hb.td("Recommend Manual Review"),
            hb.td(row['event_type']),
            hb.td(row['gene_product_type']),
            hb.td(row['fusion_splicing_pattern']),
        ]
        rows.append(hb.tr(cells))
    return rows

def build_dbsnp_url(dbsnp):
    url='https:/www.ncbi.nlm.nih.gov/snp/'+dbsnp
    return '<a href="{0}">{1}</a>'.format(url, dbsnp)


def make_germline_table_rows(row_fields):
    rows = []
    for row in row_fields:
        mutation = process_gene_context(row)
        #additional_string = process_additional(row)
        ab_class, ab = process_ab(row['ab_counts'])
        cells = [
            hb.td(row['gene'], italic=True),
            hb.td(mutation),
            '<td {0}>{1}</td>'.format(ab_class, ab),
            hb.td(row['clinvar'].replace("_"," ")),
            hb.td(build_dbsnp_url(row['dbsnp']))
            #hb.td(additional_string)
        ]
        rows.append(hb.tr(cells))
    return rows

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


def make_immune_table_rows(row_fields, ploidy, expression):
    rows = []
    for row in row_fields:
        if row['gene'] in ['PDCD1','CD274','TIGIT','CTLA4']:
            cn_class, cn = process_cn(row['copy_number'], ploidy)
            expression_class = get_expression_percentile_class(expression[row['gene']]['cohort_perc'])
            expression_percentile = make_ordinal(expression[row['gene']]['cohort_perc'])
            cells = [
                hb.td(row['gene'], italic=True),
                '<td {0}>{1}</td>'.format(expression_class, expression_percentile),
                '<td {0}>{1}</td>'.format(cn_class, cn),

            ]
            rows.append(hb.tr(cells))
    return rows

def make_hla_table_rows(row_fields, ploidy):
    all_cells = ''
    for row in row_fields:
        if row['gene'] in ['HLA-A','HLA-B','HLA-C']:
            ab_class, ab = process_ab(row['ab_counts'])
            all_cells = ''.join((
                all_cells,
                hb.td(row['gene'], italic=True),
               '<td {0}>{1}</td>'.format(ab_class, ab)
            ))
    return all_cells


def make_ordinal(n):
    '''
    Convert an integer into its ordinal representation::

        make_ordinal(0)   => '0th'
        make_ordinal(3)   => '3rd'
        make_ordinal(122) => '122nd'
        make_ordinal(213) => '213th'
    '''
    n = int(n)
    if 11 <= (n % 100) <= 13:
        suffix = 'th'
    else:
        suffix = ['th', 'st', 'nd', 'rd', 'th'][min(n % 10, 4)]
    return str(n) + suffix

def make_signatures_string(signature_dict):
    signature_strings = []
    signature_dict = dict(sorted(signature_dict.items(), key=lambda item: item[1], reverse=True))
    for this_signature in signature_dict:
        this_signature_value = round(signature_dict[this_signature] * 100)
        if this_signature_value > 10 :
            if  this_signature == 'SBS2026':
                this_signature = 'SBS26'
            this_signature_split = this_signature.split("SBS")
            this_signature = this_signature_split[1]
            this_signature_string = '<mark class="sig{0}">{0}</mark> ({1}%)'.format(this_signature, this_signature_value)
            signature_strings.append(this_signature_string)
    return(', '.join(signature_strings))
    

def make_somatic_slide_rows(row_fields, ploidy, reportable_genes):
    rows = []
    for this_gene in reportable_genes:
        theses_mutations = []
        loh = ''
        cn_class = ''
        cn = ''
        for row in row_fields:
            if row['gene'] == this_gene:
                mutation = process_gene_context_for_slide(row)
                if mutation != "":
                    theses_mutations.append(mutation)
                loh = process_ab_slide(row['ab_counts'])
                cn_class, cn = process_cn(row['copy_number'], ploidy)
        theses_mutations = list(set(theses_mutations))
        if loh != '':
            theses_mutations.append(loh)
        theses_mutations = ", ".join(theses_mutations)

        cells = [
            hb.td(this_gene, italic=True),
            hb.td(theses_mutations),
            '<td {0}>{1}</td>'.format(cn_class, cn)
        ]
        if cn != '':
            rows.append(hb.tr(cells))
    return rows

def make_somatic_table_rows(row_fields, ploidy):
    rows = []
    for row in row_fields:
        mutation = process_gene_context(row)
        additional_string = process_additional(row)
        cn_class, cn = process_cn(row['copy_number'], ploidy)
        ab_class, ab = process_ab(row['ab_counts'])
        cells = [
            hb.td(row['gene'], italic=True),
            hb.td(row['gene_chr'], italic=False),
            hb.td(mutation),
            '<td {0}>{1}</td>'.format(cn_class, cn),
            '<td {0}>{1}</td>'.format(ab_class, ab),
            hb.td(additional_string)
        ]
        rows.append(hb.tr(cells))
    return rows

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
    ab = ".".join(a_b_list)
    ab_class = ''
    if a_b_list[0] == "0":
        ab_class = 'class="cn_loh"'
    if len(ab_split) > 1:
        ab = "".join((ab,'*'))
    return(ab_class, ab)

def process_ab_slide(ab_full):
    ab_split = ab_full.split('|')
    ab = 'NA'
    loh = ''
    for this_ab in ab_split:
        if this_ab == 'NA':
            next
        else:
            ab = this_ab
            break
    a_b_list = ab.split('.')
    if a_b_list[0] == "0":
        loh = 'LOH'
    if len(ab_split) > 1 and loh == 'LOH':
        loh = "".join((loh,'*'))
    return(loh)

def process_additional(row):
    additional_string = []
    if row['dbsnp'] != 'NA':
        additional_string.append(row['dbsnp']) 
    if row['rarity'] != 'NA':
        additional_string.append(row['rarity']) 
    if row['cosmic_census_flag'] == 'cosmic_mutation':
        additional_string.append("COSMIC census hit") 
    return('; '.join(additional_string))

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
        cn = round(float(cn), 1)
        cn_class = get_cn_cell_class(cn, ploidy)
    except ValueError:
        cn_class = ''
    if len(cn_split) > 1:
        cn = "".join((str(cn),'*'))
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

def get_expression_percentile_class(percentile):
    expression_percentile_class = ''
    if float(percentile) < 20 and float(percentile) > 10:
        expression_percentile_class = 'class="cn_loh"'
    elif float(percentile) < 10  :
        expression_percentile_class = 'class="cn_loss"'
    elif float(percentile) > 80 and float(percentile) < 90:
        expression_percentile_class = 'class="cn_gain"'
    elif float(percentile) > 90:
        expression_percentile_class = 'class="cn_very_gain"'
    return(expression_percentile_class)

def process_gene_context(row):

    if row['mutation_type'] == 'NA':
        mutation_string = ''
    else:
        mutation_string = row['mutation_type']
        # if row['mutation_class'] in ("somatic snv", "somatic indel") and row['nuc_context'] != '':
        
        #     nuc_context = row['nuc_context'].split("|")
        #     nuc_context = nuc_context[0].split(":")
        #     mutation_string = ''.join((row['mutation_type']," (",nuc_context[1],')'))
            
        #     if row['aa_context'] != 'NA':
        #         aa_context = row['aa_context'].split("|")
        #         aa_context = aa_context[0].split(":")
        #         mutation_string = ''.join((row['mutation_type']," (",nuc_context[1],';',aa_context[1],')'))
        #         if len(mutation_string) > 40:
        #             mutation_string = ''.join((row['mutation_type']," (",nuc_context[1],')'))

        if row['mutation_class'] in ( "somatic snv", "somatic indel", "germline snp", "germline indel") and row['nuc_context'] != '':
 
            nuc_context = row['nuc_context']

            mutation_string = ' '.join((row['mutation_type'],"(",nuc_context,')'))
            
            if row['aa_context'] != 'NA':
                aa_context = row['aa_context']

                mutation_string = ' '.join((row['mutation_type'],"(",nuc_context,';',aa_context,')'))
                if len(mutation_string) > 40:
                    mutation_string = ' '.join((row['mutation_type'],"(",nuc_context,'; large charge)'))

        elif row['mutation_class'] in ("somatic sv", "somatic cnv"):
            if row['position'] != 'NA':
                mutation_string = ' '.join((row['mutation_type'],"(",row['position'],')'))

    return(mutation_string)

def process_gene_context_for_slide(row):

    if row['mutation_type'] == 'NA':
        mutation_string = ''
    else:
        mutation_string = row['mutation_type']
        # if row['mutation_class'] in ("somatic snv", "somatic indel") and row['nuc_context'] != '':
        
        #     if row['aa_context'] != 'NA':
        #         aa_context = row['aa_context'].split("|")
        #         aa_context = aa_context[0].split(":")
        #         mutation_string = aa_context[1]
        #         if len(mutation_string) > 40:
        #             mutation_string = nuc_context[1]
        #     else:
        #         nuc_context = row['nuc_context'].split("|")
        #         nuc_context = nuc_context[0].split(":")
        #         mutation_string = nuc_context[1]
        if row['mutation_class'] in ("germline snp", "germline indel", "somatic snv", "somatic indel") and row['nuc_context'] != '':
            
            nuc_context = row['nuc_context']

            if row['aa_context'] != 'NA':
                mutation_string = row['aa_context']
                if len(mutation_string) > 40:
                    mutation_string = nuc_context[1]
            else:
                mutation_string = nuc_context

        elif row['mutation_class'] in ("somatic sv", "somatic cnv"):
            if row['position'] != 'NA':
                mutation_string = ''.join((row['mutation_type']))

    return(mutation_string)


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

