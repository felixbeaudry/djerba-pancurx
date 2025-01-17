"""Methods to generate SNV/indel html"""

import djerba.core.constants as core_constants
import djerba.plugins.pancurx.constants as phe
from djerba.util.html import html_builder as hb

def k_comma_format(value):
    value = int(value)
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
             '<td {0}>{1}</td>'.format('style="width:25%"', ''.join(("<strong>",row[phe.REPORTING_NAME],"<strong>"))),
             '<td {0}>{1}</td>'.format('style="width:20%"', str(row[phe.VALUE]).rstrip(".0"))
        ]
        cells = ''.join(cells)
        if row_count % 2 == 0:
            row_join = ''.join((row_hold, cells))
            rows.append(hb.tr(row_join))
        else: 
            row_hold = cells
        row_count += 1
    return rows


def make_fusion_table_rows(row_fields, ploidy):
    rows = []
    already_printerd = []
    for row in row_fields:
        cn_class, cn = process_cn(row['copy_number'], ploidy)
        cells = [
            hb.td(row['gene'], italic=True),
            hb.td(row['gene_chr']),
            '<td {0}>{1}</td>'.format(cn_class, cn),
            hb.td(row['input_tpm']),
            hb.td(row['cohort_perc']),
            hb.td(row['variant'])

        ]
        
        if ' '.join((row['gene'], row['variant'])) not in already_printerd:
            rows.append(hb.tr(cells))
            already_printerd.append(' '.join((row['gene'], row['variant'])))
    return rows


def make_fusion_slide_rows(row_fields, ploidy):
    rows = []
    already_printed = []
    for row in row_fields:
        cn_class, cn = process_cn(row['copy_number'], ploidy)
        cells = [
            hb.td(row['gene'], italic=True),
            '<td {0}>{1}</td>'.format(cn_class, cn),
            hb.td(row['cohort_perc']),
            hb.td(row['variant'])

        ]
        if ' '.join((row['gene'], row['variant'])) not in already_printed and \
            (cn_class != '' or row['variant'] not in ['downregulated', 'upregulated']):
            rows.append(hb.tr(cells))
            already_printed.append(' '.join((row['gene'], row['variant'])))
    return rows

def build_dbsnp_url(dbsnp, chr, start):
    if dbsnp == 'Novel':
        chr = chr.split("chr")[1]
        spdi = ''.join(( chr,':g.', start))
        url = 'https:/www.ncbi.nlm.nih.gov/snp/?term={0}%3A{1}'.format(chr, start)
        return_string = '<a href="{0}">{1}</a>'.format(url, spdi)
    else:
        url='https:/www.ncbi.nlm.nih.gov/snp/'+dbsnp
        return_string = '<a href="{0}">{1}</a>'.format(url, dbsnp)
    return return_string

def build_cosmic_url(sbs):

    url='https://cancer.sanger.ac.uk/signatures/sbs/'+sbs
    return_string = '<a href="{0}">{1}</a>'.format(url, sbs)
    return return_string


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
            hb.td(build_dbsnp_url(row['dbsnp'], row['chr'], row['start']))
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

def make_sensitivity_rows(row_fields, this_class):
    rows = []
    for row in row_fields:
        hit_class = ''
        if row['drug_class'] == this_class:


            ic50_cmax = None
            ic50_is_smaller = None
            if row['ic50'] != 'NA' and row['cmax'] != 'NA':
                if float(row['ic50']) < float(row['cmax']):
                    ic50_cmax = ''.join((str(round(float(row['ic50']), 2)), ' (< ', str(round(float(row['cmax']), 2)),')'))
                    ic50_is_smaller = True
                    ic50_class='class="cn_gain"'
                else:
                    ic50_cmax = ''.join((str(round(float(row['ic50']), 2)), ' (> ', str(round(float(row['cmax']), 2)),')'))
                    ic50_is_smaller = False
                    ic50_class=''
            else:
                ic50_class=''
            if float(row['auc_ptile']) < 30:
                auc_percentile_is_below = True
                auc_class= 'class="cn_gain"'
            else:
                auc_percentile_is_below = False
                auc_class=''
            if auc_percentile_is_below and ic50_is_smaller:
                sensitivity_class = 'class="cn_gain"'
                pdo_is_sensitive = 'True'
            else:
                sensitivity_class = ''
                pdo_is_sensitive = 'False'
            if row['ic50'] != 'NA':
                ic50 = round(float(row['ic50']), 2)
            else:
                ic50 = 'None'
            cells = [
                hb.td(row['drug']),
                hb.td(row['mechanism']),
                '<td {0} style="text-align:center">{1}</td>'.format(auc_class, row['auc_ptile']),
                '<td  style="text-align:center">{0}</td>'.format(ic50 ),            
                ]
            if this_class == 'adopt':
              cells = [
                hb.td(row['drug']),
                hb.td(row['mechanism']),
                '<td {0} style="text-align:center">{1}</td>'.format(auc_class, row['auc_ptile']),
                '<td {0} style="text-align:center">{1}</td>'.format(ic50_class, ic50_cmax),
                '<td {0} style="text-align:center">{1}</td>'.format(sensitivity_class, pdo_is_sensitive)
            ]       
            rows.append(hb.tr(cells))
    return rows

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
    signature_dict = dict(sorted(signature_dict.items(), key=lambda item: item[1]['proportion'], reverse=True))
    for this_signature in signature_dict:
        this_signature_value = round(signature_dict[this_signature]['proportion'] * 100)
        if this_signature_value > 10 :
            this_signature_split = this_signature.split(".")
            this_signature = this_signature_split[1]
            this_signature_string = '{0} ({1}%)'.format(this_signature, this_signature_value)
            ## uncomment below to highlight signatures in same color as sigs plot
            #this_signature_string = '<mark class="sig{0}">{0}</mark> ({1}%)'.format(this_signature, this_signature_value)
            signature_strings.append(this_signature_string)
    return(', '.join(signature_strings))
    

def make_germline_slide_rows(row_fields,  reportable_genes):
    rows = []
    for this_gene in reportable_genes:
        theses_mutations = []
        loh = ''
        clinvar = ''
        for row in row_fields:
            if row['gene'] == this_gene:
                if row['mutation_type'] == 'splice':
                    theses_mutations = ['splice']
                else:
                    mutation = process_gene_context_for_slide(row)
                    if mutation != "":
                        theses_mutations.append(mutation)
                if row['clinvar'] in ['Uncertain_significance', 'Conflicting_interpretations_of_pathogenicity', 'Conflicting_interpretations_of_pathogenicity,_risk_factor']:
                    clinvar = 'VUS'
                elif row['clinvar'] == 'NA' or row['clinvar'] == 'Benign':
                    clinvar = row['dbsnp']
                else:
                    clinvar = row['clinvar'].replace("_"," ")
                loh = process_ab_slide(row['ab_counts'])
        theses_mutations = list(set(theses_mutations))
        if loh != '' :
            loh = ''.join(('s',loh))
            theses_mutations.append(loh)
            
        theses_mutations = ", ".join(theses_mutations)


        cells = [
            hb.td(this_gene, italic=True),
            hb.td(theses_mutations),
            hb.td(clinvar),
        ]
        if theses_mutations != '':
            rows.append(hb.tr(cells))
    return rows

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
        if loh != '' and "".join(theses_mutations) != 'homozygous deletion':
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

def make_somatic_table_rows(row_fields, ploidy, tier_mode, all_genes=False):
    rows = []
    for row in row_fields:
        
        mutation = process_gene_context(row)
        additional_string = process_additional(row)
        cn_class, cn = process_cn(row['copy_number'], ploidy)
        ab_class, ab = process_ab(row['ab_counts'])
        if row['tumour_freq'] == 'NA':
            tumour_freq = ''
        else:
            tumour_freq = round(float(row['tumour_freq']), 2)
        if all_genes==False and row['tier'] == 'discovery' :
            gene_name = ''.join((row['gene'],'<sup>&#9872;</sup>'))
        else:
            gene_name = row['gene']
        cells = [
            hb.td(gene_name, italic=True),
            hb.td(row['gene_chr']),
            '<td {0}>{1}</td>'.format(cn_class, cn),
            '<td {0}>{1}</td>'.format(ab_class, ab),
            hb.td(tumour_freq),
            hb.td(mutation),
           
        ]
        if all_genes:
            cells.append( hb.td(additional_string))
        if tier_mode == 'discovery':
            rows.append(hb.tr(cells))
        elif row['tier'] in ['driver', 'action']:
            rows.append(hb.tr(cells))
    return rows

def make_signature_table_rows(row_fields):
    rows = []
    for row in row_fields:
        if row_fields[row]['proportion'] > 0.3:
            sbs_class = 'class="cn_very_gain"'
        elif row_fields[row]['proportion'] > 0.1:
            sbs_class = 'class="cn_gain"'
        else:
            sbs_class = ''
        cells = [
            hb.td(build_cosmic_url(row_fields[row]['sbs'])),
            hb.td(row_fields[row]['aetiology']),
             '<td {0}>{1}</td>'.format(sbs_class, row_fields[row]['proportion']),
         ]
        rows.append(hb.tr(cells))
    return rows

def make_immune_cell_rows(row_fields):
    rows = []

    for row in row_fields:
        sbs_class = ''
        ab_class = ''
        if row['percentile'] > 95:
            sbs_class = 'class="cn_gain"' 
        if row['cibersort'] == 0:
              ab_class = 'class="cn_loh"'
        cells = [
            hb.td(row['cell_type']),

             '<td {0}>{1}</td>'.format(ab_class, row['cibersort']),
             '<td {0}>{1}</td>'.format(sbs_class, row['percentile']),

            hb.td(row['description']),
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
    ab = "|".join(a_b_list)
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
    if a_b_list[0] == "0" :
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

def process_cn(cn_full, ploidy=None):
    cn_split = cn_full.split('|')
    cn = 'NA'
    for this_cn in cn_split:
        if this_cn == 'NA':
            next
        else:
            cn = this_cn
            break
    try:
        cn = round(float(cn), 2)
        if ploidy != None:
            cn_class = get_cn_cell_class(cn, ploidy)
        else:
            cn_class = ''
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
        if 'nuc_context' not in row:
            row['nuc_context'] = 'NA'
        if 'aa_context' not in row:
            row['aa_context'] = 'NA'
        if row['mutation_class'] in ( "somatic snv", "somatic indel", "germline snp", "germline indel") and row['nuc_context'] != '':
 
            nuc_context = row['nuc_context']

            mutation_string = ' '.join((row['mutation_type'],"(",nuc_context,')'))
            
            if row['aa_context'] != 'NA':
                aa_context = row['aa_context']

                mutation_string = ' '.join((row['mutation_type'],"(",nuc_context,';',aa_context,')'))
                if len(mutation_string) > 80:
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
        if 'nuc_context' not in row:
            row['nuc_context'] = 'NA'
        if 'aa_context' not in row:
            row['aa_context'] = 'NA'

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


def add_curve_plots(row_fields):
    plots = []
    for row in row_fields:
        plot = '<tr><td ><img style="width: 80%; object-fit: contain" src="{0}"/></td></tr>'.format(row)
        plots.append(plot)
    return plots


