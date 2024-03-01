"""Methods to generate SNV/indel html"""

import djerba.core.constants as core_constants
import djerba.plugins.pancurx.classification.constants as phe
from djerba.util.html import html_builder as hb


def make_table_rows(multifactor_marker_results, type):
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


