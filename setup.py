#! /usr/bin/env python3

"""
Setup script for Djerba
"""

from setuptools import setup, find_packages

with open('src/lib/djerba/version.py') as version_file:
    exec(version_file.read()) # sets __version__
package_root = 'src/lib'

# list of wildcards, intended to capture ancillary files for plugins/helpers/mergers
# TODO make this neater and/or introduce stronger naming conventions
install_wildcards = [
    '*.json',
    '*.html',
    '*.txt',
    '*.r',
    '*.R',
    'data/*',
    'html/*',
    'resources/*',
    'R/*',
    'Rscripts/*'
]

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='djerba',
    version=__version__,
    scripts=[
        'src/bin/djerba.py',
        'src/bin/generate_ini.py',
        'src/bin/validate_plugin_json.py'
    ],
    packages=find_packages(where=package_root),
    package_dir={'' : package_root},
    package_data={
        'djerba': [
            'data/20200818-oncoKBcancerGeneList.tsv',
            'data/blurb_template.txt',
            'data/hgnc_complete_set.txt',
            'data/cosmic_census_20161115.tsv',
            'data/HUGO_synonyms_171213.txt',
            'data/dataset.xml',
            'data/ensembl_conversion_hg38.txt',
        ],
        'djerba.core': install_wildcards,
        'djerba.helpers.file_helper': install_wildcards,
        'djerba.helpers.sample_helper': install_wildcards,
        'djerba.plugins.summary': install_wildcards,
        'djerba.plugins.pancurx': install_wildcards,
        'djerba.plugins.pancurx.all_genes': install_wildcards,
        'djerba.plugins.pancurx.blurb': install_wildcards,
        'djerba.plugins.pancurx.germline': install_wildcards,
        'djerba.plugins.pancurx.slide': install_wildcards,
        'djerba.plugins.pancurx.summary': install_wildcards,
        'djerba.plugins.pancurx.appendix': install_wildcards,
        'djerba.plugins.pancurx.classification': install_wildcards,
        'djerba.plugins.pancurx.fusions': install_wildcards,
        'djerba.plugins.pancurx.somatic': install_wildcards,
        'djerba.plugins.pancurx.immune': install_wildcards,
    },
    install_requires=[
        'configparse',
        'email_validator',
        'jsonschema',
        'mako',
        'markdown',
        'pandas',
        'pdfkit',
        'pyinstaller',
        'PyPDF2',
        'requests',
        'statsmodels',
    ],
    python_requires='>=3.10.6',
    authors="Iain Bancarz, Felix Beaudry",
    author_email="ibancarz [at] oicr [dot] on [dot] ca",
    description="Create reports from metadata and workflow output",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/oicr-gsi/djerba",
    keywords=['cancer', 'bioinformatics'],
    license='GPL 3.0',
)
