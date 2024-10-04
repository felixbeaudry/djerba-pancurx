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
            'data/LBR.germline.genes.txt',
            'data/LBR.summary.csv',
            'data/PCX.germline.genes.txt',
            'data/cosmic_census_20161115.tsv',
            'data/HUGO_synonyms_171213.txt',
            'data/LBR.immune.genes.txt',
            'data/LM22.txt',
            'data/PCX.summary.csv',
            'data/dataset.xml',
            'data/LBR.all_rna.txt',
            'data/LBR.immune_rna.txt',
            'data/MANE.GRCh38.v1.3.summary.txt',
            'data/ensembl_conversion_hg38.txt',
            'data/LBR.genes.txt',
            'data/LBR.signatures.txt',
            'data/PCX.genes.txt'
        ],
        'djerba.core': install_wildcards,
        'djerba.helpers.file_helper': install_wildcards,
        'djerba.helpers.sample_helper': install_wildcards,
        'djerba.plugins.all_genes': install_wildcards,
        'djerba.plugins.blurb': install_wildcards,
        'djerba.plugins.germline': install_wildcards,
        'djerba.plugins.slide': install_wildcards,
        'djerba.plugins.summary': install_wildcards,
        'djerba.plugins.appendix': install_wildcards,
        'djerba.plugins.classification': install_wildcards,
        'djerba.plugins.fusions': install_wildcards,
        'djerba.plugins.somatic': install_wildcards,
        'djerba.plugins.pancurx': install_wildcards,
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
    author="Iain Bancarz",
    author_email="ibancarz [at] oicr [dot] on [dot] ca",
    description="Create reports from metadata and workflow output",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/oicr-gsi/djerba",
    keywords=['cancer', 'bioinformatics'],
    license='GPL 3.0',
)
