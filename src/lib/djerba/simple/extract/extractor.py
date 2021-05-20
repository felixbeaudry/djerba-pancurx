"""Extract and pre-process data, so it can be read into a clinical report JSON document"""

import json
import os
import pandas as pd
import djerba.simple.constants as constants

class extractor:
    """
    Extract the clinical report data; replaces 4-singleSample.sh
    Input: INI config from 3-configureSingleSample.sh
    Output: Directory of .txt and .json files for downstream processing
    """

    SAMPLE_INFO_KEY = 'sample_info'
    SAMPLE_PARAMS_FILENAME = 'sample_params.json'
    MAF_PARAMS_FILENAME = 'maf_params.json'
    
    def __init__(self, config, bedPath, outDir):
        # config is a ConfigParser object with required parameters (eg. from INI file)
        # INI section header is required by Python configparser, but not written by upstream script
        self.config = config
        self.outDir = outDir
        self.bedPath = bedPath # .bed file for MAF calculation; TODO check readability?
        self.configPaths = []

    def _write_json(self, config, fileName):
        outPath = os.path.join(self.outDir, fileName)
        with open(outPath, 'w') as out:
            out.write(json.dumps(config, sort_keys=True, indent=4))
        return outPath

    def getConfigPaths(self):
        """JSON configuration paths to create reader objects and build the report"""
        return self.configPaths

    def run(self):
        """Run all extractions and write output"""
        self.configPaths.append(self.writeMafParams())
        self.configPaths.append(self.writeIniParams())

    def writeIniParams(self):
        """
        Take parameters directly from the config file, and write as JSON for later use
        Output approximates data_clinical.txt in CGI-Tools, but only has fields for final JSON output
        """
        sampleParams = {}
        sampleParams['PATIENT_ID'] = self.config[constants.CONFIG_HEADER]['patientid'].strip('"')
        stringKeys = [
            'SAMPLE_TYPE',
            'CANCER_TYPE',
            'CANCER_TYPE_DETAILED',
            'CANCER_TYPE_DESCRIPTION',
            'DATE_SAMPLE_RECEIVED',
            'CLOSEST_TCGA',
            'SAMPLE_ANATOMICAL_SITE',
            'SAMPLE_PRIMARY_OR_METASTASIS',
            'SEX'
        ]
        floatKeys = [
            'MEAN_COVERAGE',
            'PCT_v7_ABOVE_80x',
            'SEQUENZA_PURITY_FRACTION',
            'SEQUENZA_PLOIDY'
        ]
        # TODO if value is empty, should we replace with NA? Or raise an error?
        # TODO can other values be used? Is 'patient'=='SAMPLE_ID'?
        for key in stringKeys:
            sampleParams[key] = self.config[constants.CONFIG_HEADER][key].strip('"')
        for key in floatKeys:
            sampleParams[key] = float(self.config[constants.CONFIG_HEADER][key])
        config = {
            constants.READER_CLASS_KEY: 'json_reader',
            self.SAMPLE_INFO_KEY: sampleParams
        }
        return self._write_json(config, self.SAMPLE_PARAMS_FILENAME)

    def writeMafParams(self):
        """Read the MAF file, extract relevant parameters, and write as JSON"""
        maf_path = self.config[constants.CONFIG_HEADER][constants.MAFFILE]
        tmb = maf_extractor(maf_path, self.bedPath).find_tmb()
        config = {
            constants.READER_CLASS_KEY: 'json_reader',
            self.SAMPLE_INFO_KEY: {
                constants.TMB_PER_MB_KEY: tmb
            }
        }
        return self._write_json(config, self.MAF_PARAMS_FILENAME)

class maf_extractor:

    def __init__(self, maf_path, bed_path):
        bed_cols = ['chrom', 'start', 'end']
        self.maf = pd.read_csv(maf_path, sep='\t', skiprows=1)
        self.bed = pd.read_csv(bed_path, sep='\t', skiprows=2, header=None, names=bed_cols)

    def find_tmb(self):
        target_space = sum(self.bed['end'] - self.bed['start']) / 1000000.0
        keep = ['Missense_Mutation', 'Frame_Shift_Del', 'In_Frame_Del', 'Frame_Shift_Ins',
                'In_Frame_Ins', 'Splice_Site', 'Translation_Start_Site', 'Nonsense_Mutation',
                'Nonstop_Mutation']
        tmb = len(self.maf.loc[self.maf["Variant_Classification"].isin(keep)]) / target_space
        return tmb