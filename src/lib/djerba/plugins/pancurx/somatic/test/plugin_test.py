#! /usr/bin/env python3

"""Test of the sample info plugin"""

import os
import unittest
import tempfile
import string

from djerba.util.validator import path_validator
from djerba.plugins.plugin_tester import PluginTester
from djerba.core.workspace import workspace
from djerba.util.environment import directory_finder

class TestSampleInfo(PluginTester):

    INI_NAME = 'test.ini'

    def setUp(self):
        self.path_validator = path_validator()
        self.maxDiff = None
        self.tmp = tempfile.TemporaryDirectory(prefix='djerba_')
        self.tmp_dir = self.tmp.name
        self.sup_dir = directory_finder().get_test_dir()

    def test(self):
        test_source_dir = os.path.realpath(os.path.dirname(__file__))
        with open(os.path.join(test_source_dir, self.INI_NAME)) as in_file:
            template_str = in_file.read()
        template = string.Template(template_str)
        ini_str = template.substitute({'DJERBA_TEST_DATA': self.sup_dir})
        input_dir = os.path.join(self.get_tmp_dir(), 'input')
        os.mkdir(input_dir)
        with open(os.path.join(input_dir, self.INI_NAME), 'w') as ini_file:
            ini_file.write(ini_str)
        params = {
            self.INI: os.path.join(input_dir, self.INI_NAME),
            self.JSON: 'somatic.json',
            self.MD5: '265107a5c86fd60254e8b459c863df6d'
        }
        self.run_basic_test(test_source_dir, params)
    def redact_json_data(self, data):
        """replaces empty method from testing.tools"""
        for key in ['whole_genome_plot_path','histbox_indel','histbox_sv','histbox_snv','indel_vaf_plot','snv_vaf_plot','sv_bin_plot','indel_bin_plot','snv_context_plot','celluloid_plot','oncoslide_cnv_plot']:
            del data['results'][key]
        return data        

if __name__ == '__main__':
    unittest.main()


