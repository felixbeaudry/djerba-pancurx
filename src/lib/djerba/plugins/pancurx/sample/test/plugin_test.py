#! /usr/bin/env python3

"""Test of the patient info plugin"""

import json
import logging
import os
import unittest
import djerba.core.constants as cc
from djerba.core.workspace import workspace
from djerba.plugins.plugin_tester import PluginTester
from djerba.plugins.patient_info.plugin import main as patient_info_plugin

class TestPatientInfo(PluginTester):

    def test(self):
        test_source_dir = os.path.realpath(os.path.dirname(__file__))
        params = {
            self.INI: 'patient_info.ini',
            self.JSON: 'patient_info.json',
            self.MD5: 'd41ee59fa3b464a665e8517fe04c12b9'
        }
        self.run_basic_test(test_source_dir, params)

if __name__ == '__main__':
    unittest.main()

