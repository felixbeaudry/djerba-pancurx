#! /usr/bin/env python3

"""Test of the provenance helper"""

import logging
import os
import time
import unittest
from configparser import ConfigParser
import djerba.core.constants as core_constants
from djerba.core.loaders import helper_loader
from djerba.core.workspace import workspace
from djerba.util.testing.tools import TestBase
import djerba.helpers.sample_helper.helper as phoenix_helper
from djerba.util.environment import directory_finder

class TestProvenance(TestBase):

    HELPER_NAME = 'sample_helper'
    CORE = 'core'

    PATH_INFO_MD5 = 'b975845036294eb340e2157429857012'

    def testGetProvenance(self):
        ws = workspace(self.tmp_dir)
        loader = helper_loader(logging.WARNING)
        helper_main = loader.load(self.HELPER_NAME, ws)

        config = helper_main.get_expected_config()
        config.add_section(self.CORE)
        config.set(self.HELPER_NAME, 'donor', 'BTC0056')
        config.set(self.HELPER_NAME, 'sample', 'BTC_0056_Lv_P_526')
        config.set(self.HELPER_NAME, 'normal', 'BTC_0056_Ly_R')
        config = helper_main.configure(config)

        path_info_path = os.path.join(self.tmp_dir, core_constants.DEFAULT_SAMPLE_INFO)
        self.assertTrue(os.path.exists(path_info_path))
        self.assertEqual(self.getMD5(path_info_path), self.PATH_INFO_MD5)

if __name__ == '__main__':
    unittest.main() 

