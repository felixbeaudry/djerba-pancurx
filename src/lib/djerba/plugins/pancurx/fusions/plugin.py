"""Plugin to generate the "Gremline" report section"""

import logging
from time import strptime
import csv
import os
import re
import collections
import djerba.core.constants as core_constants
from djerba.plugins.base import plugin_base
from djerba.util.render_mako import mako_renderer
import djerba.plugins.pancurx.constants as phe
import djerba.plugins.pancurx.tools as tools

class main(plugin_base):

    PRIORITY = 600
    PLUGIN_VERSION = '1.0.0'
    TEMPLATE_NAME = 'fusions_template.html'


    def configure(self, config):
        config = self.apply_defaults(config)
        wrapper = self.get_config_wrapper(config)

        wrapper = tools.fill_file_if_null(self, wrapper, phe.MAVIS_FUSIONS_PATH, phe.MAVIS_FUSIONS_PATH, core_constants.DEFAULT_PATH_INFO)
        
        return wrapper.get_config()

    def extract(self, config):
        wrapper = self.get_config_wrapper(config)
        attributes = wrapper.get_my_attributes()
        self.check_attributes_known(attributes)
        data = self.get_starting_plugin_data(wrapper, self.PLUGIN_VERSION)
        fusions, fusion_count = tools.parse_fusions(self, wrapper.get_my_string(phe.MAVIS_FUSIONS_PATH))
        data[core_constants.RESULTS]['fusions'] = fusions
        data[core_constants.RESULTS]['fusion_count'] = fusion_count
        return data

    def render(self, data):
        renderer = mako_renderer(self.get_module_dir())
        return renderer.render_name(self.TEMPLATE_NAME, data)

    def specify_params(self):
        discovered = [
            phe.MAVIS_FUSIONS_PATH

        ]
        for key in discovered:
            self.add_ini_discovered(key)
        self.set_ini_default(core_constants.ATTRIBUTES, 'research')
        self.set_priority_defaults(self.PRIORITY)

