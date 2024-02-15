"""Plugin to generate the "Patient & Physician Info" report section"""

import logging
from time import strptime

import djerba.core.constants as core_constants
from djerba.plugins.base import plugin_base
from djerba.util.render_mako import mako_renderer

class main(plugin_base):

    PRIORITY = 200
    PLUGIN_VERSION = '1.0.0'
    TEMPLATE_NAME = 'sample_template.html'

    TUMOUR_COVERAGE = "tumour_coverage"
    NORMAL_COVERAGE = "normal_coverage"
    CELLULARITY = "cellularity"
    PLOIDY = "ploidy"
    ALEXANDROV_CLASS = "alexandrov_class"
    COLLISSON_CLASS = "collisson_class"
    WADDELL_CLASS = "waddell_class"
    MOFFITT_CLASS = "moffitt_class"
    INFERRED_SEX = "inferred_sex"
    HLA_TYPES = "hla_types"

    def configure(self, config):
        config = self.apply_defaults(config)
        wrapper = self.get_config_wrapper(config)
        return wrapper.get_config()

    def extract(self, config):
        wrapper = self.get_config_wrapper(config)
        attributes = wrapper.get_my_attributes()
        self.check_attributes_known(attributes)
        data = self.get_starting_plugin_data(wrapper, self.PLUGIN_VERSION)
        results_keys = [
            self.TUMOUR_COVERAGE,
            self.NORMAL_COVERAGE,
            self.CELLULARITY,
            self.PLOIDY,
            self.ALEXANDROV_CLASS,
            self.COLLISSON_CLASS,
            self.WADDELL_CLASS,
            self.MOFFITT_CLASS,
            self.INFERRED_SEX,
            self.HLA_TYPES
        ]
        data[core_constants.RESULTS] = {k: wrapper.get_my_string(k) for k in results_keys}
        return data

    def render(self, data):
        renderer = mako_renderer(self.get_module_dir())
        return renderer.render_name(self.TEMPLATE_NAME, data)

    def specify_params(self):
        discovered = [
            self.TUMOUR_COVERAGE,
            self.NORMAL_COVERAGE,
            self.CELLULARITY,
            self.PLOIDY,
            self.ALEXANDROV_CLASS,
            self.COLLISSON_CLASS,
            self.WADDELL_CLASS,
            self.MOFFITT_CLASS,
            self.INFERRED_SEX,
            self.HLA_TYPES
        ]
        for key in discovered:
            self.add_ini_discovered(key)
        self.set_ini_default(core_constants.ATTRIBUTES, 'clinical')
        self.set_priority_defaults(self.PRIORITY)
