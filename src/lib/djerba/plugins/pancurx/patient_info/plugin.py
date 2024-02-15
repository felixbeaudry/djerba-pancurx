"""Plugin to generate the "Patient & Physician Info" report section"""

import logging
from time import strptime

import djerba.core.constants as core_constants
from djerba.plugins.base import plugin_base
from djerba.util.render_mako import mako_renderer

class main(plugin_base):

    PRIORITY = 100
    PLUGIN_VERSION = '1.0.0'
    REPORT_VERSION = "report_version"
    TUMOUR_ID = "tumour_id"
    DONOR_ID = "donor_id"
    NORMAL_ID = "normal_id"
    REPORT_DATE = "report_date"
    EXTERNAL_IDS = "external_ids"
    TUMOUR_TISSUE = "tumour_tissue"
    NORMAL_TISSUE = "normal_tissue"

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
            self.REPORT_VERSION,
            self.TUMOUR_ID,
            self.DONOR_ID,
            self.NORMAL_ID,
            self.REPORT_DATE,
            self.EXTERNAL_IDS,
            self.TUMOUR_TISSUE,
            self.NORMAL_TISSUE
        ]
        data[core_constants.RESULTS] = {k: wrapper.get_my_string(k) for k in results_keys}
        return data

    def render(self, data):
        renderer = mako_renderer(self.get_module_dir())
        return renderer.render_name(self.TEMPLATE_NAME, data)

    def specify_params(self):
        discovered = [
            self.REPORT_VERSION,
            self.TUMOUR_ID,
            self.DONOR_ID,
            self.NORMAL_ID,
            self.REPORT_DATE,
            self.EXTERNAL_IDS,
            self.TUMOUR_TISSUE,
            self.NORMAL_TISSUE
        ]
        for key in discovered:
            self.add_ini_discovered(key)
        self.set_ini_default(core_constants.ATTRIBUTES, 'clinical')
        self.set_priority_defaults(self.PRIORITY)
        self.set_ini_default('render_priority', 30)
