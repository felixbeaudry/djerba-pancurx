"""
Plugin to generate the Results Summary report section

"""

import logging
from time import strftime
import os

from djerba.plugins.base import plugin_base, DjerbaPluginError
from djerba.util.render_mako import mako_renderer
import djerba.core.constants as core_constants
from djerba.core.workspace import workspace
import djerba.plugins.pancurx.tools as tools
import djerba.plugins.pancurx.constants as phe

class main(plugin_base):

    PRIORITY = 50
    PLUGIN_VERSION = '0.1'
    MAKO_TEMPLATE_NAME = 'summary_report_template.html'
    SUMMARY_TEMPLATE_FILE = 'summary_template.txt'

    SUMMARY_TEXT = 'summary_text'

    def configure(self, config):
        config = self.apply_defaults(config)
        wrapper = self.get_config_wrapper(config)
        #wrapper = tools.fill_file_if_null(self, wrapper, phe.SUMMARY_FILE, phe.SUMMARY_FILE, core_constants.DEFAULT_SAMPLE_INFO)
        wrapper = tools.try_two_null_files(self, wrapper, phe.SUMMARY_FILE, phe.SUMMARY_FILE, core_constants.DEFAULT_SAMPLE_INFO, 'results_summary.txt')
        return wrapper.get_config()

    def extract(self, config):
        wrapper = self.get_config_wrapper(config)
        summary_path = wrapper.get_my_string(phe.SUMMARY_FILE)
        with open(summary_path) as in_file:
            summary_text = in_file.read()
        self.logger.debug('Read summary from {0}: "{1}"'.format(summary_path, summary_text))
        data = self.get_starting_plugin_data(wrapper, self.PLUGIN_VERSION)
        data[core_constants.RESULTS][self.SUMMARY_TEXT] = summary_text
        tools.copy_if_not_exists(wrapper.get_my_string(phe.SUMMARY_FILE), os.path.join(self.workspace.print_location(), 'results_summary.txt'))
        return data

    def specify_params(self):
        discovered = [
            phe.SUMMARY_FILE
        ]
        for key in discovered:
            self.add_ini_discovered(key)
        self.set_ini_default(core_constants.ATTRIBUTES, 'research')
        self.set_priority_defaults(self.PRIORITY)

    def render(self, data):
        renderer = mako_renderer(self.get_module_dir())
        return renderer.render_name(self.MAKO_TEMPLATE_NAME, data)
