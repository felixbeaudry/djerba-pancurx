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

        wrapper = tools.try_two_null_files(self, wrapper, phe.BLURB_FILE, phe.BLURB_FILE, core_constants.DEFAULT_PATH_INFO, 'results_summary.txt')
        return wrapper.get_config()

    def extract(self, config):
        wrapper = self.get_config_wrapper(config)
        data = self.get_starting_plugin_data(wrapper, self.PLUGIN_VERSION)

        summary_path = wrapper.get_my_string(phe.BLURB_FILE)
        tools.copy_if_not_exists(wrapper.get_my_string(phe.BLURB_FILE), os.path.join(self.workspace.print_location(), 'results_summary.txt'))

        with open(summary_path) as in_file:
            summary_text = in_file.read()

        self.logger.debug('Read summary from {0}: "{1}"'.format(summary_path, summary_text))
        data[core_constants.RESULTS][self.SUMMARY_TEXT] = summary_text

        return data

    def specify_params(self):
        discovered = [
            phe.BLURB_FILE
        ]
        for key in discovered:
            self.add_ini_discovered(key)
        self.set_ini_default(core_constants.ATTRIBUTES, 'research')
        self.set_priority_defaults(self.PRIORITY)

    def render(self, data):
        renderer = mako_renderer(self.get_module_dir())
        return renderer.render_name(self.MAKO_TEMPLATE_NAME, data)
