"""Plugin to generate the Appendix report section"""

import logging
from time import strftime

import djerba.core.constants as core_constants
from djerba.plugins.base import plugin_base
from djerba.util.render_mako import mako_renderer
import djerba.plugins.pancurx.constants as phe
import djerba.plugins.pancurx.tools as tools

class main(plugin_base):

    PRIORITY = 800
    PLUGIN_VERSION = '1.0.0'
    TEMPLATE_NAME = 'appendix_template.html'

    def configure(self, config):
        config = self.apply_defaults(config)
        wrapper = self.get_config_wrapper(config)

        for this_chromosome_plot in phe.CHROMOSOME_PLOTS.values():
            wrapper = tools.fill_categorized_file_if_null(self, wrapper, this_chromosome_plot, this_chromosome_plot, core_constants.DEFAULT_PATH_INFO, 'whole_chromosome_plots')

        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'tumour', phe.TUMOUR_COVERAGE_PATH, core_constants.DEFAULT_PATH_INFO, 'coveragePaths')
        wrapper = tools.fill_categorized_file_if_null(self, wrapper, 'normal', phe.NORMAL_COVERAGE_PATH, core_constants.DEFAULT_PATH_INFO, 'coveragePaths')
        wrapper = tools.fill_file_if_null(self, wrapper, 'template_type', 'template_type', core_constants.DEFAULT_SAMPLE_INFO)


        return wrapper.get_config()

    def extract(self, config):
        wrapper = self.get_config_wrapper(config)
        attributes = wrapper.get_my_attributes()
        self.check_attributes_known(attributes)
        data = self.get_starting_plugin_data(wrapper, self.PLUGIN_VERSION)
        data[core_constants.RESULTS][phe.TUMOUR_BAM_PATH] = wrapper.get_my_string(phe.TUMOUR_BAM_PATH)
        data[core_constants.RESULTS][phe.NORMAL_BAM_PATH] = wrapper.get_my_string(phe.NORMAL_BAM_PATH)
        data[core_constants.RESULTS]['template_type'] = '_'.join((wrapper.get_my_string('template_type'), self.TEMPLATE_NAME))

        for this_chromosome_plot in phe.CHROMOSOME_PLOTS.values():
            data[core_constants.RESULTS][this_chromosome_plot] = tools.convert_plot(self,wrapper.get_my_string(this_chromosome_plot), this_chromosome_plot)
        return(data)

    def render(self, data):
        renderer = mako_renderer(self.get_module_dir())
        template_name = data[core_constants.RESULTS]['template_type'] 
        return renderer.render_name(template_name, data)

    def specify_params(self):
        self.set_ini_default(core_constants.ATTRIBUTES, 'research')
        self.set_priority_defaults(self.PRIORITY)
        discovered = [
            phe.TUMOUR_BAM_PATH,
            phe.NORMAL_BAM_PATH,
            'template_type'
        ]
        for this_chromosome_plot in phe.CHROMOSOME_PLOTS.values():
            discovered.append(this_chromosome_plot)
        for key in discovered:
            self.add_ini_discovered(key)
