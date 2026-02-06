#!/usr/bin/env python

"""
Check that the installation of python is C enabled.
"""

from pypeit.scripts import scriptbase

class ChkPlugins(scriptbase.ScriptBase):

    @classmethod
    def main(cls, args):

        from pypeit.display import required_plugins, plugins_available
        from pypeit import log
        from pypeit import PypeItError

        # Initialize the log
        cls.init_log(args)

        success, report = plugins_available(return_report=True)
        if not success:
            raise PypeItError(report)
        log.info('All required plugins found: {0}'.format(', '.join(required_plugins)))


