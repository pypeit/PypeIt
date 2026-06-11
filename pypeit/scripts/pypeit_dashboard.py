"""
main executable for the pypeit dashboard
"""

from pypeit.scripts import scriptbase

# NOTE: I copied this from the run_pypeit.py in pypeit/scripts and am incrementally changing it
class RunDashboard(scriptbase.ScriptBase):

    @classmethod
    def name(cls):
        """
        Return the name of the executable.
        """
        return 'run_dashboard'

    @classmethod
    def usage(cls):
        """
        print dashboard usage description
        """
        # from pypeit import __version__

        descr = 'Dashboard: run pypeit with a GUI\n'
        # descr += f'Version {__version__}\n\n'
        # import textwrap
        # from pypeit.spectrographs import available_spectrographs
        # spclist = ', '.join(available_spectrographs)
        # spcl = textwrap.wrap(spclist, width=70)
        # descr += 'Available spectrographs include:\n'
        # for ispcl in spcl:
        #     descr += f'    {ispcl}\n'
        return descr

    @classmethod
    def get_parser(cls, width=None):
        import argparse

        parser = super().get_parser(description=cls.usage(),
                                    width=width, formatter=argparse.RawDescriptionHelpFormatter,
                                    default_log_file=True)
        parser.add_argument('-f','--pypeit_file', type=str, default=None,
                            help='PypeIt reduction file (must have .pypeit extension)')

        parser.add_argument('-r', '--redux_path', default=None,
                            help='Path to directory for the reduction.  Only advised for testing')

        return parser

    @classmethod
    def main(cls, args):

        from pypeit.dashboard.dashboard import main 


        cls.init_log(args)
        # make it run the dashboard with the args
        main()

        return 0

