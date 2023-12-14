"""
Built HTML for PYPIT QA

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


class QAHtml(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Script to build HTML files for PYPIT QA.',
                                    width=width)

        parser.add_argument("pypeit_file", type=str, help="PYPIT file")
        parser.add_argument("type", type=str, help="QA Type (MF, exp, all)")
        parser.add_argument("--qapath", type=str, default='QA/',
                            help="Path the QA folder including QA/)")
        return parser

    # TODO: unit_test and path aren't used, right?
    @staticmethod
    def main(args, unit_test=False, path=''):
        """Builds the HTML files.

        Args:
            args (`argparse.Namespace`_):
                Parsed command-line arguments
            unit_test (:obj:`bool`, optional):
                Method being executed for a unit test
            path (:obj:`str`, optional):
                Mainly for running the unit test
        """

        from pypeit.core import qa

        # Flags
        flg_MF, flg_exp = False, False
        if args.type == 'MF':
            flg_MF = True
        elif args.type == 'exp':
            flg_exp = True
        elif args.type == 'all':
            flg_exp, flg_MF = True, True

        if flg_MF:
            qa.gen_mf_html(args.pypeit_file, args.qapath)

        # Exposures
        if flg_exp:
            qa.gen_exp_html()


