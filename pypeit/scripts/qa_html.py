"""
Built HTML for PYPIT QA
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

    @staticmethod
    def main(args, unit_test=False, path=''):
        """ Builds the HTML files
        path : str, optional
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

        # Master Frame
        if flg_MF:
            qa.gen_mf_html(args.pypeit_file, args.qapath)

        # Exposures
        if flg_exp:
            qa.gen_exp_html()


