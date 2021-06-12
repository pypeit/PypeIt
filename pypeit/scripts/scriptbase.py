"""
Implements base classes for use with pypeit scripts
"""

from IPython import embed
import argparse

# A trick from stackoverflow to allow multi-line output in the help:
#https://stackoverflow.com/questions/3853722/python-argparse-how-to-insert-newline-in-the-help-text
class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            # TODO: We shouldn't be ignoring the terminal width here,
            # but doing anything fancier than splitlines() gets complicated quickly. I
            # think we should be careful with using this formatter to
            # make lines no longer than about 60 characters.
            # import textwrap
            # lines = np.concatenate([textwrap.wrap(t, width) if len(t) > 0 else [' ']
            #                            for t in text[2:].split('\n')]).tolist()
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)


class ScriptBase:
    @classmethod
    def entry_point(cls):
        cls.main(cls.parse_args())

    @classmethod
    def parse_args(cls, options=None):
        parser = cls.get_parser()
        return parser.parse_args() if options is None else parser.parse_args(options)

    # Base classes should override this
    @staticmethod
    def main(args):
        pass

    # Base classes should override this.  Ideally they should use this
    # base-class method to intantiate the ArgumentParser object and then fill in
    # the relevant parser arguments
    @classmethod
    def get_parser(cls, description=None, width=None,
                   formatter=argparse.RawDescriptionHelpFormatter):
        return argparse.ArgumentParser(description=description,
                                       formatter_class=lambda prog: formatter(prog, width=width))



