"""
Implements base classes for use with ``PypeIt`` scripts.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

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
        """
        Defines the main script entry point.
        """
        cls.main(cls.parse_args())

    # TODO: Combining classmethod and property works in python 3.9 and later
    # only: https://docs.python.org/3.9/library/functions.html#classmethod
    # Order matters.  In python 3.9, it would be:
    #
    # @classmethod
    # @property
    #
    # Because we're not requiring python 3.9 yet, we have to leave this as a
    # classmethod only:
    @classmethod
    def name(cls):
        """
        Provide the name of the script.  By default, this is the name of the
        module with "pypeit" prepended.
        """
        return f"pypeit_{cls.__module__.split('.')[-1]}"

    @classmethod
    def parse_args(cls, options=None):
        """
        Parse the command-line arguments.
        """
        parser = cls.get_parser()
        return parser.parse_args() if options is None else parser.parse_args(options)

    # Base classes should override this
    @staticmethod
    def main(args):
        """
        Execute the script.
        """
        pass

    # Base classes should override this.  Ideally they should use this
    # base-class method to intantiate the ArgumentParser object and then fill in
    # the relevant parser arguments
    @classmethod
    def get_parser(cls, description=None, width=None,
                   formatter=argparse.ArgumentDefaultsHelpFormatter):
        """
        Construct the command-line argument parser.

        Args:
            description (:obj:`str`, optional):
                A short description of the purpose of the script.
            width (:obj:`int`, optional):
                Restrict the width of the formatted help output to be no longer
                than this number of characters, if possible given the help
                formatter.  If None, the width is the same as the terminal
                width.
            formatter (`argparse.HelpFormatter`_):
                Class used to format the help output.

        Returns:
            `argparse.ArgumentParser`_: Command-line interpreter.
        """
        return argparse.ArgumentParser(description=description,
                                       formatter_class=lambda prog: formatter(prog, width=width))



