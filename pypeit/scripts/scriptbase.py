"""
Implements base classes for use with ``PypeIt`` scripts.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import argparse
import datetime
from functools import reduce
import logging
from pathlib import Path
import textwrap

from IPython import embed

from pypeit import log
from pypeit import PypeItError

class SmartFormatter(argparse.HelpFormatter):
    r"""
    Enable a combination of both fixed-format and wrappable lines to be
    formatted for the help statements for command-line arguments used with
    `argparse.ArgumentParser`_.

    Borrows from
    https://stackoverflow.com/questions/3853722/python-argparse-how-to-insert-newline-in-the-help-text

    Help strings that use this formatter *must* begin with "R|".  If not, the
    help string is parsed by the base class.

    When parsed by this formatter, the leading "R|" characters are stripped and
    the lines to be printed are parsed using `str.splitlines`_.  Each resulting
    line is wrapped using `textwrap.wrap`_, unless it begins with the characters
    "F|", which forces the line to remain unaltered (except for stripping the
    leading characters).

    For example, if you add an argument like this:

    .. code-block:: python

        parser.add_argument('-t', '--tell_file', type=str,
                            help='R|Configuration file to change default telluric parameters.  '
                                 'Note that the parameters in this file will be overwritten if '
                                 'you set argument in your terminal.  The --tell_file option '
                                 'requires a .tell file with the following format:\n'
                                 '\n'
                                 'F|    [tellfit]\n'
                                 'F|         objmodel = qso\n'
                                 'F|         redshift = 7.6\n'
                                 'F|         bal_wv_min_max = 10825,12060\n'
                                 'OR\n'
                                 'F|    [tellfit]\n'
                                 'F|         objmodel = star\n'
                                 'F|         star_type = A0\n'
                                 'F|         star_mag = 8.\n'
                                 'OR\n'
                                 'F|    [tellfit]\n'
                                 'F|         objmodel = poly\n'
                                 'F|         polyorder = 3\n'
                                 'F|         fit_wv_min_max = 9000.,9500.\n'
                                 '\n')

    The result will be (depending on the width of your console):

    .. code-block:: console

        -t TELL_FILE, --tell_file TELL_FILE
                          Configuration file to change default telluric
                          parameters.  Note that the parameters in this file
                          will be overwritten if you set argument in your
                          terminal.  The --tell_file option requires a .tell
                          file with the following format:

                              [tellfit]
                                   objmodel = qso
                                   redshift = 7.6
                                   bal_wv_min_max = 10825,12060
                          OR
                              [tellfit]
                                   objmodel = star
                                   star_type = A0
                                   star_mag = 8.
                          OR
                              [tellfit]
                                   objmodel = poly
                                   polyorder = 3
                                   fit_wv_min_max = 9000.,9500.
    """
    def _split_lines(self, text, width):
        """
        Split the provided text into width constrained lines.

        See the class description for formatting instructions.
        """
        if text.startswith('R|'):
            lines = text[2:].splitlines()
            for i in range(len(lines)):
                if lines[i].startswith('F|'):
                    lines[i] = [lines[i][2:]]
                elif len(lines[i]) == 0:
                    lines[i] = [' ']
                else:
                    lines[i] = textwrap.wrap(lines[i], width)
            return reduce(list.__add__, lines)
        return super()._split_lines(text, width)


class ScriptBase:
    """
    Provides a base class for all scripts.
    """
    @classmethod
    def entry_point(cls):
        """
        Defines the main script entry point.
        """
        cls.main(cls.parse_args())

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
        cls._fill_parser_cwd(parser)
        return parser.parse_args() if options is None else parser.parse_args(options)

    @staticmethod
    def _fill_parser_cwd(parser):
        """
        Replace the default of any action that is exactly ``'current working
        directory'`` with the value of ``Path.cwd()``.

        The ``parser`` is edited *in place*.

        Args:
            parser (argparse.ArgumentParser):
                The argument parsing object to edit.
        """
        for action in parser._actions:
            if action.default == 'current working directory':
                action.default = str(Path.cwd())

    # Base classes should override this
    @classmethod
    def main(cls, args):
        """
        Execute the script.
        """
        pass

    @classmethod
    def get_parser(cls, description=None, width=None,
                   formatter=argparse.ArgumentDefaultsHelpFormatter,
                   include_log_options=True,
                   default_log_file=False):
        """
        Construct the command-line argument parser.

        Derived classes should override this.  Ideally they should use this
        base-class method to instantiate the ArgumentParser object and then fill
        in the relevant parser arguments

        .. warning::

            *Any* argument that defaults to the
            string ``'current working directory'`` will be replaced by the
            result of ``Path.cwd()`` when the script is executed.  This means
            help dialogs will include this replacement, and parsing of the
            command line will use ``Path.cwd()`` as the default.  This
            functionality is largely to allow for PypeIt's automated
            documentation of script help dialogs without the "current working"
            directory being that of the developer that most recently compiled
            the docs.

        Parameters
        ----------
        description : :obj:`str`, optional
            A short description of the purpose of the script.
        width : :obj:`int`, optional
            Restrict the width of the formatted help output to be no longer than
            this number of characters, if possible given the help formatter.  If
            None, the width is the same as the terminal width.
        formatter : `argparse.HelpFormatter`_
            Class used to format the help output.
        include_log_options : :obj:`bool`, optional
            Include options that define the logging level(s) and log file.
        default_log_file : :obj:`bool`, optional
            If true, script will use the default log file name if none is
            provided.  Ignored if ``include_log_options`` is False.

        Returns
        -------
        `argparse.ArgumentParser`_
            Command-line interpreter.
        """
        parser = argparse.ArgumentParser(
            description=description, formatter_class=lambda prog: formatter(prog, width=width)
        )
        if not include_log_options:
            return parser
        # Add the logging options
        parser.add_argument(
            '-v', '--verbosity', type=int, default=2,
            help='Verbosity level, which must be 0, 1, or 2.  Level 0 includes warning and error '
                 'messages, level 1 adds informational messages, and level 2 adds debugging '
                 'messages and the calling sequence.'
        )
        parser.add_argument(
            '--log_file', type=str, default='default' if default_log_file else None,
            help='Name for the log file.  If set to "default", a default name is used.  If None, '
                 'a log file is not produced.'
        )
        parser.add_argument(
            '--log_level', type=int, default=None,
            help='Verbosity level for the log file.  If a log file is produce and this is None, '
                 'the file log will match the console stream log.'
        )
        return parser

    @classmethod
    def init_log(cls, args):
        """
        Initialize the logger provided the command-line arguments.
        """
        level = log.convert_verbosity_to_logging_level(args.verbosity)
        log_file_level = None if args.log_level is None else \
            log.convert_verbosity_to_logging_level(args.log_level)
        if args.log_file == 'default':
            _log_file = cls.default_log_file()
        elif args.log_file in ['None', None]:
            _log_file = None
        else:
            _log_file = args.log_file
        log.init(level=level,
                 log_file=_log_file,
                 log_file_level=log_file_level)

    @classmethod
    def default_log_file(cls):
        """
        Set the default name for the log file.
        """
        # Create a UT timestamp (to the minute) for the log filename
        timestamp = datetime.datetime.now(datetime.UTC).strftime("%Y%m%d-%H%M")
        return f'{cls.name()}_{timestamp}.log'

    @staticmethod
    def expandpath(path_pattern):
        """
        Expand a path pattern with wildcards into all matching files.

        Thanks to:
        https://stackoverflow.com/questions/51108256/how-to-take-a-pathname-string-with-wildcards-and-resolve-the-glob-with-pathlib

        Args:
            path_pattern (:obj:`str`):
                Search pattern for files on disk.  Wildcards can occur anywhere
                in the path.
        """
        p = Path(path_pattern).expanduser()
        parts = p.parts[p.is_absolute():]
        return Path(p.root).glob(str(Path(*parts)))



