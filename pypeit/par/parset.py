"""
Define a utility base class used to hold parameters.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import os
from pathlib import Path
import shutil
import textwrap

from astropy.io import fits
from IPython import embed
import numpy
import numpy as np

from pypeit import log
from pypeit.par import util


# TODO: Include a "mutable" attribute that decides if a parameter can be
# changed?


class ParSet:
    """
    Generic base class to handle and manipulate a list of operational
    parameters.  A glorified dictionary that constrains and types its
    components.

    .. todo::
        - Write a test for equality?

    Args:
        pars (list):
            A list of keywords for a list of parameter values.
        values (:obj:`list`, optional):
            Initialize the parameters to these values.  If not provided,
            all parameters are initialized to `None` or the provided
            default.
        defaults (:obj:`list`, optional):
            For any parameters not provided in the *values* list, use
            these default values.  If not provided, no defaults are
            assumed.
        options (:obj:`list`, optional):
            Force the parameters to be one of a list of options.  Each
            element in the list can be a list itself.  If not provided,
            all parameters are allowed to take on any value within the
            allowed data type.
        dtypes (:obj:`list`, optional):
            Force the parameter to be one of a list of data types.  Each
            element in the list can be a list itself.  If not provided,
            all parameters are allowed to have any data type.
        can_call (:obj:`list`, optional): Flag that the parameters are
            callable operations.  Default is False.
        descr (:obj:`list`, optional):
            A list of parameter descriptions.  Empty strings by default.
        cfg_section (:obj:`str`, optional): 
            The top-level designation for a configuration section
            written based on the contents of this parameter set.
        cfg_comment (:obj:`str`, optional): 
            Comment to be placed at the top-level of the configuration
            section written based on the contents of this parameter set.

    Raises:
        TypeError:
            Raised if the input parameters are not lists or if the input
            keys are not strings.
        ValueError:
            Raised if any of the optional arguments do not have the same
            length as the input list of parameter keys.

    Attributes:
        npar (int):
            Number of parameters
        data (dict):
            Dictionary with the parameter values
        default (dict):
            Dictionary with the default values
        options (dict):
            Dictionary with the allowed options for the parameter values
        dtype (dict):
            Dictionary with the allowed data types for the parameters
        can_call (dict):
            Dictionary with the callable flags
        descr (dict):
            Dictionary with the description of each parameter.
        cfg_section (str): 
            The top-level designation for a configuration section
            written based on the contents of this parameter set.
        cfg_comment (str): 
            Comment to be placed at the top-level of the configuration
            section written based on the contents of this parameter set.
    """

    # Set prefix for writing parameters to a fits header to a class
    # attribute.
    prefix = 'PAR'
    """
    Class Prefix for header keywords when writing the parset to an
    `astropy.io.fits.Header`_ object.
    """

    def __init__(self, pars, values=None, defaults=None, options=None, dtypes=None, can_call=None,
                 descr=None, cfg_section=None, cfg_comment=None):
        # Check that the list of input parameters is a list of strings
        if not isinstance(pars, list):
            raise TypeError('Input parameter keys must be provided as a list.')
        for key in pars:
            if not isinstance(key, str):
                raise TypeError('Input parameter keys must be strings.')

        # Get the length of the parameter list and make sure the list
        # has unique values
        self.npar = len(pars)
        if len(numpy.unique(numpy.array(pars))) != self.npar:
            raise ValueError('All input parameter keys must be unique.')

        # Check that the other lists, if provided, have the correct type
        # and length
        if values is not None and (not isinstance(values, list) or len(values) != self.npar):
            raise ValueError('Values must be a list with the same length as the keys list.')
        if defaults is not None and (not isinstance(defaults, list) or len(defaults) != self.npar):
            raise ValueError('Defaults must be a list with the same length as the keys list.')
        if options is not None and (not isinstance(options, list) or len(options) != self.npar):
            raise ValueError('Options must be a list with the same length as the keys list.')
        if dtypes is not None and (not isinstance(dtypes, list) or len(dtypes) != self.npar):
            raise ValueError('Data types list must have the same length as the keys list.')
        if can_call is not None and (not isinstance(can_call, list) or len(can_call) != self.npar):
            raise ValueError('List of callable flags must have the same length as keys list.')
        if descr is not None and (not isinstance(descr, list) or len(descr) != self.npar):
            raise ValueError('List of parameter descriptions must have the same length as '
                             'keys list.')

        # Set up dummy lists for no input
        _values = [None]*self.npar if values is None else values
        _defaults = [None]*self.npar if defaults is None else defaults
        _options = [None]*self.npar if options is None else options
        _dtypes = [None]*self.npar if dtypes is None else dtypes
        _can_call = [False]*self.npar if can_call is None else can_call
        _descr = ['']*self.npar if descr is None else descr

        # Set the defaults
        self.default = dict([ (p, d) for p, d in zip(pars, _defaults) ])

        # Set the valid options
        self.options = dict([ (p, [o]) if o is not None and not isinstance(o, list) else (p, o) \
                                       for p, o in zip(pars, _options) ])
        # Set the valid types
        self.dtype = dict([ (p, [t]) if t is not None and not isinstance(t, list) else (p, t) \
                                     for p, t in zip(pars, _dtypes) ])

        # Set the calling flags
        self.can_call = dict([ (p, t) for p, t in zip(pars, _can_call) ])

        # Set the calling flags
        self.descr = dict([ (p, t) for p, t in zip(pars, _descr) ])

        # Set the data dictionary using the overloaded
        # __setitem__function so that value checking is performed
        self.data = {}
        for p, d, v, t in zip(pars, _defaults, _values, _dtypes):
            # Check if 'None' is an allowed option
            none_allowed = False
            if type(t) is list:
                if type(None) in t:
                    none_allowed = True
            if v is None and not none_allowed:
                self.__setitem__(p, d)
                continue
            self.__setitem__(p, v)

        # Save the configuration file section details
        self.cfg_section = cfg_section
        self.cfg_comment = cfg_comment

    def __getitem__(self, key):
        """
        Return the value of the designated key.

        Args:
            key (str):
                Key for new parameter
        """
        return self.data[key]

    def __setitem__(self, key, value):
        """
        Set the value for a key.

        Args:
            key (str):
                Key for new parameter
            value (*dtype*) : Parameter value, must have a type in the
                list provided (*dtype*), if the list is provided

        Raises:
            ValueError: Raised if the parameter value is not among the
                allowed options (:attr:`options`).
            TypeError: Raised if the parameter value does not have an
                allowed data type (:attr:`dtype`) or if the provided
                value is not a callable object, but is expected to be by
                :attr:`can_call`.
        """
        if value is None:
            self.data[key] = value
            return

        if isinstance(value, list):
            is_parset_or_dict = [ isinstance(v, (ParSet, dict)) for v in value ]
            if numpy.any(is_parset_or_dict) and not numpy.all(is_parset_or_dict):
                log.warning(
                    "List includes a mix of ParSet and dicts with other types.  "
                    "Displaying and writing the ParSet will not be correct!"
                )

        if self.options[key] is not None:
            if isinstance(value, list):
                # `value` can be a list of items, all of which must be
                # one of the valid options.
                indx = numpy.isin(value, self.options[key], invert=True)
                if numpy.any(indx):
                    raise ValueError('Input value for {0} invalid'.format(key)
                                     + '; {0}'.format(numpy.atleast_1d(value)[indx]) 
                                     + ' are not valid options.\n'
                                     + 'Options are: {0}'.format(self.options[key]))
            elif value not in self.options[key]:
                raise ValueError('Input value for {0} invalid: {1}.\nOptions are: {2}'.format(
                                                                    key, value, self.options[key]))
        if self.dtype[key] is not None \
                and not any([ isinstance(value, d) for d in self.dtype[key]]):
            raise TypeError('Input value for {0} has incorrect type: {1}.'.format(key, value) +
                            '\nValid types are: {0}'.format(self.dtype[key]))

        if self.can_call[key] and not callable(value):
            raise TypeError('{0} is not a callable object.'.format(value))

        self.data[key] = value


    def __len__(self):
        """Return the number of parameters."""
        return self.npar
        

    def __iter__(self):
        """Return an iterable to the parameter values."""
        return iter(self.data.values())


    def __repr__(self):
        """Return a string representation of the parameters."""
        return self._output_string(header=self.cfg_section)

    
    def _output_string(self, header=None, value_only=False):
        """
        Constructs the short-format table strings for the
        :func:`__repr__` method.
        

        Args:
            header (:obj:`str`, optional):
                String header to provide for the table.  This is
                typically the name of the configuration section.
            value_only (:obj:`bool`, optional):
                By default, the table includes the parameter key, its
                current value, the default value, its data type, and if
                the value can be a callable function.  If
                `value_only=True`, only the parameter key and current
                value are returned.

        Returns:
            str: Single long string with the parameter table for the
            :func:`__repr__` method.
        """
        additional_par_strings = []
        ncol = 2 if value_only else 5
        data_table = numpy.empty((self.npar+1, ncol), dtype=object)
        data_table[0,:] = ['Parameter', 'Value'] if value_only \
                            else ['Parameter', 'Value', 'Default', 'Type', 'Callable']
        for i, k in enumerate(self.keys()):
            data_table[i+1,0] = k
            if isinstance(self.data[k], ParSet):
                _header = k if header is None else '{0}:{1}'.format(header, k)
                additional_par_strings += [ self.data[k]._output_string(header=_header,
                                                                        value_only=value_only) ]
                data_table[i+1,1] = 'see below'
                if not value_only:
                    data_table[i+1,2] = 'see below'
            else:
                data_table[i+1,1] = ParSet._data_string(self.data[k])
                if not value_only:
                    data_table[i+1,2] = ParSet._data_string(self.default[k])
            if value_only:
                continue

            data_table[i+1,3] = 'Undefined' if self.dtype[k] is None \
                                    else ', '.join([t.__name__ for t in self.dtype[k]])
# TODO: Now treating None's differently. None's shouldn't be in a list.
# Keep this code around for now in case we find a failure mode.
#                                    else ', '.join(['Undefined' if t is None else t.__name__ 
#                                                    for t in self.dtype[k]])
            data_table[i+1,4] = self.can_call[k].__repr__()

        output = [ParSet._data_table_string(data_table)]
        if header is not None:
            output = [header] + output
        if len(additional_par_strings) > 0:
            output += additional_par_strings
        return '\n'.join(output)


    @staticmethod
    def _data_table_string(data_table, delimeter='print'):
        """
        Provided the array of data, format it with equally spaced
        columns and add a header (first row) and contents delimeter.

        Args:
            data_table (`numpy.ndarray`_):
                Array of string representations of the data to print.

        Returns:
            str: Single long string with the data table.
        """
        nrows, ncols = data_table.shape
        col_width = [ numpy.amax([ len(dij) for dij in dj]) for dj in data_table.T ]
        row_string = ['']*(nrows+1) if delimeter == 'print' else ['']*(nrows+3)
        start = 2 if delimeter == 'print' else 3
        for i in range(start,nrows+start-1):
            row_string[i] = '  '.join([ data_table[1+i-start,j].ljust(col_width[j]) 
                                                                        for j in range(ncols)])
        if delimeter == 'print':
            # Heading row
            row_string[0] = '  '.join([ data_table[0,j].ljust(col_width[j]) for j in range(ncols)])
            # Delimiter
            row_string[1] = '-'*len(row_string[0])
            return '\n'.join(row_string)+'\n'

        # For an rst table
        row_string[0] = '  '.join([ '='*col_width[j] for j in range(ncols)])
        row_string[1] = '  '.join([ data_table[0,j].ljust(col_width[j]) for j in range(ncols)])
        row_string[2] = row_string[0]
        row_string[-1] = row_string[0]
        return '\n'.join(row_string)+'\n'

    @staticmethod
    def _data_string(data, use_repr=False, verbatim=False, check_dir=False):
        """
        Convert a single datum into a string
        
        Simply return strings, recursively convert the elements of any
        objects with a :attr:`__len__` attribute, and use the object's
        own :attr:`__repr__` attribute for all other objects.

        Args:
            data (object):
                The object to stringify.
            use_repr (:obj:`bool`, optional):
                Use the objects :attr:`__repr__` method; otherwise, use
                a direct string conversion.
            verbatim (:obj:`bool`, optional):
                Use quotes around the provided string to indicate that
                the string should be representated in a verbatim (fixed
                width) font.
            check_dir (:obj:`bool`, optional):
                If ``data`` is a string, check if it matches the current working
                directory and replace it with a generic string if it does.
        
        Returns:
            str: A string representation of the provided ``data``.
        """
        if isinstance(data, str):
            _data = '$PWD' if check_dir and data == os.getcwd() else data
            if verbatim:
                return '..' if len(_data) == 0 else '``' + _data + '``'
            return _data
        if isinstance(data, list):
            # When the list is empty, return an empty string, which config_lines will append a "," to.
            # This allows ConfigObj to interpret it as an empty list, instead of string, when re-reading the
            # configuration lines into a ParSet
            
            return '' if len(data) == 0 \
                        else ', '.join([ ParSet._data_string(d, use_repr=use_repr,
                                                             verbatim=verbatim) for d in data ])
        return data.__repr__() if use_repr else str(data)

    def _wrap_print(self, head, output, tcols):
        """
        Wrap the contents of an output string for a fixed terminal
        width.  Used for the long-format :func:`info` method.

        Args:
            head (str):
                The inline header for the output.  Can be an empty
                string, but cannot be None.
            output (str):
                The main body of the text to write.
            tcols (int):
                The allowed width for the output.
        """
        tail = ' '*len(head)
        if tcols is not None:
            lines = textwrap.wrap('{0}'.format(output), tcols-len(head))
            if len(lines) == 0:
                print('{0}None'.format(head))
            else:
                _head = [ head ] + [ tail ]*(len(lines)-1)
                print('\n'.join([ h+l for h,l in zip(_head, lines)]))
        else:
            print(head+'{0}'.format(output))


    def _types_list(self, key):
        """Return the string names for the specified data types."""
        return [t.__name__ for t in self.dtype[key]]


    @staticmethod
    def config_lines(par, section_name=None, section_comment=None, section_level=0,
                     exclude_defaults=False, include_descr=True):
        """
        Recursively generate the lines of a configuration file based on
        the provided :class:`ParSet` or :obj:`dict` (see ``par``).

        Args:
            section_name (:obj:`str`, optional):
                Name to give to the top-level of the configuration
                output.
            section_comment (:obj:`str`, optional):
                Description to provide for the top-level configuration
                output.
            section_level (:obj:`int`, optional):
                The level for the configuration output.  Sets the
                indentation level and the number of square brackets
                assigned to the section name.
            exclude_defaults (:obj:`bool`, optional):
                Do not include any parameters that are identical to the
                defaults.
            include_descr (:obj:`bool`, optional):
                Include the descriptions of each parameter as comments.

        Returns:
            list: The list of the lines to write to a configuration
            file.
        """
        # Get the list of parameters that are ParSets
        parset_keys = [ k for k in par.keys() if isinstance(par[k], (ParSet, dict)) ]
        n_parsets = len(parset_keys)

        # Set the top-level comment and section name
        section_indent = ' '*4*section_level
        component_indent = section_indent + ' '*4
        lines = [] if section_comment is None \
                            else ParSet._config_comment(section_comment, section_indent)
        lines += [ section_indent + '['*(section_level+1) + section_name
                   + ']'*(section_level+1) ]

        min_lines = len(lines)

        # Add all the parameters that are not ParSets
        for k in par.keys():
            # Skip it if this element is a ParSet
            if n_parsets > 0 and k in parset_keys:
                continue

            # If the value is a list, determine if all the elements of
            # the list are also dictionaries or ParSets
            if isinstance(par[k], list) and len(par[k]) > 0:
                is_parset_or_dict = [ isinstance(v, (ParSet, dict)) for v in par[k] ]
                if numpy.all(is_parset_or_dict):
                    ndig = int(numpy.log10(len(par[k])))+1
                    for i, v in enumerate(par[k]):
                        indx = str(i+1).zfill(ndig)
                        # Try to add the section comment
                        section_comment = None
                        if include_descr:
                            try:
                                section_comment = par.descr[k] + ': ' + indx
                            except:
                                pass
                        lines += ParSet.config_lines(v, section_name=k+indx,
                                                     section_comment=section_comment,
                                                     section_level=section_level+1,
                                                     exclude_defaults=exclude_defaults,
                                                     include_descr=include_descr)
                    continue

            # Working with a single element
            # Try to add the description for this parameter
            try:
                if par.descr[k] is not None and include_descr:
                    lines += ParSet._config_comment(par.descr[k], component_indent)
            except:
                pass
            if not exclude_defaults or par[k] != par.default[k]:
                argvalue = ParSet._data_string(par[k])
                if isinstance(par[k], list):
                    argvalue += ','
                lines += [ component_indent + k + ' = ' + argvalue ]
        # Then add the items that are ParSets as subsections
        for k in parset_keys:
            section_comment = None
            if include_descr:
                try:
                    section_comment = par.descr[k]
                except:
                    pass
            lines += ParSet.config_lines(par[k], section_name=k, section_comment=section_comment,
                                         section_level=section_level+1,
                                         exclude_defaults=exclude_defaults,
                                         include_descr=include_descr)
        return lines if len(lines) > min_lines else []


    @staticmethod
    def _config_comment(comment, indent, full_width=72):
        """
        Create the list of lines for the description of a given
        parameter in the configuration file.

        Args:
            comment (str):
                The description of the configuration parameter.
            indent (str):
                The string used to indent the text.
            full_width (:obj:`int`, optional):
                The full width allowed for each output string in the
                returned list.

        Returns:
            list: List of the strings to write to the output
            configuration file.
        """
        head = indent + '# '
        lines = textwrap.wrap('{0}'.format(comment), full_width-len(head))
        return [ head + l for l in lines ]
   

    def info(self, basekey=None):
        """
        A long-form version of __repr__ that includes the parameter descriptions.
        """
        # Try to get the width of the available space to print
        try:
            tr, tcols = numpy.array(os.popen('stty size', 'r').read().split()).astype(int)
            tcols -= int(tcols*0.1)
        except:
            tr = None
            tcols = None

        for k in self.data.keys():
            if isinstance(self.data[k], ParSet):
                self.data[k].info(basekey=k)
                continue
            print('{0}'.format(k) if basekey is None else '{0}:{1}'.format(basekey,k))
            self._wrap_print('        Value: ', self.data[k], tcols)
            self._wrap_print('      Default: ', self.default[k], tcols)
            self._wrap_print('      Options: ', 'None' if self.options[k] is None
                                                else ', '.join(self.options[k]), tcols)
            self._wrap_print('  Valid Types: ', 'None' if self.dtype[k] is None
                                                else ', '.join(self._types_list(k)), tcols)
            self._wrap_print('     Callable: ', self.can_call[k], tcols)
            self._wrap_print('  Description: ', self.descr[k], tcols)
            print(' ')


    def keys(self):
        """Return the list of parameter set keys."""
        return list(self.data.keys())

    
    def add(self, key, value, default=None, options=None, dtype=None, can_call=None, descr=None):
        """
        Add a new parameter.

        Args:
            key (:obj:`str`):
                Key for new parameter
            value (:obj:`dtype`):
                Parameter value, must have a type in the list provided
                by `dtype`, if the list is provided
            default (:obj:`dtype`, optional):
                Define a default value for the parameter, must have a
                type in the list provided by `dtype`, if the list
                is provided.  No default if not provided.
            options (:obj:`list`, optional):
                List of discrete values that the parameter is allowed to
                have.  Allowed to be anything if not provided.
            dtype (:obj:`list`, optional):
                List of allowed data types that the parameter can have.
                Allowed to be anything if not provided.
            can_call (:obj:`bool`, optional):
                Flag that the parameters are callable operations.
                Default is False.

        Raises:
            ValueError:
                Raised if the keyword alread exists.
        """
        if key in self.data.keys():
            raise ValueError('Keyword {0} already exists and cannot be added!')
        self.npar += 1
        self.default[key] = None if default is None else default
        self.options[key] = [options] if options is not None and not isinstance(options, list) \
                                      else options
        self.dtype[key] = [dtype] if dtype is not None and not isinstance(dtype, list) else dtype
        self.can_call[key] = False if can_call is None else can_call
        self.descr[key] = None if descr is None else descr
        try:
            self.__setitem__(key, value)
        except:
            # Delete the added components
            del self.default[key]
            del self.options[key]
            del self.dtype[key]
            del self.can_call[key]
            del self.descr[key]
            # Re-raise the exception
            raise

    def to_config(self, cfg_file=None, section_name=None, section_comment=None, section_level=0,
                  append=False, quiet=False, exclude_defaults=False, include_descr=True):
        """
        Write/Append the parameter set to a configuration file.

        Args:
            cfg_file (:obj:`str`, optional):
                The name of the file to write/append to.  If None (default), the
                function will just return the list of strings that would have
                been written to the file.  These lines can be used to construct
                a `configobj`_ instance.
            section_name (:obj:`str`, optional):
                The top-level name for the config section.  This must be
                provided if :attr:`cfg_section` is None or any of the
                parameters are not also :class:`ParSet` instances themselves.
            section_comment (:obj:`str`, optional):
                The top-level comment for the config section based on
                this :class:`ParSet`.
            section_level (:obj:`int`, optional):
                The top level of this :class:`ParSet`.  Used for recursive output
                of nested :class:`ParSet` instances.
            append (:obj:`bool`, optional):
                Append this configuration output of this :class:`ParSet` to the
                file.  False by default.  If not appending and the file
                exists, the file is automatically overwritten.
            quiet (:obj:`bool`, optional):
                Suppress all standard output from the function.
            exclude_defaults (:obj:`bool`, optional):
                Do not include any parameters that are identical to the
                defaults.
            include_descr (:obj:`bool`, optional):
                Include the descriptions of each parameter as comments.

        Raises:
            ValueError:
                Raised if there are types other than :class:`ParSet` in the
                parameter list, :attr:`cfg_section` is None, and no
                section_name argument was provided.
        """
        if cfg_file is not None and os.path.isfile(cfg_file) and not append and not quiet:
            log.warning("Selected configuration file already exists and will be overwritten!")

        config_output = []
        if numpy.all([ isinstance(d, ParSet) or d is None for d in self.data.values() ]):
            # All the elements are ParSets themselves, so just iterate
            # through each one
            for k in self.keys():
                if self.data[k] is None:
                    continue
                section_comment = self.descr[k] if include_descr else None
                config_output += ParSet.config_lines(self.data[k], section_name=k,
                                                     section_comment=section_comment,
                                                     section_level=section_level,
                                                     exclude_defaults=exclude_defaults,
                                                     include_descr=include_descr)
#                config_output += ['']
        else:
            # Cannot write the parameters as a configuration file
            # without a top-level configuration section
            if section_name is None and self.cfg_section is None:
                raise ValueError('No top-level section name available for configuration!')

            _section_name = self.cfg_section if section_name is None else section_name
            _section_comment = self.cfg_comment if section_comment is None else section_comment
            config_output += ParSet.config_lines(self, section_name=_section_name,
                                                 section_comment=_section_comment,
                                                 section_level=section_level,
                                                 exclude_defaults=exclude_defaults,
                                                 include_descr=include_descr)

        if cfg_file is None:
            # Only return the list of lines for the output file.  Useful
            # if you want to use instantly create a new ConfigObj
            # instance without having to write a file
            return config_output

        # Write the file
        with open(cfg_file, 'a' if append else 'w') as f:
            f.write('\n'.join(config_output))

    @staticmethod
    def _rst_class_name(p):
        return ':class:`~' +  type(p).__module__ + '.' + type(p).__name__ + '`'

    def to_rst_table(self, parsets_listed=[]):
        """
        Construct a reStructuredText table describing the parameter set.

        Args:
            parsets_listed (:obj:`list`, optional):
                A list of 
        
        Returns:
            list: Returns a list of lines that can be written to an
            ``*.rst`` file.
        """
        new_parsets = []
        data_table = numpy.empty((self.npar+1, 5), dtype=object)
        data_table[0,:] = ['Key', 'Type', 'Options', 'Default', 'Description']
        sorted_keys = numpy.sort(self.keys())
        for i,k in enumerate(sorted_keys):
            data_table[i+1,0] = ParSet._data_string(k, use_repr=False, verbatim=True)
            if isinstance(self.data[k], ParSet):
                if type(self.data[k]).__name__ not in parsets_listed:
                    new_parsets += [k]
                parsets_listed += [ type(self.data[k]).__name__ ]
                data_table[i+1,1] = ParSet._rst_class_name(self.data[k])
                data_table[i+1,3] = '`{0} Keywords`_'.format(type(self.data[k]).__name__)
            else: 
                data_table[i+1,1] = ', '.join([t.__name__ for t in self.dtype[k]])
                data_table[i+1,3] = '..' if self.default[k] is None \
                                    else ParSet._data_string(self.default[k], use_repr=False,
                                                             verbatim=True, check_dir=True)

            data_table[i+1,2] = '..' if self.options[k] is None \
                                    else ParSet._data_string(self.options[k], use_repr=False,
                                                             verbatim=True)
            data_table[i+1,4] = '..' if self.descr[k] is None \
                                    else ParSet._data_string(self.descr[k])

        output = [ f'.. _{type(self).__name__.lower()}:']
        output += [ '' ]
        output += [ f'{type(self).__name__} Keywords']
        output += [ '-'*len(output[2]) ]
        output += [ '' ]
        output += ['Class Instantiation: ' + ParSet._rst_class_name(self)]
        output += ['']
        output += [ParSet._data_table_string(data_table, delimeter='rst')]
        output += ['']
        for k in new_parsets:
            output += ['----']
            output += ['']
            output += self.data[k].to_rst_table(parsets_listed=parsets_listed)
        return output

    def validate_keys(self, required=None, can_be_None=None):
        if required is None and can_be_None is None:
            # No validation rules, so implicitly valid
            return

        if required is not None:
            not_defined = numpy.array([ k not in self.keys() for k in required ])
            if numpy.any(not_defined):
                raise ValueError('Required keys were not defined: {0}'.format(
                                    numpy.asarray(required)[not_defined].tolist()))

        if can_be_None is not None:
            should_not_be_None = numpy.array([ self.data[k] is None and k not in can_be_None 
                                                                    for k in self.keys()])
            if numpy.any(should_not_be_None):
                raise ValueError('These keys should not be None: {0}'.format(
                                    numpy.asarray(self.keys())[should_not_be_None].tolist()))

    def to_header(self, hdr=None, prefix=None, quiet=False):
        """
        Write the parameters to a fits header.

        Any element that has a value of None or is a :class:`ParSet` itself is
        *not* written to the header.

        Args:
            hdr (`astropy.io.fits.Header`, optional):
                Header object for the parameters. **If provided, the
                header is not copied and directly modified.** If
                None, the baseline header is empty.
            prefix (:obj:`str`, optional):
                Prefix to use for the header keywords, which
                overwrites the string defined for the class. If None,
                uses the default for the class.
            quiet (:obj:`bool`, optional):
                Suppress print statements.

        Returns:
            `astropy.io.fits.Header`_: Header with the parameter
            data included.
        """
        if hdr is None:
            hdr = fits.Header()
        if prefix is None:
            prefix = self.prefix
        ndig = int(numpy.log10(self.npar))+1 
        for i, (key, value) in enumerate(self.data.items()):
            if value is None:
                # Don't write Nones
                continue
            if isinstance(value, ParSet):
                if not quiet:
                    log.warning(
                        "ParSets within ParSets are not written to headers!  "
                        f"Skipping {key}."
                    )
                continue
            if isinstance(value, list):
                _value = '[' + ', '.join([str(v) for v in value]) + ']'
            elif isinstance(value, tuple):
                _value = '(' + ', '.join([str(v) for v in value]) + ')'
            else :
                _value = value 
            hdr['{0}K{1}'.format(prefix, str(i+1).zfill(ndig))] \
                    = (key, '{0}: Key'.format(self.__class__.__name__))
            hdr['{0}{1}'.format(prefix, str(i+1).zfill(ndig))] \
                    = (_value, '{0}: Value'.format(self.__class__.__name__))
        return hdr

    @classmethod
    def from_header(cls, hdr, prefix=None):
        """
        Instantiate the :class:`ParSet` using data parsed from a fits header.

        This is a simple wrapper for :func:`ParSet.parse_par_from_hdr` and
        :func:`ParSet.from_dict`.

        .. warning::

            The to/from header methods in the :class:`ParSet` base
            class **only** saves the parameter keys and values, not
            its other attributes (e.g., options, dtypes, etc). In
            essentially all use cases, the :class:`ParSet` should be
            used as a base class where the :func:`from_dict` method
            is overwritten such that these higher-level attributes of
            the derived class are maintained in the header I/O. See,
            e.g., :mod:`pypeit.par.pypeitpar`.

        Args:
            hdr (`astropy.io.fits.Header`):
                Header object with the parameters.
            prefix (:obj:`str`, optional):
                Prefix of the relevant header keywords, which
                overwrites the string defined for the class. If None,
                uses the default for the class.
        """
        if prefix is None:
            prefix = cls.prefix
        return cls.from_dict(util.recursive_dict_evaluate(cls.parse_par_from_hdr(hdr, prefix)))

    @staticmethod
    def parse_par_from_hdr(hdr, prefix):
        """
        Parse the dictionary of parameters written to a header

        Args:
            hdr (`astropy.io.fits.Header`):
                Header object to parse.
            prefix (:obj:`str`):
                The prefix used for the header keywords.
        
        Returns:
            dict: A dictionary with the parameter keywords and
            values.
        """
        # TODO: I don't like having to find the number of parameters first...
        # Find the numbers of parameters
        npar = 0
        for k, v in hdr.items():
            if k[:len(prefix)] == prefix:
                if k[len(prefix)] == 'K':
                    continue
                try:
                    i = int(k[len(prefix):])
                except ValueError:
                    continue
                if npar < i:
                    npar = i
        ndig = int(numpy.log10(npar))+1 

        par = {}
        for k, v in hdr.items():
            # Check if this header keyword starts with the required
            # prefix
            if k[:len(prefix)] == prefix:
                if k[len(prefix)] == 'K':
                    # Skip keys
                    continue
                try:
                    # Try to convert the keyword without the prefix into
                    # an integer.
                    i = int(k[len(prefix):])
                except ValueError:
                    # Assume the value is some other random keyword
                    # that starts with the prefix but doesn't contain a
                    # relevant parameter
                    continue
                # Assume we've found a parameter. Use the associated
                # header card to set the keyword and add both to the
                # output dictionary.
                par[hdr['{0}K{1}'.format(prefix,str(i).zfill(ndig))]] = v
        return par

    def to_dict(self):
        """
        Return a dictionary with the contents of the parameter set.

        This method recursively handles nexted :class:`ParSet` objects.

        Returns
        -------
        dict
            The contents in dictionary form.
        """
        return {k: v.to_dict() if isinstance(v, ParSet) else v for k, v in self.data.items()}


    @classmethod
    def from_dict(cls, cfg):
        """
        Instantiate a :class:`ParSet` from a dictionary.

        This simply constructs a :class:`ParSet` that behaves exactly
        like a dictionary. That is, no constraints are put on the
        types, options, etc., of the parameters.

        .. warning::

            Use of this method of instantiating a :class:`ParSet`
            **only** sets the parameter keys and values, not
            its other attributes (e.g., options, dtypes, etc). In
            essentially all use cases, the :class:`ParSet` should be
            used as a base class where the :func:`from_dict` method
            is overwritten such that these higher-level attributes of
            the derived class are maintained in the header I/O. See,
            e.g., :mod:`pypeit.par.pypeitpar`.

        Args:
            cfg (dict-like):
                Dictionary-like object used to set the parameters and
                parameter values for the :class:`ParSet`.
        """
        pars, values = map(lambda x : list(x), zip(*cfg.items()))
        return cls(pars, values=values)


def set_parameter_definition(dtype=None, default=None, options=None, descr=None):
    """
    Define a parameter for a :class:`ParSet`.

    This should be used by the ``parameters`` attribute of :class:`ParSet`
    subclasses to ensure each parameter has all its necessary components.

    Parameters
    ----------
    dtype : type, list, optional
        A single type or list of types that are allowed for the parameter.
        Parameters cannot be dictionaries; instead create a new :class:`ParSet`
        subclass for that parameter.  If a parameter is a :class:`ParSet`, it
        *cannot have any other type* and it must be a single :class:`ParSet`
        subclass.
    default : object, optional
        The default value for the parameter.
    options : object, list, optional
        A list of valid options for the parameter.
    descr : str, optional
        A description of the parameter.

    Returns
    -------
    dict
        A dictionary containing the parameter specifications.

    Raises
    ------
    ValueError
        Raised if the parameter definition does not adhere to the rules outlined
        in the ``dtype`` argument.
    """
    _dtype = None if dtype is None else np.atleast_1d(dtype).tolist()
    if _dtype is not None:
        # Parameter types are not allowed to be dictionaries
        if any(d is dict for d in _dtype):
            raise ValueError(
                'The dtype of a parameter cannot be dict! Use a ParSet subclass instead.'
            )
        # If a parameter can be a ParSet, it is not allowed to be any other type
        if len(_dtype) > 1 and any(issubclass(d, ParSet) for d in _dtype):
            raise ValueError(
                'If a parameter can be a ParSet, it is not allowed to be any other type and must '
                'be a specific ParSet subclass!'
            )

    return {
        'dtype': _dtype,
        'default': default,
        'options': None if options is None else np.atleast_1d(options).tolist(),
        'descr': descr
    }


# TODO: We might need to define the types that allowed for list parameters.
# They should be single element objects (ints, floats, strings), NOT more
# complex things like dicts or ParSets.  Nested ParSets are allowed.
class NewParSet:
    """
    """

    default_key = None
    """
    The default key to use for this parameter set in a configuration file.
    """

    default_comment = None
    """
    The default commment to use for this parameter set in a configuration file.
    """

    card_prefix = None
    """
    The header card prefix to use when writing this parameter set to a FITS
    header.
    """

    parameters = None
    """
    Definition of the parameters used by this parameter set.
    """

    def __init__(self, **kwargs):
        if self.parameters is None:
            raise AttributeError(
                f'CODING ERROR: The parameters attribute for {self.__class__.__name__} has not '
                'been defined!'
            )

        # The keys of self.parameters define the allowed keywords.
        allowed_keys = self.keys()
        badkeys = np.array([key for key in kwargs.keys() if key not in allowed_keys])
        if np.any(badkeys):
            raise TypeError(
                f'One or more unrecognized parameters for {self.__class__.__name__}: {badkeys}'
            )
        
        # The number of parameters is set by the parameters attribute
        self.npar = len(self.parameters)
        # Instantiate the data dictionary with the keys provided by the
        # parameters attribute
        self.data = dict.fromkeys(allowed_keys)
        # First set all the parameters to their default values; this ensures
        # that the defaults in the parameters attribute adhere to their
        # definition, just as any user-defined value should.
        for key in self.data.keys():
            try:
                self.__setitem__(key, self.parameters[key]['default'])
            except TypeError as e:
                raise TypeError(
                    f'CODING ERROR: The object type for the default value of {key} in '
                    f'{self.__class__.__name__} does not adhere to its definition!'
                ) from e
            except ValueError as e:
                raise ValueError(
                    f'CODING ERROR: The default value for {key} in {self.__class__.__name__} is '
                    'not one of the allowed options!'
                ) from e

        # Now set all the user-defined values.
        for key in kwargs.keys():
            self.__setitem__(key, kwargs[key])

        # Finally, we validate
        self.validate()

    def __getitem__(self, key):
        """
        Get a parameter value.

        Parameters
        ----------
        key : str
            Keyword of the parameter
        """
        return self.data[key]

    def __setitem__(self, key, value):
        """
        Set a parameter value.

        The parameter type and value must adhere to the :attr:`parameters`
        definition; however, a value can *always* be set to None.

        Parameters
        ----------
        key : str
            Key for new parameter
        value
            New parameter value to set

        Raises
        ------
        ValueError
            Raised if the parameter value is not among the allowed options.
        TypeError
            Raised if the parameter value is not an instance of an allowed data
            type.
        """
        # Check that the key is valid
        if key not in self.parameters.keys():
            raise KeyError(f'{key} is not a valid parameter for {self.__class__.__name__}!')

        # Set the value to None
        if value is None:
            self.data[key] = value
            return

        # Check that the value has an allowed data type
        if self.parameters[key]['dtype'] is not None \
                and not any([ isinstance(value, d) for d in self.parameters[key]['dtype']]):
            raise TypeError(
                f'Unable to set {key} in {self.__class__.__name__} to an object with type '
                f'{type(value)}.'
            )

        # Disallow elements of a list to be ParSets or dicts 
        if (
            list in self.parameters[key]['dtype']
            and isinstance(value, list)
            and any(isinstance(v, (NewParSet, dict)) for v in value)
        ):
            raise TypeError('Elements of a list should never be a dict or ParSet!')

        # Check that the value is one among a set of allowed options.
        if self.parameters[key]['options'] is not None:
            _value = [value] if not isinstance(value, list) else value
            for v in _value:
                if v not in self.parameters[key]['options']:
                    raise ValueError(
                        f'Unable to set {key} in {self.__class__.__name__} to {v}.  '
                        f'Allowed options are: {self.parameters[key]["options"]}.'
                    )

        # Otherwise, set the value
        self.data[key] = value

    def __len__(self):
        """Return the number of parameters."""
        return self.npar
        
    def __iter__(self):
        """Return an iterable to the parameter values."""
        return iter(self.data.values())

    def __repr__(self):
        """Return a string representation of the parameters."""
        return self._output_string() #header=self.cfg_section)

    def _output_string(self, header=None, value_only=False):
        """
        Constructs the short-format table strings for the
        :func:`__repr__` method.
        
        Parameters
        ----------
        header : :obj:`str`, optional
            String header to provide for the table.  This is typically the name
            of the configuration section.
        value_only : :obj:`bool`, optional
            By default, the table includes the parameter key, its current value,
            the default value, and its data type.  If `value_only=True`, only
            the parameter key and current value are returned.

        Returns
        -------
        str
            Single long string with the parameter table for the :func:`__repr__`
            method.
        """
        # Instantiate the table
        # TODO: Make this an astropy table?
        if value_only:
            data_table = np.empty((self.npar+1, 2), dtype=object)
            data_table[0] = ['Parameter', 'Value']
        else:
            data_table = np.empty((self.npar+1, 4), dtype=object)
            data_table[0] = ['Parameter', 'Value', 'Default', 'Type']

        additional_par_strings = []
        for i, key in enumerate(self.keys()):
            data_table[i+1,0] = key
            if isinstance(self.data[key], NewParSet):
                _header = key if header is None else '{0}:{1}'.format(header, key)
                additional_par_strings += [
                    self.data[key]._output_string(header=_header, value_only=value_only)
                ]
                data_table[i+1,1] = 'see below'
                if not value_only:
                    data_table[i+1,2] = 'see below'
            else:
                data_table[i+1,1] = NewParSet._data_string(self.data[key])
                if not value_only:
                    data_table[i+1,2] = NewParSet._data_string(self.parameters[key]['default'])
            if value_only:
                continue

            data_table[i+1,3] = (
                'Undefined' if self.parameters[key]['dtype'] is None
                else ', '.join([t.__name__ for t in self.parameters[key]['dtype']])
            )

        output = [NewParSet._data_table_string(data_table)]
        if header is not None:
            output = [header] + output
        if len(additional_par_strings) > 0:
            output += additional_par_strings
        return '\n'.join(output)

    @staticmethod
    def _data_table_string(data_table, delimeter='print'):
        """
        Provided the array of data, format it with equally spaced
        columns and add a header (first row) and contents delimeter.

        Parameters
        ----------
        data_table : :class:`numpy.ndarray`
            Array of string representations of the data to print.
        delimeter : :obj:`str`, optional
            The type of delimeter used to separate the header row from the rest
            of the table.  If 'print', the delimeter is a row of dashes below
            the header row.  If anything else, the delimeter is a row of equal
            signs above and below the header row; this is used to construct
            reStructuredText (rst) tables.

        Returns
        -------
        str
            Single long string with the data table.
        """
        nrows, ncols = data_table.shape
        col_width = [np.max([len(dij) for dij in dj]) for dj in data_table.T]
        row_string = ['']*(nrows+1) if delimeter == 'print' else ['']*(nrows+3)
        start = 2 if delimeter == 'print' else 3
        for i in range(start,nrows+start-1):
            row_string[i] = '  '.join([
                data_table[1+i-start,j].ljust(col_width[j]) for j in range(ncols)
            ])

        if delimeter == 'print':
            # Heading row
            row_string[0] = '  '.join([data_table[0,j].ljust(col_width[j]) for j in range(ncols)])
            # Delimiter
            row_string[1] = '-'*len(row_string[0])
            return '\n'.join(row_string)+'\n'

        # For an rst table
        row_string[0] = '  '.join(['='*col_width[j] for j in range(ncols)])
        row_string[1] = '  '.join([data_table[0,j].ljust(col_width[j]) for j in range(ncols)])
        row_string[2] = row_string[0]
        row_string[-1] = row_string[0]
        return '\n'.join(row_string)+'\n'

    @staticmethod
    def _data_string(data, use_repr=False, verbatim=False, check_dir=False):
        """
        Convert a single datum into a string
        
        Simply return strings, recursively convert the elements of any
        objects with a :attr:`__len__` attribute, and use the object's
        own :attr:`__repr__` attribute for all other objects.

        Parameters
        ----------
        data : object
            The object to stringify.
        use_repr : :obj:`bool`, optional
            Use the objects :attr:`__repr__` method; otherwise, use a direct
            string conversion.
        verbatim : :obj:`bool`, optional
            Use quotes around the provided string to indicate that the string
            should be representated in a verbatim (fixed width) font.
        check_dir : :obj:`bool`, optional
            If ``data`` is a string, check if it matches the current working
            directory and replace it with a generic string if it does.
        
        Returns
        -------
        str
            A string representation of the provided ``data``.
        """
        if isinstance(data, str):
            _data = '$PWD' if check_dir and data == str(Path().cwd()) else data
            if verbatim:
                return '..' if len(_data) == 0 else f'``{_data}``'
            return _data

        if isinstance(data, list):
            # When the list is empty, return an empty string, which config_lines
            # will append a "," to.  This allows ConfigObj to interpret it as an
            # empty list, instead of string, when re-reading the configuration
            # lines into a ParSet
            return (
                '' if len(data) == 0
                else ', '.join([
                    NewParSet._data_string(d, use_repr=use_repr, verbatim=verbatim) for d in data
                ])
            )

        return data.__repr__() if use_repr else str(data)

    def _wrap_print(self, head, output, tcols):
        """
        Wrap the contents of an output string for a fixed terminal
        width.  Used for the long-format :func:`info` method.

        Parameters
        ----------
        head : :obj:`str`
            The inline header for the output.  Can be an empty string, but
            cannot be None.
        output : :obj:`str`
            The main body of the text to write.
        tcols : :obj:`int`
            The allowed width for the output.
        """
        tail = ' '*len(head)
        if tcols is not None:
            lines = textwrap.wrap('{0}'.format(output), tcols-len(head))
            if len(lines) == 0:
                print('{0}None'.format(head))
            else:
                _head = [ head ] + [ tail ]*(len(lines)-1)
                print('\n'.join([ h+l for h,l in zip(_head, lines)]))
        else:
            print(head+'{0}'.format(output))

    def _types_list(self, key):
        """Return the string names for the specified data types."""
        return [t.__name__ for t in self.parameters[key]['dtype']]

    @staticmethod
    def _config_comment(comment, indent, full_width=72):
        """
        Create a list of lines for the description of a given parameter.

        Parameters
        ----------
        comment : str
            The description of the configuration parameter.
        indent : str
            The string used to indent the text.
        full_width : int, optional
            The full width allowed for each output string in the returned list.

        Returns
        -------
        list
            List of the strings to write to the output configuration file.
        """
        head = indent + '# '
        lines = textwrap.wrap(f'{comment}', full_width-len(head))
        return [head + l for l in lines]

    @classmethod
    def _rst_class_name(cls):
        return ':class:`~' +  cls.__module__ + '.' + cls.__name__ + '`'

    def keys(self):
        """Return the list of parameter set keys."""
        return list(self.parameters.keys())

    def validate(self):
        """
        Validate the parameter set.
        """
        pass

    def to_rst_table(self, parsets_listed=[]):
        """
        Construct a reStructuredText table describing the parameter set.

        This works recursively for nested :class:`ParSet` instances.
        
        Parameters
        ----------
        parsets_listed : :obj:`list`, optional
            For nested :class:`ParSet` instances, this is the list of
            :class:`ParSet` subclass names that already have already a table in
            the string list (so that they're not repeated).
        
        Returns
        -------
        list
            A list of lines that can be written to an ``*.rst`` file.
        """
        new_parsets = []
        data_table = np.empty((self.npar+1, 5), dtype=object)
        data_table[0,:] = ['Key', 'Type', 'Options', 'Default', 'Description']
        sorted_keys = np.sort(self.keys())
        for i, key in enumerate(sorted_keys):
            data_table[i+1,0] = NewParSet._data_string(key, use_repr=False, verbatim=True)
            if isinstance(self.data[key], NewParSet):
                if type(self.data[key]).__name__ not in parsets_listed:
                    new_parsets += [key]
                parsets_listed += [ type(self.data[key]).__name__ ]
                data_table[i+1,1] = type(self.data[key])._rst_class_name()
                data_table[i+1,3] = f'`{type(self.data[key]).__name__} Keywords`_'
            else: 
                data_table[i+1,1] = ', '.join([t.__name__ for t in self.dtype[key]])
                data_table[i+1,3] = (
                    '..' if self.parameters[key]['default'] is None
                    else NewParSet._data_string(
                        self.parameters[key]['default'], use_repr=False, verbatim=True,
                        check_dir=True
                    )
                )

            data_table[i+1,2] = (
                '..' if self.parameters[key]['options'] is None
                else NewParSet._data_string(
                    self.parameters[key]['options'], use_repr=False, verbatim=True
                )
            )
            data_table[i+1,4] = (
                '..' if self.parameters[key]['descr'] is None
                else NewParSet._data_string(self.parameters[key]['descr'])
            )

        output = [ f'.. _{self.__class__.__name__.lower()}:']
        output += [ '' ]
        output += [ f'{self.__class__.__name__} Keywords']
        output += [ '-'*len(output[2]) ]
        output += [ '' ]
        output += ['Class Instantiation: ' + self.__class__._rst_class_name()]
        output += ['']
        output += [NewParSet._data_table_string(data_table, delimeter='rst')]
        output += ['']
        for k in new_parsets:
            output += ['----']
            output += ['']
            output += self.data[k].to_rst_table(parsets_listed=parsets_listed)
        return output

    def info(self, basekey=None):
        """
        A long-form version of __repr__ that includes the parameter descriptions.
        """
        tcols = int(0.9*shutil.get_terminal_size(fallback=(80, 25)).columns)
        for key in self.parameters.keys():
            if isinstance(self.data[key], NewParSet):
                self.data[key].info(basekey=k)
                continue
            print(f'{key}' if basekey is None else f'{basekey}:{key}')
            self._wrap_print('        Value: ', self.data[key], tcols)
            self._wrap_print('      Default: ', self.parameters[key]['default'], tcols)
            self._wrap_print(
                '      Options: ',
                'None' if self.parameters[key]['options'] is None
                else ', '.join(self.parameters[key]['options']),
                tcols
            )
            self._wrap_print(
                '  Valid Types: ',
                'None' if self.parameters[key]['dtype'] is None
                else ', '.join(self._types_list(key)),
                tcols
            )
            self._wrap_print('  Description: ', self.parameters[key]['descr'], tcols)
            print(' ')

    def config_lines(self, section_name=None, section_comment=None, section_level=0,
                     exclude_defaults=False, include_descr=True):
        """
        Recursively generate the lines of a configuration file that can be used
        to instantiate this parameter set.

        Parameters
        ----------
        section_name : :obj:`str`, optional
            Name to give to the top-level of the configuration output.
        section_comment : :obj:`str`, optional
            Description to provide for the top-level configuration output.
        section_level : :obj:`int`, optional
            The level for the configuration output.  Sets the indentation level
            and the number of square brackets assigned to the section name.
        exclude_defaults : :obj:`bool`, optional
            Do not include any parameters that are identical to the defaults.
        include_descr : :obj:`bool`, optional
            Include the descriptions of each parameter as comments.

        Returns
        -------
        list
            The list of the lines to write to a configuration file.
        """
        # Get the list of parameters that are ParSets
        parset_keys = [key for key in self.keys() if isinstance(self.data[key], ParSet)]
        n_parsets = len(parset_keys)

        # Set the top-level comment and section name
        section_indent = ' '*4*section_level
        component_indent = section_indent + ' '*4
        lines = (
            [] if section_comment is None
            else self._config_comment(section_comment, section_indent)
        )
        lines += [section_indent + '['*(section_level+1) + section_name + ']'*(section_level+1)]

        min_lines = len(lines)

        # First add all the parameters that are not ParSets
        for key in self.keys():
            # Skip it if this element is a ParSet
            if n_parsets > 0 and key in parset_keys:
                continue

            # Add a comment with the parameter description
            if self.parameters[key]['descr'] is not None and include_descr:
                lines += self._config_comment(self.parameters[key]['descr'], component_indent)

            # Add the parameter and its value
            if not exclude_defaults or self.data[key] != self.parameters[key]['default']:
                argvalue = self._data_string(self.data[key])
                if isinstance(self.data[key], list):
                    argvalue += ','
                lines += [component_indent + key + ' = ' + argvalue]

        # Then add the items that are ParSets as subsections
        for key in parset_keys:
            section_comment = None
            if include_descr:
                section_comment = self.parameters[key]['descr']
            lines += self.data[key].config_lines(
                section_name=key, section_comment=section_comment, section_level=section_level+1,
                exclude_defaults=exclude_defaults, include_descr=include_descr
            )
        return lines if len(lines) > min_lines else []

    def to_config(self, cfg_file=None, section_name=None, section_comment=None, section_level=0,
                  exclude_defaults=False, include_descr=True):
        """
        Write the parameter set to a configuration file.

        Parameters
        ----------
        cfg_file : :obj:`str`, optional
            The name of the file to write.  If None, the function returns the
            list of strings that would have been written to the file.  These
            lines can be used to construct a `configobj`_ instance.
        section_name : :obj:`str`, optional
            The top-level keyword used for this section of the configuration
            file.  If None, use :attr:`default_key`.  If :attr:`default_key` is
            None and ``section_level`` is 0, this *must* be provided.
        section_comment : :obj:`str`, optional
            The top-level comment used for this section of the configuration
            file.  If None, use :attr:`default_comment`.
        section_level : :obj:`int`, optional
            The top level of this :class:`ParSet`.  Used for recursive output
            of nested :class:`ParSet` instances.
        exclude_defaults : :obj:`bool`, optional
            Do not include any parameters that are identical to the default
            values.
        include_descr : :obj:`bool`, optional
            Include the descriptions of each parameter as comments.

        Returns
        -------
        list
            The list of lines to write to a configuration file.  This is
            returned regardless of whether or not a file is written.

        Raises
        ------
        ValueError
            Raised if :attr:`default_key` is None, ``section_level`` is 0, and
            ``section_name`` is None.
        """
        _cfg_file = cfg_file if cfg_file is None else Path(cfg_file).absolute()
        if _cfg_file is not None and _cfg_file.is_file():
            log.warning("Selected configuration file already exists and will be overwritten!")

        config_output = []
        if all(isinstance(d, ParSet) or d is None for d in self.data.values()):
            # All the elements are ParSets themselves, so just iterate
            # through each one
            for key in self.keys():
                if self.data[key] is None:
                    continue
                comment = self.parameters[key]['descr'] if include_descr else None
                config_output += self.data[key].config_lines(
                    section_name=key, section_comment=comment, section_level=section_level,
                    exclude_defaults=exclude_defaults, include_descr=include_descr
                )
        else:
            # Cannot write the parameters as a configuration file
            # without a top-level configuration section
            if self.default_key is None and section_level == 0 and section_name is None:
                raise ValueError('A top-level section name must be provided.')

            _section_name = self.default_key if section_name is None else section_name
            _section_comment = self.default_comment if section_comment is None else section_comment
            config_output += self.config_lines(
                section_name=_section_name, section_comment=_section_comment,
                section_level=section_level, exclude_defaults=exclude_defaults,
                include_descr=include_descr
            )

        if _cfg_file is not None:
            # Write the file
            with open(_cfg_file, 'w') as f:
                f.write('\n'.join(config_output))

        return config_output

    def to_dict(self):
        """
        Return a dictionary with the contents of the parameter set.

        This method recursively handles nexted :class:`ParSet` subclasses.

        Returns
        -------
        dict
            The contents in dictionary form.
        """
        return {key: v.to_dict() if isinstance(v, ParSet) else v for key, v in self.data.items()}

    @classmethod
    def from_dict(cls, cfg):
        """
        Instantiate from a dictionary.

        This recursively handles elements of the dictionary that are expected to
        be (subclasses of) :class:`ParSet` objects by calling their
        :func:`from_dict` methods.

        For subclasses that do not have any nested :class:`ParSet` objects, note
        that ``cls(**cfg)`` and ``cls.from_dict(cfg)`` are identical.

        Parameters
        ----------
        cfg : dict
            Dictionary to use to instantiate the parameter set.
        """
        pars, values = map(lambda x : list(x), zip(*cfg.items()))
        for i, key in enumerate(pars):
            if key not in cls.parameters.keys():
                raise KeyError(f'{key} is not a valid {cls.__name__} parameter!')
            if isinstance(values[i], dict) and issubclass(cls.parameters[key]['dtype'][0], ParSet):
                values[i] = cls.parameters[key]['dtype'][0].from_dict(values[i])
        return cls(**dict(zip(pars, values)))
    
    @classmethod
    def class_header_card(cls):
        """
        Return the header card that specifies the class of the parameter set.

        Returns
        -------
        str
            The header card that specifies the class of the parameter set.

        Raises
        ------
        ValueError
            Raised if :attr:`card_prefix` is None for this subclass of
            :class:`ParSet`.
        """
        if cls.card_prefix is None:
            raise ValueError(
                f'CODING ERROR: The card_prefix attribute for {cls.__name__} has not been defined!'
            )
        return f'{cls.card_prefix.upper()}C'
    
    @classmethod
    def dict_header_card(cls):
        """
        Return the header card that specifies the dictionary version of the
        parameter set.

        Returns
        -------
        str
            The header card that specifies the dictionary version of the
            parameter set.

        Raises
        ------
        ValueError
            Raised if :attr:`card_prefix` is None for this subclass of
            :class:`ParSet`.
        """
        if cls.card_prefix is None:
            raise ValueError(
                f'CODING ERROR: The card_prefix attribute for {cls.__name__} has not been defined!'
            )
        return f'{cls.card_prefix.upper()}D'
    
    def to_header(self, hdr=None):
        """
        Include the parameter set in a FITS header.

        This adds two header cards:

            - ``self.class_header_card()``: This is the full name of the class
              for the :class:`ParSet` subclass, including the module.

            - ``self.dict_header_card()``: This is a string representation of
              the dictionary version of the parameter set, as given by
              ``str(self.to_dict())``.

        Parameters
        ----------
        hdr : :class:`astropy.io.fits.Header`, optional
            The header to which to add the parameter set.  If None, a new
            header is instantiated.

        Returns
        -------
        :class:`astropy.io.fits.Header`
            The header with the parameter set included.
        """
        if hdr is None:
            hdr = fits.Header()
        hdr[self.class_header_card()] = f'{self.__class__.__module__}.{self.__class__.__name__}'
        hdr[self.dict_header_card()] = str(self.to_dict())
        return hdr
    
    @classmethod
    def from_header(cls, hdr):
        """
        Instantiate a :class:`ParSet` from a FITS header.

        This looks for the header cards with the provided prefix and uses
        the information in those cards to instantiate the parameter set.

        Parameters
        ----------
        hdr : :class:`astropy.io.fits.Header`
            The header from which to read the parameter set.

        Raises
        ------
        ValueError
            Raised if :attr:`card_prefix` is None for this subclass of
            :class:`ParSet`, or if the class header card does not match the
            expected class.
        """
        if cls.card_prefix is None:
            raise ValueError(
                f'CODING ERROR: The card_prefix attribute for {cls.__name__} has not been defined!'
            )
        expected_cls = hdr.get(cls.class_header_card())
        if expected_cls is None:
            raise ValueError(
                'Header does not include card specifying the ParSet class: '
                f'{cls.class_header_card()}.'
            )
        if expected_cls != f'{cls.__module__}.{cls.__name__}':
            raise ValueError(
                f'Expected {cls.__module__}.{cls.__name__} in header card {cls.class_header_card()}, '
                f'but found {expected_cls}.'
            )
        contents = hdr.get(cls.dict_header_card())
        if contents is None:
            raise ValueError(
                'Header does not include card specifying the ParSet contents: '
                f'{cls.dict_header_card()}.'
            )
        cfg = util.ast_literal_eval(contents)
        if not isinstance(cfg, dict):
            raise ValueError(
                f'Unable to parse the contents of header card {cls.dict_header_card()} as a '
                f'dictionary: {contents}.'
            )
        return cls.from_dict(cfg)
