"""
Define a utility base class used to hold parameters.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import os
import textwrap

from IPython import embed

import numpy

from astropy.io import fits

from pypeit import msgs
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
                msgs.warn(
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
            data_table (:obj:`numpy.ndarray`):
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
    def _data_string(data, use_repr=True, verbatim=False):
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
        
        Returns:
            str: A string representation of the provided ``data``.
        """
        if isinstance(data, str):
            if verbatim:
                return '..' if len(data) == 0 else '``' + data + '``'
            return data
        if isinstance(data, list):
            # TODO: When the list is empty, should the return include the
            # brackets?
            return '[]' if len(data) == 0 \
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
        the provided ParSet or dict (par).

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
                lines += [ component_indent + k + ' = ' + ParSet._data_string(par[k]) ]

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
                The name of the file to write/append to.  If None
                (default), the function will just return the list of
                strings that would have been written to the file.  These
                lines can be used to construct a :class:`ConfigObj`
                instance.
            section_name (:obj:`str`, optional):
                The top-level name for the config section.  This must be
                provided if :attr:`cfg_section` is None or any of the
                parameters are not also ParSet instances themselves.
            section_comment (:obj:`str`, optional):
                The top-level comment for the config section based on
                this ParSet.
            section_level (:obj:`int`, optional):
                The top level of this ParSet.  Used for recursive output
                of nested ParSets.
            append (:obj:`bool`, optional):
                Append this configuration output of this ParSet to the
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
                Raised if there are types other than ParSet in the
                parameter list, :attr:`cfg_section` is None, and no
                section_name argument was provided.
        """
        if cfg_file is not None and os.path.isfile(cfg_file) and not append and not quiet:
            msgs.warn("Selected configuration file already exists and will be overwritten!")

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
                                                             verbatim=True)

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

        Any element that has a value of None or is a ParSet itself is
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
                    msgs.warn(
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
        Instantiate the ParSet using data parsed from a fits header.

        This is a simple wrapper for
        :func:`ParSet.parse_par_from_hdr` and
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


class ParDatabase(object):
    """
    NOTE: This isn't used in pypeit yet...

    Class used as a list of ParSets in a glorified structured numpy
    array.

    Very similar to yanny when converted to a numpy array.

    Can be initialized using a list of ParSet objects, or an SDSS
    parameter file.

    .. todo::

        - Check that the data types are the same for all ParSet objects
          in the list
        - Better handle the NaN values when converting None to a float
          type
        - Add from_par_file classmethod?
    
    """
    def __init__(self, inp):
        """
        nsets - number of parameter sets
        npar - number of parameters in each parameter set
        data - parameter set data
        options - allowed options for values
        dtype - allowed datatypes
        can_call - parameter is a callable function
        """
        _inp = [inp] if isinstance(inp, ParSet) else inp
        if not isinstance(_inp, list):
            raise TypeError('Input must be a list.')
        for i in _inp:
            if not isinstance(i, ParSet):
                raise TypeError('Input must be a list of ParSet objects.')
        self.npar = _inp[0].npar
        self.nsets = len(_inp)
        keys = _inp[0].keys()
        for i in range(1,self.nsets):
            if _inp[i].npar != self.npar:
                raise ValueError('Not all ParSet objects have the same number of parameters.')
            if _inp[i].keys() != keys:
                raise ValueError('Not all ParSet objects have the same keys.')
            # Other checks?

        record_dtypes = self._set_dtypes(_inp, 0)

        data = []
        for i in range(self.nsets):
            data += [ tuple([_inp[i][k] for k in keys]) ]

        # WARNING: None values are converted to nan if data type is
        # float
        self.data = numpy.array(data, dtype=record_dtypes ).view(numpy.recarray)
        self.options = inp[0].options.copy()
        self.dtype = inp[0].dtype.copy()
        self.can_call = inp[0].can_call.copy()
   

    def __getitem__(self, key):
        """
        Return the value of the designated key.

        Args:
            key (str) : Key for new parameter
        """
        return self.data[key]


    @staticmethod
    def _set_dtypes(inp, i):
        keys = inp[i].keys()
        dtypes = []
        for k in keys:
            if inp[i].dtype[k] is None:
                dtypes += [(k,object)]
                continue
            # inp.dtype is always a list
            if any([t in inp[i].dtype[k] for t in [int , float]]) \
                and any([t in inp[i].dtype[k] for t in [list, numpy.ndarray]]):
                msgs.warn(
                    "Parameter set has elements that can be either individual "
                    f"ints/floats or lists/arrays.  Database column {k} will have type "
                    "'object'."
                )
                dtypes += [(k,object)]
            elif len(list({int, float} - set(inp[i].dtype[k]))) == 0:
                dtypes += [(k,float)]
            elif len(list({list, numpy.ndarray} - set(inp[i].dtype[k]))) == 0 \
                    or inp[i].dtype[k] == numpy.ndarray:
                _inp = numpy.asarray(inp[i][k])
                dtypes += [(k,_inp.dtype,_inp.shape)]
            elif isinstance(inp[i][k], str):
                if any([ _inp[k] is None for _inp in inp]):
                    dtypes += [(k, object)]
                else:
                    dtypes += [(k,'<U{0:d}'.format(max([ len(_inp[k]) for _inp in inp])))]
            else:
                dtypes += [(k,type(inp[i][k]))]
        return dtypes


    def append(self, pdb):
        if not isinstance(pdb, ParDatabase):
            raise TypeError('Can only append ParDatabase object.')

        try:
            self.data = numpy.append(self.data, pdb.data)
        except TypeError as e:
            raise TypeError('Could not append data:: {0}'.format(e))
            



