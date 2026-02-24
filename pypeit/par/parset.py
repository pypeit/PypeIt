"""
Define a utility base class used to hold parameters.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from pathlib import Path
import shutil
import textwrap

from astropy.io import fits
from IPython import embed
import numpy as np

from pypeit import log
from pypeit.par import util


def tuple_force(par):
    """
    Cast object as tuple.

    Args:
        par (object):
            Object to cast as a tuple.  Can be None; if so, the returned value
            is also None (*not* an empty tuple).  If this is already a tuple,
            this is the returned value.  If the input is a list with one tuple,
            the returned value is just the single tuple in the list (i.e, this
            does not convert the result to a tuple of one tuple).

    Returns:
        :obj:`tuple`: Casted result.
    """
    # Already has correct type
    if par is None or isinstance(par, tuple):
        return par

    # If the value is a list of one tuple, return the tuple.
    # TODO: This is a hack, and we should probably revisit how this is done.
    # The issue is that pypeit.par.util.eval_tuple always returns a list of
    # tuples, something that's required for allowing lists of detector mosaics.
    # But elements of TelluricPar are forced to be tuples.  When constructing
    # the parameters to use in a given run, the sequence of merging the
    # defaults, configuration-specific, and user-provided parameters leads to
    # converting these TelluricPar parameters into multiply nested tuples.  This
    # hook avoids that.
    if isinstance(par, list) and len(par) == 1 and isinstance(par[0], tuple):
        return par[0]

    return tuple(par)


# TODO: May need to disallow parameters being either a tuple or list.  Because
# of issues with ConfigObj, these may need to be mutually exclusive.
def set_parameter_definition(dtype=None, default=None, options=None, descr=None):
    """
    Define a parameter for a :class:`~pypeit.par.parset.ParSet`.

    This should be used by the ``parameters`` attribute of
    :class:`~pypeit.par.parset.ParSet` subclasses to ensure each parameter
    has all its necessary components.

    Parameters
    ----------
    dtype : type, list, optional
        A single type or list of types that are allowed for the parameter.
        Parameters cannot be dictionaries; instead create a new
        :class:`~pypeit.par.parset.ParSet` subclass for that parameter.  If a
        parameter is a :class:`~pypeit.par.parset.ParSet`, it
        *cannot have any other type* and it must be a single
        :class:`~pypeit.par.parset.ParSet` subclass.
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
class ParSet:
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
            try:
                self.__setitem__(key, kwargs[key])
            except:
                embed(header='setting item in init')
                exit()

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
            and any(isinstance(v, (ParSet, dict)) for v in value)
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
            if isinstance(self.data[key], ParSet):
                _header = key if header is None else '{0}:{1}'.format(header, key)
                additional_par_strings += [
                    self.data[key]._output_string(header=_header, value_only=value_only)
                ]
                data_table[i+1,1] = 'see below'
                if not value_only:
                    data_table[i+1,2] = 'see below'
            else:
                data_table[i+1,1] = ParSet._data_string(self.data[key])
                if not value_only:
                    data_table[i+1,2] = ParSet._data_string(self.parameters[key]['default'])
            if value_only:
                continue

            data_table[i+1,3] = (
                'Undefined' if self.parameters[key]['dtype'] is None
                else ', '.join([t.__name__ for t in self.parameters[key]['dtype']])
            )

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
                    ParSet._data_string(d, use_repr=use_repr, verbatim=verbatim) for d in data
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

    def validate_keys(self, required=None, can_be_None=None):
        """
        Validate that a set of keys are present and that they have values that
        are not None.

        Parameters
        ----------
        required : :obj:`list`, optional
            A list of strings that provide the set of required keys
        can_be_None : :obj:`list`, optional
            A list of strings with the keywords that are allowed to be None

        Raises
        ------
        ValueError
            Raised if any of the required keywords are not present or if any of
            the keyword values are None if they are not allowed to be.
        """
        if required is None and can_be_None is None:
            # No validation rules, so implicitly valid
            return

        if required is not None:
            not_defined = np.array([key not in self.keys() for key in required])
            if np.any(not_defined):
                raise ValueError(
                    f'Required keys were not defined: {np.asarray(required)[not_defined].tolist()}'
                )

        if can_be_None is not None:
            should_not_be_None = np.array([
                self.data[key] is None and key not in can_be_None for key in self.keys()
            ])
            if np.any(should_not_be_None):
                raise ValueError(
                    'The following keys cannot be None: '
                    f'{np.asarray(self.keys())[should_not_be_None].tolist()}'
                )

    def to_rst_table(self, parsets_listed=[]):
        """
        Construct a reStructuredText table describing the parameter set.

        This works recursively for nested :class:`~pypeit.par.parset.ParSet` instances.
        
        Parameters
        ----------
        parsets_listed : :obj:`list`, optional
            For nested :class:`~pypeit.par.parset.ParSet` instances, this is the list of
            :class:`~pypeit.par.parset.ParSet` subclass names that already have already a table in
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
            data_table[i+1,0] = ParSet._data_string(key, use_repr=False, verbatim=True)
            if isinstance(self.data[key], ParSet):
                if type(self.data[key]).__name__ not in parsets_listed:
                    new_parsets += [key]
                parsets_listed += [ type(self.data[key]).__name__ ]
                data_table[i+1,1] = type(self.data[key])._rst_class_name()
                data_table[i+1,3] = f'`{type(self.data[key]).__name__} Keywords`_'
            else: 
                data_table[i+1,1] = ', '.join([t.__name__ for t in self.parameters[key]['dtype']])
                data_table[i+1,3] = (
                    '..' if self.parameters[key]['default'] is None
                    else ParSet._data_string(
                        self.parameters[key]['default'], use_repr=False, verbatim=True,
                        check_dir=True
                    )
                )

            data_table[i+1,2] = (
                '..' if self.parameters[key]['options'] is None
                else ParSet._data_string(
                    self.parameters[key]['options'], use_repr=False, verbatim=True
                )
            )
            data_table[i+1,4] = (
                '..' if self.parameters[key]['descr'] is None
                else ParSet._data_string(self.parameters[key]['descr'])
            )

        output = [ f'.. _{self.__class__.__name__.lower()}:']
        output += [ '' ]
        output += [ f'{self.__class__.__name__} Keywords']
        output += [ '-'*len(output[2]) ]
        output += [ '' ]
        output += ['Class Instantiation: ' + self.__class__._rst_class_name()]
        output += ['']
        output += [ParSet._data_table_string(data_table, delimeter='rst')]
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
            if isinstance(self.data[key], ParSet):
                self.data[key].info(basekey=key)
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
            The top level of this :class:`~pypeit.par.parset.ParSet`.  Used for recursive output
            of nested :class:`~pypeit.par.parset.ParSet` instances.
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

        This method recursively handles nexted :class:`~pypeit.par.parset.ParSet` subclasses.

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
        be (subclasses of) :class:`~pypeit.par.parset.ParSet` objects by
        calling their :func:`from_dict` methods.

        For subclasses that do not have any nested
        :class:`~pypeit.par.parset.ParSet` objects, note that ``cls(**cfg)``
        and ``cls.from_dict(cfg)`` are identical.

        Parameters
        ----------
        cfg : dict
            Dictionary to use to instantiate the parameter set.
        """
        pars, values = map(lambda x : list(x), zip(*cfg.items()))
        for i, key in enumerate(pars):
            if key not in cls.parameters.keys():
                embed()
                exit()
                raise KeyError(f'{key} is not a valid {cls.__name__} parameter!')
            if isinstance(values[i], dict) and issubclass(cls.parameters[key]['dtype'][0], ParSet):
                values[i] = cls.parameters[key]['dtype'][0].from_dict(values[i])
            if len(cls.parameters[key]['dtype']) == 1 and cls.parameters[key]['dtype'][0] is tuple:
                values[i] = tuple_force(values[i])
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
            :class:`~pypeit.par.parset.ParSet`.
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
            :class:`~pypeit.par.parset.ParSet`.
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
              for the :class:`~pypeit.par.parset.ParSet` subclass, including the module.

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
        Instantiate a :class:`~pypeit.par.parset.ParSet` from a FITS header.

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
            :class:`~pypeit.par.parset.ParSet`, or if the class header card does not match the
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
