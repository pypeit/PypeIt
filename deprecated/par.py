
class TraceSlitsPar(ParSet):
    """
    The parameter set used to hold arguments for tracing the slit
    positions along the dispersion axis.
    
    For a table with the current keywords, defaults, and descriptions,
    see :ref:`pypeitpar`.
    """
    prefix = 'TSP'  # Prefix for writing parameters to a header is a class attribute
    def __init__(self, function=None, medrep=None, number=None, trim=None, maxgap=None,
                 maxshift=None, pad=None, sigdetect=None, min_slit_width=None, add_slits=None,
                 rm_slits=None, diffpolyorder=None, single=None, sobel_mode=None, pcaextrap=None,
                 smash_range=None, trace_npoly=None, mask_frac_thresh=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[1:]])      # "1:" to skip 'self'

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set

        defaults['function'] = 'legendre'
        options['function'] = TraceSlitsPar.valid_functions()
        dtypes['function'] = str
        descr['function'] = 'Function use to trace the slit center.  ' \
                            'Options are: {0}'.format(', '.join(options['function']))

        defaults['medrep'] = 0
        dtypes['medrep'] = int
        descr['medrep'] = 'Median-smoothing iterations to perform on sqrt(trace) image before ' \
                          'applying to Sobel filter, which detects slit/order edges.'

        # TODO: Never used?
        # Force number to be an integer
        if values['number'] == 'auto':
            values['number'] = -1
        defaults['number'] = -1
        dtypes['number'] = int
        descr['number'] = 'Manually set the number of slits to identify (>=1). \'auto\' or -1 ' \
                          'will automatically identify the number of slits.'

        # Force trim to be a tuple
        if pars['trim'] is not None and not isinstance(pars['trim'], tuple):
            try:
                pars['trim'] = tuple(pars['trim'])
            except:
                raise TypeError('Could not convert provided trim to a tuple.')
        defaults['trim'] = (0,0)
        dtypes['trim'] = tuple
        descr['trim'] = 'How much to trim off each edge of each slit.  Each number should be 0 ' \
                        'or positive'

        dtypes['maxgap'] = int
        descr['maxgap'] = 'Maximum number of pixels to allow for the gap between slits.  Use ' \
                          'None if the neighbouring slits are far apart or of similar ' \
                          'illumination.'

        defaults['maxshift'] = 0.15
        dtypes['maxshift'] = [int, float]
        descr['maxshift'] = 'Maximum shift in trace crude. Use a larger number for more curved ' \
                            'slits/orders.'

        defaults['pad'] = 0
        dtypes['pad'] = int
        descr['pad'] = 'Integer number of pixels to consider beyond the slit edges.'

        defaults['sigdetect'] = 20.0
        dtypes['sigdetect'] = [int, float]
        descr['sigdetect'] = 'Sigma detection threshold for edge detection'

        defaults['mask_frac_thresh'] = 0.6
        dtypes['mask_frac_thresh'] = float
        descr['mask_frac_thresh'] = 'Minimum fraction of the slit edge that was *not* masked ' \
                                    'to use in initial PCA.'

        defaults['smash_range'] = [0., 1.]
        dtypes['smash_range'] = list
        descr['smash_range'] = 'Range of the slit in the spectral direction (in fractional ' \
                               'units) to smash when searching for slit edges.  If the ' \
                               'spectrum covers only a portion of the image, use that range.'

        defaults['trace_npoly'] = 5
        dtypes['trace_npoly'] = int
        descr['trace_npoly'] = 'Order of legendre polynomial fits to slit/order boundary traces.'

        # TODO: slit *width* is a misnomer.  Should be slit *length*
        defaults['min_slit_width'] = 6.0  # arcseconds!
        dtypes['min_slit_width'] = float
        descr['min_slit_width'] = 'If a slit spans less than this number of arcseconds over ' \
                                   'the spatial direction of the detector, it will be ignored.' \
                                   '  Use this option to prevent the alignment (box) slits ' \
                                   'from multislit reductions, which typically cannot be ' \
                                   'reduced without a significant struggle.'

        defaults['diffpolyorder'] = 2
        dtypes['diffpolyorder'] = int
        descr['diffpolyorder'] = 'Order of the 2D function used to fit the 2d solution for the ' \
                                 'spatial size of all orders.'

        # TO BE DEPRECATED?
        defaults['single'] = []
        dtypes['single'] = list
        descr['single'] = 'Add a single, user-defined slit based on its location on each ' \
                          'detector.  Syntax is a list of values, 2 per detector, that define ' \
                          'the slit according to column values.  The second value (for the ' \
                          'right edge) must be greater than 0 to be applied.  LRISr example: ' \
                          'setting single = -1, -1, 7, 295 means the code will skip the ' \
                          'user-definition for the first detector but adds one for the second. ' \
                          ' None means no user-level slits defined.'

        dtypes['add_slits'] = [str, list]
        descr['add_slits'] = 'Add one or more user-defined slits.  The syntax to define a ' \
                             'slit to add is: \'det:spec:spat_left:spat_right\' where ' \
                             'det=detector, spec=spectral pixel, spat_left=spatial pixel of ' \
                             'left slit boundary, and spat_righ=spatial pixel of right slit ' \
                             'boundary.  For example, \'2:2000:2121:2322,3:2000:1201:1500\' ' \
                             'will add a slit to detector 2 passing through spec=2000 ' \
                             'extending spatially from 2121 to 2322 and another on detector 3 ' \
                             'at spec=2000 extending from 1201 to 1500.'

        dtypes['rm_slits'] = [str, list]
        descr['rm_slits'] = 'Remove one or more user-specified slits.  The syntax used to ' \
                            'define a slit to remove is: \'det:spec:spat\' where det=detector, ' \
                            'spec=spectral pixel, spat=spatial pixel.  For example, ' \
                            '\'2:2000:2121,3:2000:1500\' will remove the slit on detector 2 ' \
                            'that contains pixel (spat,spec)=(2000,2121) and on detector 3 ' \
                            'that contains pixel (2000,2121).'

        defaults['sobel_mode'] = 'nearest'
        options['sobel_mode'] = TraceSlitsPar.valid_sobel_modes()
        dtypes['sobel_mode'] = str
        descr['sobel_mode'] = 'Mode for Sobel filtering.  Default is \'nearest\' but the ' \
                              'developers find \'constant\' works best for DEIMOS.'

#        # DEPRECATED
#        defaults['pcatype'] = 'pixel'
#        options['pcatype'] = TraceSlitsPar.valid_pca_types()
#        dtypes['pcatype'] = str
#        descr['pcatype'] = 'Select to perform the PCA using the pixel position (pcatype=pixel) ' \
#                           'or by spectral order (pcatype=order).  Pixel positions can be used ' \
#                           'for multi-object spectroscopy where the gap between slits is ' \
#                           'irregular.  Order is used for echelle spectroscopy or for slits ' \
#                           'with separations that are a smooth function of the slit number.'
#
#        # DEPRECATED
#        defaults['pcapar'] = [3, 2, 1, 0]
#        dtypes['pcapar'] = list
#        descr['pcapar'] = 'Order of the polynomials to be used to fit the principle ' \
#                          'components.  The list length must be equal to or less than ' \
#                          'polyorder+1. TODO: Provide more explanation'

        defaults['pcaextrap'] = [0, 0]
        dtypes['pcaextrap'] = list
        descr['pcaextrap'] = 'The number of extra orders to predict in the negative (first ' \
                             'number) and positive (second number) direction.  Must be two ' \
                             'numbers in the list and they must be integers.'

        # Instantiate the parameter set
        super(TraceSlitsPar, self).__init__(list(pars.keys()),
                                            values=list(pars.values()),
                                            defaults=list(defaults.values()),
                                            options=list(options.values()),
                                            dtypes=list(dtypes.values()),
                                            descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = cfg.keys()
        parkeys = ['function', 'medrep', 'number', 'trim', 'maxgap', 'maxshift', 'pad',
                   'sigdetect', 'min_slit_width', 'diffpolyorder', 'sobel_mode', 
                   'pcaextrap', 'add_slits', 'rm_slits', 'smash_range', 'trace_npoly',
                   'mask_frac_thresh']
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    @staticmethod
    def valid_functions():
        """
        Return the list of valid functions to use for slit tracing.
        """
        return [ 'polynomial', 'legendre', 'chebyshev' ]

    @staticmethod
    def valid_sobel_modes():
        """Return the valid sobel modes."""
        return [ 'nearest', 'constant' ]

    @staticmethod
    def valid_pca_types():
        """
        Return the valid PCA types.
        """
        return ['pixel', 'order']

    def validate(self):
        if self.data['number'] == 0:
            raise ValueError('Number of slits must be -1 for automatic identification or '
                             'greater than 0')
        if len(self.data['pcaextrap']) != 2:
            raise ValueError('PCA extrapolation parameters must be a list with two values.')
        for e in self.data['pcaextrap']:
            if not isinstance(e, int):
                raise ValueError('PCA extrapolation values must be integers.')


