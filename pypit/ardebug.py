def init():
    """
    Returns
    -------
    debug : dict
        default debug dict
    """
    debug = dict(develop=False,
                 arc=False,
                 obj_profile=False,
                 slit_profile=False,
                 sky_sub=False,
                 trace=False,
                 wave=False,
                 tilts=False,
                 flexure=False,
                 no_qa=False,
                 trace_obj=False,
                 )
    return debug
