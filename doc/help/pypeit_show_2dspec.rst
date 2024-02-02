.. code-block:: console

    $ pypeit_show_2dspec -h
    usage: pypeit_show_2dspec [-h] [--list] [--det DET] [--spat_id SPAT_ID]
                              [--maskID MASKID] [--showmask [SHOWMASK ...]]
                              [--removetrace] [--embed] [--ignore_extract_mask]
                              [--channels CHANNELS] [--prefix PREFIX] [--no_clear]
                              [-v VERBOSITY] [--try_old]
                              file
    
    Display sky subtracted, spec2d image in a ginga viewer.
    
    positional arguments:
      file                  Path to a PypeIt spec2d file
    
    options:
      -h, --help            show this help message and exit
      --list                List the extensions only? (default: False)
      --det DET             Detector name or number. If a number, the name is
                            constructed assuming the reduction is for a single
                            detector. If a string, it must match the name of the
                            detector object (e.g., DET01 for a detector, MSC01 for a
                            mosaic). (default: 1)
      --spat_id SPAT_ID     Restrict plotting to this slit (PypeIt ID notation)
                            (default: None)
      --maskID MASKID       Restrict plotting to this maskID (default: None)
      --showmask [SHOWMASK ...]
                            Include a channel showing the mask. If no arguments are
                            provided, the mask bit values are provided directly. You
                            can also specify one or more mask flags used to
                            construct an image identifying which pixels are flagged
                            by any of these issues. E.g., to show pixels flagged by
                            the instrument specific bad-pixel mask or cosmic arrays,
                            use --showmask BPM CR . See
                            https://pypeit.readthedocs.io/en/release/out_masks.html
                            for the list of flags. (default: None)
      --removetrace         Do not overplot traces in the skysub, sky_resid, and
                            resid channels (default: False)
      --embed               Upon completion embed in ipython shell (default: False)
      --ignore_extract_mask
                            Ignore the extraction mask (default: False)
      --channels CHANNELS   Only show a subset of the channels (0-indexed), e.g. 1,3
                            (default: None)
      --prefix PREFIX       Channel name prefix [lets you display more than one set]
                            (default: )
      --no_clear            Do *not* clear all existing tabs (default: True)
      -v VERBOSITY, --verbosity VERBOSITY
                            Verbosity level between 0 [none] and 2 [all] (default:
                            1)
      --try_old             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    