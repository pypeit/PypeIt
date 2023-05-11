.. code-block:: console

    $ pypeit_chk_noise_2dspec -h
    usage: pypeit_chk_noise_2dspec [-h] [--det DET] [--z [Z ...]]
                                   [--maskdef_id MASKDEF_ID] [--pypeit_id PYPEIT_ID]
                                   [--pad PAD] [--aspect_ratio ASPECT_RATIO]
                                   [--wavemin WAVEMIN] [--wavemax WAVEMAX]
                                   [--mode MODE] [--list]
                                   [files ...]
    
    Examine the noise in a PypeIt slit/order
    
    positional arguments:
      files                 PypeIt spec2d file(s) (default: None)
    
    optional arguments:
      -h, --help            show this help message and exit
      --det DET             Detector name or number. If a number, the name is
                            constructed assuming the reduction is for a single
                            detector. If a string, it must match the name of the
                            detector object (e.g., DET01 for a detector, MSC01 for a
                            mosaic). (default: 1)
      --z [Z ...]           Object redshift (default: None)
      --maskdef_id MASKDEF_ID
                            MASKDEF_ID of the slit that you want to plot. If
                            maskdef_id is not provided, nor a pypeit_id, all the 2D
                            spectra in the file(s) will be plotted. (default: None)
      --pypeit_id PYPEIT_ID
                            PypeIt ID of the slit that you want to plot. If
                            pypeit_id is not provided, nor a maskdef_id, all the 2D
                            spectra in the file(s) will be plotted. (default: None)
      --pad PAD             Padding for the selected slit. Negative value will trim.
                            (default: -5)
      --aspect_ratio ASPECT_RATIO
                            Aspect ratio when plotting the spec2d (default: 3)
      --wavemin WAVEMIN     Wavelength min. This is for selecting a region of the
                            spectrum to analyze. (default: None)
      --wavemax WAVEMAX     Wavelength max.This is for selecting a region of the
                            spectrum to analyze. (default: None)
      --mode MODE           Options are: plot, save, printDo you want to save to
                            disk or open a plot in a mpl window. If you choose save,
                            a folder called spec2d*_noisecheck will be created and
                            all the relevant plot will be placed there. If you
                            choose print, check noise value are printed in the
                            terminal (default: plot)
      --list                List the extensions only? (default: False)
    