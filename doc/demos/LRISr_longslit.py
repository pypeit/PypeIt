{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# Fun with the PypitSetup Class [v1.2]"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 1,
   "metadata": {
    "collapsed": true
   },
   "outputs": [],
   "source": [
    "# import\n",
    "from importlib import reload\n",
    "import os\n",
    "import glob\n",
    "import numpy as np\n",
    "\n",
    "from astropy.table import Table\n",
    "\n",
    "from pypit import pypitsetup"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## To play along, you need the Development suite and the $PYPIT_DEV environmental variable pointed at it"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "'/data/Projects/Python/PYPIT-development-suite/'"
      ]
     },
     "execution_count": 2,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "os.getenv('PYPIT_DEV')"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Kuldge for settings"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "metadata": {},
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marparse.py 75 load_file()\u001b[0m - Loading default settings\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marparse.py 87 load_file()\u001b[0m - Loading base settings from settings.baseargflag\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marparse.py 1393 run_ncpus()\u001b[0m - Setting 19 CPUs\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marparse.py 75 load_file()\u001b[0m - Loading default settings\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marparse.py 75 load_file()\u001b[0m - Loading default settings\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marparse.py 87 load_file()\u001b[0m - Loading base settings from settings.basespect\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marparse.py 75 load_file()\u001b[0m - Loading default settings\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marparse.py 75 load_file()\u001b[0m - Loading default settings\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marparse.py 87 load_file()\u001b[0m - Loading base settings from settings.basespect\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marparse.py 75 load_file()\u001b[0m - Loading default settings\n"
     ]
    }
   ],
   "source": [
    "from pypit import arparse as settings\n",
    "settings.dummy_settings()\n",
    "settings.argflag['run']['spectrograph'] = 'shane_kast_blue'\n",
    "settings.argflag['reduce']['masters']['setup'] = 'C_01_aa'\n",
    "#\n",
    "# Load default spectrograph settings\n",
    "spect = settings.get_spect_class(('ARMLSD', 'shane_kast_blue', 'pypit'))# '.'.join(redname.split('.')[:-1])))\n",
    "lines = spect.load_file(base=True)  # Base spectrograph settings\n",
    "spect.set_paramlist(lines)\n",
    "lines = spect.load_file()  # Instrument specific\n",
    "spect.set_paramlist(lines)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Build the fitstbl"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "27"
      ]
     },
     "execution_count": 4,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "kast_blue_files = glob.glob(os.getenv('PYPIT_DEV')+'RAW_DATA/Shane_Kast_blue/600_4310_d55/b*')\n",
    "len(kast_blue_files)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "metadata": {
    "collapsed": true
   },
   "outputs": [],
   "source": [
    "# Init\n",
    "setupc = pypitsetup.PypitSetup(settings.argflag, settings.spect)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "metadata": {},
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b27.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b27.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b3.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b3.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b1.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b1.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b10.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b10.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b13.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b13.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b11.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b11.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b22.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b22.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b17.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b17.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b7.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b7.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b24.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b24.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b16.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b16.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b12.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b12.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b29.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b29.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b2.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b2.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b9.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b9.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b14.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b14.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b4.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b4.fits.gz\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b8.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b8.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b5.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b5.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b28.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b28.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b18.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b18.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b15.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b15.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b23.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b23.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b21.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b21.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b19.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b19.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b6.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b6.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b20.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b20.fits.gz\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 197 load_headers()\u001b[0m - Checking spectrograph settings for required header information\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 206 load_headers()\u001b[0m - Headers loaded for 27 files successfully\n"
     ]
    }
   ],
   "source": [
    "fitstbl = setupc.build_fitstbl(kast_blue_files)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<i>Table length=5</i>\n",
       "<table id=\"table139792909909464\" class=\"table-striped table-bordered table-condensed\">\n",
       "<thead><tr><th>directory</th><th>filename</th><th>utc</th><th>target</th><th>idname</th><th>time</th><th>date</th><th>equinox</th><th>ra</th><th>dec</th><th>airmass</th><th>naxis0</th><th>naxis1</th><th>binning</th><th>exptime</th><th>filter1</th><th>filter2</th><th>hatch</th><th>shutopen</th><th>shutclose</th><th>decker</th><th>lamps</th><th>slitwid</th><th>slitlen</th><th>detrot</th><th>dichroic</th><th>dispname</th><th>dispangle</th><th>lampname01</th><th>lampstat01</th><th>lampname02</th><th>lampstat02</th><th>lampname03</th><th>lampstat03</th><th>lampname04</th><th>lampstat04</th><th>lampname05</th><th>lampstat05</th><th>lampname06</th><th>lampstat06</th><th>lampname07</th><th>lampstat07</th><th>lampname08</th><th>lampstat08</th><th>lampname09</th><th>lampstat09</th><th>lampname10</th><th>lampstat10</th><th>lampname11</th><th>lampstat11</th><th>lampname12</th><th>lampstat12</th><th>lampname13</th><th>lampstat13</th><th>lampname14</th><th>lampstat14</th><th>lampname15</th><th>lampstat15</th><th>lampname16</th><th>lampstat16</th></tr></thead>\n",
       "<thead><tr><th>str84</th><th>str11</th><th>str4</th><th>str10</th><th>str6</th><th>float64</th><th>str22</th><th>str4</th><th>str10</th><th>str10</th><th>float64</th><th>int64</th><th>int64</th><th>str4</th><th>int64</th><th>str4</th><th>str4</th><th>str4</th><th>str4</th><th>str4</th><th>str10</th><th>str4</th><th>str4</th><th>str4</th><th>str4</th><th>str3</th><th>str8</th><th>str4</th><th>str4</th><th>str3</th><th>str3</th><th>str3</th><th>str6</th><th>str3</th><th>str6</th><th>str3</th><th>str8</th><th>str3</th><th>str6</th><th>str3</th><th>str6</th><th>str3</th><th>str6</th><th>str3</th><th>str8</th><th>str3</th><th>str8</th><th>str3</th><th>str2</th><th>str3</th><th>str5</th><th>str3</th><th>str4</th><th>str3</th><th>str6</th><th>str3</th><th>str4</th><th>str3</th><th>str5</th><th>str3</th></tr></thead>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b1.fits.gz</td><td>None</td><td>Arcs</td><td>OBJECT</td><td>397801.59944444447</td><td>2015-05-20T01:35:58.10</td><td>None</td><td>09:21:46.0</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>30</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>0.5 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>off</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>on</td><td>Hg-Cd</td><td>on</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b2.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.7936111111</td><td>2015-05-20T01:47:37.23</td><td>None</td><td>09:33:26.9</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>30</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b3.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.8230555556</td><td>2015-05-20T01:49:23.29</td><td>None</td><td>09:35:28.3</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b4.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.83194444445</td><td>2015-05-20T01:49:55.93</td><td>None</td><td>09:36:01.1</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b5.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.8411111111</td><td>2015-05-20T01:50:28.57</td><td>None</td><td>09:36:33.7</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "</table>"
      ],
      "text/plain": [
       "<Table length=5>\n",
       "                                     directory                                       ...\n",
       "                                       str84                                         ...\n",
       "------------------------------------------------------------------------------------ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ..."
      ]
     },
     "execution_count": 7,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "fitstbl[0:5]"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Write"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 8,
   "metadata": {
    "collapsed": true
   },
   "outputs": [],
   "source": [
    "setupc.write_fitstbl('fitstbl.fits')"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Load"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "metadata": {},
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34mpypitsetup.py 233 load_fitstbl()\u001b[0m - Loaded fitstbl from fitstbl.fits\n"
     ]
    },
    {
     "data": {
      "text/html": [
       "<i>Table length=27</i>\n",
       "<table id=\"table139793334734976\" class=\"table-striped table-bordered table-condensed\">\n",
       "<thead><tr><th>directory</th><th>filename</th><th>utc</th><th>target</th><th>idname</th><th>time</th><th>date</th><th>equinox</th><th>ra</th><th>dec</th><th>airmass</th><th>naxis0</th><th>naxis1</th><th>binning</th><th>exptime</th><th>filter1</th><th>filter2</th><th>hatch</th><th>shutopen</th><th>shutclose</th><th>decker</th><th>lamps</th><th>slitwid</th><th>slitlen</th><th>detrot</th><th>dichroic</th><th>dispname</th><th>dispangle</th><th>lampname01</th><th>lampstat01</th><th>lampname02</th><th>lampstat02</th><th>lampname03</th><th>lampstat03</th><th>lampname04</th><th>lampstat04</th><th>lampname05</th><th>lampstat05</th><th>lampname06</th><th>lampstat06</th><th>lampname07</th><th>lampstat07</th><th>lampname08</th><th>lampstat08</th><th>lampname09</th><th>lampstat09</th><th>lampname10</th><th>lampstat10</th><th>lampname11</th><th>lampstat11</th><th>lampname12</th><th>lampstat12</th><th>lampname13</th><th>lampstat13</th><th>lampname14</th><th>lampstat14</th><th>lampname15</th><th>lampstat15</th><th>lampname16</th><th>lampstat16</th></tr></thead>\n",
       "<thead><tr><th>bytes84</th><th>bytes11</th><th>bytes4</th><th>bytes10</th><th>bytes6</th><th>float64</th><th>bytes22</th><th>bytes4</th><th>bytes10</th><th>bytes10</th><th>float64</th><th>int64</th><th>int64</th><th>bytes4</th><th>int64</th><th>bytes4</th><th>bytes4</th><th>bytes4</th><th>bytes4</th><th>bytes4</th><th>bytes10</th><th>bytes4</th><th>bytes4</th><th>bytes4</th><th>bytes4</th><th>bytes3</th><th>bytes8</th><th>bytes4</th><th>bytes4</th><th>bytes3</th><th>bytes3</th><th>bytes3</th><th>bytes6</th><th>bytes3</th><th>bytes6</th><th>bytes3</th><th>bytes8</th><th>bytes3</th><th>bytes6</th><th>bytes3</th><th>bytes6</th><th>bytes3</th><th>bytes6</th><th>bytes3</th><th>bytes8</th><th>bytes3</th><th>bytes8</th><th>bytes3</th><th>bytes2</th><th>bytes3</th><th>bytes5</th><th>bytes3</th><th>bytes4</th><th>bytes3</th><th>bytes6</th><th>bytes3</th><th>bytes4</th><th>bytes3</th><th>bytes5</th><th>bytes3</th></tr></thead>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b1.fits.gz</td><td>None</td><td>Arcs</td><td>OBJECT</td><td>397801.59944444447</td><td>2015-05-20T01:35:58.10</td><td>None</td><td>09:21:46.0</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>30</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>0.5 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>off</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>on</td><td>Hg-Cd</td><td>on</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b2.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.7936111111</td><td>2015-05-20T01:47:37.23</td><td>None</td><td>09:33:26.9</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>30</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b3.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.8230555556</td><td>2015-05-20T01:49:23.29</td><td>None</td><td>09:35:28.3</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b4.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.83194444445</td><td>2015-05-20T01:49:55.93</td><td>None</td><td>09:36:01.1</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b5.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.8411111111</td><td>2015-05-20T01:50:28.57</td><td>None</td><td>09:36:33.7</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b6.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.85027777776</td><td>2015-05-20T01:51:01.10</td><td>None</td><td>09:37:06.5</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b7.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.8591666667</td><td>2015-05-20T01:51:33.56</td><td>None</td><td>09:37:39.1</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b8.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.86833333335</td><td>2015-05-20T01:52:06.03</td><td>None</td><td>09:38:11.7</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b9.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.8772222222</td><td>2015-05-20T01:52:38.46</td><td>None</td><td>09:38:44.3</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b10.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.8861111111</td><td>2015-05-20T01:53:10.95</td><td>None</td><td>09:39:16.9</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b18.fits.gz</td><td>None</td><td>Bias</td><td>DARK</td><td>397803.72</td><td>2015-05-20T03:43:12.82</td><td>None</td><td>11:30:32.9</td><td>37:06:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>0</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>off</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b19.fits.gz</td><td>None</td><td>Bias</td><td>DARK</td><td>397803.725</td><td>2015-05-20T03:43:30.27</td><td>None</td><td>11:30:50.8</td><td>37:11:10.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>0</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>off</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b20.fits.gz</td><td>None</td><td>Bias</td><td>DARK</td><td>397803.7297222222</td><td>2015-05-20T03:43:47.75</td><td>None</td><td>11:31:08.2</td><td>37:15:16.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>0</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>off</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b21.fits.gz</td><td>None</td><td>Bias</td><td>DARK</td><td>397803.7347222222</td><td>2015-05-20T03:44:05.14</td><td>None</td><td>11:31:25.7</td><td>37:19:25.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>0</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>off</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b22.fits.gz</td><td>None</td><td>Bias</td><td>DARK</td><td>397803.7394444444</td><td>2015-05-20T03:44:22.50</td><td>None</td><td>11:31:43.2</td><td>37:23:33.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>0</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>off</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b23.fits.gz</td><td>None</td><td>Bias</td><td>DARK</td><td>397803.7444444444</td><td>2015-05-20T03:44:40.04</td><td>None</td><td>11:32:00.4</td><td>37:25:21.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>0</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>off</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b24.fits.gz</td><td>None</td><td>Feige 66</td><td>OBJECT</td><td>397804.21277777775</td><td>2015-05-20T04:12:46.96</td><td>None</td><td>12:37:54.8</td><td>24:59:47.0</td><td>1.039999961853</td><td>350</td><td>2112</td><td>None</td><td>30</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>off</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b27.fits.gz</td><td>None</td><td>J1217p3905</td><td>OBJECT</td><td>397804.95916666667</td><td>2015-05-20T04:57:33.56</td><td>None</td><td>12:17:36.7</td><td>39:00:40.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>1200</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>off</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b28.fits.gz</td><td>None</td><td>J1217p3905</td><td>OBJECT</td><td>397805.3002777778</td><td>2015-05-20T05:18:01.47</td><td>None</td><td>12:17:37.0</td><td>39:00:40.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>1200</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>off</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b29.fits.gz</td><td>None</td><td>J1217p3905</td><td>OBJECT</td><td>397805.6383333333</td><td>2015-05-20T05:38:18.97</td><td>None</td><td>12:17:37.2</td><td>39:00:41.0</td><td>1.009999990463</td><td>350</td><td>2112</td><td>None</td><td>1200</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>off</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "</table>"
      ],
      "text/plain": [
       "<Table length=27>\n",
       "                                     directory                                       ...\n",
       "                                      bytes84                                        ...\n",
       "------------------------------------------------------------------------------------ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "                                                                                 ... ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ..."
      ]
     },
     "execution_count": 9,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "setupc.load_fitstbl('fitstbl.fits')"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Image type\n",
    "    Classifies the images\n",
    "    Adds image type columns to the fitstbl"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 10,
   "metadata": {},
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 85 type_data()\u001b[0m - Typing files\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 144 type_data()\u001b[0m - Making forced file identification changes\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marsort.py 145 type_data()\u001b[0m - Note that the image will have *only* the specified type\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 177 type_data()\u001b[0m - Typing completed!\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34mpypitsetup.py 211 type_data()\u001b[0m - Adding file type information to the fitstbl\n"
     ]
    }
   ],
   "source": [
    "filetypes = setupc.type_data()"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Show"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 11,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<i>Table length=27</i>\n",
       "<table id=\"table139792907041368\" class=\"table-striped table-bordered table-condensed\">\n",
       "<thead><tr><th>directory</th><th>filename</th><th>utc</th><th>target</th><th>idname</th><th>time</th><th>date</th><th>equinox</th><th>ra</th><th>dec</th><th>airmass</th><th>naxis0</th><th>naxis1</th><th>binning</th><th>exptime</th><th>filter1</th><th>filter2</th><th>hatch</th><th>shutopen</th><th>shutclose</th><th>decker</th><th>lamps</th><th>slitwid</th><th>slitlen</th><th>detrot</th><th>dichroic</th><th>dispname</th><th>dispangle</th><th>lampname01</th><th>lampstat01</th><th>lampname02</th><th>lampstat02</th><th>lampname03</th><th>lampstat03</th><th>lampname04</th><th>lampstat04</th><th>lampname05</th><th>lampstat05</th><th>lampname06</th><th>lampstat06</th><th>lampname07</th><th>lampstat07</th><th>lampname08</th><th>lampstat08</th><th>lampname09</th><th>lampstat09</th><th>lampname10</th><th>lampstat10</th><th>lampname11</th><th>lampstat11</th><th>lampname12</th><th>lampstat12</th><th>lampname13</th><th>lampstat13</th><th>lampname14</th><th>lampstat14</th><th>lampname15</th><th>lampstat15</th><th>lampname16</th><th>lampstat16</th><th>arc</th><th>bias</th><th>dark</th><th>pinhole</th><th>pixelflat</th><th>science</th><th>standard</th><th>trace</th><th>unknown</th></tr></thead>\n",
       "<thead><tr><th>bytes84</th><th>bytes11</th><th>bytes4</th><th>bytes10</th><th>bytes6</th><th>float64</th><th>bytes22</th><th>bytes4</th><th>bytes10</th><th>bytes10</th><th>float64</th><th>int64</th><th>int64</th><th>bytes4</th><th>int64</th><th>bytes4</th><th>bytes4</th><th>bytes4</th><th>bytes4</th><th>bytes4</th><th>bytes10</th><th>bytes4</th><th>bytes4</th><th>bytes4</th><th>bytes4</th><th>bytes3</th><th>bytes8</th><th>bytes4</th><th>bytes4</th><th>bytes3</th><th>bytes3</th><th>bytes3</th><th>bytes6</th><th>bytes3</th><th>bytes6</th><th>bytes3</th><th>bytes8</th><th>bytes3</th><th>bytes6</th><th>bytes3</th><th>bytes6</th><th>bytes3</th><th>bytes6</th><th>bytes3</th><th>bytes8</th><th>bytes3</th><th>bytes8</th><th>bytes3</th><th>bytes2</th><th>bytes3</th><th>bytes5</th><th>bytes3</th><th>bytes4</th><th>bytes3</th><th>bytes6</th><th>bytes3</th><th>bytes4</th><th>bytes3</th><th>bytes5</th><th>bytes3</th><th>bool</th><th>bool</th><th>bool</th><th>bool</th><th>bool</th><th>bool</th><th>bool</th><th>bool</th><th>bool</th></tr></thead>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b1.fits.gz</td><td>None</td><td>Arcs</td><td>OBJECT</td><td>397801.59944444447</td><td>2015-05-20T01:35:58.10</td><td>None</td><td>09:21:46.0</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>30</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>0.5 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>off</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>on</td><td>Hg-Cd</td><td>on</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td><td>True</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b2.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.7936111111</td><td>2015-05-20T01:47:37.23</td><td>None</td><td>09:33:26.9</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>30</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td><td>False</td><td>False</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>True</td><td>False</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b3.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.8230555556</td><td>2015-05-20T01:49:23.29</td><td>None</td><td>09:35:28.3</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td><td>False</td><td>False</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>True</td><td>False</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b4.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.83194444445</td><td>2015-05-20T01:49:55.93</td><td>None</td><td>09:36:01.1</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td><td>False</td><td>False</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>True</td><td>False</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b5.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.8411111111</td><td>2015-05-20T01:50:28.57</td><td>None</td><td>09:36:33.7</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td><td>False</td><td>False</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>True</td><td>False</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b6.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.85027777776</td><td>2015-05-20T01:51:01.10</td><td>None</td><td>09:37:06.5</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td><td>False</td><td>False</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>True</td><td>False</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b7.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.8591666667</td><td>2015-05-20T01:51:33.56</td><td>None</td><td>09:37:39.1</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td><td>False</td><td>False</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>True</td><td>False</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b8.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.86833333335</td><td>2015-05-20T01:52:06.03</td><td>None</td><td>09:38:11.7</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td><td>False</td><td>False</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>True</td><td>False</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b9.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.8772222222</td><td>2015-05-20T01:52:38.46</td><td>None</td><td>09:38:44.3</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td><td>False</td><td>False</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>True</td><td>False</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b10.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.8861111111</td><td>2015-05-20T01:53:10.95</td><td>None</td><td>09:39:16.9</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td><td>False</td><td>False</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>True</td><td>False</td></tr>\n",
       "<tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b18.fits.gz</td><td>None</td><td>Bias</td><td>DARK</td><td>397803.72</td><td>2015-05-20T03:43:12.82</td><td>None</td><td>11:30:32.9</td><td>37:06:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>0</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>off</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td><td>False</td><td>True</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b19.fits.gz</td><td>None</td><td>Bias</td><td>DARK</td><td>397803.725</td><td>2015-05-20T03:43:30.27</td><td>None</td><td>11:30:50.8</td><td>37:11:10.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>0</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>off</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td><td>False</td><td>True</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b20.fits.gz</td><td>None</td><td>Bias</td><td>DARK</td><td>397803.7297222222</td><td>2015-05-20T03:43:47.75</td><td>None</td><td>11:31:08.2</td><td>37:15:16.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>0</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>off</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td><td>False</td><td>True</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b21.fits.gz</td><td>None</td><td>Bias</td><td>DARK</td><td>397803.7347222222</td><td>2015-05-20T03:44:05.14</td><td>None</td><td>11:31:25.7</td><td>37:19:25.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>0</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>off</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td><td>False</td><td>True</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b22.fits.gz</td><td>None</td><td>Bias</td><td>DARK</td><td>397803.7394444444</td><td>2015-05-20T03:44:22.50</td><td>None</td><td>11:31:43.2</td><td>37:23:33.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>0</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>off</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td><td>False</td><td>True</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b23.fits.gz</td><td>None</td><td>Bias</td><td>DARK</td><td>397803.7444444444</td><td>2015-05-20T03:44:40.04</td><td>None</td><td>11:32:00.4</td><td>37:25:21.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>0</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>off</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td><td>False</td><td>True</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b24.fits.gz</td><td>None</td><td>Feige 66</td><td>OBJECT</td><td>397804.21277777775</td><td>2015-05-20T04:12:46.96</td><td>None</td><td>12:37:54.8</td><td>24:59:47.0</td><td>1.039999961853</td><td>350</td><td>2112</td><td>None</td><td>30</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>off</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b27.fits.gz</td><td>None</td><td>J1217p3905</td><td>OBJECT</td><td>397804.95916666667</td><td>2015-05-20T04:57:33.56</td><td>None</td><td>12:17:36.7</td><td>39:00:40.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>1200</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>off</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>False</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b28.fits.gz</td><td>None</td><td>J1217p3905</td><td>OBJECT</td><td>397805.3002777778</td><td>2015-05-20T05:18:01.47</td><td>None</td><td>12:17:37.0</td><td>39:00:40.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>1200</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>off</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>False</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b29.fits.gz</td><td>None</td><td>J1217p3905</td><td>OBJECT</td><td>397805.6383333333</td><td>2015-05-20T05:38:18.97</td><td>None</td><td>12:17:37.2</td><td>39:00:41.0</td><td>1.009999990463</td><td>350</td><td>2112</td><td>None</td><td>1200</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>off</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>False</td></tr>\n",
       "</table>"
      ],
      "text/plain": [
       "<Table length=27>\n",
       "                                     directory                                       ...\n",
       "                                      bytes84                                        ...\n",
       "------------------------------------------------------------------------------------ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "                                                                                 ... ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ..."
      ]
     },
     "execution_count": 11,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "setupc.fitstbl"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Do some logic (for fun)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 12,
   "metadata": {
    "collapsed": true
   },
   "outputs": [],
   "source": [
    "both = filetypes['pixelflat'] & filetypes['trace']"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 13,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<i>Table length=12</i>\n",
       "<table id=\"table139792907422800\" class=\"table-striped table-bordered table-condensed\">\n",
       "<thead><tr><th>directory</th><th>filename</th><th>utc</th><th>target</th><th>idname</th><th>time</th><th>date</th><th>equinox</th><th>ra</th><th>dec</th><th>airmass</th><th>naxis0</th><th>naxis1</th><th>binning</th><th>exptime</th><th>filter1</th><th>filter2</th><th>hatch</th><th>shutopen</th><th>shutclose</th><th>decker</th><th>lamps</th><th>slitwid</th><th>slitlen</th><th>detrot</th><th>dichroic</th><th>dispname</th><th>dispangle</th><th>lampname01</th><th>lampstat01</th><th>lampname02</th><th>lampstat02</th><th>lampname03</th><th>lampstat03</th><th>lampname04</th><th>lampstat04</th><th>lampname05</th><th>lampstat05</th><th>lampname06</th><th>lampstat06</th><th>lampname07</th><th>lampstat07</th><th>lampname08</th><th>lampstat08</th><th>lampname09</th><th>lampstat09</th><th>lampname10</th><th>lampstat10</th><th>lampname11</th><th>lampstat11</th><th>lampname12</th><th>lampstat12</th><th>lampname13</th><th>lampstat13</th><th>lampname14</th><th>lampstat14</th><th>lampname15</th><th>lampstat15</th><th>lampname16</th><th>lampstat16</th></tr></thead>\n",
       "<thead><tr><th>str84</th><th>str11</th><th>str4</th><th>str10</th><th>str6</th><th>float64</th><th>str22</th><th>str4</th><th>str10</th><th>str10</th><th>float64</th><th>int64</th><th>int64</th><th>str4</th><th>int64</th><th>str4</th><th>str4</th><th>str4</th><th>str4</th><th>str4</th><th>str10</th><th>str4</th><th>str4</th><th>str4</th><th>str4</th><th>str3</th><th>str8</th><th>str4</th><th>str4</th><th>str3</th><th>str3</th><th>str3</th><th>str6</th><th>str3</th><th>str6</th><th>str3</th><th>str8</th><th>str3</th><th>str6</th><th>str3</th><th>str6</th><th>str3</th><th>str6</th><th>str3</th><th>str8</th><th>str3</th><th>str8</th><th>str3</th><th>str2</th><th>str3</th><th>str5</th><th>str3</th><th>str4</th><th>str3</th><th>str6</th><th>str3</th><th>str4</th><th>str3</th><th>str5</th><th>str3</th></tr></thead>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b2.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.7936111111</td><td>2015-05-20T01:47:37.23</td><td>None</td><td>09:33:26.9</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>30</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b3.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.8230555556</td><td>2015-05-20T01:49:23.29</td><td>None</td><td>09:35:28.3</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b4.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.83194444445</td><td>2015-05-20T01:49:55.93</td><td>None</td><td>09:36:01.1</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b5.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.8411111111</td><td>2015-05-20T01:50:28.57</td><td>None</td><td>09:36:33.7</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b6.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.85027777776</td><td>2015-05-20T01:51:01.10</td><td>None</td><td>09:37:06.5</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b7.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.8591666667</td><td>2015-05-20T01:51:33.56</td><td>None</td><td>09:37:39.1</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b8.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.86833333335</td><td>2015-05-20T01:52:06.03</td><td>None</td><td>09:38:11.7</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b9.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.8772222222</td><td>2015-05-20T01:52:38.46</td><td>None</td><td>09:38:44.3</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b10.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.8861111111</td><td>2015-05-20T01:53:10.95</td><td>None</td><td>09:39:16.9</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b11.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.8952777778</td><td>2015-05-20T01:53:43.42</td><td>None</td><td>09:39:49.2</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b12.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.9041666667</td><td>2015-05-20T01:54:15.96</td><td>None</td><td>09:40:21.8</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "<tr><td>/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/</td><td>b13.fits.gz</td><td>None</td><td>Dome Flat</td><td>OBJECT</td><td>397801.91333333333</td><td>2015-05-20T01:54:48.52</td><td>None</td><td>09:40:54.7</td><td>37:25:56.0</td><td>1.0</td><td>350</td><td>2112</td><td>None</td><td>15</td><td>None</td><td>None</td><td>None</td><td>None</td><td>None</td><td>2.0 arcsec</td><td>None</td><td>None</td><td>None</td><td>None</td><td>d55</td><td>600/4310</td><td>None</td><td>Blue</td><td>off</td><td>Red</td><td>off</td><td>Spare3</td><td>off</td><td>Spare4</td><td>off</td><td>Sup_Blue</td><td>on</td><td>Spare1</td><td>off</td><td>Spare2</td><td>off</td><td>Spare3</td><td>off</td><td>Spare_Ar</td><td>off</td><td>Dim_Neon</td><td>off</td><td>He</td><td>off</td><td>Hg-Cd</td><td>off</td><td>Hg-A</td><td>off</td><td>Spare9</td><td>off</td><td>Neon</td><td>off</td><td>Laser</td><td>off</td></tr>\n",
       "</table>"
      ],
      "text/plain": [
       "<Table length=12>\n",
       "                                     directory                                       ...\n",
       "                                       str84                                         ...\n",
       "------------------------------------------------------------------------------------ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ...\n",
       "/data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/ ..."
      ]
     },
     "execution_count": 13,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "fitstbl[both]"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Match to science"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 14,
   "metadata": {},
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 850 match_to_science()\u001b[0m - Matching calibrations to Science frames\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 859 match_to_science()\u001b[0m - =================================================\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 861 match_to_science()\u001b[0m - Matching calibrations to J1217p3905: b27.fits.gz\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 1 arc frame for J1217p3905 (1 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 10 bias frame for J1217p3905 (5 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 874 match_to_science()\u001b[0m -   Dark frames not required.  Not matching..\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 888 match_to_science()\u001b[0m -    No pinhole frames are required.  Not matching..\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 12 pixelflat frame for J1217p3905 (5 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 1 standard frame for J1217p3905 (1 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 12 trace frame for J1217p3905 (3 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 859 match_to_science()\u001b[0m - =================================================\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 861 match_to_science()\u001b[0m - Matching calibrations to J1217p3905: b28.fits.gz\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 1 arc frame for J1217p3905 (1 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 10 bias frame for J1217p3905 (5 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 874 match_to_science()\u001b[0m -   Dark frames not required.  Not matching..\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 888 match_to_science()\u001b[0m -    No pinhole frames are required.  Not matching..\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 12 pixelflat frame for J1217p3905 (5 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 1 standard frame for J1217p3905 (1 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 12 trace frame for J1217p3905 (3 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 859 match_to_science()\u001b[0m - =================================================\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 861 match_to_science()\u001b[0m - Matching calibrations to J1217p3905: b29.fits.gz\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 1 arc frame for J1217p3905 (1 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 10 bias frame for J1217p3905 (5 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 874 match_to_science()\u001b[0m -   Dark frames not required.  Not matching..\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 888 match_to_science()\u001b[0m -    No pinhole frames are required.  Not matching..\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 12 pixelflat frame for J1217p3905 (5 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 1 standard frame for J1217p3905 (1 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 12 trace frame for J1217p3905 (3 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 933 match_to_science()\u001b[0m - Science frames successfully matched to calibration frames\n"
     ]
    }
   ],
   "source": [
    "fitstbl = setupc.match_to_science()"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Run it all -- as if in PYPIT run mode"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 15,
   "metadata": {
    "collapsed": true
   },
   "outputs": [],
   "source": [
    "setupc = pypitsetup.PypitSetup(settings.argflag, settings.spect)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 16,
   "metadata": {},
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b27.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b27.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b3.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b3.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b1.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b1.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b10.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b10.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b13.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b13.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b11.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b11.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b22.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b22.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b17.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b17.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b7.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b7.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b24.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b24.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b16.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b16.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b12.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b12.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b29.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b29.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b2.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b2.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b9.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b9.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b14.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b14.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b4.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b4.fits.gz\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b8.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b8.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b5.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b5.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b28.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b28.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b18.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b18.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b15.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b15.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b23.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b23.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b21.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b21.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b19.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b19.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b6.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b6.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b20.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b20.fits.gz\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 197 load_headers()\u001b[0m - Checking spectrograph settings for required header information\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 206 load_headers()\u001b[0m - Headers loaded for 27 files successfully\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 85 type_data()\u001b[0m - Typing files\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 144 type_data()\u001b[0m - Making forced file identification changes\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marsort.py 145 type_data()\u001b[0m - Note that the image will have *only* the specified type\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 177 type_data()\u001b[0m - Typing completed!\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34mpypitsetup.py 211 type_data()\u001b[0m - Adding file type information to the fitstbl\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 850 match_to_science()\u001b[0m - Matching calibrations to Science frames\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 859 match_to_science()\u001b[0m - =================================================\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 861 match_to_science()\u001b[0m - Matching calibrations to J1217p3905: b27.fits.gz\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 1 arc frame for J1217p3905 (1 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 10 bias frame for J1217p3905 (5 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 874 match_to_science()\u001b[0m -   Dark frames not required.  Not matching..\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 888 match_to_science()\u001b[0m -    No pinhole frames are required.  Not matching..\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 12 pixelflat frame for J1217p3905 (5 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 1 standard frame for J1217p3905 (1 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 12 trace frame for J1217p3905 (3 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 859 match_to_science()\u001b[0m - =================================================\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 861 match_to_science()\u001b[0m - Matching calibrations to J1217p3905: b28.fits.gz\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 1 arc frame for J1217p3905 (1 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 10 bias frame for J1217p3905 (5 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 874 match_to_science()\u001b[0m -   Dark frames not required.  Not matching..\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 888 match_to_science()\u001b[0m -    No pinhole frames are required.  Not matching..\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 12 pixelflat frame for J1217p3905 (5 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 1 standard frame for J1217p3905 (1 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 12 trace frame for J1217p3905 (3 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 859 match_to_science()\u001b[0m - =================================================\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 861 match_to_science()\u001b[0m - Matching calibrations to J1217p3905: b29.fits.gz\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 1 arc frame for J1217p3905 (1 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 10 bias frame for J1217p3905 (5 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 874 match_to_science()\u001b[0m -   Dark frames not required.  Not matching..\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 888 match_to_science()\u001b[0m -    No pinhole frames are required.  Not matching..\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 12 pixelflat frame for J1217p3905 (5 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 1 standard frame for J1217p3905 (1 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 12 trace frame for J1217p3905 (3 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 933 match_to_science()\u001b[0m - Science frames successfully matched to calibration frames\n"
     ]
    }
   ],
   "source": [
    "code, fitstbl, setup_dict = setupc.run(kast_blue_files)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Now as calcheck"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 17,
   "metadata": {},
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b27.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b27.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b3.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b3.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b1.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b1.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b10.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b10.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b13.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b13.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b11.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b11.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b22.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b22.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b17.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b17.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b7.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b7.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b24.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b24.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b16.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b16.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b12.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b12.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b29.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b29.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b2.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b2.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b9.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b9.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b14.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b14.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b4.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b4.fits.gz\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b8.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b8.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b5.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b5.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b28.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b28.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b18.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b18.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b15.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b15.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b23.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b23.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b21.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b21.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b19.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b19.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b6.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b6.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b20.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b20.fits.gz\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 197 load_headers()\u001b[0m - Checking spectrograph settings for required header information\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 206 load_headers()\u001b[0m - Headers loaded for 27 files successfully\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 85 type_data()\u001b[0m - Typing files\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 144 type_data()\u001b[0m - Making forced file identification changes\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marsort.py 145 type_data()\u001b[0m - Note that the image will have *only* the specified type\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 177 type_data()\u001b[0m - Typing completed!\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34mpypitsetup.py 211 type_data()\u001b[0m - Adding file type information to the fitstbl\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 850 match_to_science()\u001b[0m - Matching calibrations to Science frames\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 859 match_to_science()\u001b[0m - =================================================\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 861 match_to_science()\u001b[0m - Matching calibrations to J1217p3905: b27.fits.gz\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 1 arc frame for J1217p3905 (1 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 10 bias frame for J1217p3905 (5 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 874 match_to_science()\u001b[0m -   Dark frames not required.  Not matching..\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 888 match_to_science()\u001b[0m -    No pinhole frames are required.  Not matching..\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 12 pixelflat frame for J1217p3905 (5 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 1 standard frame for J1217p3905 (1 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 12 trace frame for J1217p3905 (3 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 859 match_to_science()\u001b[0m - =================================================\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 861 match_to_science()\u001b[0m - Matching calibrations to J1217p3905: b28.fits.gz\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 1 arc frame for J1217p3905 (1 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 10 bias frame for J1217p3905 (5 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 874 match_to_science()\u001b[0m -   Dark frames not required.  Not matching..\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 888 match_to_science()\u001b[0m -    No pinhole frames are required.  Not matching..\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 12 pixelflat frame for J1217p3905 (5 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 1 standard frame for J1217p3905 (1 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 12 trace frame for J1217p3905 (3 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 859 match_to_science()\u001b[0m - =================================================\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 861 match_to_science()\u001b[0m - Matching calibrations to J1217p3905: b29.fits.gz\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 1 arc frame for J1217p3905 (1 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 10 bias frame for J1217p3905 (5 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 874 match_to_science()\u001b[0m -   Dark frames not required.  Not matching..\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 888 match_to_science()\u001b[0m -    No pinhole frames are required.  Not matching..\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 12 pixelflat frame for J1217p3905 (5 required)\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 1 standard frame for J1217p3905 (1 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 12 trace frame for J1217p3905 (3 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 933 match_to_science()\u001b[0m - Science frames successfully matched to calibration frames\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34mpypitsetup.py 309 run()\u001b[0m - Inspect the .calib file: tmp.calib\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34mpypitsetup.py 310 run()\u001b[0m - *********************************************************\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34mpypitsetup.py 311 run()\u001b[0m - Calibration check complete and successful!\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34mpypitsetup.py 312 run()\u001b[0m - Set 'run calcheck False' to continue with data reduction\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34mpypitsetup.py 313 run()\u001b[0m - *********************************************************\n"
     ]
    }
   ],
   "source": [
    "settings.argflag['run']['calcheck'] = True\n",
    "setupc = pypitsetup.PypitSetup(settings.argflag, settings.spect)\n",
    "code, fitstbl, setup_dict = setupc.run(kast_blue_files)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 18,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "('calcheck', None)"
      ]
     },
     "execution_count": 18,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "code, setup_dict"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Now as setup"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 19,
   "metadata": {},
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b27.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b27.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b3.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b3.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b1.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b1.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b10.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b10.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b13.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b13.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b11.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b11.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b22.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b22.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b17.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b17.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b7.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b7.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b24.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b24.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b16.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b16.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b12.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b12.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b29.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b29.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b2.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b2.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b9.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b9.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b14.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b14.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b4.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b4.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b8.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b8.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b5.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b5.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b28.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b28.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b18.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b18.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b15.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b15.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b23.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b23.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b21.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b21.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b19.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b19.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b6.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b6.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 136 load_headers()\u001b[0m - UTC is not listed as a header keyword in file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b20.fits.gz\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marload.py 159 load_headers()\u001b[0m - BINNING keyword not in header. Setting to None\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 194 load_headers()\u001b[0m - Successfully loaded headers for file:\n",
      "             /data/Projects/Python/PYPIT-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55/b20.fits.gz\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 197 load_headers()\u001b[0m - Checking spectrograph settings for required header information\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marload.py 206 load_headers()\u001b[0m - Headers loaded for 27 files successfully\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 85 type_data()\u001b[0m - Typing files\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 144 type_data()\u001b[0m - Making forced file identification changes\n",
      "\u001b[1;31m[WARNING] ::\u001b[0m \u001b[1;34marsort.py 145 type_data()\u001b[0m - Note that the image will have *only* the specified type\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 177 type_data()\u001b[0m - Typing completed!\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34mpypitsetup.py 211 type_data()\u001b[0m - Adding file type information to the fitstbl\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 850 match_to_science()\u001b[0m - Matching calibrations to Science frames\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 859 match_to_science()\u001b[0m - =================================================\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 861 match_to_science()\u001b[0m - Matching calibrations to J1217p3905: b27.fits.gz\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 1 arc frame for J1217p3905 (1 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 10 bias frame for J1217p3905 (0 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 874 match_to_science()\u001b[0m -   Dark frames not required.  Not matching..\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 895 match_to_science()\u001b[0m - No matching criteria for pinhole frames with this instrument\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 0 pinhole frame for J1217p3905 (0 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 12 pixelflat frame for J1217p3905 (0 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 1 standard frame for J1217p3905 (0 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 12 trace frame for J1217p3905 (0 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 859 match_to_science()\u001b[0m - =================================================\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 861 match_to_science()\u001b[0m - Matching calibrations to J1217p3905: b28.fits.gz\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 1 arc frame for J1217p3905 (1 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 10 bias frame for J1217p3905 (0 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 874 match_to_science()\u001b[0m -   Dark frames not required.  Not matching..\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 895 match_to_science()\u001b[0m - No matching criteria for pinhole frames with this instrument\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 0 pinhole frame for J1217p3905 (0 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 12 pixelflat frame for J1217p3905 (0 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 1 standard frame for J1217p3905 (0 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 12 trace frame for J1217p3905 (0 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 859 match_to_science()\u001b[0m - =================================================\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 861 match_to_science()\u001b[0m - Matching calibrations to J1217p3905: b29.fits.gz\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 1 arc frame for J1217p3905 (1 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 10 bias frame for J1217p3905 (0 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 874 match_to_science()\u001b[0m -   Dark frames not required.  Not matching..\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 895 match_to_science()\u001b[0m - No matching criteria for pinhole frames with this instrument\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 0 pinhole frame for J1217p3905 (0 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 12 pixelflat frame for J1217p3905 (0 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 1 standard frame for J1217p3905 (0 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 911 match_to_science()\u001b[0m -   Found 12 trace frame for J1217p3905 (0 required)\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34marsort.py 933 match_to_science()\u001b[0m - Science frames successfully matched to calibration frames\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34mpypitsetup.py 109 build_group_dict()\u001b[0m - Wrote group dict to tmp.sorted\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34mpypitsetup.py 323 run()\u001b[0m - Setup is complete.\n",
      "\u001b[1;32m[INFO]    ::\u001b[0m \u001b[1;34mpypitsetup.py 324 run()\u001b[0m - Inspect the .setups file\n"
     ]
    }
   ],
   "source": [
    "settings.argflag['run']['calcheck'] = False\n",
    "settings.argflag['run']['setup'] = True\n",
    "#\n",
    "setupc = pypitsetup.PypitSetup(settings.argflag, settings.spect)\n",
    "code, fitstbl, setup_dict = setupc.run(kast_blue_files)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 20,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "'setup'"
      ]
     },
     "execution_count": 20,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "code"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Check the sci_idx column"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 21,
   "metadata": {
    "collapsed": true
   },
   "outputs": [],
   "source": [
    "cols = ['filename', 'target']+setupc.ftypes+['failures', 'sci_ID']"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 22,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<i>Table length=27</i>\n",
       "<table id=\"table139792908426768\" class=\"table-striped table-bordered table-condensed\">\n",
       "<thead><tr><th>filename</th><th>target</th><th>arc</th><th>bias</th><th>dark</th><th>pinhole</th><th>pixelflat</th><th>science</th><th>standard</th><th>trace</th><th>unknown</th><th>failures</th><th>sci_ID</th></tr></thead>\n",
       "<thead><tr><th>str11</th><th>str10</th><th>bool</th><th>bool</th><th>bool</th><th>bool</th><th>bool</th><th>bool</th><th>bool</th><th>bool</th><th>bool</th><th>bool</th><th>int64</th></tr></thead>\n",
       "<tr><td>b1.fits.gz</td><td>Arcs</td><td>True</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>7</td></tr>\n",
       "<tr><td>b2.fits.gz</td><td>Dome Flat</td><td>False</td><td>False</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>7</td></tr>\n",
       "<tr><td>b3.fits.gz</td><td>Dome Flat</td><td>False</td><td>False</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>7</td></tr>\n",
       "<tr><td>b4.fits.gz</td><td>Dome Flat</td><td>False</td><td>False</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>7</td></tr>\n",
       "<tr><td>b5.fits.gz</td><td>Dome Flat</td><td>False</td><td>False</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>7</td></tr>\n",
       "<tr><td>b6.fits.gz</td><td>Dome Flat</td><td>False</td><td>False</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>7</td></tr>\n",
       "<tr><td>b7.fits.gz</td><td>Dome Flat</td><td>False</td><td>False</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>7</td></tr>\n",
       "<tr><td>b8.fits.gz</td><td>Dome Flat</td><td>False</td><td>False</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>7</td></tr>\n",
       "<tr><td>b9.fits.gz</td><td>Dome Flat</td><td>False</td><td>False</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>7</td></tr>\n",
       "<tr><td>b10.fits.gz</td><td>Dome Flat</td><td>False</td><td>False</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>7</td></tr>\n",
       "<tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>\n",
       "<tr><td>b18.fits.gz</td><td>Bias</td><td>False</td><td>True</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>7</td></tr>\n",
       "<tr><td>b19.fits.gz</td><td>Bias</td><td>False</td><td>True</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>7</td></tr>\n",
       "<tr><td>b20.fits.gz</td><td>Bias</td><td>False</td><td>True</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>7</td></tr>\n",
       "<tr><td>b21.fits.gz</td><td>Bias</td><td>False</td><td>True</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>7</td></tr>\n",
       "<tr><td>b22.fits.gz</td><td>Bias</td><td>False</td><td>True</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>7</td></tr>\n",
       "<tr><td>b23.fits.gz</td><td>Bias</td><td>False</td><td>True</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>7</td></tr>\n",
       "<tr><td>b24.fits.gz</td><td>Feige 66</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>False</td><td>7</td></tr>\n",
       "<tr><td>b27.fits.gz</td><td>J1217p3905</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>False</td><td>False</td><td>1</td></tr>\n",
       "<tr><td>b28.fits.gz</td><td>J1217p3905</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>False</td><td>False</td><td>2</td></tr>\n",
       "<tr><td>b29.fits.gz</td><td>J1217p3905</td><td>False</td><td>False</td><td>False</td><td>False</td><td>False</td><td>True</td><td>False</td><td>False</td><td>False</td><td>False</td><td>4</td></tr>\n",
       "</table>"
      ],
      "text/plain": [
       "<Table length=27>\n",
       "  filename    target    arc   bias  dark ... trace unknown failures sci_ID\n",
       "   str11      str10     bool  bool  bool ...  bool   bool    bool   int64 \n",
       "----------- ---------- ----- ----- ----- ... ----- ------- -------- ------\n",
       " b1.fits.gz       Arcs  True False False ... False   False    False      7\n",
       " b2.fits.gz  Dome Flat False False False ...  True   False    False      7\n",
       " b3.fits.gz  Dome Flat False False False ...  True   False    False      7\n",
       " b4.fits.gz  Dome Flat False False False ...  True   False    False      7\n",
       " b5.fits.gz  Dome Flat False False False ...  True   False    False      7\n",
       " b6.fits.gz  Dome Flat False False False ...  True   False    False      7\n",
       " b7.fits.gz  Dome Flat False False False ...  True   False    False      7\n",
       " b8.fits.gz  Dome Flat False False False ...  True   False    False      7\n",
       " b9.fits.gz  Dome Flat False False False ...  True   False    False      7\n",
       "b10.fits.gz  Dome Flat False False False ...  True   False    False      7\n",
       "        ...        ...   ...   ...   ... ...   ...     ...      ...    ...\n",
       "b18.fits.gz       Bias False  True False ... False   False    False      7\n",
       "b19.fits.gz       Bias False  True False ... False   False    False      7\n",
       "b20.fits.gz       Bias False  True False ... False   False    False      7\n",
       "b21.fits.gz       Bias False  True False ... False   False    False      7\n",
       "b22.fits.gz       Bias False  True False ... False   False    False      7\n",
       "b23.fits.gz       Bias False  True False ... False   False    False      7\n",
       "b24.fits.gz   Feige 66 False False False ... False   False    False      7\n",
       "b27.fits.gz J1217p3905 False False False ... False   False    False      1\n",
       "b28.fits.gz J1217p3905 False False False ... False   False    False      2\n",
       "b29.fits.gz J1217p3905 False False False ... False   False    False      4"
      ]
     },
     "execution_count": 22,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "setupc.fitstbl[cols]"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "----"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Development"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 23,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "{'combine': {'method': 'mean',\n",
       "  'reject': {'cosmics': 20.0,\n",
       "   'level': [3.0, 3.0],\n",
       "   'lowhigh': [0, 0],\n",
       "   'replace': 'median'},\n",
       "  'satpix': 'reject'},\n",
       " 'useframe': 'bias'}"
      ]
     },
     "execution_count": 23,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "settings.argflag['bias']"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 24,
   "metadata": {
    "collapsed": true
   },
   "outputs": [],
   "source": [
    "ftbl = setupc.fitstbl"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 25,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "['b27.fits.gz', 'b28.fits.gz', 'b29.fits.gz']"
      ]
     },
     "execution_count": 25,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "ftbl['filename'][ftbl['science']].tolist()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "collapsed": true
   },
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.6.5"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
