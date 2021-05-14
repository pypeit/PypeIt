import numpy as np
from astropy import table
import astropy.constants as cons
import astropy.units as u


def blackbody_func(a, teff, outname):
    waves = np.arange(3000.0, 25000.0, 0.1) * u.AA
    # Setup the units
    teff *= u.K
    a *= 1.0E-23
    # Calculate the function
    flam = ((a*2*cons.h*cons.c**2)/waves**5)/(np.exp((cons.h*cons.c/(waves*cons.k_B*teff)).to(u.m/u.m).value)-1.0)
    flam = flam.to(u.erg / u.s / u.cm ** 2 / u.AA).value / 1.0E-17
    bb_tab = table.Table([waves.value, flam], names=("WAVELENGTH", "FLUX"))
    bb_tab.write(outname, format='fits')
    print("Wrote file: {0:s}".format(outname))


# Load the table, generate the blackbody, and save the output to a file
star_tbl = table.Table.read("blackbody_info.txt", comment='#', format='ascii')
for ss, star in enumerate(star_tbl['File']):
    blackbody_func(star_tbl['a_x10m23'][ss], star_tbl['T_K'][ss], star_tbl['File'][ss])
