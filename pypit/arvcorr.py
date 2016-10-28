import numpy as np
from pypit import armsgs
from pypit import arparse

# Logging and settings
msgs = armsgs.get_logger()
spect = arparse.get_spect().__dict__['_spect']

jyear = 365.25       # Julian Year
j2000 = 2000.0       # J2000 epoch
jd2000 = 2451545.0   # J2000 Julian Date
sjyear = 31557600.0  # convert seconds to Julian year
uttol = 4.0          # Tolerance for UT diffs in header values (in minutes)


def coord(ao, bo, ap, bp, a1, b1):
    x = np.cos(a1)*np.cos(b1)
    y = np.sin(a1)*np.cos(b1)
    z = np.sin(b1)
    xp = np.cos(ap)*np.cos(bp)
    yp = np.sin(ap)*np.cos(bp)
    zp = np.sin(bp)

    # Rotate the origin about z
    sao = np.sin(ao)
    cao = np.cos(ao)
    sbo = np.sin(bo)
    cbo = np.cos(bo)
    temp = -xp*sao + yp*cao
    xp = xp*cao + yp*sao
    yp = temp
    temp = -x*sao + y*cao
    x = x*cao + y*sao
    y = temp

    # Rotate the origin about y
    temp = -xp*sbo + zp*cbo
    # xp = xp*cbo + zp*sbo
    zp = temp
    temp = -x*sbo + z*cbo
    x = x*cbo + z*sbo
    z = temp

    # Rotate pole around x
    sbp = zp
    cbp = yp
    temp = y*cbp + z*sbp
    y = y*sbp - z*cbp
    z = temp
  
    # Final angular coordinates
    l = np.arctan2(y, x)
    b = np.arcsin(z)
    return l, b


def date_to_epoch(year, month, day, ut):
    jd = date_to_jd(year, month, day, ut)
    return j2000+(jd-jd2000)/jyear


def date_to_jd(year, month, day, ut):
    y = year
    if (month < 1) or (month > 12):
        msgs.error("Invalid month ({0:d}) when converting date to JD".format(month))
    elif month > 2:
        m = month+1
    else:
        m = month + 13
        y -= 1
    jd = float(int(jyear*y) + int(30.6001*m) + day + 1720995)
    if day+31*(m+12*y) >= 588829:
        d = int(y/100)
        m = int(y/400)
        jd += float(2-d+m)
    jd += float(int(ut*360000.0+0.5))/(360000.0*24.0)-0.5
    return jd


def epoch_to_jd(epoch):
    jd = jd2000 + (epoch - j2000)*jyear
    return jd


def jd_to_date(j):
    ja = int(j+0.5)
    ut = 24.0*(j+0.5-float(ja))
    if ja >= 2299161:
        jb = int((float(ja-1867216)-0.25)/36524.25)
        ja = ja+1+jb-int(jb/4)
    jb = ja + 1524
    jc = int(6680.0+(float(jb-2439870)-122.1)/jyear)
    jd = 365*jc + int(jc/4)
    je = int(float(jb-jd)/30.6001)
    day = jb - jd - int(30.6001*je)
    month = je-1
    if month > 12:
        month -= 12
    year = jc-4715
    if month > 2:
        year -= 1
    if year < 0:
        year -= 1
    return year, month, day, ut


def mst(epoch, lon):
    # Determine JD and UT, and T (JD in centuries from J2000.0)
    jd = epoch_to_jd(epoch)
    ut = (jd-int(jd)-0.5)*24.0
    t = (jd-jd2000)/(100.0*jyear)

    # The GMST at 0 UT in seconds is a power series in T
    st = 24110.54841+t*(8640184.812866+t*(0.093104-t*6.2e-6))

    # Correct for longitude and convert to standard hours
    st = (st/3600.0 + ut - lon/15.0)-(int((st/3600.0 + ut - lon/15.0)/24.0))*24.0
    if st < 0.0:
        st += 24.0
    return st


def radec_to_decdeg(ra, dec):
    if ":" in ra:
        raspl = ra.split(":")
    else:
        raspl = ra.split(" ")
    if ":" in dec:
        dcspl = dec.split(":")
    else:
        dcspl = dec.split(" ")
    raout = 15.0*(float(raspl[0]) + float(raspl[1])/60.0 + float(raspl[2])/3600.0)
    if "-" in dcspl[0]:
        dcout = float(dcspl[0]) - float(dcspl[1])/60.0 - float(dcspl[2])/3600.0
    else:
        dcout = float(dcspl[0]) + float(dcspl[1])/60.0 + float(dcspl[2])/3600.0
    return raout/15.0, dcout


def rotate(epoch):
    # The rotation matrix coefficients are polynomials in time measured in
    # Julian centuries from the standard epoch. The coefficients are in
    # radians

    t = (epoch_to_jd(epoch)-jd2000)/(100.0*jyear)
    a = t*(0.6406161+t*(0.0000839+t*0.0000050))*np.pi/180.0
    b = t*(0.6406161+t*(0.0003041+t*0.0000051))*np.pi/180.0
    c = t*(0.5567530-t*(0.0001185+t*0.0000116))*np.pi/180.0

    # Compute the cosines and sines of these angles
    ca = np.cos(a)
    sa = np.sin(a)
    cb = np.cos(b)
    sb = np.sin(b)
    cc = np.cos(c)
    sc = np.sin(c)

    # Now compute and then return the rotation matrix
    p = np.zeros((3, 3))
    p[0, 0] = ca*cb*cc-sa*sb
    p[1, 0] = -sa*cb*cc-ca*sb
    p[2, 0] = -cb*sc
    p[0, 1] = ca*sb*cc+sa*cb
    p[1, 1] = -sa*sb*cc+ca*cb
    p[2, 1] = -sb*sc
    p[0, 2] = ca*sc
    p[1, 2] = -sa*sc
    p[2, 2] = cc
    return p


def precess(ra1, dec1, epoch1, epoch2):
    r0 = np.zeros(3)
    r1 = np.zeros(3)

    if epoch1 == 0.0 or epoch1 == epoch2:
        return ra1, dec1

    # Rectangular equitorial coordinates (direction cosines)
    ra2 = ra1 * 15.0 * np.pi/180.0
    dec2 = dec1 * np.pi/180.0
    r0[0] = np.cos(ra2) * np.cos(dec2)
    r0[1] = np.sin(ra2) * np.cos(dec2)
    r0[2] = np.sin(dec2)

    if epoch1 != j2000:
        p = rotate(epoch1)
        r1[0] = p[0, 0]*r0[0] + p[0, 1]*r0[1] + p[0, 2]*r0[2]
        r1[1] = p[1, 0]*r0[0] + p[1, 1]*r0[1] + p[1, 2]*r0[2]
        r1[2] = p[2, 0]*r0[0] + p[2, 1]*r0[1] + p[2, 2]*r0[2]
        r0[0] = r1[0]
        r0[1] = r1[1]
        r0[2] = r1[2]

    if epoch2 != j2000:
        p = rotate(epoch2)
        r1[0] = p[0, 0]*r0[0] + p[1, 0]*r0[1] + p[2, 0]*r0[2]
        r1[1] = p[0, 1]*r0[0] + p[1, 1]*r0[1] + p[2, 1]*r0[2]
        r1[2] = p[0, 2]*r0[0] + p[1, 2]*r0[1] + p[2, 2]*r0[2]
        r0[0] = r1[0]
        r0[1] = r1[1]
        r0[2] = r1[2]

    # Convert from radians to hours and degrees
    ra2 = np.arctan2(r0[1], r0[0]) / (15.0*(np.pi/180.0))
    dec2 = np.arcsin(r0[2]) / (np.pi/180.0)
    if ra2 < 0.0:
        ra2 += 24.0
    return ra2, dec2


def helio_corr(slf, idx):
    """
    Correct the wavelength calibration solution to the heliocentric reference frame
    """
    # Get the modified Julian day, and convert to the Julian day
    mjd = float(slf._fitsdict["time"][idx])
    jd = mjd + 2400000.5
    # Get the exposure time (in seconds)
    exptime = float(slf._fitsdict["exptime"][idx])
    # Get the equinox
    equinox = float(slf._fitsdict["equinox"][idx])
    # Get the coordinates of the target (convert to hours for RA, and degrees for DEC)
    tra = slf._fitsdict["ra"][idx]
    tdec = slf._fitsdict["dec"][idx]
    ra, dec = radec_to_decdeg(tra, tdec)
    hdr_jd = jd
    hdr_exptime = exptime
    hdr_ra = ra
    hdr_dec = dec
    hdr_equ = equinox
    hdr_lat = spect['mosaic']['latitude']
    hdr_lon = spect['mosaic']['longitude']
    hdr_alt = spect['mosaic']['elevation']
    vhel = vhelio(hdr_jd, hdr_exptime, hdr_ra, hdr_dec, hdr_equ, hdr_lat, hdr_lon, hdr_alt)
    msgs.info("Heliocentric velocity correction = {0:+.4f} km/s for file:".format(vhel) + msgs.newline() +
              slf._fitsdict["filename"][idx])
    w = np.where(slf._waveids == -999999.9)
    slf._waveids += slf._waveids*vhel/299792.458
    slf._waveids[w] = -999999.9
    return slf._waveids


def vhelio(hdr_jd, hdr_exptime, hdr_ra, hdr_dec, hdr_equ, hdr_lat, hdr_lon, hdr_alt):
    year, month, day, ut = jd_to_date(hdr_jd)
    msgs.work("Check the year, month, date and UT just calculated match the header")

    # Assume the MJD has the most accurate start-of-exspoure date+time information
    hdr_year = year
    hdr_month = month
    hdr_day = day
    hdr_ut = ut
    hdr_epoch = date_to_epoch(hdr_year, hdr_month, hdr_day, hdr_ut)

    # Modify the JD of the observation to be approximately in the middle of the exposure
    msgs.work("Make sure exposure time is in seconds")
    hdr_epoch += 0.5*hdr_exptime/sjyear

    # Precess the object coordinates (J2000 in most cases) to the epoch in
    # which the observations were carried out
    ra, dec = precess(hdr_ra, hdr_dec, hdr_equ, hdr_epoch)

    # Calculate the Earth-moon barycentric velocity relative to Sun (annual)
    vorb = vorbit(ra, dec, hdr_epoch)

    # Calculate the Earth's velocity around Earth-Moon barycentre (lunar)
    vbary = vbarymoon(ra, dec, hdr_epoch)

    # Calculate the observer's motion relative to the Earth's centre (diurnal)
    vrot = vrotate(ra, dec, hdr_epoch, hdr_lat, hdr_lon, hdr_alt)

    # Calculate heliocentric velocity
    vhel = vorb + vbary + vrot

    return vhel


def vbarymoon(ra, dec, epoch):
    # T is the number of Julian centuries since J1900
    t = (epoch_to_jd(epoch)-2415020.0)/36525.0

    # OBLQ is the mean obliquity of the ecliptic
    # OMEGA is the longitude of the mean ascending node
    # LLONG is the mean lunar longitude (should be 13.1763965268)
    # LPERI is the mean lunar longitude of perigee
    # INCLIN is the inclination of the lunar orbit to the ecliptic
    # EM is the eccentricity of the lunar orbit (dimensionless)
    # All quantities except the eccentricity are in degrees

    oblq = 23.452294-t*(0.0130125+t*(0.00000164-t*0.000000503))
    omega = 259.183275-t*(1934.142008+t*(0.002078+t*0.000002))
    llong = 270.434164+t*(481267.88315+t*(-0.001133+t*0.0000019))-omega
    lperi = 334.329556+t*(4069.034029-t*(0.010325+t*0.000012))-omega
    em = 0.054900489
    inclin = 5.1453964
    emsq = em*em
    emcu = emsq*em

    # Determine true longitude.  Compute mean anomaly, convert to true
    # anomaly (approximate formula), and convert back to longitude. The mean
    # anomaly is only approximate because LPERI should be the true rather
    # than the mean longitude of lunar perigee
    lperi *= np.pi/180.0
    llong *= np.pi/180.0
    anom = llong-lperi
    anom += (2.0*em-0.25*emcu)*np.sin(anom)+1.25*emsq*np.sin(2.0*anom)+13.0/12.0*emcu*np.sin(3.0*anom)
    llong = anom+lperi

    # L and B are the ecliptic longitude and latitude of the observation. LM
    # and BM are the lunar longitude and latitude of the observation in the
    # lunar orbital plane relative to the ascending node
    r = (np.pi/180.0)*ra*15.0
    d = (np.pi/180.0)*dec
    omega *= np.pi/180.0
    oblq *= np.pi/180.0
    inclin *= np.pi/180.0

    l, b = coord(0.0, 0.0, -0.5*np.pi, 0.5*np.pi-oblq, r, d)
    lm, bm = coord(omega, 0.0, omega-0.5*np.pi, 0.5*np.pi-inclin, l, b)

    # VMOON is the component of the lunar velocity perpendicular to the
    # radius vector.  V is the projection onto the line of sight to the
    # observation of the velocity of the Earth's center with respect to the
    # Earth-Moon barycenter.  The 81.53 is the ratio of the Earth's mass to
    # the Moon's mass
    vmoon = (2.0*np.pi/27.321661)*384403.12040/(np.sqrt(1.0-emsq)*86400.0)
    vmoon *= np.cos(bm)*(np.sin(llong-lm)-em*np.sin(lperi-lm))/81.53
    return vmoon


def vorbit(ra, dec, epoch):

    # t is the number of Julian centuries since J1900
    t = (epoch_to_jd(epoch)-2415020.0)/36525.0
    # MANOM is the mean anomaly of the Earth's orbit (degrees)
    # LPERI is the mean longitude of perihelion (degrees)
    # OBLQ is the mean obliquity of the ecliptic (degrees)
    # ECCEN is the eccentricity of the Earth's orbit (dimensionless)
    manom = 358.47583+t*(35999.04975-t*(0.000150+t*0.000003))
    lperi = 101.22083+t*(1.7191733+t*(0.000453+t*0.000003))
    oblq = 23.452294-t*(0.0130125+t*(0.00000164-t*0.000000503))
    eccen = 0.01675104-t*(0.00004180+t*0.000000126)
    eccensq = eccen*eccen
    eccencu = eccensq*eccen

    # Convert to principle angles
    manom -= int(manom/360.0)*360.0
    lperi -= int(lperi/360.0)*360.0
    # Convert to radians
    r = ra*15.0*np.pi/180.0
    d = dec*np.pi/180.0
    manom *= np.pi/180.0
    lperi *= np.pi/180.0
    oblq *= np.pi/180.0

    # TANOM is the true anomaly (approximate formula) (radians)
    tanom = manom + (2.0*eccen-0.25*eccencu)*np.sin(manom)
    tanom += 1.25*eccensq*np.sin(2.0*manom)
    tanom += (13.0/12.0)*eccencu*np.sin(3.0*manom)
    # SLONG is the true longitude of the Sun seen from the Earth (radians)
    slong = lperi + tanom + np.pi

    # L and B are the longitude and latitude of the star in the orbital
    # plane of the Earth (radians)
    l, b = coord(0.0, 0.0, -np.pi/2.0, np.pi/2.0-oblq, r, d)

    # VORB is the component of the Earth's orbital velocity perpendicular to
    # the radius vector (km/s) where the Earth's semi-major axis is
    # 149598500 km and the year is 365.2564 days
    vorb = ((2.0*np.pi/365.2564)*149598500.0/np.sqrt(1.0-eccensq))/86400.0

    # V is the projection onto the line of sight to the observation of the
    # velocity of the Earth-Moon barycenter with respect to the Sun (km/s)
    vorb *= np.cos(b)*(np.sin(slong-l)-eccen*np.sin(lperi-l))
    return vorb


def vrotate(ra, dec, epoch, lat, lon, alt):
    # LATD is the latitude in radians
    latd = (np.pi/180.0)*lat

    # Reduction of geodetic latitude to geocentric latitude. Dlatd is in arcseconds
    dlatd = -(11.0*60.0+32.743000)*np.sin(2.0*latd)+1.163300*np.sin(4.0*latd)-0.002600*np.sin(6.0*latd)
    latd += (np.pi/180.0)*dlatd/3600.0

    # R is the radius vector from the Earth's center to the observer
    # (meters).  Vc is the corresponding circular velocity (meters/sidereal
    # day converted to km / sec). (sidereal day = 23.934469591229 hours (1986))
    r = 6378160.0*(0.998327073+0.00167643800*np.cos(2.0*latd)-0.00000351*np.cos(4.0*latd)+0.000000008*np.cos(6.0*latd))
    r += alt
    vc = 2.0*np.pi*(r/1000.0)/(23.934469591229*3600.0)

    # Project the velocity onto the line of sight to the star
    lmst = mst(epoch, lon)
    vc *= np.cos(latd)*np.cos((np.pi/180.0)*dec)*np.sin((np.pi/180.0)*(ra-lmst)*15.0)
    return vc
