"""Experimental support for the VISTA/4MOST LRS-B red detector."""

import os
from functools import cached_property
from pathlib import Path

import numpy as np

from pypeit import PypeItError, io, log
from pypeit.core import framematch
from pypeit.images.detector_container import DetectorContainer
from pypeit.par.pypeitpar import TelescopePar
from pypeit.spectrographs.spectrograph import Spectrograph


class VISTA4MOSTLRSBRedSpectrograph(Spectrograph):
    """Single-detector commissioning mode with an external trace calibration."""

    name = 'vista_4most_lrsb_red'
    camera = '4MOST LRS-B red'
    header_name = 'QMOST'
    telescope = TelescopePar(name='VISTA', longitude=-70.3976, latitude=-24.6157,
                             elevation=2529.6, diameter=4.1)
    pypeline = 'Fiber'
    fiber_require_ifu_metadata = False
    ndet = 1
    allowed_extensions = ['.fits', '.fits.fz']
    supported = False
    fiber_use_compact_profiles = True
    url = 'https://www.4most.eu/'

    @cached_property
    def aperture_scale(self):
        """
        Return the opt-in boxcar-to-wide-aperture calibration, if set.
        """
        raw = os.environ.get('PYPEIT_4MOST_APERTURE_SCALE')
        if raw is None:
            return None
        try:
            scale = float(raw)
        except ValueError as error:
            raise PypeItError('PYPEIT_4MOST_APERTURE_SCALE must be a positive number.') from error
        if not np.isfinite(scale) or scale <= 0:
            raise PypeItError('PYPEIT_4MOST_APERTURE_SCALE must be a finite positive number.')
        return scale

    def subheader_for_spec(
        self, row_fitstbl, raw_header, extra_header_cards=None, allow_missing=False
    ):
        """
        Record the opt-in aperture calibration in the spectrum primary header.
        """
        header = super().subheader_for_spec(
            row_fitstbl, raw_header, extra_header_cards=extra_header_cards,
            allow_missing=allow_missing
        )
        if self.aperture_scale is not None:
            header['APERSCL'] = (self.aperture_scale, '4MOST boxcar-to-wide aperture scale')
        return header

    def calibrate_fiber_aperture(self, sobjs, image):
        """
        Anchor bright, isolated optimal spectra to aperture-corrected boxcar counts.

        Parameters
        ----------
        sobjs : :class:`~pypeit.specobjs.SpecObjs`
            This exposure's extracted spectra with observed-frame wavelengths.
        image : `numpy.ndarray`_
            Sky-subtracted detector image in electrons, before fiber throughput division.

        The boxcar and optimal medians use identical valid rows from 8440 to
        8860 Angstrom. The grey normalization preserves spectral shape. This
        calibration only rescales the conditional uncertainties and does not
        represent uncertainty in the shared normalization.
        """
        if self.aperture_scale is None:
            return

        positions, core_counts, anchors = [], [], []
        for sobj in sobjs:
            if sobj.OPT_WAVE is None or sobj.OPT_MASK is None or sobj.BOX_MASK is None:
                positions.append(np.median(sobj.TRACE_SPAT))
                core_counts.append(np.inf)
                anchors.append(None)
                continue
            band = (sobj.OPT_WAVE >= 8440) & (sobj.OPT_WAVE < 8860) & sobj.OPT_MASK
            rows = np.flatnonzero(band)
            trace = sobj.TRACE_SPAT[rows]
            positions.append(np.median(trace) if rows.size > 0 else np.median(sobj.TRACE_SPAT))
            if rows.size < 300:
                core_counts.append(np.inf)
                anchors.append(None)
                continue
            cols = np.rint(trace).astype(int)[:, None] + np.arange(-2, 3)
            if cols.min() < 0 or cols.max() >= image.shape[1]:
                core_counts.append(np.inf)
                anchors.append(None)
                continue
            dx = cols - trace[:, None]
            weights = np.maximum(0., np.minimum(dx + .5, 2.) - np.maximum(dx - .5, -2.))
            core = np.sum(weights * image[rows[:, None], cols], axis=1)
            finite = np.isfinite(core)
            core_median = np.median(core[finite]) if np.count_nonzero(finite) >= 300 else np.inf
            core_counts.append(core_median)

            box, opt = sobj.BOX_COUNTS[rows], sobj.OPT_COUNTS[rows]
            common = sobj.BOX_MASK[rows] & np.isfinite(box) & np.isfinite(opt)
            anchors.append((np.median(box[common]), np.median(opt[common]))
                           if np.count_nonzero(common) >= 300 else None)

        positions = np.asarray(positions)
        core_counts = np.asarray(core_counts)
        count = 0
        for i, (sobj, anchor) in enumerate(zip(sobjs, anchors)):
            name = sobj.MASKDEF_OBJNAME
            if (name is None or str(name).startswith(('SKY', 'UNUSED-')) or anchor is None
                    or not np.isfinite(core_counts[i]) or core_counts[i] < 900):
                continue
            distance = abs(positions - positions[i])
            neighbors = (distance > 0) & (distance < 8.)
            if np.any(core_counts[neighbors] > .10 * core_counts[i]):
                continue
            box, opt = anchor
            if box <= 0 or opt <= 0:
                continue
            factor = self.aperture_scale * box / opt
            for field in ('OPT_COUNTS', 'OPT_COUNTS_SKY', 'OPT_COUNTS_SIG', 'OPT_COUNTS_SIG_DET'):
                value = getattr(sobj, field)
                if value is not None:
                    setattr(sobj, field, value * factor)
            for field in ('OPT_COUNTS_IVAR', 'OPT_COUNTS_NIVAR'):
                value = getattr(sobj, field)
                if value is not None:
                    setattr(sobj, field, value / factor**2)
            log.info(f'4MOST aperture calibration fiber {sobj.MASKDEF_ID}: {factor:.6f}')
            count += 1
        log.info(f'Applied the 4MOST optimal aperture calibration to {count} fibers.')

    def init_meta(self):
        self.meta = {
            key: dict(ext=0, card=card, required=False) for key, card in
            {'ra': 'RA', 'dec': 'DEC', 'target': 'OBJECT', 'mjd': 'MJD-OBS',
             'exptime': 'EXPTIME', 'airmass': 'ESO TEL AIRM START',
             'decker': 'ESO INS PATH', 'dispname': 'ESO INS PATH',
             'idname': 'ESO DPR TYPE', 'instrument': 'INSTRUME'}.items()
        }
        self.meta['binning'] = dict(ext=0, card=None, default='1,1')

    def configuration_keys(self):
        return ['dispname', 'binning']

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        good = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        kinds = {
            'science': ['OBJECT'],
            'bias': ['BIAS'],
            'arc': ['WAVE,DAY,SIMUARC'],
            'tilt': ['WAVE,DAY,SIMUARC'],
            'trace': ['FLAT,DAY'],
            'pixelflat': ['FLAT,DAY'],
        }
        match = np.isin(fitstbl['idname'], kinds.get(ftype, []))
        return good & (fitstbl['decker'] == 'CCDLRB') & match

    @classmethod
    def default_pypeit_par(cls):
        par = super().default_pypeit_par()
        par.reset_all_processimages_par(
            use_overscan=True, overscan_method='savgol', overscan_par=[0, 1], use_darkimage=False,
            use_illumflat=False, use_specillum=False, use_pixelflat=False
        )
        for kind in ('biasframe', 'arcframe', 'tiltframe', 'traceframe', 'pixelflatframe'):
            par['calibrations'][kind]['exprng'] = [None, None]
        par['calibrations']['biasframe']['process']['combine'] = 'median'
        par['calibrations']['wavelengths']['lamps'] = ['ThAr']
        par['calibrations']['wavelengths']['refframe'] = 'observed'
        par['reduce']['skysub']['joint_fit'] = True
        par['reduce']['skysub']['joint_fit_use_sci'] = False
        par['reduce']['skysub']['bspline_spacing'] = 0.7
        par['reduce']['findobj']['skip_second_find'] = True
        par['reduce']['findobj']['skip_final_global'] = True
        par['flexure']['spec_method'] = 'skip'
        par['scienceframe']['process']['spat_flexure_correct'] = False
        return par

    def get_detector_par(self, det, hdu=None):
        if det != 1:
            raise PypeItError('4MOST LRS-B red provides only detector 1.')
        gain = np.ones(4)
        noise = np.full(4, 3.0)
        if hdu is not None:
            header = hdu['RED_DATA'].header
            gain = np.array([header[f'ESO DET OUT{i} CONAD'] for i in range(1, 5)])
            noise = np.array([header[f'ESO DET OUT{i} RON'] for i in range(1, 5)])
        return DetectorContainer(det=1, binning='1,1', dataext=1, specaxis=1,
                                 specflip=True, spatflip=False, platescale=0.30,
                                 darkcurr=0.0, saturation=55000.0, nonlinear=0.82,
                                 mincounts=-1.e10, numamplifiers=4, gain=gain, ronoise=noise)

    def get_rawimage(self, raw_file, det):
        hdu = io.fits_open(raw_file)
        try:
            primary = hdu[0].header
            detector = self.get_detector_par(det, hdu)
            image = hdu['RED_DATA'].data
            setup = tuple(primary[key] for key in ('ESO INS PATH', 'ESO DET BINX',
                                                    'ESO DET BINY', 'ESO DET READ CURNAME'))
            if setup != ('CCDLRB', 1, 1, 'Science Mode') or image.shape != (6288, 6272):
                raise PypeItError(f'{raw_file} has unsupported 4MOST geometry or binning.')
            data = np.zeros(image.shape, dtype=np.uint8)
            overscan = np.zeros_like(data)
            for amp, (left, bottom) in enumerate(((True, True), (False, True),
                                                  (False, False), (True, False)), 1):
                width, height = 3072, 3080
                y = 32 if bottom else image.shape[0] - 32 - height
                x = 32 if left else image.shape[1] - 32 - width
                data[y:y+height, x:x+width] = amp
                scan = slice(x+width+4, x+width+28) if left else slice(x-28, x-4)
                overscan[y:y+height, scan] = amp
            if primary['ESO DPR CATG'] == 'SCIENCE':
                assignments = hdu['FIBINFO'].data
                self.assignments = {}
                for row in assignments[assignments['FIB_ROOT'] == 3]:
                    fibre = int(row['FIB_ID'])
                    name = {2: f'UNUSED-{fibre}', 3: 'SKY'}.get(row['FIB_USE'], str(row['OBJECT']))
                    self.assignments[int(row['SLIT_POS'])] = fibre, name
            return detector, image.astype(np.float32), hdu, primary['EXPTIME'], data, overscan
        except Exception:
            hdu.close()
            raise

    @cached_property
    def trace_calibration(self):
        trace_file = os.environ.get('PYPEIT_4MOST_TRACE')
        if trace_file is None:
            raise PypeItError('Set PYPEIT_4MOST_TRACE to the commissioning trace NPZ.')
        with np.load(Path(trace_file)) as saved:
            return {key: saved[key] for key in saved.files}

    def get_block_slit_edges(self, trace_image, det):
        cal = self.trace_calibration
        return cal['left'], cal['right']

    def bpm(self, filename, det, shape=None, msbias=None):
        mask = super().bpm(filename, det, shape=shape, msbias=msbias)
        if 'PYPEIT_4MOST_TRACE' in os.environ:
            mask[~self.trace_calibration['valid_rows'], :] = 1
        return mask

    def get_fiber_blocks(self, det):
        cal = self.trace_calibration
        positions = cal['trace'][cal['trace'].shape[0]//2]
        return [{'fiber_positions': positions[cal['block'] == i]}
                for i in range(cal['left'].shape[1])]

    def identify_fibers_in_block(self, det, block_idx, detected_positions):
        cal = self.trace_calibration
        slit = cal['slit_position'][cal['block'] == block_idx]
        assignments = self.assignments
        ids, names = zip(*(assignments[int(position)] for position in slit))
        return {'fiber_id': np.asarray(ids), 'fiber_name': np.asarray(names)}

    def get_fiber_traces(self, det, fiber_ids):
        cal = self.trace_calibration
        by_id = {int(fibre): slit for slit, (fibre, _) in self.assignments.items()}
        idx = {int(slit): i for i, slit in enumerate(cal['slit_position'])}
        return cal['trace'][:, [idx[by_id[fibre]] for fibre in fiber_ids]].T
