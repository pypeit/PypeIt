"""Experimental support for the VISTA/4MOST LRS-B red detector."""

import os
from functools import cached_property
from pathlib import Path

import numpy as np

from pypeit import PypeItError, io
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
