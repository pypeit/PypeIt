""" States to monitor the progress of the reduction """
from pydantic import BaseModel, Field
from typing import List, Optional, Dict, Literal

# Hopefully this isn't circular
import os
import io
import json

from IPython import embed

# Calibration state
class BaseCalibState(BaseModel):
    calib_id: int # Calibration ID
    det: int | List[int]  # Detector number or mosaic tuple
    step: str
    required: bool = False
    input_files: Optional[List[str]] = None
    output_file: Optional[str] = None
    qa_files: Optional[List[str]] = None
    status: Literal["complete", "fail", "undone", "running", "success"] = "undone"

class BiasCalibState(BaseCalibState):
    step: Literal["bias"] = "bias"
    # Metrics
    mean: Optional[float] = None
    std: Optional[float] = None

class WvCalibSlit(BaseModel):
    status: Literal["success", "fail", "undone", ] = "undone"
    # Metrics
    rms: Optional[float] = None

class WvCalibState(BaseCalibState):
    step: Literal["wv_calib"] = "wv_calib"
    slits: Optional[Dict[int, WvCalibSlit]] = Field(default_factory=dict)

class SlitEdges(BaseModel):
    status: Literal["success", "fail", "undone", ] = "undone"
    # Metrics
    center: Optional[float] = None
    slitord_id: Optional[int] = None

class SlitEdgesState(BaseCalibState):
    step: Literal["slits"] = "slits"
    nslits: Optional[int] = None
    slits: Optional[Dict[int, SlitEdges]] = Field(default_factory=dict)

class TiltsSlit(BaseModel):
    status: Literal["success", "fail", "undone", ] = "undone"
    # Metrics
    rms: Optional[float] = None

class TiltsState(BaseCalibState):
    step: Literal["tilts"] = "tilts"
    slits: Optional[Dict[int, TiltsSlit]] = Field(default_factory=dict)

class FlatsState(BaseCalibState):
    step: Literal["flats"] = "flats"
    types: Optional[List[str]] = []

class DarkCalibState(BaseCalibState):
    step: Literal["dark"] = "dark"

class ArcCalibState(BaseCalibState):
    step: Literal["arc"] = "arc"

class TiltImgCalibState(BaseCalibState):
    step: Literal["tiltimg"] = "tiltimg"

class ScattLightCalibState(BaseCalibState):
    step: Literal["scattlight"] = "scattlight"

class AlignCalibState(BaseCalibState):
    step: Literal["align"] = "align"

calib_classes = {
    'bias': BiasCalibState,
    'dark': DarkCalibState,
    'arc': ArcCalibState,
    'tiltimg': TiltImgCalibState,
    'wv_calib': WvCalibState,
    'tilts': TiltsState,
    'scattlight': ScattLightCalibState,
    'flats': FlatsState,
    'slits': SlitEdgesState,
    'align': AlignCalibState,
}

slit_classes = {
    'wv_calib': WvCalibSlit,
    'tilts': TiltsSlit,
    'slits': SlitEdges
}

class RunPypeItState(BaseModel):
    """ State of the PypeIt run """

    # Required
    pypeit_file: str
    current_step: str
    current_det: int
    current_calibID: int

    # Optional
    previous_step: str = 'none'

    # Calibrations
    bias: Optional[List[BiasCalibState]] = Field(default_factory=list)
    dark: Optional[List[DarkCalibState]] = Field(default_factory=list)
    arc: Optional[List[ArcCalibState]] = Field(default_factory=list)
    tiltimg: Optional[List[TiltImgCalibState]] = Field(default_factory=list)
    slits: Optional[List[SlitEdgesState]] = Field(default_factory=list)
    wv_calib: Optional[List[WvCalibState]] = Field(default_factory=list)
    tilts: Optional[List[TiltsState]] = Field(default_factory=list)
    scattlight: Optional[List[ScattLightCalibState]] = Field(default_factory=list)
    flats: Optional[List[FlatsState]] = Field(default_factory=list)
    align: Optional[List[AlignCalibState]] = Field(default_factory=list)
    path: Optional[str] = None

    @property
    def outfile(self):
        outfile = self.pypeit_file.replace('.pypeit', '_state.json') if self.path is None else self.path
        return outfile

    # Load existing state 
    def load(self, path:str=None):
        if not os.path.isfile(self.outfile):
            return self
        print("Loading existing state from {:s}".format(self.outfile))
        with open(self.outfile, 'rt') as fh:
            update_dict = json.load(fh)
        # Return
        return RunPypeItState.model_validate(update_dict)
        

    def update_calib(self, step:str, calib_id: int, det: str, key:str, value,
                     slit:str=None):
        # Current step
        if self.current_step != step:
            self.previous_step = self.current_step
        self.current_step = step
        # Select items so far
        if step not in calib_classes:
            return
        # Grab the entry
        self_items = getattr(self, step)
        found_it = False
        # Grab the tiem
        for index, item in enumerate(self_items):
            # TODO -- if det is a tuple, this will probably fail
            if item.calib_id == calib_id and item.det == det:
                found_it = True
                break

        # Create it?
        if not found_it:
            item = calib_classes[step](calib_id=calib_id, det=det)
            self_items.append(item)
            index = -1

        # Set
        if slit is None:
            if isinstance(getattr(self_items[index],key), list):
                getattr(self_items[index],key).append(value)
            else:
                setattr(self_items[index], key, value)
        else:
            if slit not in self_items[index].slits.keys():
                self_items[index].slits[slit] = slit_classes[step]()
            setattr(self_items[index].slits[slit], key, value)

    def write(self):
        json_string = self.model_dump_json(exclude_none=True, indent=4, round_trip=True)
        # Write
        with io.open(self.outfile, 'w', encoding='utf-8') as f:
            #f.write(json.dumps(obj, sort_keys=True, indent=4,
            #                   separators=(',', ': '), **kwargs))
            f.write(json_string)

    def print_status(self):
        """
        Print a pretty-formatted summary of the calibration status
        to stdout.
        """
        # Collect all unique (calib_id, det) pairs across all steps
        pairs = set()
        for step in calib_classes:
            for item in getattr(self, step):
                pairs.add((item.calib_id, item.det))

        print(f'PypeIt Reduction Status: {os.path.basename(self.pypeit_file)}')
        print('=' * 70)

        if len(pairs) == 0:
            print('  No calibration state entries found.')
            return

        # Sort by calib_id then det
        for calib_id, det in sorted(pairs):
            print(f'\n  Calibration Group: {calib_id}, Detector: {det}')
            print(f'  {"Step":<14s} {"Required":<10s} {"Status":<10s} {"Output File"}')
            print(f'  {"----":<14s} {"--------":<10s} {"------":<10s} {"-----------"}')
            # Loop over steps in the order defined by calib_classes
            for step in calib_classes:
                items = getattr(self, step)
                # Find the entry for this (calib_id, det)
                entry = None
                for item in items:
                    if item.calib_id == calib_id and item.det == det:
                        entry = item
                        break
                if entry is None:
                    print(f'  {step:<14s} {"--":<10s} {"--":<10s} --')
                else:
                    req_str = str(entry.required)
                    stat_str = entry.status
                    outfile = os.path.basename(entry.output_file) \
                        if entry.output_file is not None else '--'
                    print(f'  {step:<14s} {req_str:<10s} {stat_str:<10s} {outfile}')
        print()
