from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional

from ginga.canvas.types.basic import Polygon

from pypeit.slittrace import SlitTraceSet


@dataclass
class QLViewState:
    raw_filepath: Optional[str] = None
    reduced_filepath: Optional[str] = None
    redux_path: str = "/hqdrpdata/DRP_TESTING/outputs"
    slittracesets: Optional[Dict[str, SlitTraceSet]] = None
    active_slit_key: Optional[str] = None
    slit_polys: Dict[str, Polygon] = field(default_factory=dict)
    manual_extract_str: Optional[str] = None
    ab_partner_filepath: Optional[str] = None
