from __future__ import annotations

from typing import Dict, Optional

import numpy as np
from ginga.canvas.types.basic import Polygon, Text

from pypeit.slittrace import SlitTraceSet


class SlitOverlay:
    def __init__(self, dc) -> None:
        self.dc = dc
        self.slit_polys: Dict[str, Polygon] = {}
        self.slit_labels: Dict[str, Text] = {}
        self._labels_visible: bool = False

    def build(
        self,
        slittraceset_dict: Dict[str, SlitTraceSet],
        canvas,
        slit_names: Optional[Dict[str, str]] = None,
        show_labels: bool = False,
    ) -> Dict[str, Polygon]:
        """Add slit polygons (and optionally labels) to *canvas*.

        Parameters
        ----------
        slittraceset_dict : dict
            Mapping of mosaic index → SlitTraceSet.
        canvas : DrawingCanvas
            Canvas already attached to the viewer. Caller must clear it first.
        slit_names : dict, optional
            Mapping of slit key → OBJNAME from mask design.
        show_labels : bool
            Whether labels are visible immediately after build.

        Returns
        -------
        dict
            Mapping of slit key → Polygon.
        """
        self.slit_polys.clear()
        self.slit_labels.clear()
        self._labels_visible = show_labels

        for msc_idx, slittrace in slittraceset_dict.items():
            spatial_ids = slittrace.spat_id
            left_init = slittrace.left_init.T
            right_init = slittrace.right_init.T
            sampling = 200
            y_values_left = np.arange(slittrace.nspec)[::sampling]
            y_values_right = np.arange(slittrace.nspec)[::-sampling]
            x_offset = (int(msc_idx) - 1) * slittrace.nspat

            # Label position: left edge of slit, 85% down (near bottom)
            label_row = int(0.85 * slittrace.nspec)

            for idx, spat_id in enumerate(spatial_ids):
                slit_key = f"S{spat_id}"
                x_vals = np.concatenate(
                    (left_init[idx][::sampling], right_init[idx][::-sampling]), axis=0
                ) + x_offset
                y_vals = np.concatenate((y_values_left, y_values_right), axis=0)
                poly = Polygon(
                    list(zip(x_vals.tolist(), y_vals.tolist())),
                    color="green",
                    linewidth=2,
                    fill=True,
                    fillcolor="green",
                    fillalpha=0.05,
                )
                canvas.add(poly, redraw=False)
                self.slit_polys[slit_key] = poly

                label_text = slit_key
                if slit_names and slit_key in slit_names:
                    label_text = f"{slit_key} ({slit_names[slit_key]})"
                label_x = float(left_init[idx][label_row]) + x_offset
                label_y = float(label_row)
                label = Text(
                    label_x, label_y, label_text,
                    color="yellow",
                    fontsize=10,
                    rot_deg=90,
                )
                self.slit_labels[slit_key] = label
                if show_labels:
                    canvas.add(label, redraw=False)

        canvas.update_canvas(whence=3)
        return dict(self.slit_polys)

    def set_labels_visible(self, visible: bool, canvas) -> None:
        if visible == self._labels_visible:
            return
        self._labels_visible = visible
        for label in self.slit_labels.values():
            if visible:
                if not canvas.has_object(label):
                    canvas.add(label, redraw=False)
            else:
                if canvas.has_object(label):
                    canvas.delete_object(label, redraw=False)
        canvas.update_canvas(whence=3)

    def activate(self, slit_key: str, canvas) -> None:
        poly = self.slit_polys.get(slit_key)
        if poly is None:
            return
        poly.color = "blue"
        canvas.update_canvas(whence=3)

    def deactivate(self, slit_key: str, canvas) -> None:
        poly = self.slit_polys.get(slit_key)
        if poly is None:
            return
        poly.color = "green"
        canvas.update_canvas(whence=3)
