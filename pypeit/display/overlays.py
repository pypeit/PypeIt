"""
Fast ginga canvas objects used to overplot PypeIt data products.

The objects in this module decide what to draw each time they are rendered,
based on the region of the image currently in view:

    - Only the curves (and parts of curves) inside the view are drawn, with a
      limited number of vertices per curve.  Ginga's rendering cost scales
      with the number of vertices drawn, so this keeps redraws fast for
      detectors with thousands of spectral pixels and hundreds of traces.

    - For very large numbers of curves (e.g., the hundreds of fibers of a
      fiber-fed IFU), the curves are not drawn when the view is zoomed out so
      far that they could not be distinguished; a note asking the user to
      zoom in is shown instead.

.. warning::

    Importing this module imports ginga's widget modules, which requires
    ginga's widget toolkit to have been chosen.  Import it only after the
    viewer has been created (see :mod:`pypeit.display.viewer`).

.. include:: ../include/links.rst
"""
from ginga import toolkit
if toolkit.toolkit == 'choose':
    toolkit.choose()

import numpy as np

from ginga.canvas.types.layer import CompoundObject
from ginga.canvas.types.basic import Path, Text, SquareBox


def _view_key(viewer):
    """Return a hashable description of what the viewer currently shows."""
    x1, y1, x2, y2 = viewer.get_datarect()
    sx, sy = viewer.get_scale_xy()
    wd, ht = viewer.get_window_size()
    return (round(x1, 1), round(y1, 1), round(x2, 1), round(y2, 1), round(sx, 6),
            round(sy, 6), wd, ht)


class _ViewDependentOverlay(CompoundObject):
    """
    Base class for compound objects whose children depend on the view.

    Subclasses implement :meth:`build`, which returns the list of ginga canvas
    objects to draw for a given viewer.  The objects are cached for each
    viewer (the same canvas can be drawn by more than one viewer, e.g. the
    channel viewer and ginga's Pan plugin) and rebuilt only when the view
    changes.
    """
    def __init__(self):
        super().__init__()
        self.kind = 'pypeit_overlay'
        self._cache = {}

    def build(self, viewer):
        """Return the objects to draw in ``viewer``; implemented by subclasses."""
        raise NotImplementedError

    def state_key(self):
        """
        Return a hashable description of any state, other than the view,
        that changes what is drawn (e.g., the curve under the cursor).
        """
        return None

    def draw(self, viewer):
        """
        Draw the overlay.

        Parameters
        ----------
        viewer : ``ginga.ImageView.ImageViewBase``
            The ginga viewer doing the rendering.
        """
        if viewer.get_image() is None:
            return
        key = (_view_key(viewer), self.state_key())
        cached = self._cache.get(id(viewer))
        if cached is None or cached[0] != key:
            objs = self.build(viewer)
            for obj in objs:
                obj.initialize(None, viewer, self.logger)
            self._cache[id(viewer)] = (key, objs)
        else:
            objs = cached[1]
        for obj in objs:
            obj.draw(viewer)


class TraceOverlay(_ViewDependentOverlay):
    r"""
    Draw a set of curves (object traces, slit edges, ...) and their labels.

    Parameters
    ----------
    x : `numpy.ndarray`_
        Horizontal (spatial) image coordinates of the curves.  Shape must be
        :math:`(N_{\rm spec},)` or :math:`(N_{\rm spec}, N_{\rm curve})`.
    y : `numpy.ndarray`_, optional
        Vertical (spectral) image coordinates of the curves.  Shape must be
        :math:`(N_{\rm spec},)` or the same as ``x``.  If None, the
        coordinates are the pixel indices along the spectral axis.
    colors : :obj:`str`, array-like, optional
        Color for all curves or one color per curve.
    labels : array-like, optional
        Text label for each curve.  If None, no labels are drawn.
    label_colors : :obj:`str`, array-like, optional
        Color for all labels or one per curve.  If None, use ``colors``.
    fontsize : :obj:`float`, optional
        Font size for the labels.
    rot_deg : :obj:`float`, optional
        Rotation of the labels in degrees.
    label_offset : :obj:`tuple`, optional
        Offset of the labels with respect to their curve, in *screen* pixels
        (horizontal, vertical).
    draw_lines : :obj:`bool`, optional
        Draw the curves.  If False, only the labels are drawn (e.g., to
        label slits at their center).
    linewidth : :obj:`int`, optional
        Line width in screen pixels.
    max_vertices : :obj:`int`, optional
        Maximum number of vertices drawn for the visible part of each curve.
    min_stride : :obj:`int`, optional
        Minimum step between the curve points drawn.
    lod_ncurve : :obj:`int`, optional
        If the number of curves is larger than this, curves are only drawn
        when the view is zoomed in enough that adjacent curves are at least
        ``min_spacing`` screen pixels apart.
    min_spacing : :obj:`float`, optional
        Minimum separation in screen pixels between adjacent curves needed
        to draw them when there are more than ``lod_ncurve`` curves.
    what : :obj:`str`, optional
        Name for the curves used in the note shown when they are not drawn;
        e.g. "traces".
    note_index : :obj:`int`, optional
        Index of the line in the top-left corner of the viewer used for the
        note shown when the curves are not drawn.
    owner : ``ginga.ImageView.ImageViewBase``, optional
        The viewer of the channel with this overlay.  The labels and the
        note are only drawn in this viewer, not in other viewers that draw
        the same canvas (e.g., ginga's Pan plugin).  If None, they are drawn
        in all viewers.

    Notes
    -----
    Labels that would overlap a label already drawn (to its left, for
    vertical labels) are not drawn.  The name of any curve can be seen by
    hovering over it (see :meth:`enable_hover`).
    """
    def __init__(self, x, y=None, colors='orange', labels=None, label_colors=None,
                 fontsize=12., rot_deg=90., label_offset=(0., 0.), draw_lines=True, linewidth=1,
                 max_vertices=100, min_stride=1, lod_ncurve=200, min_spacing=4.,
                 what='traces', note_index=0, owner=None):
        super().__init__()
        self.owner = owner
        _x = np.asarray(x, dtype=float)
        self.x = _x.reshape(-1, 1) if _x.ndim == 1 else _x
        nspec, ncurve = self.x.shape
        if y is None:
            self.y = np.broadcast_to(np.arange(nspec, dtype=float)[:, None], self.x.shape)
        else:
            _y = np.asarray(y, dtype=float)
            self.y = np.broadcast_to(_y.reshape(-1, 1) if _y.ndim == 1 else _y, self.x.shape)
        self.ncurve = ncurve
        self.colors = [colors]*ncurve if colors is None or isinstance(colors, str) \
                        else list(colors)
        self.labels = None if labels is None else [str(lbl) for lbl in labels]
        if label_colors is None:
            self.label_colors = self.colors
        elif isinstance(label_colors, str):
            self.label_colors = [label_colors]*ncurve
        else:
            self.label_colors = list(label_colors)
        self.fontsize = fontsize
        self.rot_deg = rot_deg
        self.label_offset = label_offset
        self.draw_lines = draw_lines
        self.linewidth = linewidth
        self.max_vertices = max_vertices
        self.min_stride = max(int(min_stride), 1)
        self.lod_ncurve = lod_ncurve
        self.min_spacing = min_spacing
        self.what = what
        self.note_index = note_index
        # Index of the curve under the cursor (see enable_hover)
        self.hover = None

        # Typical separation between adjacent curves, in image pixels, used
        # to decide if the curves can be distinguished at the current zoom.
        mid = nspec // 2
        xmid = np.sort(self.x[mid][np.isfinite(self.x[mid])])
        self.spacing = np.median(np.diff(xmid)) if xmid.size > 1 else np.inf

    def build(self, viewer):
        """
        Construct the ginga objects to draw for the current view.

        Parameters
        ----------
        viewer : ``ginga.ImageView.ImageViewBase``
            The ginga viewer doing the rendering.

        Returns
        -------
        :obj:`list`
            The ginga canvas objects to draw.
        """
        x1, y1, x2, y2 = viewer.get_datarect()
        sx, sy = viewer.get_scale_xy()

        if self.ncurve > self.lod_ncurve and self.spacing * sx < self.min_spacing:
            if not self._is_owner(viewer):
                return []
            return [Text(10, 20 + 20*self.note_index,
                         f'Zoom in to see the {self.ncurve} {self.what}', coord='window',
                         color=self.colors[0], fontsize=12, bgcolor='black', bgalpha=0.6)]

        # Only draw the labels if there is room for them between adjacent
        # curves; otherwise they overlap and are unreadable.
        draw_labels = self.labels is not None and self._is_owner(viewer) \
                        and (self.ncurve <= self.lod_ncurve
                             or self.spacing * sx >= 1.2 * self.fontsize)

        # Expand the view by a margin so that curves entering the view are
        # drawn continuously up to the edge
        mx = 0.05*(x2 - x1) + 2.
        my = 0.05*(y2 - y1) + 2.
        inview = (self.x >= x1 - mx) & (self.x <= x2 + mx) \
                    & (self.y >= y1 - my) & (self.y <= y2 + my)
        visible = np.where(np.any(inview, axis=0))[0]
        if visible.size == 0:
            return []
        first = np.argmax(inview[:, visible], axis=0)
        last = inview.shape[0] - 1 - np.argmax(inview[::-1, visible], axis=0)

        objs = []
        labels = []
        hover_span = None
        dx = self.label_offset[0] / sx
        dy = self.label_offset[1] / sy
        for k, i0, i1 in zip(visible, first, last):
            i0 = max(i0 - 1, 0)
            i1 = min(i1 + 1, self.x.shape[0] - 1)
            if self.draw_lines and self.colors[k] is not None:
                step = max(self.min_stride, int(np.ceil((i1 - i0 + 1) / self.max_vertices)))
                indx = np.arange(i0, i1 + 1, step)
                if indx[-1] != i1:
                    indx = np.append(indx, i1)
                pts = np.column_stack([self.x[indx, k], self.y[indx, k]])
                pts = pts[np.all(np.isfinite(pts), axis=1)]
                if len(pts) > 1:
                    objs.append(Path(pts, color=self.colors[k], linewidth=self.linewidth))
            if draw_labels and self.label_colors[k] is not None:
                im = (i0 + i1) // 2
                if np.isfinite(self.x[im, k]) and np.isfinite(self.y[im, k]):
                    labels.append((self.x[im, k] + dx, self.y[im, k] + dy, k))
            if k == self.hover:
                hover_span = (i0, i1)
        objs += self._declutter(labels, sx, sy)
        if self.hover is not None and self.labels is not None and hover_span is not None \
                and self._is_owner(viewer):
            objs += self._hover_objects(viewer, *hover_span)
        return objs


    def _is_owner(self, viewer):
        """Check if labels and notes should be drawn in this viewer."""
        return self.owner is None or viewer is self.owner

    def _declutter(self, labels, sx, sy):
        """
        Construct the label objects, skipping labels that would overlap.

        Labels are considered from left to right (vertical labels) or bottom
        to top (horizontal labels), and a label is skipped if it is closer
        than the label height to the previous label drawn.

        Parameters
        ----------
        labels : :obj:`list`
            List of (x, y, curve index) tuples with the label positions in
            image coordinates.
        sx, sy : :obj:`float`
            Image to screen scale factors.

        Returns
        -------
        :obj:`list`
            The ginga text objects.
        """
        if len(labels) == 0:
            return []
        vertical = abs(np.sin(np.radians(self.rot_deg))) > 0.5
        coo = np.array([lbl[0] * sx if vertical else lbl[1] * sy for lbl in labels])
        mindist = 1.2 * self.fontsize
        objs = []
        last = -np.inf
        for i in np.argsort(coo):
            if coo[i] - last < mindist:
                continue
            last = coo[i]
            x, y, k = labels[i]
            objs.append(Text(x, y, self.labels[k], color=self.label_colors[k],
                             fontsize=self.fontsize, rot_deg=self.rot_deg))
        return objs

    def state_key(self):
        """Return the index of the curve under the cursor."""
        return self.hover

    def enable_hover(self, viewer):
        """
        Highlight and label the curve under the cursor.

        The curve closest to the cursor (within half the typical separation
        between curves) is drawn in cyan with its label, even if the labels
        are otherwise not drawn.

        Parameters
        ----------
        viewer : ``ginga.ImageView.ImageViewBase``
            The ginga viewer showing the canvas with this overlay.
        """
        if self.labels is not None:
            viewer.add_callback('cursor-changed', self._cursor_cb)

    def _cursor_cb(self, viewer, button, data_x, data_y):
        """Cursor callback: update the curve under the cursor."""
        canvas = viewer.get_canvas()
        if self.tag is None or not canvas.has_tag(self.tag):
            # The overlay has been removed from the canvas
            viewer.remove_callback('cursor-changed', self._cursor_cb)
            return
        hover = self.nearest(viewer, data_x, data_y)
        if hover != self.hover:
            self.hover = hover
            canvas.update_canvas(whence=3)

    def nearest(self, viewer, data_x, data_y):
        """
        Find the curve closest to an image position.

        Parameters
        ----------
        viewer : ``ginga.ImageView.ImageViewBase``
            The ginga viewer.
        data_x, data_y : :obj:`float`
            Image coordinates.

        Returns
        -------
        :obj:`int`
            Index of the closest curve, or None if no curve is close enough.
        """
        if data_x is None or data_y is None or not np.isfinite(data_x) \
                or not np.isfinite(data_y):
            return None
        i = np.argmin(np.absolute(self.y[:, 0] - data_y))
        if not (min(self.y[0, 0], self.y[-1, 0]) - 1 <= data_y
                    <= max(self.y[0, 0], self.y[-1, 0]) + 1):
            return None
        dx = np.absolute(self.x[i] - data_x)
        if not np.any(np.isfinite(dx)):
            return None
        k = np.nanargmin(dx)
        sx = viewer.get_scale_xy()[0]
        tol = max(0.5 * self.spacing if np.isfinite(self.spacing) else 0., 5. / sx)
        return int(k) if dx[k] <= tol else None

    def _hover_objects(self, viewer, i0, i1):
        """Objects highlighting the curve under the cursor."""
        k = self.hover
        sx, sy = viewer.get_scale_xy()
        step = max(1, int(np.ceil((i1 - i0 + 1) / self.max_vertices)))
        indx = np.arange(i0, i1 + 1, step)
        pts = np.column_stack([self.x[indx, k], self.y[indx, k]])
        pts = pts[np.all(np.isfinite(pts), axis=1)]
        objs = [] if len(pts) < 2 else [Path(pts, color='cyan', linewidth=2)]
        im = (i0 + i1) // 2
        if np.isfinite(self.x[im, k]) and np.isfinite(self.y[im, k]):
            objs.append(Text(self.x[im, k] + 6 / sx, self.y[im, k], self.labels[k],
                             color='cyan', fontsize=max(self.fontsize, 14.), bgcolor='black',
                             bgalpha=0.7))
        return objs


class MarkerOverlay(_ViewDependentOverlay):
    """
    Draw a (possibly large) set of square markers.

    Only the markers inside the view are drawn, and at most one marker is
    drawn for each patch of the screen about the size of a marker, so that
    the number of markers drawn does not grow when zooming out.  At most
    ``max_markers`` markers are drawn.  The markers are drawn directly with
    the viewer's renderer (not as individual ginga canvas objects), which is
    much faster for thousands of markers.

    Parameters
    ----------
    x : array-like
        Horizontal image coordinates of the markers.
    y : array-like
        Vertical image coordinates of the markers.
    color : :obj:`str`, optional
        Marker color.
    size : :obj:`float`, optional
        Half-size of the markers in image pixels.
    max_markers : :obj:`int`, optional
        Maximum number of markers drawn.
    min_screen_size : :obj:`float`, optional
        Minimum size in screen pixels of the patch of the screen in which only
        one marker is drawn.
    owner : ``ginga.ImageView.ImageViewBase``, optional
        The viewer of the channel with this overlay.  The markers are only
        drawn in this viewer, not in other viewers that draw the same canvas
        (e.g., ginga's Pan plugin).  If None, they are drawn in all viewers.
    """
    def __init__(self, x, y, color='cyan', size=1., max_markers=2000, min_screen_size=3.,
                 owner=None):
        super().__init__()
        self.x = np.asarray(x, dtype=float).ravel()
        self.y = np.asarray(y, dtype=float).ravel()
        self.color = color
        self.size = size
        self.max_markers = max_markers
        self.min_screen_size = min_screen_size
        self.owner = owner
        # Shape used to set the line style for the renderer
        self._style = SquareBox(0., 0., size, color=color, fill=False)

    def select(self, viewer):
        """
        Select the markers to draw for the current view.

        Parameters
        ----------
        viewer : ``ginga.ImageView.ImageViewBase``
            The ginga viewer doing the rendering.

        Returns
        -------
        `numpy.ndarray`_
            Indices of the markers to draw.
        """
        x1, y1, x2, y2 = viewer.get_datarect()
        sx, sy = viewer.get_scale_xy()
        indx = np.where((self.x >= x1) & (self.x <= x2) & (self.y >= y1) & (self.y <= y2))[0]
        if indx.size == 0:
            return indx
        # Keep one marker per screen patch
        patch = max(2. * self.size * min(sx, sy), self.min_screen_size)
        cells = np.column_stack([np.floor((self.x[indx] - x1) * sx / patch),
                                 np.floor((self.y[indx] - y1) * sy / patch)]).astype(np.int64)
        _, first = np.unique(cells, axis=0, return_index=True)
        indx = indx[np.sort(first)]
        if indx.size > self.max_markers:
            indx = indx[::int(np.ceil(indx.size / self.max_markers))]
        return indx

    def build(self, viewer):
        """
        Construct the screen coordinates of the markers to draw.

        Parameters
        ----------
        viewer : ``ginga.ImageView.ImageViewBase``
            The ginga viewer doing the rendering.

        Returns
        -------
        `numpy.ndarray`_
            Array with shape :math:`(N_{\\rm marker}, 4, 2)` with the corners of
            each marker in the renderer's (native) coordinates.
        """
        indx = self.select(viewer)
        if indx.size == 0:
            return np.empty((0, 4, 2))
        s = self.size
        offsets = np.array([[-s, -s], [s, -s], [s, s], [-s, s]])
        corners = np.column_stack([self.x[indx], self.y[indx]])[:, None, :] + offsets[None]
        native = np.asarray(viewer.tform['data_to_native'].to_(corners.reshape(-1, 2)))
        return native[:, :2].reshape(-1, 4, 2)

    def draw(self, viewer):
        """
        Draw the markers.

        Parameters
        ----------
        viewer : ``ginga.ImageView.ImageViewBase``
            The ginga viewer doing the rendering.
        """
        if viewer.get_image() is None or (self.owner is not None and viewer is not self.owner):
            return
        key = _view_key(viewer)
        cached = self._cache.get(id(viewer))
        if cached is None or cached[0] != key:
            corners = self.build(viewer)
            self._cache[id(viewer)] = (key, corners)
        else:
            corners = cached[1]
        if len(corners) == 0:
            return
        cr = viewer.renderer.setup_cr(self._style)
        cr.initialize_from_shape(self._style, line=True, fill=False)
        for cpoints in corners:
            cr.draw_polygon(cpoints, line=cr.line, fill=None)


class PolylineOverlay(_ViewDependentOverlay):
    """
    Draw a (possibly large) set of short polylines (e.g., the fitted tilts of
    the arc lines).

    Only the lines that overlap the view are drawn, each with at most
    ``max_vertices`` vertices.  The lines are drawn directly with the viewer's
    renderer (not as individual ginga canvas objects), which is much faster for
    hundreds of lines.

    Parameters
    ----------
    lines : :obj:`list`
        List of (x, y) tuples with the image coordinates of the vertices of
        each line.
    color : :obj:`str`, optional
        Line color.
    linewidth : :obj:`int`, optional
        Line width in screen pixels.
    max_vertices : :obj:`int`, optional
        Maximum number of vertices drawn for each line.
    owner : ``ginga.ImageView.ImageViewBase``, optional
        The viewer of the channel with this overlay.  The lines are only drawn
        in this viewer, not in other viewers that draw the same canvas (e.g.,
        ginga's Pan plugin).  If None, they are drawn in all viewers.
    """
    def __init__(self, lines, color='blue', linewidth=1, max_vertices=20, owner=None):
        super().__init__()
        self.lines = []
        for x, y in lines:
            x = np.asarray(x, dtype=float).ravel()
            y = np.asarray(y, dtype=float).ravel()
            gpm = np.isfinite(x) & np.isfinite(y)
            if np.sum(gpm) < 2:
                continue
            x, y = x[gpm], y[gpm]
            if x.size > max_vertices:
                indx = np.unique(np.linspace(0, x.size - 1, max_vertices).astype(int))
                x, y = x[indx], y[indx]
            self.lines.append(np.column_stack([x, y]))
        # Bounding boxes used to select the lines in view
        self.bbox = np.array([[ln[:, 0].min(), ln[:, 1].min(), ln[:, 0].max(), ln[:, 1].max()]
                              for ln in self.lines]).reshape(-1, 4)
        self.owner = owner
        # Shape used to set the line style for the renderer
        self._style = Path(np.zeros((2, 2)), color=color, linewidth=linewidth)

    def select(self, viewer):
        """
        Select the lines that overlap the view.

        Parameters
        ----------
        viewer : ``ginga.ImageView.ImageViewBase``
            The ginga viewer doing the rendering.

        Returns
        -------
        `numpy.ndarray`_
            Indices of the lines to draw.
        """
        x1, y1, x2, y2 = viewer.get_datarect()
        b = self.bbox
        return np.where((b[:, 2] >= x1) & (b[:, 0] <= x2) & (b[:, 3] >= y1) & (b[:, 1] <= y2))[0]

    def draw(self, viewer):
        """
        Draw the lines.

        Parameters
        ----------
        viewer : ``ginga.ImageView.ImageViewBase``
            The ginga viewer doing the rendering.
        """
        if viewer.get_image() is None or len(self.lines) == 0 \
                or (self.owner is not None and viewer is not self.owner):
            return
        key = _view_key(viewer)
        cached = self._cache.get(id(viewer))
        if cached is None or cached[0] != key:
            tform = viewer.tform['data_to_native']
            cpoints = [np.asarray(tform.to_(self.lines[i]))[:, :2] for i in self.select(viewer)]
            self._cache[id(viewer)] = (key, cpoints)
        else:
            cpoints = cached[1]
        if len(cpoints) == 0:
            return
        cr = viewer.renderer.setup_cr(self._style)
        cr.initialize_from_shape(self._style, line=True)
        for cpts in cpoints:
            cr.draw_path(cpts, line=cr.line)
