
    - ``calib`` assigns each frame to one or more :ref:`calibration-groups`.
      Calibration frames with the same calibration group number will be used to
      reduce a science frame with that calibration group number.  Importantly,
      calibration frames (e.g., biases) can be part of multiple calibrations
      groups, but each science frame must be assigned to *only one* calibration
      group.  Calibration groups should form a running sequence from
      :math:`1...N` for :math:`N` calibration groups, where :math:`N\leq 63`.
      The value can also be set to ``all`` meaning the frame is part of *all*
      calibration groups.

    - ``comb_id`` represents a combination ID assigned to each science frame.
      Frames with the same value of ``comb_id`` will be combined. Note that this
      is an unweighted co-add (and hence may not be necessarily be "optimal" in
      terms of S/N ratio).  The ``comb_id`` must be a single integer, but the
      integers can be anything.  Science frames that are combined together can
      have the same ``calib`` value if they use the same set of calibrations.
      For the calibration frames, ``comb_id`` is irrelevant and its value should
      be set to ``-1``.

    - ``bkg_id`` represents a combination ID assigned to each frame that will be
      *used as a background image*.  Frames with the same value ``bkg_id`` will
      be combined.  Note that this is an unweighted co-add (and hence may not be
      necessarily be "optimal" in terms of S/N ratio).  The ``bkd_id`` must be a
      single integer and the integers can be anything.  However, the integer
      must match one of the provided ``comb_id`` values, and this pairs the
      background images with the associated science images.
