.. include:: include/links.rst

.. _bspline:

================
B-Spline Fitting
================

PypeIt uses weighted least-squares B-spline fitting throughout its reduction
pipeline, primarily for flat-field modelling and sky subtraction.  This page
describes the underlying algorithm, the public API of
:class:`~pypeit.bspline.refactor.BSpline` and
:class:`~pypeit.bspline.refactor.BSpline2D`, and how to migrate from the
legacy :class:`~pypeit.bspline.bspline.bspline` class and
:func:`~pypeit.core.fitting.bspline_profile` function.

.. note::

    The implementation of these B-spline routines draws heavily from legacy
    software produced for the Sloan Digital Sky Survey in IDL [idlutils]_ and
    then rewritten in python [pydl]_.  This most recent version is a refactor
    facilitated by
    `Claude Code from Anthropic <https://www.anthropic.com/claude-code>`__.

.. _bspline-algorithm:

Algorithm
=========

1D B-Spline
-----------

A B-spline of order :math:`n` over the knot vector
:math:`t_1 \le t_2 \le \cdots \le t_T` is a piecewise polynomial of degree
:math:`n - 1` with :math:`K = T - n` basis functions :math:`B_k(x)`.  Each
:math:`B_k` has compact support on :math:`[t_k, t_{k+n}]` and is evaluated
via the de Boor--Cox recurrence [Cox1972]_ [deBoor1972]_.  For :math:`N`
observations :math:`(x_i, y_i)` with inverse-variance weights
:math:`w_i = \varepsilon_i^{-2}`, the model is

.. math::

    \hat{y}(x) = \sum_{k=1}^{K} c_k\, B_k(x).

The weighted least-squares problem

.. math::

    \min_{\mathbf{c}}\; \sum_{i=1}^{N} w_i \bigl[y_i - \hat{y}(x_i)\bigr]^2

leads to the normal equations

.. math::

    \bigl(\mathbf{A}^\top \mathbf{W} \mathbf{A}\bigr)\,\mathbf{c}
    \;=\; \mathbf{A}^\top \mathbf{W}\,\mathbf{y},

where :math:`A_{ik} = B_k(x_i)` is the design matrix and
:math:`\mathbf{W} = \operatorname{diag}(w_i)`.  Because each :math:`B_k` is
non-zero on at most :math:`n` consecutive knot spans,
:math:`\mathbf{A}^\top \mathbf{W} \mathbf{A}` is symmetric positive
semi-definite and banded with half-bandwidth :math:`n`.  The system is solved
via :func:`scipy.linalg.cholesky_banded` and
:func:`scipy.linalg.cho_solve_banded`.

Quasi-2D Extension
------------------

:class:`~pypeit.bspline.refactor.BSpline2D` extends the 1D model by
introducing a polynomial modulation in a second variable :math:`u`:

.. math::

    \hat{y}(x, u) = \sum_{k=1}^{K} \sum_{j=1}^{p} c_{kj}\,
                    B_k(x)\, P_j(u),

where :math:`P_j` are orthogonal polynomial basis functions (Legendre,
Chebyshev, or simple power-law) evaluated at :math:`u` after normalisation to
:math:`[-1, 1]`, and :math:`p` is the number of polynomial terms.  The
combined design matrix is

.. math::

    \tilde{A}_{i,\,(k-1)p+j} = B_k(x_i)\, P_j(u_i).

Because the column grouping preserves the banded structure of
:math:`\mathbf{A}` in the knot index, the normal equations remain banded with
half-bandwidth :math:`np` and are again solved via
:func:`scipy.linalg.cholesky_banded` and
:func:`scipy.linalg.cho_solve_banded`.

Iterative Sigma-Clipping
------------------------

:func:`~pypeit.bspline.refactor.bspline_profile_refactor` wraps either class
in an iterative outlier-rejection loop.  After each fit the normalised
residuals

.. math::

    r_i = \frac{y_i - \hat{y}(x_i)}{\varepsilon_i}

are tested against upper and lower clipping thresholds.  Points that exceed
either threshold are assigned zero weight for the next iteration.  The loop
repeats until the mask converges (no new rejections) or the maximum number of
iterations is reached.

.. _bspline-usage:

Usage
=====

.. _bspline-1d:

BSpline (1D)
------------

Construct a :class:`~pypeit.bspline.refactor.BSpline` with a
:class:`~pypeit.bspline.refactor.Knots` specification and call
:meth:`~pypeit.bspline.refactor.BSpline.fit`.  The knot spacing is controlled
by the ``spacing``, ``count``, or ``stride`` arguments to
:class:`~pypeit.bspline.refactor.Knots`.

.. code-block:: python

    import numpy as np
    from pypeit.bspline.refactor import BSpline, Knots

    rng  = np.random.default_rng(0)
    x    = np.sort(rng.uniform(0, 10, 500))
    y    = np.sin(x) + rng.normal(0, 0.05, 500)
    ivar = np.full(500, 400.0)

    bspl       = BSpline(x=x, knots=Knots(spacing=0.5), nord=4)
    err, yfit  = bspl.fit(x, y, ivar=ivar)

    # Evaluate at arbitrary positions
    x_new          = np.linspace(0, 10, 200)
    y_new, gpm_new = bspl.value(x_new)

The return code ``err`` is 0 on success, -1 if breakpoints were masked
(requiring another call), and -2 on failure.

.. _bspline-2d:

BSpline2D (quasi-2D)
--------------------

Pass ``basis`` (a polynomial family name) and ``basis_x`` (the second
variable array) to :meth:`~pypeit.bspline.refactor.BSpline2D.fit`.

.. code-block:: python

    import numpy as np
    from pypeit.bspline.refactor import BSpline2D, Knots

    rng     = np.random.default_rng(1)
    x       = np.sort(rng.uniform(0, 10, 500))
    basis_x = rng.uniform(-1, 1, 500)
    y       = (1 + 0.4 * basis_x) * np.sin(x) + rng.normal(0, 0.05, 500)
    ivar    = np.full(500, 400.0)

    bspl      = BSpline2D(x=x, knots=Knots(spacing=0.5), nord=4)
    err, yfit = bspl.fit(
        x, y, ivar=ivar,
        basis='legendre', basis_x=basis_x, npoly=2,
        xmin=-1.0, xmax=1.0,
    )

    # Evaluate at new positions — pass basis_x to recompute the polynomial basis
    x_new    = np.linspace(0, 10, 200)
    u_new    = np.linspace(-1, 1, 200)
    y_new, _ = bspl.value(x_new, basis_x=u_new)

When the fit was performed with a **pre-built array** rather than a string
family name, :meth:`~pypeit.bspline.refactor.BSpline2D.value` cannot
reconstruct the polynomial basis from ``basis_x`` alone (no family name is
stored).  In that case a corresponding evaluation basis of shape
``(M, npoly)`` must be passed explicitly as the ``basis`` argument:

.. code-block:: python

    from pypeit.core.basis import flegendre

    P_train = flegendre(basis_x, npoly)   # used during fit
    bspl.fit(x, y, ivar=ivar, basis=P_train)

    # Evaluate at new positions — supply a matching evaluation basis
    u_new  = np.linspace(-1, 1, 200)
    P_eval = flegendre(u_new, npoly)      # shape (200, npoly)
    y_new, _ = bspl.value(x_new, basis=P_eval)

.. _bspline-profile-refactor:

bspline_profile_refactor
------------------------

:func:`~pypeit.bspline.refactor.bspline_profile_refactor` provides the
iterative sigma-clipping loop around either class.  The fitting class is
selected automatically: omitting ``basis`` (or passing ``None``) uses the 1D
:class:`~pypeit.bspline.refactor.BSpline`; providing ``basis`` uses the
quasi-2D :class:`~pypeit.bspline.refactor.BSpline2D`.

**1D fit:**

.. code-block:: python

    from pypeit.bspline.refactor import bspline_profile_refactor

    bspl, outmask, yfit, reduced_chi, exit_status = bspline_profile_refactor(
        x, y, ivar=ivar,
        nord=4, upper=5.0, lower=5.0,
        kwargs_knots={'spacing': 0.5},
    )

**Quasi-2D fit with a string basis:**

.. code-block:: python

    bspl, outmask, yfit, reduced_chi, exit_status = bspline_profile_refactor(
        x, y, ivar=ivar, gpm=gpm,
        nord=4, basis='legendre', basis_x=basis_x, npoly=3,
        upper=5.0, lower=5.0,
        kwargs_knots={'spacing': 1.0},
        kwargs_reject={'groupbadpix': True, 'maxrej': 5},
    )

**Quasi-2D fit with a pre-built basis array:**

.. code-block:: python

    from pypeit.core.basis import flegendre

    profile_basis = flegendre(basis_x, npoly)   # shape (N, npoly)
    bspl, outmask, yfit, reduced_chi, exit_status = bspline_profile_refactor(
        x, y, ivar=ivar,
        nord=4, basis=profile_basis,
        upper=5.0, lower=5.0,
        kwargs_knots={'spacing': 1.0},
    )

The five return values are:

- ``bspl`` — the fitted :class:`~pypeit.bspline.refactor.BSpline` or
  :class:`~pypeit.bspline.refactor.BSpline2D` object.
- ``outmask`` — boolean array indicating the final good-pixel mask.
- ``yfit`` — best-fit model values at the input ``x``.
- ``reduced_chi`` — reduced :math:`\chi^2` of the final fit.
- ``exit_status`` — integer convergence code: 0 (converged normally),
  1 (maximum iterations reached), 2 (too few good points during
  iteration), 3 (degenerate or singular fit), 4 (fewer good points than
  ``nord`` on entry).

.. _bspline-migration:

Migration from Legacy Code
==========================

The sections below map the legacy
:class:`~pypeit.bspline.bspline.bspline` /
:func:`~pypeit.core.fitting.bspline_profile` API to the new classes and
function.

Class Instantiation
-------------------

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Legacy
     - New
   * - ``bspl = bspline.bspline(x, bkspace=s)``
     - ``bspl = BSpline(x=x, knots=Knots(spacing=s))``
   * - ``bspl = bspline.bspline(x, nbkpts=n)``
     - ``bspl = BSpline(x=x, knots=Knots(count=n))``
   * - ``bspl = bspline.bspline(x, everyn=n)``
     - ``bspl = BSpline(x=x, knots=Knots(stride=n))``
   * - ``bspl = bspline.bspline(x, fullbkpt=t)``
     - ``bspl = BSpline(knots=Knots(full=t))``

Fitting
-------

The legacy ``bspline.fit`` accepts the raw second variable as ``x2``.  The
new :meth:`~pypeit.bspline.refactor.BSpline2D.fit` accepts either a string
family name (and ``basis_x``) or a pre-built polynomial basis matrix.

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Legacy
     - New
   * - ``bspl.fit(x, y, invvar)``
     - ``bspl = BSpline(); bspl.fit(x, y, ivar=invvar)``
   * - ``bspl.fit(x, y, invvar, x2=x2)``
     - ``bspl = BSpline2D(); bspl.fit(x, y, ivar=invvar, basis='legendre', basis_x=x2)``

Evaluation
----------

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Legacy
     - New
   * - ``yfit, mask = bspl.value(x)``
     - ``yfit, mask = bspl.value(x)``
   * - ``yfit, mask = bspl.value(x, x2=x2)``
     - ``yfit, mask = bspl.value(x, basis_x=x2)``

Profile Fitting Function
------------------------

The legacy :func:`~pypeit.core.fitting.bspline_profile` requires positional
arguments for both the inverse variance and the polynomial basis (``invvar`` and
``profile_basis``); it also uses ``kwargs_bspline={'bkspace': s}`` for knot
spacing.  The new :func:`~pypeit.bspline.refactor.bspline_profile_refactor`
makes both optional (omit it for a 1D fit), renames ``invvar`` to ``ivar``
and ``ingpm`` to ``gpm``, and uses ``kwargs_knots={'spacing': s}``.

.. code-block:: python

    # Legacy — 2D case
    bspl, outmask, yfit, rchi, status = bspline_profile(
        xdata, ydata, invvar, profile_basis,
        ingpm=gpm, upper=5, lower=5, nord=4,
        kwargs_bspline={'bkspace': s},
        kwargs_reject=kw,
    )

    # New — 2D case
    bspl, outmask, yfit, rchi, status = bspline_profile_refactor(
        xdata, ydata, ivar=invvar, gpm=gpm,
        basis=profile_basis,
        upper=5, lower=5, nord=4,
        kwargs_knots={'spacing': s},
        kwargs_reject=kw,
    )

    # New — 1D case (legacy passed profile_basis=np.ones((N, 1)))
    bspl, outmask, yfit, rchi, status = bspline_profile_refactor(
        xdata, ydata, ivar=invvar, gpm=gpm,
        upper=5, lower=5, nord=4,
        kwargs_knots={'spacing': s},
        kwargs_reject=kw,
    )

The ``quiet`` parameter of :func:`~pypeit.core.fitting.bspline_profile` has
no equivalent; :func:`~pypeit.bspline.refactor.bspline_profile_refactor`
never prints to standard output.

API Reference
-------------

- :class:`pypeit.bspline.refactor.BSpline`
- :class:`pypeit.bspline.refactor.BSpline2D`
- :class:`pypeit.bspline.refactor.Knots`
- :func:`pypeit.bspline.refactor.bspline_profile_refactor`

References
==========

.. [Cox1972] Cox, M. G. (1972). "The numerical evaluation of B-splines."
   *Journal of the Institute of Mathematics and its Applications*, 10,
   134--149.

.. [deBoor1972] de Boor, C. (1972). "On calculating with B-splines."
   *Journal of Approximation Theory*, 6(1), 50--62.
   https://doi.org/10.1016/0021-9045(72)90080-9

.. [idlutils] https://www.sdss4.org/dr17/software/idlutils/

.. [pydl] https://doi.org/10.5281/zenodo.1095150
