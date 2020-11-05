********
MetaData
********

.. index:: Metadata

Overview
========

PypeIt will slurp from the raw data frames a set of
meta data for performing the reduction.  These are
held in an astropy Table within the PypeItMetaData class.

The following doc is mainly for developers but we begin
with a point or two for the user.

User
====

The meta data is primarily ingested via code but the
user will over-ride any such data via the PypeIt file.
So tread carefully.

When perform a reduction run (with run_pypeit), the code
will require that all required meta data be ingested.
There are cases, however, when the header has been mucked
(e.g. Keck).  It is best to fix the file and proceed, but
if you are feeling adventursome, you may set the
`ignore_image_headers` parameter in the `rdx` block, e.g.::

    [rdx]
       ignore_bad_headers = True

Buyer beware in this case..

Developers
==========

Data Model
++++++++++

The data model for the meta data is specified in a series
of methods at the bottom of the metadata.py module.  Here
are a few examples::

    core_meta['target'] = dict(dtype=str, comment='Name of the target')
    core_meta['binning'] = dict(dtype=str, comment='(spatial,spectral) binning')

The key specifies the name, there is a data type, and a comment.
If the meta data is a `float` and is intended to be used for
defining a configuration, then an `rtol` key should be present, e.g.::

    additional_meta['dispangle'] = dict(dtype=float, comment='Angle of the disperser', rtol=0.)

The core meta are required for the code to run.  The additional meta
are used primarily for configurations.

Ingesting Meta
++++++++++++++

Meta data is read from the `Spectrograph` classes.
The main method is `get_meta_value()` in the parent class.

The children specify how a specific meta data item is to
be ingested.  The standard method is via the header and
one specifies the extension of the header array (typically 0)
and the header card, e.g.::

        meta['dec'] = dict(ext=0, card='DEC')
        meta['target'] = dict(ext=0, card='OBJECT')

If an item of required meta data is not provided by the data file
and has a set value, one sets it by `default`::

        meta['binning'] = dict(ext=0, card=None, default='1,1')

If the meta is only required for a set of frametypes (e.g. VLT which does
not provide RA/DEC for calibrations), you may specify which frametypes
*require* it with::

        meta['dec'] = dict(ext=0, card='DEC', required_ftypes=['science'])

For meta data that depends on more than one header card or has some
other complexity, the `compound_meta()` method is used.  Here is
an example from Keck/DEIMOS::

    def compound_meta(self, headarr, meta_key):
        """

        Args:
            headarr: list
            meta_key: str

        Returns:
            value

        """
        if meta_key == 'binning':
            binspatial, binspec = parse.parse_binning(headarr[0]['BINNING'])
            binning = parse.binning2string(binspatial, binspec)
            return binning
        elif meta_key == 'dispangle':
            if headarr[0]['GRATEPOS'] == 3:
                return headarr[0]['G3TLTWAV']
            elif headarr[0]['GRATEPOS'] == 4:
                return headarr[0]['G4TLTWAV']
            else:
                debugger.set_trace()
        else:
            msgs.error("Not ready for this compound meta")



