.. highlight:: rest

***************
Reduction Files
***************

Several file are generated in preparing to run the
full reduction and when running PYPIT.  These
are distinsuished by extension:


=========== ===========  ====== ===========================================
Extension   File Name    Format Description
=========== ===========  ====== ===========================================
.pypit      PYPIT        ASCII  Primary file that guides the data reduction
 ..          ..                 We *recommend* one per instrument configuration.
 ..          ..                 See :doc:`pypit_file` for further details.
.setup      setup        YAML   Lists the various setups of the instrument.
 ..          ..                 See :doc:`setups` for further details
.group      group        YAML   Lists the various setup groups
 ..          ..                 See :ref:`groupings` for further details
.lst        list         ASCII  Lists all of the raw files to be analyzed
 ..          ..                 In particular, frametype is given
.xml        xml list     xml    Also lists all of the raw files to be analyzed
.spect      spect        ASCII  A record of every spectrograph setting used
 ..          ..                 during the data reduction.
.settings   settings     ASCII  A record of every PYPIT reduction setting used
 ..          ..                 during the data reduction.
.log        log          ASCII  Log file of the reduction messages
=========== ===========  ====== ===========================================

