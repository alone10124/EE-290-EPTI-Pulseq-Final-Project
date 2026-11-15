"""mrtools: building blocks for designing, simulating and reconstructing MRI pulse sequences (EE 290 project).

* ``mrtools.config``: reads a run's TOML config (``configs/<run>.toml``) and the file-naming rules
* ``mrtools.hardware``: ``HardwareProfile`` (scanner limits, rasters, receive coils; also builds the KomaMRI scanner)
* ``mrtools.seq``: PyPulseq building blocks (EPI readout, fat-sat, slice ordering, TRID labels), hardware checks
* ``mrtools.sim``: KomaMRI phantom with ground truth, simulator wrapper
* ``mrtools.recon``: ESPIRiT coil maps and coil combination

Shared code lives here and changes rarely. Each sequence's own code lives in its scripts:
``sequences/scripts/<name>.py`` (builds the .seq) and ``recon/<name>.py`` (reconstructs it). Scanner values live in
the config files, not here.

The earlier EPTI port (incl. SMS pulses, ramp-sampling regridding, PNS and forbidden-ESP checks) is on the ``epti``
branch.
"""
