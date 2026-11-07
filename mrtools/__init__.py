"""mrtools: building blocks for designing, simulating and reconstructing MRI pulse sequences (EE 290 project).

* ``mrtools.hardware``: scanner profiles (GE MR750w)
* ``mrtools.seq``: PyPulseq building blocks (EPI readout, fat-sat, slice ordering, SMS pulses and CAIPI), GE checks
* ``mrtools.sim``: KomaMRI phantom with ground truth, simulator wrapper
* ``mrtools.recon``: ESPIRiT coil maps and coil combination

Shared code lives here and changes rarely. Each sequence's own code lives in its scripts:
``sequences/scripts/<name>.py`` (builds the .seq) and ``recon/<name>.py`` (reconstructs it).

The earlier EPTI port lives on the ``epti`` branch; docs/LESSONS_LEARNED.md summarizes what carried over.
"""
