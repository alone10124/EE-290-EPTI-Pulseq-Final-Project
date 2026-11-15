"""KomaMRI simulation: a numerical brain phantom with known ground truth, and a wrapper around the simulator.

Importing this package sets the juliacall environment variables KomaMRI needs, but doesn't start Julia; that
happens on the first call that needs it. They only take effect if set before Julia starts (the first
``import juliacall``/``komamripy``), so import ``mrtools.sim`` first; variables you set yourself win:

* ``PYTHON_JULIACALL_THREADS=auto``: Julia (and so KomaMRI's simulation loop) uses all CPU cores, not one.
* ``PYTHON_JULIACALL_HANDLE_SIGNALS=yes``: required with more than one Julia thread, or Julia can crash.
"""

import os

os.environ.setdefault("PYTHON_JULIACALL_THREADS", "auto")
os.environ.setdefault("PYTHON_JULIACALL_HANDLE_SIGNALS", "yes")
