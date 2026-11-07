"""Long-time position statistics of a droplet in a harmonic trap.

Usage: python trap.py [memory] [bounces]

Around memory 20 the droplet settles on a precessing orbit; at
memory 50 and above it wanders chaotically and the histogram builds
up a ring pattern on the Faraday wavelength.
"""
import sys
from itertools import islice

import numpy as np
import matplotlib.pyplot as plt

from pilot_wave import field, spring, walk
from render import statistics

STIFFNESS = 0.002
SKIP = 5000

memory = float(sys.argv[1]) if len(sys.argv) > 1 else 50
bounces = int(sys.argv[2]) if len(sys.argv) > 2 else 100_000
f = field(memory, walkers=1)
walkers = walk(f, [[0.3, 0]], [[0, 0]], spring(STIFFNESS))
trajectory = np.array([p[0] for p in islice(walkers, bounces)])
statistics(trajectory[SKIP:], half=3)
plt.show()
