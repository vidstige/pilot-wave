"""Two droplets started close together.

Usage: python pair.py [memory] [distance]

A small tangential nudge sets them orbiting. At memory 8 they lock
into orbit at a distance of 0.73 or 1.69 wavelengths depending on
where they start; try distances 0.7 and 1.6. Started from rest they
instead bind as a static pair at 0.58 wavelengths.
"""
import sys

import matplotlib.pyplot as plt

from pilot_wave import field, free, walk
from render import animate

memory = float(sys.argv[1]) if len(sys.argv) > 1 else 8
distance = float(sys.argv[2]) if len(sys.argv) > 2 else 0.7
f = field(memory, walkers=2)
start = [[-distance / 2, 0], [distance / 2, 0]]
nudge = [[0, 0.02], [0, -0.02]]
walkers = walk(f, start, nudge, free)
animation = animate(f, walkers, half=2.5, follow=True)
plt.show()
