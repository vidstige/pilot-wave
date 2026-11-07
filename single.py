"""A single droplet started at rest. Usage: python single.py [memory]"""
import sys

import matplotlib.pyplot as plt

from pilot_wave import field, free, walk
from render import animate

memory = float(sys.argv[1]) if len(sys.argv) > 1 else 20
f = field(memory, walkers=1)
walkers = walk(f, [[0, 0]], [[0, 0]], free)
animation = animate(f, walkers, half=3, follow=True)
plt.show()
