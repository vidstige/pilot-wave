"""Matplotlib views of the wave field and walker trajectories."""
from itertools import islice

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from pilot_wave import Field, Positions, height

GRID = 120
TRAIL = 400
BOUNCES_PER_FRAME = 4
COLORS = ["white", "orange"]


def height_grid(f: Field, center, half: float) -> np.ndarray:
    axis = np.linspace(-half, half, GRID)
    x, y = np.meshgrid(axis + center[0], axis + center[1])
    points = np.column_stack([x.ravel(), y.ravel()])
    return height(f, points).reshape(GRID, GRID)


def animate(f: Field, walkers, half: float, follow: bool):
    """Live view of the field with walkers and their trails.

    The field colormap is re-scaled every frame so the pattern stays
    visible as the memory fills up. With follow, the view tracks the
    walkers' centroid.
    """
    positions = next(walkers)
    trails = [[p.copy()] for p in positions]
    fig, ax = plt.subplots(figsize=(6, 6))
    image = ax.imshow(
        np.zeros((GRID, GRID)), cmap="RdBu_r", origin="lower")
    lines = [ax.plot([], [], color=c, lw=0.8)[0] for c in COLORS]
    dots = ax.scatter(
        positions[:, 0], positions[:, 1], c=COLORS[:len(positions)],
        edgecolors="black", zorder=3)
    ax.set_xticks([])
    ax.set_yticks([])

    def frame(_):
        for positions in islice(walkers, BOUNCES_PER_FRAME):
            for trail, p in zip(trails, positions):
                trail.append(p.copy())
        center = positions.mean(axis=0) if follow else np.zeros(2)
        h = height_grid(f, center, half)
        image.set_data(h)
        image.set_clim(-abs(h).max(), abs(h).max())
        extent = (center[0] - half, center[0] + half,
                  center[1] - half, center[1] + half)
        image.set_extent(extent)
        ax.set_xlim(extent[:2])
        ax.set_ylim(extent[2:])
        for line, trail in zip(lines, trails):
            del trail[:-TRAIL]
            line.set_data(*np.array(trail).T)
        dots.set_offsets(positions)
        return [image, dots, *lines]

    return FuncAnimation(fig, frame, interval=30, cache_frame_data=False)


def statistics(trajectory: Positions, half: float):
    """Trajectory tail, 2D position histogram and radial density."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    tail = trajectory[-3000:]
    axes[0].plot(tail[:, 0], tail[:, 1], lw=0.5)
    axes[0].set_aspect("equal")
    axes[0].set_title("last 3000 bounces")
    box = [[-half, half], [-half, half]]
    axes[1].hist2d(
        trajectory[:, 0], trajectory[:, 1], bins=100, range=box)
    axes[1].set_aspect("equal")
    axes[1].set_title("position histogram")
    r = np.linalg.norm(trajectory, axis=1)
    counts, edges = np.histogram(r, bins=60, range=(0, half))
    area = np.pi * np.diff(edges**2)
    axes[2].plot(edges[:-1], counts / area, drawstyle="steps-post")
    axes[2].set_xlabel("radius (Faraday wavelengths)")
    axes[2].set_title("radial density")
    fig.tight_layout()
    return fig
