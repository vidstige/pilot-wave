"""README figures from the runner's output: frame montages, trails and
screen plots. Run the experiments into out/<name> first (see README)."""
from __future__ import annotations

import os
import re
import sys
from typing import Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

# Python 3.8 has no builtin generics at runtime.
Box = Tuple[float, float, float, float]

OUT = "out"
IMAGES = "images"
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
MUTED = "#52514e"
HITS = "#9a9890"
FLUX = "#2a78d6"
THEORY = "#eb6834"
CLASSICAL = "#1baf7a"
DPI = 100
# Figures larger than this in bytes get a 256-colour palette.
LIMIT = 400_000
# Domain side of the 1024 and 256 grids at dx = 0.1.
LARGE = 102.4
SMALL = 25.6
# Distance from the slits to the screen, for the single slit's zeros.
SCREEN = 29.8
SINGLE_SLIT_WIDTH = 2.0
# The entangled pair's line, the widths of its centre and separation,
# particle 2's wall, and the line view's colours of the two particles.
LONG = 51.2
EPR_CENTRE = 1.0
EPR_GAP = 0.6
EPR_WALL = 5.0
ONE = "#ff8c32"
TWO = "#46aaff"

plt.rcParams.update({
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "axes.edgecolor": MUTED,
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "font.size": 10,
})


def crop(image: np.ndarray, extent: float, box: Box) -> np.ndarray:
    """Cut the world box (x0, x1, y0, y1) out of an image of the domain."""
    n = image.shape[0]
    px = lambda v: int(round((v + extent / 2) / extent * n))
    x0, x1, y0, y1 = box
    return image[n - px(y1):n - px(y0), px(x0):px(x1)]


def show(ax, path: str, extent: float, box: Box, title: str):
    ax.imshow(crop(plt.imread(path), extent, box), extent=box)
    ax.set_title(title, fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])


def frames(name: str, run: str, indices: list[int], extent: float, box: Box) -> plt.Figure:
    x0, x1, y0, y1 = box
    width = 3.3 * (x1 - x0) / (y1 - y0) + 0.3
    figure, axes = plt.subplots(1, len(indices) + 1, figsize=(width * (len(indices) + 1), 3.6))
    for ax, i in zip(axes, indices):
        show(ax, f"{OUT}/{run}/{name}_{i:04}.png", extent, box, f"t ≈ {i / 2:g}")
    show(axes[-1], f"{OUT}/{run}/{name}_trajectories.png", extent, box, "trails")
    figure.tight_layout()
    return figure


def screen(run: str, name: str) -> dict[str, np.ndarray]:
    table = np.genfromtxt(f"{OUT}/{run}/{name}_screen.csv", delimiter=",", names=True)
    return {key: table[key] for key in table.dtype.names}


def screen_plot(ax, name: str, run: str, theory: str, reach: float):
    """Arrivals of one run against the quantum current and the theory."""
    quantum = screen(name, name)
    x = quantum["x"]
    kind = "classical" if run.endswith("classical") else "quantum"
    ax.fill_between(x, screen(run, name)["hits"], step="mid", lw=0, alpha=0.6,
                    color=CLASSICAL if kind == "classical" else HITS,
                    label=f"{kind} arrivals")
    ax.plot(x, quantum["flux"], color=FLUX, lw=2, label="probability current")
    ax.plot(x, quantum["theory"], color=THEORY, lw=1.5, ls="--", label=theory)
    ax.set_xlim(-reach, reach)
    ax.set_ylim(0, None)
    ax.set_xlabel("x on the screen (wavelengths)")
    ax.set_ylabel("probability per bin")
    ax.legend(frameon=False, fontsize=9)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def comparison(name: str, theory: str, box: Box, reach: float) -> plt.Figure:
    figure = plt.figure(figsize=(10, 12))
    grid = figure.add_gridspec(3, 2, height_ratios=[1.3, 1, 1])
    runs = [name, f"{name}-classical"]
    for column, run in enumerate(runs):
        title = run.replace(name, "quantum").replace("quantum-", "")
        show(figure.add_subplot(grid[0, column]), f"{OUT}/{run}/{name}_trajectories.png",
             LARGE, box, f"{title} trails over the time-integrated density")
    for row, run in enumerate(runs, 1):
        screen_plot(figure.add_subplot(grid[row, :]), name, run, theory, reach)
    figure.tight_layout()
    return figure


def single_slit() -> plt.Figure:
    figure, (left, right) = plt.subplots(1, 2, figsize=(12, 5), width_ratios=[1, 1.4])
    show(left, f"{OUT}/single-slit/single-slit_trajectories.png", LARGE, (-25, 25, -10, 40),
         "trails over the time-integrated density")
    screen_plot(right, "single-slit", "single-slit", "Fraunhofer", 40)
    zero = SCREEN * np.tan(np.arcsin(1 / SINGLE_SLIT_WIDTH))
    for x in (-zero, zero):
        right.axvline(x, color=MUTED, lw=0.8, ls=":")
    figure.tight_layout()
    return figure


def harmonic() -> plt.Figure:
    figure, axes = plt.subplots(1, 5, figsize=(15, 3.3))
    box = (-7, 7, -7, 7)
    for ax, i in zip(axes, [6, 12, 19, 25]):
        show(ax, f"{OUT}/harmonic/harmonic_{i:04}.png", SMALL, box, f"t ≈ {i / 2:g}")
    show(axes[-1], f"{OUT}/harmonic-few/harmonic_trajectories.png", SMALL, box, "12 trails")
    figure.tight_layout()
    return figure


def ground() -> plt.Figure:
    figure, axes = plt.subplots(1, 3, figsize=(9, 3.3))
    box = (-4, 4, -4, 4)
    for ax, i in zip(axes, [1, 50]):
        show(ax, f"{OUT}/ground/ground_{i:04}.png", SMALL, box, f"t ≈ {i / 2:g}")
    show(axes[-1], f"{OUT}/ground-few/ground_trajectories.png", SMALL, box, "60 trails")
    figure.tight_layout()
    return figure


def vortex() -> plt.Figure:
    """The trap's m = +1, m = -1 and p orbital states, each with 40
    trails over two periods."""
    figure, axes = plt.subplots(1, 3, figsize=(9, 3.3))
    box = (-4, 4, -4, 4)
    titles = ["m = +1: counterclockwise", "m = -1: clockwise", "p orbital: at rest"]
    for ax, name, title in zip(axes, ["vortex", "antivortex", "orbital"], titles):
        show(ax, f"{OUT}/{name}-few/{name}_trajectories.png", SMALL, box, title)
    figure.tight_layout()
    return figure


def numbers(path: str) -> list[float]:
    """Every decimal number in a run's stdout summary, in order."""
    return [float(w.strip("(),")) for w in open(path).read().split()
            if w.strip("(),").replace(".", "", 1).isdigit() and "." in w]


def barrier() -> plt.Figure:
    figure, axes = plt.subplots(1, 4, figsize=(14, 4), width_ratios=[1, 1, 1, 2.2])
    box = (-12, 12, -24, 18)
    for ax, i in zip(axes, [5, 7, 9]):
        show(ax, f"{OUT}/barrier-few/barrier_{i:04}.png", LARGE, box, f"t ≈ {i / 2:g}")
    particles, field, analytic = numbers(f"{OUT}/barrier/summary.txt")[:3]
    ax = axes[-1]
    # The packet starts at y = -20 with width 2; the barrier is ahead.
    y = np.linspace(-27, -13, 400)
    density = np.exp(-(y + 20) ** 2 / 8) / np.sqrt(8 * np.pi)
    cumulative = np.cumsum(density[::-1])[::-1] * (y[1] - y[0])
    split = y[np.argmin(abs(cumulative - particles))]
    ax.plot(y, density, color=INK, lw=1.5)
    ax.fill_between(y, density, where=y >= split, color=HITS, alpha=0.6, lw=0,
                    label=f"particles that got through: {particles:.3f}")
    for value, color, label in [(field, FLUX, "field"), (analytic, THEORY, "analytic")]:
        at = y[np.argmin(abs(cumulative - value))]
        ax.axvline(at, color=color, lw=1.5, ls="--",
                   label=f"front {value:.3f} of the packet ({label})")
    ax.set_xlabel("starting y (the barrier is at y = 0)")
    ax.set_ylabel("starting density")
    ax.set_title("the transmitted particles are the front of the packet", fontsize=10)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    figure.tight_layout()
    return figure


def pair_frame(ax, path: str, reach: tuple[float, float], title: str, extent: float = SMALL):
    """Configuration space over the line view, both cut to x in reach."""
    image = plt.imread(path)
    n = image.shape[1]
    field, line = image[:n], image[n:]
    px = lambda v: int(round((v + extent / 2) / extent * n))
    a, b = map(px, reach)
    ax.imshow(np.vstack([field[n - b:n - a, a:b], line[:, a:b]]))
    ax.set_title(title, fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])


def trails(ax, run: str, reach: tuple[float, float], title: str):
    x0, x1 = reach
    show(ax, f"{OUT}/{run}/{run}_trajectories.png", SMALL, (x0, x1, x0, x1), title)
    ax.set_xlabel("x₁")
    ax.set_ylabel("x₂")


def fermions() -> plt.Figure:
    figure = plt.figure(figsize=(11, 8.4), layout="constrained")
    top, bottom = figure.subfigures(2, 1, height_ratios=[1.25, 1])
    reach = (-8, 8)
    for ax, i in zip(top.subplots(1, 4), [0, 3, 4, 6]):
        pair_frame(ax, f"{OUT}/fermions/fermions_{i:04}.png", reach, f"fermions, t = {i / 2:g}")
    for ax, run in zip(bottom.subplots(1, 3), ["distinguishable", "bosons", "fermions"]):
        trails(ax, run, reach, f"{run}: 300 runs")
    return figure


def separations(path: str) -> np.ndarray:
    """Rows of time, mean and spread of the separation in a summary."""
    rows = [line.split("t =")[1] for line in open(path) if " t = " in line]
    return np.array([[float(w.strip(",")) for w in r.split() if w[0].isdigit()] for r in rows])


def molecule() -> plt.Figure:
    figure = plt.figure(figsize=(11, 11.5), layout="constrained")
    top, middle, bottom = figure.subfigures(3, 1, height_ratios=[1.2, 1.2, 0.8])
    for ax, i in zip(top.subplots(1, 4), [0, 12, 15, 26]):
        path = f"{OUT}/molecule/molecule_{i:04}.png"
        pair_frame(ax, path, (-9, 6), f"molecule, t = {i / 2:g}")
    for ax, i in zip(middle.subplots(1, 4), [4, 6, 8, 12]):
        path = f"{OUT}/molecule-break/molecule-break_{i:04}.png"
        pair_frame(ax, path, (-12.8, 12.8), f"kicked harder, t = {i / 2:g}")
    ax = bottom.subplots()
    for run, color in [("molecule", FLUX), ("molecule-break", THEORY)]:
        t, mean, spread = separations(f"{OUT}/{run}-many/summary.txt").T
        ax.plot(t, mean, color=color, lw=2, label=f"{run}: mean |x₁ − x₂|")
        ax.fill_between(t, mean - spread, mean + spread, color=color, alpha=0.2, lw=0)
    ax.set_xlabel("t")
    ax.set_ylabel("separation over 20000 runs")
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
def flux_trap() -> plt.Figure:
    figure, axes = plt.subplots(1, 3, figsize=(9, 3.3))
    box = (-4, 4, -4, 4)
    for ax, i in zip(axes, [1, 50]):
        show(ax, f"{OUT}/flux-trap/flux-trap_{i:04}.png", SMALL, box, f"t ≈ {i / 2:g}")
    show(axes[-1], f"{OUT}/flux-trap-few/flux-trap_trajectories.png", SMALL, box, "60 trails")
    figure.tight_layout()
    return figure


def flux_slit() -> plt.Figure:
    figure = plt.figure(figsize=(10, 7.5))
    grid = figure.add_gridspec(2, 2, height_ratios=[1.3, 1])
    box = (-25, 25, -10, 40)
    show(figure.add_subplot(grid[0, 0]), f"{OUT}/flux-slit-few/flux-slit_0012.png", LARGE,
         (-15, 15, -5, 25), "t ≈ 6, the cut along the dark middle")
    show(figure.add_subplot(grid[0, 1]), f"{OUT}/flux-slit/flux-slit_trajectories.png", LARGE,
         box, "trails over the time-integrated density")
    ax = figure.add_subplot(grid[1, :])
    screen_plot(ax, "flux-slit", "flux-slit", "Fraunhofer, half a fringe over", 40)
    plain = screen("double-slit", "double-slit")
    ax.plot(plain["x"], plain["flux"], color=MUTED, lw=1, ls=":",
            label="current without the flux line")
    ax.legend(frameon=False, fontsize=9)
    figure.tight_layout()
    return figure


def table(path: str) -> dict[str, np.ndarray]:
    data = np.genfromtxt(path, delimiter=",", names=True)
    return {key: data[key] for key in data.dtype.names}


def correlations(path: str) -> np.ndarray:
    """Rows of time and the correlation over the runs in a summary."""
    pattern = re.compile(r"t = +([\d.]+), correlation (-?[\d.]+) over the runs")
    return np.array([[float(g) for g in m.groups()] for m in map(pattern.search, open(path)) if m])


def closed_form(t: np.ndarray) -> np.ndarray:
    """Correlation of x₁ and x₂ for the free Gaussian pair."""
    centre = EPR_CENTRE ** 2 + t ** 2 / (16 * EPR_CENTRE ** 2)
    half_gap = (EPR_GAP ** 2 + t ** 2 / EPR_GAP ** 2) / 4
    return (centre - half_gap) / (centre + half_gap)


def tidy(ax, xlabel: str, ylabel: str):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(frameon=False, fontsize=8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def paths(ax, run: str, follow: int):
    """x₁ and x₂ of a few runs over time, with the wall and without."""
    data = table(f"{OUT}/{run}/{run}_paired.csv")
    for k in range(follow):
        row = data["run"] == k
        t = data["t"][row]
        for axis, sub, color in [("1", "₁", ONE), ("2", "₂", TWO)]:
            free = dict(label=f"x{sub} without the wall") if k == 0 else {}
            wall = dict(label=f"x{sub} with the wall") if k == 0 else {}
            ax.plot(t, data[f"x{axis}_free"][row], color=color, lw=1, ls="--", alpha=0.6, **free)
            ax.plot(t, data[f"x{axis}_wall"][row], color=color, lw=1.4, **wall)
    ax.axhline(EPR_WALL, color=MUTED, lw=0.8, ls="--")
    ax.set_title(f"{follow} runs, particle 2's wall at x = {EPR_WALL:g}", fontsize=10)
    tidy(ax, "t", "position")


def marginals(ax, run: str):
    """The densities of x₁ and x₂ at the end, with the wall and without."""
    field = table(f"{OUT}/{run}/{run}_marginals.csv")
    runs = table(f"{OUT}/{run}/{run}_histograms.csv")
    for axis, color in [("one", ONE), ("two", TWO)]:
        n = axis.replace("one", "x₁").replace("two", "x₂")
        ax.plot(field["x"], field[f"{axis}_free"], color=color, lw=2.5, alpha=0.4,
                label=f"{n} without the wall")
        ax.plot(field["x"], field[f"{axis}_wall"], color=color, lw=1, label=f"{n} with it")
        ax.step(runs["x"], runs[f"{axis}_wall"], where="mid", color=color, lw=0.6, alpha=0.7)
    ax.set_xlim(-22, 22)
    ax.set_title("densities at t = 3.5, runs as steps", fontsize=10)
    tidy(ax, "x", "density")


def entangled() -> plt.Figure:
    figure = plt.figure(figsize=(12, 8.6), layout="constrained")
    top, bottom = figure.subfigures(2, 1, height_ratios=[1.25, 1])
    low, high = -18, 8
    for ax, i in zip(top.subplots(1, 4), [0, 3, 5, 7]):
        path = f"{OUT}/entangled-wall-few/entangled-wall_{i:04}.png"
        pair_frame(ax, path, (low, high), f"with the wall, t = {i / 2:g}", LONG)
        # Particle 2's wall: x₂ = EPR_WALL across configuration space.
        height = ax.get_images()[0].get_array().shape[1]
        ax.axhline((high - EPR_WALL) / (high - low) * height, color=TWO, lw=0.8, ls="--")
    left, middle, right = bottom.subplots(1, 3)
    paths(left, "entangled-wall", 8)
    marginals(middle, "entangled-wall")
    t = np.linspace(0, 3.5, 100)
    right.plot(t, closed_form(t), color=INK, lw=1, label="closed form")
    for run, color in [("entangled", THEORY), ("unentangled", FLUX)]:
        t, rho = correlations(f"{OUT}/{run}/summary.txt").T
        right.plot(t, rho, "o", color=color, label=f"{run}, 20000 runs")
    right.axhline(0, color=MUTED, lw=0.5)
    right.set_title("correlation of x₁ and x₂, no wall", fontsize=10)
    tidy(right, "t", "correlation")
    return figure


FIGURES = {
    "double-slit": lambda: comparison("double-slit", "Fraunhofer", (-36, 36, -25, 45), 40),
    "disc": lambda: comparison("disc", "Fresnel (Kirchhoff)", (-30, 30, -25, 40), 25),
    "single-slit": single_slit,
    "barrier": barrier,
    "harmonic": harmonic,
    "ground": ground,
    "vortex": vortex,
    "box": lambda: frames("box", "box", [1, 3, 5, 8], SMALL, (-6.5, 6.5, -6.5, 6.5)),
    "free": lambda: frames("free", "free", [1, 4, 8, 12], LARGE, (-12, 12, -24, 24)),
    "fermions": fermions,
    "molecule": molecule,
    "flux-trap": flux_trap,
    "flux-slit": flux_slit,
    "entangled": entangled,
}

def save(figure: plt.Figure, path: str):
    figure.savefig(path, dpi=DPI)
    if os.path.getsize(path) > LIMIT:
        Image.open(path).convert("RGB").quantize(256).save(path, optimize=True)


if __name__ == "__main__":
    for name in sys.argv[1:] or FIGURES:
        save(FIGURES[name](), f"{IMAGES}/{name}.png")
