"""README figures from a run's output in out/<experiment>."""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = "out"
IMAGES = "images"
DPI = 100
HITS = "#9a9890"
FLUX = "#2a78d6"
THEORY = "#eb6834"
INK = "#2b2a28"

plt.rcParams.update({"figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "font.size": 10})


def crop(image: np.ndarray, extent: float, box: tuple) -> np.ndarray:
    """Cut the world box (x0, x1, y0, y1) out of an image of the domain."""
    n = image.shape[0]
    px = lambda v: int(round((v + extent / 2) / extent * n))
    x0, x1, y0, y1 = box
    return image[n - px(y1):n - px(y0), px(x0):px(x1)]


def show(ax, path: str, extent: float, box: tuple, title: str, aspect="equal"):
    ax.imshow(crop(plt.imread(path), extent, box), extent=box, aspect=aspect)
    ax.set_title(title, fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])


def frame(name: str, index: int) -> str:
    return f"{OUT}/{name}/{name}_{index:04}.png"


def trails(name: str) -> str:
    return f"{OUT}/{name}/{name}_trails.png"


def tidy(ax, xlabel: str, ylabel: str):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def screen(ax, name: str, reach: float):
    d = np.genfromtxt(f"{OUT}/{name}/{name}_screen.csv", delimiter=",", names=True)
    width = d["x"][1] - d["x"][0]
    ax.bar(d["x"], d["hits"], width=width, color=HITS, label="arrivals")
    ax.plot(d["x"], d["flux"], color=FLUX, lw=2, label="probability current")
    ax.plot(d["x"], d["theory"], color=THEORY, lw=1.5, ls="--", label="Fraunhofer")
    ax.set_xlim(-reach, reach)
    ax.set_xlabel("x on the screen (wavelengths)")
    ax.set_ylabel("probability per bin")
    ax.legend(frameon=False, fontsize=9)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def buildup(axes, name: str, counts: list, reach: float):
    """The screen filling up one arrival at a time."""
    a = np.genfromtxt(f"{OUT}/{name}/{name}_arrivals.csv", delimiter=",", names=True)
    x = a["x"][np.argsort(a["t"])]
    rng = np.random.default_rng(1)
    for ax, count in zip(axes, counts):
        shown = x[:count]
        ax.scatter(shown, rng.uniform(0, 1, len(shown)), s=2, color="white", lw=0)
        ax.set_facecolor("black")
        ax.set_xlim(-reach, reach)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_title(f"first {count} arrivals" if count < len(x) else f"all {len(x)} arrivals", fontsize=10)


def slit(name: str) -> plt.Figure:
    """Frames, trails, the build-up and the screen of a slit experiment."""
    extent, reach = 102.4, 40
    figure = plt.figure(figsize=(12, 10), layout="constrained")
    top, middle, bottom = figure.subfigures(3, 1, height_ratios=[1.6, 0.5, 1])
    axes = top.subplots(1, 3)
    show(axes[0], f"{OUT}/{name}/{name}_0008.png", extent, (-20, 20, -25, 15), "t = 4: at the slit")
    show(axes[1], f"{OUT}/{name}/{name}_0016.png", extent, (-30, 30, -25, 35), "t = 8: at the screen")
    show(axes[2], f"{OUT}/{name}/{name}_trails.png", extent, (-30, 30, -25, 35), "300 trails over the time-integrated density")
    buildup(middle.subplots(1, 4), name, [20, 200, 2000, 10**9], reach)
    screen(bottom.subplots(), name, reach)
    return figure


def sanity() -> plt.Figure:
    """Trails of the free packet, the mirror and the barrier."""
    extent = 51.2
    panels = [
        ("free", (-8, 8, -20, 16)),
        ("mirror", (-8, 8, -20, 8)),
        ("barrier", (-8, 8, -16, 10)),
    ]
    figure, axes = plt.subplots(1, 3, figsize=(8, 4.2), layout="constrained")
    for ax, (name, box) in zip(axes, panels):
        show(ax, trails(name), extent, box, f"{name}: 300 trails", aspect="auto")
    return figure


def spin() -> plt.Figure:
    """Trajectories up the channel with the spin in and out of the
    plane, and when the particles arrive."""
    figure = plt.figure(figsize=(12, 7), layout="constrained")
    left, right = figure.subfigures(1, 2, width_ratios=[1, 1.4])
    axes = left.subplots(1, 2, sharey=True)
    for ax, name in zip(axes, ["spin-in-plane", "spin-normal"]):
        r = np.genfromtxt(f"{OUT}/{name}/{name}_trails.csv", delimiter=",", names=True)
        for k in range(8):
            trail = r[r["trail"] == k]
            ax.plot(trail["x"], trail["y"], color=FLUX, lw=0.6)
        ax.axhline(0, color=THEORY, ls="--", lw=1)
        ax.set_xlim(-1.5, 1.5)
        ax.set_ylim(-40.5, 3)
        tidy(ax, "x", "y")
        ax.set_title(f"{name}: 8 trajectories", fontsize=10)
    ax = right.subplots()
    edges = np.arange(0, 120, 1.0)
    for name, color in [("spin-in-plane", HITS), ("spin-normal", FLUX)]:
        a = np.genfromtxt(f"{OUT}/{name}/{name}_arrivals.csv", delimiter=",", names=True)
        ax.hist(a["t"], edges, histtype="step", color=color, lw=1.5, label=name)
    ax.set_yscale("log")
    tidy(ax, "first arrival time at the line", "particles per unit time")
    ax.legend(frameon=False, fontsize=9)
    return figure


FIGURES = {
    "single-slit": lambda: slit("single-slit"),
    "double-slit": lambda: slit("double-slit"),
    "sanity": sanity,
    "spin": spin,
}

if __name__ == "__main__":
    for name in sys.argv[1:] or FIGURES:
        FIGURES[name]().savefig(f"{IMAGES}/{name}.png", dpi=DPI)
