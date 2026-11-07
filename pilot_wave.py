"""Stroboscopic path-memory model of walking droplets. No I/O.

Length unit is the Faraday wavelength, time unit is one bounce. Every
bounce leaves a standing J0 wave that decays over `memory` bounces;
the droplet slides down the slope of the summed wave field and
loses a fixed fraction of its momentum per bounce.

A bouncer at rest sits on top of its own wave hill. Once
WAVENUMBER**2 * FORCE * memory**2 > 2 * DRAG that hill is steep
enough that the slightest nudge grows into steady walking.
"""
from dataclasses import dataclass
from typing import Callable

import numpy as np

from bessel import j0, j1_over_x

WAVENUMBER = 2 * np.pi
FORCE = 6e-4
DRAG = 0.3
NOISE = 1e-4
HORIZON = 6  # impacts older than this many memories are forgotten

Positions = np.ndarray  # (n, 2)
Forces = np.ndarray  # (n, 2)
External = Callable[[Positions], Forces]


@dataclass
class Field:
    """The bath's memory: past impacts with decaying weights."""
    impacts: Positions
    weights: np.ndarray
    decay: float


def field(memory: float, walkers: int) -> Field:
    n = walkers * int(HORIZON * memory)
    return Field(np.zeros((n, 2)), np.zeros(n), np.exp(-1 / memory))


def height(f: Field, points: Positions) -> np.ndarray:
    r = np.linalg.norm(points[:, None] - f.impacts, axis=-1)
    return j0(WAVENUMBER * r) @ f.weights


def slope(f: Field, points: Positions) -> Forces:
    d = points[:, None] - f.impacts
    r = np.linalg.norm(d, axis=-1)
    w = f.weights * j1_over_x(WAVENUMBER * r)
    return -WAVENUMBER**2 * np.einsum("pm,pmk->pk", w, d)


def remember(f: Field, impacts: Positions) -> None:
    n = len(impacts)
    f.impacts = np.roll(f.impacts, n, axis=0)
    f.weights = np.roll(f.weights * f.decay, n)
    f.impacts[:n] = impacts
    f.weights[:n] = 1


def step(f: Field, positions: Positions, velocities: Positions,
         external: External, rng: np.random.Generator):
    push = -FORCE * slope(f, positions) + external(positions)
    kick = NOISE * rng.standard_normal(positions.shape)
    velocities = (1 - DRAG) * velocities + push + kick
    positions = positions + velocities
    remember(f, positions)
    return positions, velocities


def walk(f: Field, positions, velocities, external: External,
         seed: int = 0):
    """Yield walker positions after each bounce, forever."""
    rng = np.random.default_rng(seed)
    positions = np.array(positions, float)
    velocities = np.array(velocities, float)
    while True:
        positions, velocities = step(
            f, positions, velocities, external, rng)
        yield positions


def free(positions: Positions) -> Forces:
    return np.zeros_like(positions)


def spring(stiffness: float) -> External:
    return lambda positions: -stiffness * positions
