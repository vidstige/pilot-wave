from itertools import islice

import numpy as np
import pytest

from bessel import j0, j1
from pilot_wave import field, free, height, remember, slope, walk

H = 1e-5


def trajectory(memory, start, bounces, nudge=None):
    f = field(memory, len(start))
    nudge = np.zeros_like(start) if nudge is None else nudge
    walkers = walk(f, start, nudge, free)
    return np.array(list(islice(walkers, bounces)))


def speed(trajectory):
    return np.linalg.norm(np.diff(trajectory, axis=0), axis=-1)


def test_j0_derivative_is_minus_j1():
    x = np.linspace(-30, 30, 2001)
    dj0 = (j0(x + H) - j0(x - H)) / (2 * H)
    assert np.abs(dj0 + j1(x)).max() < 1e-6


def test_j0_vanishes_at_its_first_zero():
    assert j0(np.array([2.404825557695773])) == pytest.approx(0, abs=1e-7)


def test_slope_is_gradient_of_height():
    rng = np.random.default_rng(1)
    f = field(10, 1)
    for p in rng.normal(size=(30, 2)):
        remember(f, p[None])
    points = rng.normal(size=(5, 2))
    dx = (height(f, points + [H, 0]) - height(f, points - [H, 0])) / (2 * H)
    dy = (height(f, points + [0, H]) - height(f, points - [0, H])) / (2 * H)
    gradient = np.stack([dx, dy], axis=1)
    assert slope(f, points) == pytest.approx(gradient, abs=1e-6)


def test_short_memory_bouncer_stays_put():
    t = trajectory(3, [[0, 0]], 1500)
    assert np.linalg.norm(t[-1]) < 0.05


def test_long_memory_bouncer_walks_at_steady_speed():
    t = trajectory(10, [[0, 0]], 1500)
    v = speed(t)[-300:]
    assert v.mean() > 0.02
    assert v.std() < 0.1 * v.mean()


def test_close_pair_locks_into_orbit():
    start = [[-0.35, 0], [0.35, 0]]
    nudge = [[0, 0.02], [0, -0.02]]
    t = trajectory(8, start, 4000, nudge)
    separation = np.linalg.norm(t[:, 0] - t[:, 1], axis=1)[-1000:]
    assert separation.std() < 0.01 * separation.mean()
    relative = t[:, 1] - t[:, 0]
    angle = np.unwrap(np.arctan2(relative[:, 1], relative[:, 0]))
    assert abs(angle[-1] - angle[0]) > 4 * np.pi
