"""Bessel functions J0 and J1 for NumPy arrays.

Rational approximations from Abramowitz & Stegun 9.4, accurate to
about 1e-8, so no SciPy dependency is needed.
"""
import numpy as np

SWITCH = 8.0


def _asymptotic(x, phase, p, q):
    z = SWITCH / x
    y = z * z
    xx = x - phase
    p = np.polyval(p, y)
    q = np.polyval(q, y)
    return np.sqrt(2 / (np.pi * x)) * (np.cos(xx) * p - z * np.sin(xx) * q)


def j0(x: np.ndarray) -> np.ndarray:
    ax = np.abs(x)
    y = ax * ax
    p = np.polyval(
        [-184.9052456, 77392.33017, -11214424.18, 651619640.7,
         -13362590354.0, 57568490574.0], y)
    q = np.polyval(
        [1.0, 267.8532712, 59272.64853, 9494680.718,
         1029532985.0, 57568490411.0], y)
    big = _asymptotic(
        np.maximum(ax, SWITCH), 0.785398164,
        [0.2093887211e-6, -0.2073370639e-5, 0.2734510407e-4,
         -0.1098628627e-2, 1.0],
        [-0.934935152e-7, 0.7621095161e-6, -0.6911147651e-5,
         0.1430488765e-3, -0.1562499995e-1])
    return np.where(ax < SWITCH, p / q, big)


def j1(x: np.ndarray) -> np.ndarray:
    ax = np.abs(x)
    y = ax * ax
    p = ax * np.polyval(
        [-30.16036606, 15704.48260, -2972611.439, 242396853.1,
         -7895059235.0, 72362614232.0], y)
    q = np.polyval(
        [1.0, 376.9991397, 99447.43394, 18583304.74,
         2300535178.0, 144725228442.0], y)
    big = _asymptotic(
        np.maximum(ax, SWITCH), 2.356194491,
        [-0.240337019e-6, 0.2457520174e-5, -0.3516396496e-4,
         0.183105e-2, 1.0],
        [0.105787412e-6, -0.88228987e-6, 0.8449199096e-5,
         -0.2002690873e-3, 0.04687499995])
    return np.sign(x) * np.where(ax < SWITCH, p / q, big)


def j1_over_x(x: np.ndarray) -> np.ndarray:
    """J1(x)/x, continuous through x = 0 where it equals 1/2."""
    small = x < 1e-6
    return np.where(small, 0.5, j1(x) / np.where(small, 1, x))
