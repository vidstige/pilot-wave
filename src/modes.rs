//! Standing wave modes of a disc with a reflecting rim: Jm(k r) times
//! cos or sin of m·θ, with Jm'(k R) = 0. Jm comes from its integral
//! representation, which the trapezoid rule evaluates to machine
//! precision once the number of points exceeds m + x by a margin.
//! That is slow, so each mode tabulates Jm and Jm' along the radius
//! once and interpolates linearly after that.
use std::f64::consts::{PI, TAU};

use crate::bath::Vec2;

const QUADRATURE: usize = 128;
const SCAN_START: f64 = 0.1;
const SCAN_STEP: f64 = 0.05;
const BISECTIONS: usize = 60;
const TABLE_STEP: f64 = 0.002;

pub fn jm(m: u32, x: f64) -> f64 {
    let n = QUADRATURE as f64;
    (0..QUADRATURE)
        .map(|j| TAU * j as f64 / n)
        .map(|t| (m as f64 * t - x * t.sin()).cos())
        .sum::<f64>()
        / n
}

pub fn jm_prime(m: u32, x: f64) -> f64 {
    if m == 0 { -jm(1, x) } else { (jm(m - 1, x) - jm(m + 1, x)) / 2.0 }
}

fn bisect(m: u32, mut lo: f64, mut hi: f64) -> f64 {
    for _ in 0..BISECTIONS {
        let mid = (lo + hi) / 2.0;
        if jm_prime(m, lo).signum() == jm_prime(m, mid).signum() {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (lo + hi) / 2.0
}

/// Positive zeros of Jm' up to the given bound, in increasing order.
/// The first one lies above m, and below it Jm' sinks under rounding
/// noise, so the scan starts there.
pub fn neumann_zeros(m: u32, upto: f64) -> Vec<f64> {
    let start = SCAN_START.max(m as f64);
    let steps = ((upto - start) / SCAN_STEP).max(0.0) as usize;
    (0..steps)
        .map(|i| start + SCAN_STEP * i as f64)
        .filter(|&x| jm_prime(m, x).signum() != jm_prime(m, x + SCAN_STEP).signum())
        .map(|x| bisect(m, x, x + SCAN_STEP))
        .collect()
}

#[derive(Clone, Debug)]
pub struct Mode {
    pub m: u32,
    pub k: f64,
    pub sine: bool,
    table: Vec<(f64, f64)>,
}

fn polar(p: Vec2) -> (f64, f64) {
    ((p[0] * p[0] + p[1] * p[1]).sqrt().max(1e-9), p[1].atan2(p[0]))
}

impl Mode {
    /// A mode usable out to the given radius.
    pub fn new(m: u32, k: f64, sine: bool, reach: f64) -> Mode {
        let table = (0..=(reach / TABLE_STEP) as usize + 1)
            .map(|i| k * TABLE_STEP * i as f64)
            .map(|x| (jm(m, x), jm_prime(m, x)))
            .collect();
        Mode { m, k, sine, table }
    }

    /// Jm(k r) and Jm'(k r) by interpolation.
    fn radial(&self, r: f64) -> (f64, f64) {
        let x = r / TABLE_STEP;
        let i = (x as usize).min(self.table.len() - 2);
        let t = x - i as f64;
        let (a, b) = (self.table[i], self.table[i + 1]);
        (a.0 + (b.0 - a.0) * t, a.1 + (b.1 - a.1) * t)
    }

    fn angular(&self, theta: f64) -> f64 {
        let a = self.m as f64 * theta;
        if self.sine { a.sin() } else { a.cos() }
    }

    /// Derivative of the angular factor with respect to θ.
    fn angular_prime(&self, theta: f64) -> f64 {
        let a = self.m as f64 * theta;
        self.m as f64 * if self.sine { a.cos() } else { -a.sin() }
    }

    pub fn value(&self, p: Vec2) -> f64 {
        let (r, theta) = polar(p);
        self.radial(r).0 * self.angular(theta)
    }

    /// The value straight from the integral, slow; the reference
    /// that tests hold the table and gradient against.
    #[cfg(test)]
    pub fn exact(&self, p: Vec2) -> f64 {
        let (r, theta) = polar(p);
        jm(self.m, self.k * r) * self.angular(theta)
    }

    pub fn gradient(&self, p: Vec2) -> Vec2 {
        let (r, theta) = polar(p);
        let (j, j_prime) = self.radial(r);
        let radial = self.k * j_prime * self.angular(theta);
        let tangential = j * self.angular_prime(theta) / r;
        [
            radial * theta.cos() - tangential * theta.sin(),
            radial * theta.sin() + tangential * theta.cos(),
        ]
    }

    /// Integral of the squared mode over the disc, in closed form for
    /// Neumann zeros.
    pub fn norm(&self, radius: f64) -> f64 {
        let angular = if self.m == 0 { TAU } else { PI };
        let kr = self.k * radius;
        let ratio = self.m as f64 / kr;
        angular * radius * radius / 2.0 * (1.0 - ratio * ratio) * jm(self.m, kr).powi(2)
    }
}

/// Every mode of the disc with wavenumber in [lo, hi], tabulated out
/// to the given reach.
pub fn modes(radius: f64, lo: f64, hi: f64, reach: f64) -> Vec<Mode> {
    let mut all = Vec::new();
    for m in 0.. {
        let zeros = neumann_zeros(m, hi * radius);
        if zeros.is_empty() {
            break;
        }
        for k in zeros.iter().map(|z| z / radius).filter(|&k| k >= lo) {
            all.push(Mode::new(m, k, false, reach));
            if m > 0 {
                all.push(Mode::new(m, k, true, reach));
            }
        }
    }
    all
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bessel::{j0, j1};

    const H: f64 = 1e-5;

    #[test]
    fn integral_matches_rational_j0_and_j1() {
        for x in [0.0, 0.3, 2.5, 7.0, 15.0, 24.0] {
            assert!((jm(0, x) - j0(x)).abs() < 1e-7);
            assert!((jm(1, x) - j1(x)).abs() < 1e-7);
        }
    }

    #[test]
    fn neumann_zeros_are_where_the_derivative_vanishes() {
        let zeros = neumann_zeros(3, 20.0);
        assert_eq!(zeros.len(), 5);
        assert!((zeros[0] - 4.201189).abs() < 1e-5);
        assert!(zeros.iter().all(|&z| jm_prime(3, z).abs() < 1e-9));
    }

    #[test]
    fn high_order_zeros_start_above_the_order() {
        let zeros = neumann_zeros(20, 40.0);
        assert!(zeros[0] > 20.0 && zeros[0] < 25.0, "{zeros:?}");
        assert!(zeros.iter().all(|&z| jm_prime(20, z).abs() < 1e-9));
    }

    #[test]
    fn tabulated_values_match_the_integral() {
        let mode = Mode::new(4, 5.3, false, 2.0);
        for p in [[0.0, 0.0], [0.3, 0.1], [-0.7, 0.2], [1.5, -0.9]] {
            assert!((mode.value(p) - mode.exact(p)).abs() < 1e-4);
        }
    }

    #[test]
    fn mode_gradient_matches_finite_differences() {
        let mode = Mode::new(4, 5.3, true, 2.0);
        for p in [[0.4, 0.2], [-0.7, 0.9], [0.0, -1.1]] {
            let dx = (mode.exact([p[0] + H, p[1]]) - mode.exact([p[0] - H, p[1]])) / (2.0 * H);
            let dy = (mode.exact([p[0], p[1] + H]) - mode.exact([p[0], p[1] - H])) / (2.0 * H);
            let g = mode.gradient(p);
            assert!((g[0] - dx).abs() < 1e-4 && (g[1] - dy).abs() < 1e-4);
        }
    }

    #[test]
    fn closed_form_norm_matches_quadrature() {
        let radius = 2.0;
        let mode = Mode::new(2, neumann_zeros(2, 10.0)[1] / radius, false, 3.0);
        let n = 200;
        let cell = 2.0 * radius / n as f64;
        let numeric: f64 = (0..n * n)
            .map(|i| [(i % n) as f64 + 0.5, (i / n) as f64 + 0.5])
            .map(|c| [-radius + c[0] * cell, -radius + c[1] * cell])
            .filter(|p| p[0] * p[0] + p[1] * p[1] <= radius * radius)
            .map(|p| mode.value(p).powi(2) * cell * cell)
            .sum();
        assert!((numeric - mode.norm(radius)).abs() < 0.01 * mode.norm(radius));
    }

    #[test]
    fn modes_cover_the_band_and_pair_sines_with_cosines() {
        let band = modes(2.0, 6.0, 7.0, 2.0);
        assert!(band.iter().all(|m| m.k >= 6.0 && m.k <= 7.0));
        let cosines = band.iter().filter(|m| !m.sine).count();
        let sines = band.iter().filter(|m| m.sine).count();
        let zero_order = band.iter().filter(|m| m.m == 0).count();
        assert_eq!(cosines, sines + zero_order);
    }
}
