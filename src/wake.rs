//! A droplet's wake: its past impacts, each a standing J0 wave whose
//! weight is the decay raised to its age, so nothing has to be aged.
//!
//! Length unit is the Faraday wavelength, time unit is one bounce.
//! The J0 and J1 the wake needs are read from a table built once.
use std::f64::consts::PI;
use std::sync::OnceLock;

use crate::bath::Vec2;
use crate::bessel::{j0, j1_over_x};

pub const WAVENUMBER: f64 = 2.0 * PI;
/// Impacts older than this many memories are forgotten.
const HORIZON: f64 = 6.0;
/// Table spacing and reach in k·r; past the reach the functions run.
const TABLE_STEP: f64 = 0.01;
const TABLE_REACH: f64 = 40.0 * WAVENUMBER;

static TABLE: OnceLock<Vec<(f64, f64)>> = OnceLock::new();

/// J0(x) and J1(x)/x by table lookup.
fn kernel(x: f64) -> (f64, f64) {
    let table = TABLE.get_or_init(|| {
        (0..=(TABLE_REACH / TABLE_STEP) as usize + 1)
            .map(|i| i as f64 * TABLE_STEP)
            .map(|x| (j0(x), j1_over_x(x)))
            .collect()
    });
    if x >= TABLE_REACH {
        return (j0(x), j1_over_x(x));
    }
    let s = x / TABLE_STEP;
    let i = s as usize;
    let t = s - i as f64;
    let (a, b) = (table[i], table[i + 1]);
    (a.0 + (b.0 - a.0) * t, a.1 + (b.1 - a.1) * t)
}

pub struct Wake {
    impacts: Vec<Vec2>,
    size: usize,
    next: usize,
    decay: f64,
}

impl Wake {
    pub fn new(memory: f64) -> Wake {
        let size = (HORIZON * memory) as usize;
        Wake {
            impacts: Vec::with_capacity(size),
            size,
            next: 0,
            decay: (-1.0 / memory).exp(),
        }
    }

    /// Impacts newest first, each with its weight.
    fn weighted(&self) -> impl Iterator<Item = (Vec2, f64)> + '_ {
        let n = self.impacts.len();
        let mut weight = 1.0 / self.decay;
        (0..n).map(move |age| {
            weight *= self.decay;
            (self.impacts[(self.next + n - 1 - age) % n], weight)
        })
    }

    pub fn height(&self, p: Vec2) -> f64 {
        self.weighted()
            .map(|(q, w)| {
                let d = [p[0] - q[0], p[1] - q[1]];
                w * kernel(WAVENUMBER * (d[0] * d[0] + d[1] * d[1]).sqrt()).0
            })
            .sum()
    }

    pub fn slope(&self, p: Vec2) -> Vec2 {
        let mut g = [0.0; 2];
        for (q, w) in self.weighted() {
            let d = [p[0] - q[0], p[1] - q[1]];
            let r = (d[0] * d[0] + d[1] * d[1]).sqrt();
            let f = -WAVENUMBER * WAVENUMBER * w * kernel(WAVENUMBER * r).1;
            g[0] += f * d[0];
            g[1] += f * d[1];
        }
        g
    }

    pub fn remember(&mut self, p: Vec2) {
        if self.impacts.len() < self.size {
            self.impacts.push(p);
        } else {
            self.impacts[self.next] = p;
        }
        self.next = (self.next + 1) % self.size;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const H: f64 = 1e-5;

    /// Height straight from the J0 function, without the table.
    fn exact_height(w: &Wake, p: Vec2) -> f64 {
        w.weighted()
            .map(|(q, weight)| {
                let r = ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2)).sqrt();
                weight * j0(WAVENUMBER * r)
            })
            .sum()
    }

    fn scattered() -> Wake {
        let mut wake = Wake::new(10.0);
        for i in 0..30 {
            let t = i as f64 * 0.7;
            wake.remember([t.sin() * 1.3, t.cos() * 0.8 + t * 0.05]);
        }
        wake
    }

    #[test]
    fn table_matches_the_functions() {
        for x in [0.0, 0.004, 1.0, 7.5, 31.4, 200.0, 300.0] {
            let (a, b) = kernel(x);
            assert!((a - j0(x)).abs() < 1e-4 && (b - j1_over_x(x)).abs() < 1e-4);
        }
    }

    #[test]
    fn slope_is_gradient_of_height() {
        let w = scattered();
        for p in [[0.3, -0.2], [1.1, 0.4], [-0.7, 0.9]] {
            let exact = |x, y| exact_height(&w, [x, y]);
            let dx = (exact(p[0] + H, p[1]) - exact(p[0] - H, p[1])) / (2.0 * H);
            let dy = (exact(p[0], p[1] + H) - exact(p[0], p[1] - H)) / (2.0 * H);
            let g = w.slope(p);
            assert!((g[0] - dx).abs() < 1e-3 && (g[1] - dy).abs() < 1e-3);
        }
    }

    #[test]
    fn fresh_impact_has_unit_height_and_no_slope() {
        let mut w = Wake::new(10.0);
        w.remember([0.5, 0.5]);
        assert!((w.height([0.5, 0.5]) - 1.0).abs() < 1e-4);
        assert_eq!(w.slope([0.5, 0.5]), [0.0, 0.0]);
    }

    /// Height at the origin due to one impact there, aged by the given
    /// number of later impacts far away, with those impacts' own small
    /// contribution subtracted.
    fn aged(bounces: usize) -> f64 {
        let (mut with, mut without) = (Wake::new(10.0), Wake::new(10.0));
        with.remember([0.0, 0.0]);
        for _ in 0..bounces {
            with.remember([1000.0, 1000.0]);
            without.remember([1000.0, 1000.0]);
        }
        with.height([0.0, 0.0]) - without.height([0.0, 0.0])
    }

    #[test]
    fn impacts_fade_by_the_decay_per_bounce_and_are_forgotten_past_the_horizon() {
        assert!((aged(5) - (-0.5f64).exp()).abs() < 1e-6);
        assert!(aged(60).abs() < 1e-6);
    }
}
