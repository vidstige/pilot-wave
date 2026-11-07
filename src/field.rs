//! The bath's memory: past impacts with decaying weights. No I/O.
//!
//! Length unit is the Faraday wavelength, time unit is one bounce.
//! Every bounce leaves a standing J0 wave that decays by
//! exp(-1/memory) per bounce.
use std::f64::consts::PI;

use crate::bessel::{j0, j1_over_x};

pub type Vec2 = [f64; 2];

pub const WAVENUMBER: f64 = 2.0 * PI;
/// Impacts older than this many memories are forgotten.
const HORIZON: f64 = 6.0;

pub struct Field {
    impacts: Vec<Vec2>,
    weights: Vec<f64>,
    next: usize,
    decay: f64,
}

impl Field {
    pub fn new(memory: f64, walkers: usize) -> Field {
        let n = walkers * (HORIZON * memory) as usize;
        Field {
            impacts: vec![[0.0; 2]; n],
            weights: vec![0.0; n],
            next: 0,
            decay: (-1.0 / memory).exp(),
        }
    }

    pub fn height(&self, p: Vec2) -> f64 {
        self.impacts.iter().zip(&self.weights)
            .map(|(q, w)| w * j0(WAVENUMBER * distance(p, *q)))
            .sum()
    }

    pub fn slope(&self, p: Vec2) -> Vec2 {
        let mut g = [0.0; 2];
        for (q, w) in self.impacts.iter().zip(&self.weights) {
            let d = [p[0] - q[0], p[1] - q[1]];
            let r = (d[0] * d[0] + d[1] * d[1]).sqrt();
            let f = -WAVENUMBER * WAVENUMBER * w
                * j1_over_x(WAVENUMBER * r);
            g[0] += f * d[0];
            g[1] += f * d[1];
        }
        g
    }

    /// Age every wave by one bounce and add this bounce's impacts.
    pub fn remember(&mut self, impacts: &[Vec2]) {
        for w in &mut self.weights {
            *w *= self.decay;
        }
        for &p in impacts {
            self.impacts[self.next] = p;
            self.weights[self.next] = 1.0;
            self.next = (self.next + 1) % self.impacts.len();
        }
    }
}

pub fn distance(a: Vec2, b: Vec2) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    const H: f64 = 1e-5;

    fn scattered_field() -> Field {
        let mut f = Field::new(10.0, 1);
        for i in 0..30 {
            let t = i as f64 * 0.7;
            f.remember(&[[t.sin() * 1.3, t.cos() * 0.8 + t * 0.05]]);
        }
        f
    }

    #[test]
    fn slope_is_gradient_of_height() {
        let f = scattered_field();
        for p in [[0.3, -0.2], [1.1, 0.4], [-0.7, 0.9]] {
            let dx = (f.height([p[0] + H, p[1]]) - f.height([p[0] - H, p[1]]))
                / (2.0 * H);
            let dy = (f.height([p[0], p[1] + H]) - f.height([p[0], p[1] - H]))
                / (2.0 * H);
            let g = f.slope(p);
            assert!((g[0] - dx).abs() < 1e-6 && (g[1] - dy).abs() < 1e-6);
        }
    }

    #[test]
    fn fresh_impact_has_unit_height_and_no_slope() {
        let mut f = Field::new(10.0, 1);
        f.remember(&[[0.5, 0.5]]);
        assert!((f.height([0.5, 0.5]) - 1.0).abs() < 1e-6);
        assert_eq!(f.slope([0.5, 0.5]), [0.0, 0.0]);
    }
}
