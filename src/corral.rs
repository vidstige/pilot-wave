//! A circular corral: the bath's memory lives in the disc's standing
//! modes near the Faraday wavenumber, so the rim is built into the
//! waves. No I/O.
//!
//! An impact at p adds gain·φ(p) to each mode amplitude, where the
//! gain is a Gaussian window on the mode's wavenumber divided by the
//! mode's norm. The window is scaled so that the mean over the disc of
//! the height an impact leaves under itself is one, matching J0(0) in
//! the unbounded bath, so the same FORCE applies.
use std::f64::consts::PI;

use crate::bath::{Bath, Vec2};
use crate::field::WAVENUMBER;
use crate::modes::{modes, Mode};

/// Width of the Faraday resonance in wavenumber.
pub const BANDWIDTH: f64 = 0.5;
/// Modes further than this many bandwidths off resonance are dropped.
const REACH: f64 = 3.0;
/// Modes are tabulated a little beyond the rim so interpolation
/// stays in range there.
const OVERHANG: f64 = 1.01;

pub struct Corral {
    radius: f64,
    modes: Vec<Mode>,
    gains: Vec<f64>,
    amplitudes: Vec<f64>,
    decay: f64,
}

fn window(k: f64) -> f64 {
    (-0.5 * ((k - WAVENUMBER) / BANDWIDTH).powi(2)).exp()
}

impl Corral {
    pub fn new(radius: f64, memory: f64) -> Corral {
        let reach = REACH * BANDWIDTH;
        let band = (WAVENUMBER - reach, WAVENUMBER + reach);
        let modes = modes(radius, band.0, band.1, OVERHANG * radius);
        let windows: Vec<f64> = modes.iter().map(|m| window(m.k)).collect();
        let scale = PI * radius * radius / windows.iter().sum::<f64>();
        let gains = modes.iter().zip(&windows)
            .map(|(m, w)| scale * w / m.norm(radius))
            .collect();
        Corral {
            radius,
            amplitudes: vec![0.0; modes.len()],
            modes,
            gains,
            decay: (-1.0 / memory).exp(),
        }
    }
}

impl Bath for Corral {
    /// Zero beyond the rim, where the bath is flat.
    fn height(&self, p: Vec2) -> f64 {
        if p[0] * p[0] + p[1] * p[1] > self.radius * self.radius {
            return 0.0;
        }
        self.modes.iter().zip(&self.amplitudes)
            .map(|(m, a)| a * m.value(p))
            .sum()
    }

    fn slope(&self, p: Vec2) -> Vec2 {
        let mut g = [0.0; 2];
        for (m, a) in self.modes.iter().zip(&self.amplitudes) {
            let d = m.gradient(p);
            g[0] += a * d[0];
            g[1] += a * d[1];
        }
        g
    }

    fn remember(&mut self, impacts: &[Vec2]) {
        for a in &mut self.amplitudes {
            *a *= self.decay;
        }
        for &p in impacts {
            for ((m, gain), a) in self.modes.iter().zip(&self.gains).zip(&mut self.amplitudes) {
                *a += gain * m.value(p);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const H: f64 = 1e-5;
    const RADIUS: f64 = 2.0;

    fn disc_points(n: usize) -> Vec<Vec2> {
        let cell = 2.0 * RADIUS / n as f64;
        (0..n * n)
            .map(|i| [(i % n) as f64 + 0.5, (i / n) as f64 + 0.5])
            .map(|c| [-RADIUS + c[0] * cell, -RADIUS + c[1] * cell])
            .filter(|p| p[0] * p[0] + p[1] * p[1] < RADIUS * RADIUS)
            .collect()
    }

    fn exact_height(c: &Corral, p: Vec2) -> f64 {
        c.modes.iter().zip(&c.amplitudes).map(|(m, a)| a * m.exact(p)).sum()
    }

    fn stirred() -> Corral {
        let mut c = Corral::new(RADIUS, 10.0);
        for i in 0..20 {
            let t = i as f64 * 0.9;
            c.remember(&[[t.cos() * 0.1 * i as f64, t.sin() * 0.08 * i as f64]]);
        }
        c
    }

    #[test]
    fn mean_height_under_a_fresh_impact_is_one() {
        let points = disc_points(24);
        let mut forgetful = Corral::new(RADIUS, 1e-3);
        let mean: f64 = points.iter().map(|&p| {
            forgetful.remember(&[p]);
            forgetful.height(p)
        }).sum::<f64>() / points.len() as f64;
        assert!((mean - 1.0).abs() < 0.05, "{mean}");
    }

    #[test]
    fn slope_is_gradient_of_height() {
        let c = stirred();
        let scale = c.height([0.0; 2]).abs().max(1.0);
        for p in [[0.3, -0.2], [1.1, 0.4], [-0.7, 0.9]] {
            let exact = |x, y| exact_height(&c, [x, y]);
            let dx = (exact(p[0] + H, p[1]) - exact(p[0] - H, p[1])) / (2.0 * H);
            let dy = (exact(p[0], p[1] + H) - exact(p[0], p[1] - H)) / (2.0 * H);
            let g = c.slope(p);
            assert!((g[0] - dx).abs() < 1e-3 * scale && (g[1] - dy).abs() < 1e-3 * scale);
        }
    }

    #[test]
    fn rim_has_no_radial_slope() {
        let c = stirred();
        for i in 0..12 {
            let t = i as f64 * 0.5;
            let g = c.slope([RADIUS * t.cos(), RADIUS * t.sin()]);
            let radial = g[0] * t.cos() + g[1] * t.sin();
            assert!(radial.abs() < 1e-3 * c.height([0.0; 2]).abs().max(1.0), "{radial}");
        }
    }
}
