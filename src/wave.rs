//! The guiding field: a wavefunction on a square grid stepped by the
//! Schrödinger equation. No I/O.
//!
//! Units: ħ = m = 1, so a wavenumber k is also a velocity and the
//! kinetic energy is k²/2. Real and imaginary parts are stepped
//! alternately, the imaginary part living half a step ahead (Visscher
//! 1991), which keeps the scheme explicit and stable while
//! dt < 2 / (max kinetic + max potential). Hard walls are cells kept
//! at zero; an absorbing sponge at the edges swallows what reaches it.
use std::thread;

use crate::grid::Vec2;

pub type Index = (usize, usize);

/// Fourth order Laplacian stencil along one axis.
const STENCIL: [f32; 5] = [-1.0 / 12.0, 4.0 / 3.0, -5.0 / 2.0, 4.0 / 3.0, -1.0 / 12.0];
/// Sponge depth at each edge, in cells, and its strength.
const SPONGE: usize = 24;
const SPONGE_RATE: f32 = 40.0;

pub struct Wave {
    pub n: usize,
    pub dx: f64,
    pub re: Vec<f32>,
    pub im: Vec<f32>,
    pub potential: Vec<f32>,
    /// Per-cell survival factor applied each step: 0 in walls, below
    /// 1 in the sponge, 1 elsewhere.
    pub keep: Vec<f32>,
    /// Steps taken so far times dt.
    pub time: f64,
    scratch: Vec<f32>,
}

/// Largest stable step for the given cell size and potential.
pub fn stable_dt(dx: f64, potential_max: f64) -> f64 {
    let kinetic_max = 2.0 * 16.0 / 3.0 / (dx * dx) / 2.0;
    0.5 * 2.0 / (kinetic_max + potential_max)
}

impl Wave {
    pub fn new(n: usize, dx: f64) -> Wave {
        Wave {
            n,
            dx,
            re: vec![0.0; n * n],
            im: vec![0.0; n * n],
            potential: vec![0.0; n * n],
            keep: vec![1.0; n * n],
            time: 0.0,
            scratch: vec![0.0; n * n],
        }
    }

    /// Side length of the domain, centred on the origin.
    pub fn extent(&self) -> f64 {
        self.n as f64 * self.dx
    }

    pub fn position(&self, (i, j): Index) -> Vec2 {
        let half = self.extent() / 2.0;
        [(i as f64 + 0.5) * self.dx - half, (j as f64 + 0.5) * self.dx - half]
    }

    /// Fractional cell coordinates of a position.
    pub fn cell(&self, p: Vec2) -> (f64, f64) {
        let half = self.extent() / 2.0;
        ((p[0] + half) / self.dx - 0.5, (p[1] + half) / self.dx - 0.5)
    }

    pub fn index(&self, (i, j): Index) -> usize {
        j * self.n + i
    }

    /// Sponge factors near the edges, left at 1 for a hard wall edge.
    pub fn absorb_edges(&mut self, dt: f64) {
        let n = self.n;
        for j in 0..n {
            for i in 0..n {
                let depth = i.min(j).min(n - 1 - i).min(n - 1 - j);
                if depth < SPONGE {
                    let s = (SPONGE - depth) as f32 / SPONGE as f32;
                    let rate = SPONGE_RATE * s * s;
                    self.keep[j * n + i] = (-rate * dt as f32).exp();
                }
            }
        }
    }

    pub fn density(&self, k: usize) -> f32 {
        self.re[k] * self.re[k] + self.im[k] * self.im[k]
    }

    /// Total probability, which stays near 1 until the sponge or a wall
    /// takes some.
    pub fn norm(&self) -> f64 {
        let sum: f64 = (0..self.n * self.n).map(|k| self.density(k) as f64).sum();
        sum * self.dx * self.dx
    }

    /// Phase gradient at a cell from the phase difference between
    /// its neighbours, which is exact for a plane wave and bounded by
    /// π/(2dx) however close to a node of ψ the cell lies.
    fn flow(&self, k: usize) -> Vec2 {
        let arg = |a: usize, b: usize| {
            let (ra, ia, rb, ib) = (self.re[a], self.im[a], self.re[b], self.im[b]);
            (ia * rb - ra * ib).atan2(ra * rb + ia * ib) as f64
        };
        let n = self.n;
        [arg(k + 1, k - 1) / (2.0 * self.dx), arg(k + n, k - n) / (2.0 * self.dx)]
    }

    /// Guiding velocity ∇ arg ψ at a position, interpolated between
    /// the surrounding cells.
    pub fn velocity(&self, p: Vec2) -> Vec2 {
        let (x, y) = self.cell(p);
        let top = (self.n - 3) as f64;
        let (x, y) = (x.clamp(1.0, top), y.clamp(1.0, top));
        let (i, j) = (x as usize, y as usize);
        let (fx, fy) = (x - i as f64, y - j as f64);
        let corners = [
            (self.index((i, j)), (1.0 - fx) * (1.0 - fy)),
            (self.index((i + 1, j)), fx * (1.0 - fy)),
            (self.index((i, j + 1)), (1.0 - fx) * fy),
            (self.index((i + 1, j + 1)), fx * fy),
        ];
        corners.iter().fold([0.0; 2], |v, &(k, w)| {
            let f = self.flow(k);
            [v[0] + w * f[0], v[1] + w * f[1]]
        })
    }

    /// Hamiltonian applied to one component at every interior cell,
    /// rows split among the cores.
    fn hamiltonian(n: usize, dx: f64, potential: &[f32], f: &[f32], out: &mut [f32]) {
        let inv = 1.0 / (dx * dx) as f32;
        let threads = thread::available_parallelism().map_or(1, |t| t.get());
        let band = n.div_ceil(threads);
        let mut rows: Vec<(usize, &mut [f32])> = out.chunks_mut(n).enumerate().collect();
        thread::scope(|s| {
            for chunk in rows.chunks_mut(band) {
                s.spawn(move || {
                    for (j, row) in chunk.iter_mut() {
                        if *j < 2 || *j >= n - 2 {
                            continue;
                        }
                        for i in 2..n - 2 {
                            let k = *j * n + i;
                            let lap: f32 = STENCIL.iter().zip(-2i64..)
                                .map(|(c, o)| {
                                    let (a, b) = (k as i64 + o, k as i64 + o * n as i64);
                                    c * (f[a as usize] + f[b as usize])
                                })
                                .sum();
                            row[i] = -0.5 * lap * inv + potential[k] * f[k];
                        }
                    }
                });
            }
        });
    }

    pub fn step(&mut self, dt: f64) {
        let dt = dt as f32;
        let h = &mut self.scratch;
        Wave::hamiltonian(self.n, self.dx, &self.potential, &self.im, h);
        for (r, (h, keep)) in self.re.iter_mut().zip(h.iter().zip(&self.keep)) {
            *r = (*r + dt * h) * keep;
        }
        Wave::hamiltonian(self.n, self.dx, &self.potential, &self.re, h);
        for (i, (h, keep)) in self.im.iter_mut().zip(h.iter().zip(&self.keep)) {
            *i = (*i - dt * h) * keep;
        }
        self.time += dt as f64;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::packet::Packet;

    const DX: f64 = 0.1;

    fn packet() -> Packet {
        Packet { center: [-2.0, 0.0], width: [1.0, 1.0], wavevector: [4.0, 0.0] }
    }

    fn free_wave() -> Wave {
        let mut wave = Wave::new(128, DX);
        packet().write(&mut wave);
        wave
    }

    fn centroid(wave: &Wave) -> Vec2 {
        let mut sum = [0.0; 2];
        let mut total = 0.0;
        for j in 0..wave.n {
            for i in 0..wave.n {
                let d = wave.density(wave.index((i, j))) as f64;
                let p = wave.position((i, j));
                sum[0] += d * p[0];
                sum[1] += d * p[1];
                total += d;
            }
        }
        [sum[0] / total, sum[1] / total]
    }

    #[test]
    fn packet_starts_normalised_and_keeps_its_norm() {
        let mut wave = free_wave();
        assert!((wave.norm() - 1.0).abs() < 1e-3, "{}", wave.norm());
        let dt = stable_dt(DX, 0.0);
        for _ in 0..300 {
            wave.step(dt);
        }
        assert!((wave.norm() - 1.0).abs() < 1e-3, "{}", wave.norm());
    }

    #[test]
    fn free_packet_moves_at_its_wavevector() {
        let mut wave = free_wave();
        let dt = stable_dt(DX, 0.0);
        let steps = 400;
        let start = centroid(&wave);
        for _ in 0..steps {
            wave.step(dt);
        }
        let moved = centroid(&wave)[0] - start[0];
        let expected = 4.0 * dt * steps as f64;
        assert!((moved - expected).abs() < 0.03 * expected, "{moved} vs {expected}");
    }

    #[test]
    fn guiding_velocity_of_a_broad_packet_is_its_wavevector() {
        let mut wave = Wave::new(128, DX);
        Packet { center: [0.0; 2], width: [3.0, 3.0], wavevector: [4.0, 1.0] }.write(&mut wave);
        let v = wave.velocity([0.2, -0.3]);
        assert!((v[0] - 4.0).abs() < 1e-3 && (v[1] - 1.0).abs() < 1e-3, "{v:?}");
    }

    #[test]
    fn sponge_removes_what_reaches_the_edge() {
        let mut wave = free_wave();
        let dt = stable_dt(DX, 0.0);
        wave.absorb_edges(dt);
        for _ in 0..3000 {
            wave.step(dt);
        }
        assert!(wave.norm() < 0.05, "{}", wave.norm());
    }
}
