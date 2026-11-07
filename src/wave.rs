//! The guiding field: a wavefunction on a square grid stepped by the
//! Schrödinger equation. No I/O.
//!
//! Units: ħ = m = 1, so a wavenumber k is also a velocity and the
//! kinetic energy is k²/2. Real and imaginary parts are stepped
//! alternately, the imaginary part living half a step ahead (Visscher
//! 1991), which keeps the scheme explicit and stable while
//! dt < 2 / (max kinetic + max potential). Hard walls are cells kept
//! at zero; an absorbing sponge at the edges swallows what reaches it.
use std::sync::Mutex;
use std::thread;

use crate::grid::Vec2;

pub type Index = (usize, usize);

/// Fourth order Laplacian stencil along one axis.
const STENCIL: [f32; 5] = [-1.0 / 12.0, 4.0 / 3.0, -5.0 / 2.0, 4.0 / 3.0, -1.0 / 12.0];
/// Sponge depth at each edge, in cells, and its strength.
const SPONGE: usize = 24;
const SPONGE_RATE: f32 = 40.0;
/// Work per thread below which another thread costs more to spawn
/// than it saves: on an M1 a 512² step takes 0.125 ms on 4 threads
/// and 0.165 ms on 8.
const CELLS_PER_THREAD: usize = 1 << 16;
/// Bands per thread, so a fast core can take over bands from a slow
/// one: 0.45 ms for a 1024² step on an M1 (4 fast and 4 slow cores)
/// against 0.48 ms with one band per thread.
const BANDS_PER_THREAD: usize = 2;

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
}

/// Largest stable step for the given cell size and potential.
pub fn stable_dt(dx: f64, potential_max: f64) -> f64 {
    let kinetic_max = 2.0 * 16.0 / 3.0 / (dx * dx) / 2.0;
    0.5 * 2.0 / (kinetic_max + potential_max)
}

/// One interior row of f ← (f + c H g) keep, given the five rows of g
/// centred on it. H g is zero on the outer two cells at each end.
fn update_row(
    f: &mut [f32],
    g: [&[f32]; 5],
    potential: &[f32],
    keep: &[f32],
    c: f32,
    inv: f32,
) {
    // c H folded into weights for the symmetric stencil, with fused
    // multiply-adds: 11 flops a cell instead of 23. Measured on an M1,
    // 1024² step: 0.50 → 0.45 ms on 8 threads, 1.58 → 1.43 ms on one.
    let [w0, w1, w2, _, _] = STENCIL.map(|s| -0.5 * inv * c * s);
    let n = f.len();
    let [a, b, m, d, e] = g.map(|row| &row[..n]);
    let (potential, keep) = (&potential[..n], &keep[..n]);
    for i in 2..n - 2 {
        let outer = (a[i] + e[i]) + (m[i - 2] + m[i + 2]);
        let inner = (b[i] + d[i]) + (m[i - 1] + m[i + 1]);
        let centre = c.mul_add(potential[i], 2.0 * w2);
        let stepped = w1.mul_add(inner, centre.mul_add(m[i], f[i]));
        f[i] = w0.mul_add(outer, stepped) * keep[i];
    }
    for i in [0, 1, n - 2, n - 1] {
        f[i] *= keep[i];
    }
}

/// Rows of one component from row `first` on: copies of the rows
/// above a band, the band itself, and copies of the rows below it.
struct Rows<'a> {
    first: usize,
    parts: [&'a [f32]; 3],
}

impl<'a> Rows<'a> {
    fn row(&self, n: usize, j: usize) -> &'a [f32] {
        let mut k = (j - self.first) * n;
        for part in self.parts {
            if k < part.len() {
                return &part[k..k + n];
            }
            k -= part.len();
        }
        unreachable!("row {j} is outside the band and its halo")
    }
}

/// What a row update needs besides the rows themselves.
struct Lattice<'a> {
    n: usize,
    inv: f32,
    potential: &'a [f32],
    keep: &'a [f32],
}

impl Lattice<'_> {
    /// Row j of f ← (f + c H g) keep.
    fn update(&self, j: usize, c: f32, g: &Rows, f: &mut [f32]) {
        let n = self.n;
        let keep = &self.keep[j * n..][..n];
        if j < 2 || j >= n - 2 {
            for (f, k) in f.iter_mut().zip(keep) {
                *f *= k;
            }
        } else {
            let rows = [0, 1, 2, 3, 4].map(|o| g.row(n, j + o - 2));
            let potential = &self.potential[j * n..][..n];
            update_row(f, rows, potential, keep, c, self.inv);
        }
    }
}

/// Copies of the rows of f up to `depth` above and below rows j0..j1.
fn halo(f: &[f32], n: usize, j0: usize, j1: usize, depth: usize) -> [Vec<f32>; 2] {
    let above = j0.saturating_sub(depth) * n..j0 * n;
    let below = j1 * n..(j1 + depth).min(n) * n;
    [f[above].to_vec(), f[below].to_vec()]
}

/// Steps the band of rows starting at j0 in a single sweep, re on row
/// j then im on row j - 2, so each row is read from memory once. The
/// band needs nothing from its neighbours during the sweep: it gets
/// copies of im four rows deep around it, enough to advance its own
/// copies of the two re rows on either side that im needs.
fn advance(
    lattice: &Lattice,
    dt: f32,
    j0: usize,
    (re, im): (&mut [f32], &mut [f32]),
    [mut re_above, mut re_below]: [Vec<f32>; 2],
    [im_above, im_below]: [Vec<f32>; 2],
) {
    let n = lattice.n;
    let j1 = j0 + re.len() / n;
    let im_first = j0 - im_above.len() / n;
    let re_first = j0 - re_above.len() / n;
    let g = Rows { first: im_first, parts: [&im_above, &*im, &im_below] };
    for (k, f) in re_above.chunks_mut(n).enumerate() {
        lattice.update(re_first + k, dt, &g, f);
    }
    for (k, f) in re_below.chunks_mut(n).enumerate() {
        lattice.update(j1 + k, dt, &g, f);
    }
    for j in j0..j1 + 2 {
        if j < j1 {
            let g = Rows { first: im_first, parts: [&im_above, &*im, &im_below] };
            lattice.update(j, dt, &g, &mut re[(j - j0) * n..][..n]);
        }
        if j >= j0 + 2 {
            let g = Rows { first: re_first, parts: [&re_above, &*re, &re_below] };
            lattice.update(j - 2, -dt, &g, &mut im[(j - 2 - j0) * n..][..n]);
        }
    }
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

    /// One step, rows split into bands that the threads take in turn.
    pub fn step(&mut self, dt: f64) {
        let (n, dt) = (self.n, dt as f32);
        let lattice = Lattice {
            n,
            inv: 1.0 / (self.dx * self.dx) as f32,
            potential: &self.potential,
            keep: &self.keep,
        };
        let cores = thread::available_parallelism().map_or(1, |t| t.get());
        let threads = (n * n / CELLS_PER_THREAD).clamp(1, cores);
        let band = n.div_ceil(threads * BANDS_PER_THREAD);
        let starts = (0..n).step_by(band);
        let halos: Vec<_> = starts
            .clone()
            .map(|j0| {
                let j1 = (j0 + band).min(n);
                (halo(&self.re, n, j0, j1, 2), halo(&self.im, n, j0, j1, 4))
            })
            .collect();
        let bands = self.re.chunks_mut(band * n).zip(self.im.chunks_mut(band * n));
        let queue = Mutex::new(starts.zip(bands).zip(halos));
        let work = || loop {
            let Some(((j0, rows), (re, im))) = queue.lock().unwrap().next() else {
                break;
            };
            advance(&lattice, dt, j0, rows, re, im);
        };
        thread::scope(|s| {
            for _ in 1..threads {
                s.spawn(work);
            }
            work();
        });
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

    /// H f cell by cell, zero on the outer two cells at each edge.
    fn plain_hamiltonian(wave: &Wave, f: &[f32]) -> Vec<f32> {
        let n = wave.n;
        let inv = 1.0 / (wave.dx * wave.dx) as f32;
        let cell = |k: usize| {
            let (i, j) = (k % n, k / n);
            if i < 2 || j < 2 || i >= n - 2 || j >= n - 2 {
                return 0.0;
            }
            let lap: f32 = (0..5)
                .map(|o| STENCIL[o] * (f[k + o - 2] + f[k + o * n - 2 * n]))
                .sum();
            -0.5 * lap * inv + wave.potential[k] * f[k]
        };
        (0..n * n).map(cell).collect()
    }

    /// The scheme written out plainly, one component at a time.
    fn plain_step(wave: &mut Wave, dt: f64) {
        let dt = dt as f32;
        let h = plain_hamiltonian(wave, &wave.im);
        for k in 0..h.len() {
            wave.re[k] = (wave.re[k] + dt * h[k]) * wave.keep[k];
        }
        let h = plain_hamiltonian(wave, &wave.re);
        for k in 0..h.len() {
            wave.im[k] = (wave.im[k] - dt * h[k]) * wave.keep[k];
        }
    }

    /// A packet heading for a potential bump and a wall, with a
    /// sponge, on a grid that splits unevenly among the cores.
    fn obstacle_course() -> Wave {
        let mut wave = Wave::new(101, DX);
        packet().write(&mut wave);
        wave.absorb_edges(stable_dt(DX, 10.0));
        for j in 30..70 {
            for i in 55..60 {
                let k = wave.index((i, j));
                wave.potential[k] = 10.0;
                wave.keep[k + 10] = 0.0;
            }
        }
        wave
    }

    #[test]
    fn steps_match_the_plain_scheme() {
        let (mut wave, mut plain) = (obstacle_course(), obstacle_course());
        let dt = stable_dt(DX, 10.0);
        for _ in 0..50 {
            wave.step(dt);
            plain_step(&mut plain, dt);
        }
        let parts = [(&wave.re, &plain.re), (&wave.im, &plain.im)];
        for (a, b) in parts {
            let worst = a.iter().zip(b).map(|(x, y)| (x - y).abs()).fold(0.0, f32::max);
            assert!(worst < 1e-5, "{worst}");
        }
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
