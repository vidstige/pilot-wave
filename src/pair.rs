//! Two particles on a line as one point in the plane: ψ(x₁, x₂) on
//! the solver's grid, x₁ along i and x₂ along j. With unit masses the
//! grid's Laplacian is both particles' kinetic energy. States, their
//! marginals and imaginary-time relaxation. No I/O.
use crate::wave::{Vec2, Wave};

pub type Complex = [f64; 2];

/// The solver's fourth order Laplacian stencil along one axis.
const STENCIL: [f32; 5] = [-1.0 / 12.0, 4.0 / 3.0, -5.0 / 2.0, 4.0 / 3.0, -1.0 / 12.0];
/// Subdivisions of a cell along each axis when measuring a region.
const SUBCELLS: usize = 4;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Symmetry {
    Distinguishable,
    /// Symmetric under x₁ ↔ x₂.
    Bosons,
    /// Antisymmetric under x₁ ↔ x₂, so ψ vanishes on the diagonal.
    Fermions,
}

/// One particle's Gaussian packet on the line.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Orbital {
    pub center: f64,
    /// Standard deviation of |φ|².
    pub width: f64,
    pub wavenumber: f64,
}

impl Orbital {
    /// φ(x), unnormalised.
    fn at(&self, x: f64) -> Complex {
        let d = x - self.center;
        let a = (-d * d / (4.0 * self.width * self.width)).exp();
        let phase = self.wavenumber * x;
        [a * phase.cos(), a * phase.sin()]
    }
}

fn mul(a: Complex, b: Complex) -> Complex {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}

pub fn normalise(wave: &mut Wave) {
    let s = (1.0 / wave.norm()).sqrt() as f32;
    for v in wave.re.iter_mut().chain(&mut wave.im) {
        *v *= s;
    }
}

/// ψ = f(x₁, x₂) over the grid, zero in walls, normalised.
pub fn fill(wave: &mut Wave, f: impl Fn(f64, f64) -> Complex) {
    for j in 0..wave.n {
        for i in 0..wave.n {
            let p = wave.position((i, j));
            let k = wave.index((i, j));
            let wall = wave.keep[k] == 0.0;
            let v = if wall { [0.0; 2] } else { f(p[0], p[1]) };
            wave.re[k] = v[0] as f32;
            wave.im[k] = v[1] as f32;
        }
    }
    normalise(wave);
}

/// a(x₁) b(x₂), plus or minus b(x₁) a(x₂) for identical particles.
pub fn write(wave: &mut Wave, a: Orbital, b: Orbital, symmetry: Symmetry) {
    let sign = match symmetry {
        Symmetry::Distinguishable => 0.0,
        Symmetry::Bosons => 1.0,
        Symmetry::Fermions => -1.0,
    };
    fill(wave, |x1, x2| {
        let p = mul(a.at(x1), b.at(x2));
        let q = mul(b.at(x1), a.at(x2));
        [p[0] + sign * q[0], p[1] + sign * q[1]]
    });
}

/// Multiply by e^{ik(x₁+x₂)/2}: total momentum k shared by the pair,
/// each particle moving at k/2 on top of its relative motion.
pub fn kick(wave: &mut Wave, k: f64) {
    for j in 0..wave.n {
        for i in 0..wave.n {
            let p = wave.position((i, j));
            let phase = k * (p[0] + p[1]) / 2.0;
            let (sin, cos) = (phase.sin() as f32, phase.cos() as f32);
            let k = wave.index((i, j));
            let (re, im) = (wave.re[k], wave.im[k]);
            wave.re[k] = re * cos - im * sin;
            wave.im[k] = re * sin + im * cos;
        }
    }
}

/// Densities of x₁ and of x₂ on the line, one value per cell.
pub fn marginals(wave: &Wave) -> [Vec<f64>; 2] {
    let n = wave.n;
    let mut m = [vec![0.0; n], vec![0.0; n]];
    for j in 0..n {
        for i in 0..n {
            let d = wave.density(wave.index((i, j))) as f64 * wave.dx;
            m[0][i] += d;
            m[1][j] += d;
        }
    }
    m
}

/// Correlation coefficient of x₁ and x₂ in the field.
pub fn correlation(wave: &Wave) -> f64 {
    let n = wave.n;
    let mean = |f: &dyn Fn(Vec2) -> f64| {
        let weighted = |k: usize| wave.density(k) as f64 * f(wave.position((k % n, k / n)));
        (0..n * n).map(weighted).sum::<f64>() * wave.dx * wave.dx / wave.norm()
    };
    let (x, y) = (mean(&|p| p[0]), mean(&|p| p[1]));
    let xy = mean(&|p| p[0] * p[1]) - x * y;
    let xx = mean(&|p| p[0] * p[0]) - x * x;
    let yy = mean(&|p| p[1] * p[1]) - y * y;
    xy / (xx * yy).sqrt()
}

/// Mutual information of x₁ and x₂ in the field, in nats, the cells
/// as bins: the Kullback-Leibler divergence of the joint from the
/// product of its marginals, zero exactly when it is that product.
pub fn mutual_information(wave: &Wave) -> f64 {
    let n = wave.n;
    let [m1, m2] = marginals(wave);
    let total: f64 = m1.iter().sum::<f64>() * wave.dx;
    let term = |k: usize| {
        let p = wave.density(k) as f64 * wave.dx * wave.dx / total;
        let q = m1[k % n] * m2[k / n] * wave.dx * wave.dx / (total * total);
        if p > 0.0 { p * (p / q).ln() } else { 0.0 }
    };
    (0..n * n).map(term).sum()
}

/// Probability that the pair is in the region of configurations,
/// with the density constant over each cell as the sampler takes it.
/// Counting whole cells instead would move the edge of a band
/// |x₁ − x₂| < d by up to half a cell, 10% of a band 0.5 wide.
pub fn probability(wave: &Wave, region: impl Fn(Vec2) -> bool) -> f64 {
    let (n, s) = (wave.n, SUBCELLS);
    let offset = |a: usize| ((a as f64 + 0.5) / s as f64 - 0.5) * wave.dx;
    let share = |p: Vec2| {
        let inside = (0..s * s).filter(|q| region([p[0] + offset(q % s), p[1] + offset(q / s)]));
        inside.count() as f64 / (s * s) as f64
    };
    let at = |k: usize| wave.density(k) as f64 * share(wave.position((k % n, k / n)));
    (0..n * n).map(at).sum::<f64>() * wave.dx * wave.dx
}

/// H f with the solver's stencil, zero on the outer two cells.
fn hamiltonian(wave: &Wave, f: &[f32]) -> Vec<f32> {
    let n = wave.n;
    let inv = 1.0 / (wave.dx * wave.dx) as f32;
    let cell = |k: usize| {
        let (i, j) = (k % n, k / n);
        if i < 2 || j < 2 || i >= n - 2 || j >= n - 2 {
            return 0.0;
        }
        let along = |o: usize| f[k + o - 2] + f[k + o * n - 2 * n];
        let lap: f32 = (0..5).map(|o| STENCIL[o] * along(o)).sum();
        -0.5 * lap * inv + wave.potential[k] * f[k]
    };
    (0..n * n).map(cell).collect()
}

/// ⟨ψ|H|ψ⟩ / ⟨ψ|ψ⟩.
pub fn energy(wave: &Wave) -> f64 {
    let dot = |f: &[f32]| {
        let h = hamiltonian(wave, f);
        f.iter().zip(&h).map(|(a, b)| (a * b) as f64).sum::<f64>()
    };
    (dot(&wave.re) + dot(&wave.im)) * wave.dx * wave.dx / wave.norm()
}

/// f ← f − dτ H f, zero in walls.
fn descend(f: &mut [f32], h: &[f32], keep: &[f32], dtau: f32) {
    for ((f, h), keep) in f.iter_mut().zip(h).zip(keep) {
        *f = if *keep == 0.0 { 0.0 } else { *f - dtau * h };
    }
}

/// Imaginary-time relaxation, ψ ← ψ − dτ H ψ renormalised each
/// step. Every eigenstate fades against the lowest one ψ overlaps as
/// e^{−(E − E₀)τ}, so ψ settles into the ground state; dτ must stay
/// below 2 over the largest energy on the grid. Walls stay zero, the
/// sponge plays no part.
pub fn relax(wave: &mut Wave, dtau: f64, steps: usize) {
    let dtau = dtau as f32;
    for _ in 0..steps {
        let h = hamiltonian(wave, &wave.re);
        descend(&mut wave.re, &h, &wave.keep, dtau);
        let h = hamiltonian(wave, &wave.im);
        descend(&mut wave.im, &h, &wave.keep, dtau);
        normalise(wave);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::wave::stable_dt;

    const DX: f64 = 0.1;

    fn orbitals() -> [Orbital; 2] {
        [
            Orbital { center: -1.5, width: 0.8, wavenumber: 2.0 },
            Orbital { center: 1.0, width: 0.6, wavenumber: -1.0 },
        ]
    }

    fn pair(symmetry: Symmetry) -> Wave {
        let mut wave = Wave::new(96, DX);
        let [a, b] = orbitals();
        write(&mut wave, a, b, symmetry);
        wave
    }

    fn swapped(wave: &Wave, k: usize) -> usize {
        wave.index((k / wave.n, k % wave.n))
    }

    #[test]
    fn antisymmetric_state_vanishes_on_the_diagonal() {
        let wave = pair(Symmetry::Fermions);
        let diagonal = (0..wave.n).map(|i| wave.density(wave.index((i, i))));
        assert!(diagonal.fold(0.0, f32::max) < 1e-12 * wave.peak());
        for k in 0..wave.n * wave.n {
            assert_eq!(wave.re[k], -wave.re[swapped(&wave, k)]);
        }
    }

    #[test]
    fn symmetric_state_is_unchanged_by_exchange() {
        let wave = pair(Symmetry::Bosons);
        for k in 0..wave.n * wave.n {
            assert_eq!(wave.im[k], wave.im[swapped(&wave, k)]);
        }
        let diagonal = (0..wave.n).map(|i| wave.density(wave.index((i, i))));
        assert!(diagonal.fold(0.0, f32::max) > 1e-3 * wave.peak());
    }

    #[test]
    fn distinguishable_marginals_are_the_two_packets_and_integrate_to_one() {
        let wave = pair(Symmetry::Distinguishable);
        let x = |i: usize| wave.position((i, 0))[0];
        for (m, orbital) in marginals(&wave).iter().zip(orbitals()) {
            let total: f64 = m.iter().sum::<f64>() * DX;
            let mean: f64 = m.iter().enumerate().map(|(i, d)| x(i) * d).sum::<f64>() * DX;
            assert!((total - 1.0).abs() < 1e-4, "{total}");
            assert!((mean - orbital.center).abs() < 1e-3, "{mean}");
        }
    }

    #[test]
    fn a_product_shares_no_information_and_two_branches_share_ln_2() {
        assert!(mutual_information(&pair(Symmetry::Distinguishable)).abs() < 1e-9);
        let mut wave = Wave::new(96, DX);
        let blob = |x: f64, y: f64| (-(x * x + y * y) / 0.5).exp();
        fill(&mut wave, |x1, x2| [blob(x1 - 2.0, x2 + 2.0) + blob(x1 + 2.0, x2 - 2.0), 0.0]);
        assert!((mutual_information(&wave) - 2f64.ln()).abs() < 1e-3);
    }

    #[test]
    fn relaxation_in_a_well_lowers_the_energy_to_a_bound_state() {
        let mut wave = Wave::new(64, DX);
        let well = |r2: f64| -4.0 * (-r2 / 2.0).exp();
        for k in 0..wave.n * wave.n {
            let p = wave.position((k % wave.n, k / wave.n));
            wave.potential[k] = well(p[0] * p[0] + p[1] * p[1]) as f32;
        }
        let a = Orbital { center: 0.5, width: 0.5, wavenumber: 1.0 };
        write(&mut wave, a, a, Symmetry::Distinguishable);
        let start = energy(&wave);
        let dtau = stable_dt(DX, 4.0);
        relax(&mut wave, dtau, 3000);
        let relaxed = energy(&wave);
        relax(&mut wave, dtau, 1000);
        assert!(relaxed < start - 1.0 && relaxed < 0.0, "{start} {relaxed}");
        assert!((energy(&wave) - relaxed).abs() < 1e-3 * relaxed.abs());
        assert!((wave.norm() - 1.0).abs() < 1e-4);
    }
}
