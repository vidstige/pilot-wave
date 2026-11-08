//! An experiment as data: the stage, the source, how long it runs and
//! what is expected, in the engine's units with the conversion to a
//! physical realisation. No I/O.
//!
//! Engine units: ħ = m = 1. A beam's wavelength is 1, so its wavenumber
//! and speed are both 2π and its kinetic energy 2π². A realisation
//! names the particle and the SI length of one unit, from which the
//! units of time, energy and speed follow.
//!
//! The grid is a plane. It is either the plane one particle moves in,
//! or the configuration space (x₁, x₂) of two particles on a line:
//! with unit masses the Laplacian in (x₁, x₂) is both particles'
//! kinetic energy, so the solver steps it unchanged, and one sampled
//! point is one run of the two-particle experiment.
use std::f64::consts::{PI, TAU};

use crate::wave::{stable_dt, Vec2, Wave};

pub const SPEED: f64 = TAU;
pub const ENERGY: f64 = SPEED * SPEED / 2.0;
pub const HBAR: f64 = 1.054_571_817e-34;
pub const ELECTRON_MASS: f64 = 9.109_383_7e-31;
pub const ELECTRON_VOLT: f64 = 1.602_176_634e-19;
/// Two particles on a line count as together within this distance.
pub const TOGETHER: f64 = 0.5;

/// A rectangle [[x₀, x₁], [y₀, y₁]] of hard wall.
pub type Rect = [[f64; 2]; 2];
/// The potential at a position.
pub type Potential = fn(Vec2) -> f64;
/// Writes the initial ψ onto the field.
pub type Source = fn(&mut Wave);
/// A closed form along the screen.
pub type Curve = fn(f64) -> f64;

/// A physical particle and the SI length of one engine unit: for a
/// beam, its de Broglie wavelength.
#[derive(Clone)]
pub struct Realisation {
    pub particle: &'static str,
    /// Mass in kg and unit length in m.
    pub mass: f64,
    pub length: f64,
    pub reference: &'static str,
}

impl Realisation {
    /// SI value of one engine unit of time, energy and speed.
    pub fn time(&self) -> f64 {
        self.mass * self.length * self.length / HBAR
    }

    pub fn energy(&self) -> f64 {
        HBAR * HBAR / (self.mass * self.length * self.length)
    }

    pub fn speed(&self) -> f64 {
        HBAR / (self.mass * self.length)
    }
}

#[derive(Clone)]
pub enum Measure {
    /// Fraction of the particles that crossed the screen, and the
    /// probability that did.
    Arrived,
    Through,
    /// Fraction of the particles that first crossed the screen after
    /// the given time.
    ArrivedAfter { time: f64 },
    /// Mean time of the first crossings.
    MeanArrivalTime,
    /// Rms gap between the screen's current and the theory, both
    /// scaled to the same total, as a fraction of the current's peak.
    CurrentVsTheory,
    /// Rms gap between the arrivals and the current, likewise.
    ArrivalsVsCurrent,
    /// The current at x on the screen, as a fraction of its peak.
    CurrentAt { x: f64 },
    /// Where the current's local maximum nearest x lies.
    MaximumNear { x: f64 },
    /// Two particles on a line: fraction of the runs whose point ever
    /// crossed the diagonal x₁ = x₂, the particles passing each other.
    Crossed,
    /// The most of the runs ever within TOGETHER of each other at one
    /// time, and the most probability.
    TogetherPeak,
    ProbabilityTogetherPeak,
    /// Fraction of the runs more than d apart at the end, and the
    /// probability, counting what the sponge took as apart.
    Apart { beyond: f64 },
    ProbabilityApart { beyond: f64 },
    /// Correlation coefficient of x₁ and x₂ over the runs at the end,
    /// and in the field.
    Correlation,
    FieldCorrelation,
    /// Median change of a coordinate between each run and the same run
    /// of the twin: the experiment without its potential.
    MovedVsTwin { axis: usize },
    /// Largest gap between the field's density of a coordinate and the
    /// twin's, as a fraction of its peak.
    MarginalVsTwin { axis: usize },
}

#[derive(Clone)]
pub struct Expectation {
    pub what: &'static str,
    pub measure: Measure,
    pub expected: f64,
    pub tolerance: f64,
}

#[derive(Clone)]
pub struct Experiment {
    pub name: &'static str,
    pub summary: &'static str,
    pub realisation: Realisation,
    pub cells: usize,
    pub dx: f64,
    pub walls: Vec<Rect>,
    pub potential: Potential,
    pub source: Source,
    /// The spin direction of the pilot wave, a unit vector or zero for
    /// a spinless particle. Only the component normal to the plane
    /// steers the trajectories; a spin in the plane leaves them as
    /// spinless.
    pub spin: [f64; 3],
    pub particles: usize,
    pub duration: f64,
    /// y of the line where arrivals are counted.
    pub screen: Option<f64>,
    /// The closed form the screen's current is compared with.
    pub theory: Option<Curve>,
    pub expectations: Vec<Expectation>,
}

pub fn free(_: Vec2) -> f64 {
    0.0
}

/// A wall across the plane at height y with openings of the given
/// width at the given x, as the rectangles between the openings.
pub fn slits(y: f64, thickness: f64, centers: &[f64], width: f64) -> Vec<Rect> {
    let edges: Vec<f64> = centers.iter().flat_map(|c| [c - width / 2.0, c + width / 2.0]).collect();
    let ends = [f64::NEG_INFINITY].into_iter().chain(edges).chain([f64::INFINITY]);
    let pieces: Vec<f64> = ends.collect();
    let piece = |x: &[f64]| [[x[0], x[1]], [y - thickness / 2.0, y + thickness / 2.0]];
    pieces.chunks(2).map(piece).collect()
}

/// A wall filling everything with the coordinate along the axis
/// between `from` and `to`; along axis 0 on a line, a wall only
/// particle 1 feels.
pub fn band(axis: usize, from: f64, to: f64) -> Rect {
    let mut rect = [[f64::NEG_INFINITY, f64::INFINITY]; 2];
    rect[axis] = [from, to];
    rect
}

fn inside(rect: &Rect, p: Vec2) -> bool {
    (rect[0][0]..rect[0][1]).contains(&p[0]) && (rect[1][0]..rect[1][1]).contains(&p[1])
}

/// Fraunhofer intensity per unit x on a screen at distance d from
/// slits of the given width at the given centres. The cylindrical
/// wave falls as 1/r and meets the flat screen obliquely, which
/// together give a factor cos²θ.
pub fn fraunhofer(x: f64, d: f64, centers: &[f64], width: f64) -> f64 {
    let r = x.hypot(d);
    let (sin, cos) = (x / r, d / r);
    let beta = PI * width * sin;
    let envelope = if beta.abs() < 1e-9 { 1.0 } else { (beta.sin() / beta).powi(2) };
    let (re, im) = centers.iter().fold((0.0, 0.0), |(re, im), c| {
        let phase = TAU * c * sin;
        (re + phase.cos(), im + phase.sin())
    });
    envelope * (re * re + im * im) * cos * cos
}

impl Experiment {
    pub fn extent(&self) -> f64 {
        self.cells as f64 * self.dx
    }

    /// The largest stable step, which the potential's largest
    /// magnitude on the grid bounds.
    pub fn dt(&self) -> f64 {
        let wave = Wave::new(self.cells, self.dx);
        let at = |k: usize| (self.potential)(wave.position((k % wave.n, k / wave.n))).abs();
        stable_dt(self.dx, (0..wave.n * wave.n).map(at).fold(0.0, f64::max))
    }

    /// The field at the start: sponge at the edges, walls cut out
    /// after it so they stay whole, the potential, then the source.
    pub fn build(&self) -> Wave {
        let mut wave = Wave::new(self.cells, self.dx);
        wave.spin = self.spin[2];
        wave.absorb_edges(self.dt());
        for j in 0..wave.n {
            for i in 0..wave.n {
                let p = wave.position((i, j));
                let k = wave.index((i, j));
                if self.walls.iter().any(|w| inside(w, p)) {
                    wave.keep[k] = 0.0;
                }
                wave.potential[k] = (self.potential)(p) as f32;
            }
        }
        (self.source)(&mut wave);
        wave
    }

    /// The same experiment without its potential, for measures that
    /// compare runs with and without it.
    pub fn twin(&self) -> Experiment {
        Experiment { potential: free, expectations: Vec::new(), ..self.clone() }
    }

    pub fn measures(&self, f: impl Fn(&Measure) -> bool) -> bool {
        self.expectations.iter().any(|x| f(&x.measure))
    }

    /// The theory's value at x on the screen, zero without a theory.
    pub fn theory(&self, x: f64) -> f64 {
        self.theory.map_or(0.0, |theory| theory(x))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_slit_fraunhofer_has_its_first_zero_at_the_slit_angle() {
        let zero = (0.5f64).asin().tan() * 10.0;
        assert!(fraunhofer(zero, 10.0, &[0.0], 2.0) < 1e-12);
        assert!(fraunhofer(0.0, 10.0, &[0.0], 2.0) > 0.99);
    }

    #[test]
    fn electron_realisation_recovers_its_kinetic_energy() {
        let r = Realisation { particle: "electron", mass: ELECTRON_MASS, length: 50e-12, reference: "" };
        let ev = ENERGY * r.energy() / ELECTRON_VOLT;
        assert!((ev - 602.0).abs() < 2.0, "{ev}");
    }

    #[test]
    fn slits_leave_only_their_openings_in_the_wall() {
        let walls = slits(0.0, 0.4, &[-2.0, 2.0], 1.0);
        assert_eq!(walls.len(), 3);
        let open = |x: f64| !walls.iter().any(|w| inside(w, [x, 0.1]));
        assert!(open(-2.0) && open(2.4) && !open(0.0) && !open(-9.0) && !open(9.0));
        assert!(!walls.iter().any(|w| inside(w, [0.0, 0.3])));
    }

    #[test]
    fn a_band_spans_the_other_axis() {
        let wall = band(0, 5.0, f64::INFINITY);
        assert!(inside(&wall, [6.0, -1e9]) && !inside(&wall, [4.9, 0.0]));
    }
}
