//! An experiment as data: the stage, the source, how long it runs and
//! what is expected, in the engine's units with the conversion to a
//! physical realisation. No I/O.
//!
//! Engine units: ħ = m = 1. A beam's wavelength is 1, so its wavenumber
//! and speed are both 2π and its kinetic energy 2π². A realisation
//! names the particle and the SI length of one unit, from which the
//! units of time, energy and speed follow.
//!
//! The grid is either the plane one particle moves in, or the
//! configuration space (x₁, x₂) of two particles on a line, which the
//! solver steps the same way since with unit masses the Laplacian is
//! both particles' kinetic energy.
use std::f64::consts::{PI, TAU};

use crate::packet::Packet;
use crate::pair::{self, Orbital, Symmetry};
use crate::wave::{stable_dt, Vec2, Wave};

pub const SPEED: f64 = TAU;
pub const ENERGY: f64 = SPEED * SPEED / 2.0;
pub const HBAR: f64 = 1.054_571_817e-34;
pub const ELECTRON_MASS: f64 = 9.109_383_7e-31;
pub const ELECTRON_VOLT: f64 = 1.602_176_634e-19;
/// Two particles on a line count as together within this distance.
pub const TOGETHER: f64 = 0.5;

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

/// What the grid is: the plane one particle moves in, or the
/// configuration space of two particles on a line, x₁ along the
/// grid's x and x₂ along its y.
#[derive(Clone, Copy, PartialEq)]
pub enum Space {
    Plane,
    Line,
}

#[derive(Clone)]
pub enum Wall {
    /// A wall across the plane at height y, with openings of the
    /// given width at the given x.
    Slits { y: f64, thickness: f64, centers: Vec<f64>, width: f64 },
    /// Everything with the coordinate along the axis between `from`
    /// and `to`. On a line, a band along axis 0 is a wall only
    /// particle 1 feels.
    Band { axis: usize, from: f64, to: f64 },
}

#[derive(Clone)]
pub enum Potential {
    /// A steep ramp, slope·(p − from) beyond `from` along the axis: a
    /// soft wall, which on a line only one particle feels.
    Ramp { axis: usize, from: f64, slope: f64 },
    /// ½ω²p² along the axis, held at `cap` far out, where the state
    /// is negligible, so that the time step stays large.
    Harmonic { axis: usize, omega: f64, cap: f64 },
    /// A Gaussian well in the separation of two particles on a line,
    /// −depth·exp(−(x₁ − x₂)²/2 range²).
    Bond { depth: f64, range: f64 },
}

#[derive(Clone)]
pub enum Source {
    Packet(Packet),
    /// Two particles on a line, a(x₁)b(x₂), plus or minus b(x₁)a(x₂)
    /// for bosons or fermions.
    Pair { a: Orbital, b: Orbital, symmetry: Symmetry },
    /// Two particles flying apart from the origin, Gaussian in their
    /// centre c = (x₁ + x₂)/2 and their separation r = x₁ − x₂ with
    /// relative wavenumber k, ψ ∝ exp(−c²/4σ_c² − r²/4σ_r²) e^{−ikr}:
    /// positions correlated, momenta anticorrelated. Unentangled, it is
    /// the product of that state's two marginal packets instead.
    Epr { widths: [f64; 2], wavenumber: f64, entangled: bool },
    /// Two bound particles: a Gaussian in the centre and in the
    /// separation relaxed in imaginary time into the ground state of
    /// the bond, then kicked with total wavenumber k.
    Molecule { center: f64, widths: [f64; 2], relax: f64, kick: f64 },
    /// The ground state of a box `from` < y < `from` + length at the
    /// closed end of a channel along y with harmonic confinement of
    /// frequency ω across it; the box's far wall is removed at t = 0.
    Channel { from: f64, length: f64, omega: f64 },
}

/// The closed form the screen's current is compared with.
#[derive(Clone, Copy)]
pub enum Theory {
    Fraunhofer,
}

#[derive(Clone)]
pub enum Measure {
    /// Fraction of the particles that crossed the screen.
    Arrived,
    /// That fraction minus the probability that crossed the screen.
    ArrivedMinusThrough,
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
    /// On a line: fraction of the runs whose point ever crossed the
    /// diagonal x₁ = x₂, the two particles passing each other.
    Crossed,
    /// On a line: the most of the runs ever within TOGETHER of each
    /// other at one time, and the most probability.
    TogetherPeak,
    ProbabilityTogetherPeak,
    /// On a line: fraction of the runs more than d apart at the end,
    /// and that minus the probability apart, counting what the sponge
    /// took as apart.
    Apart { beyond: f64 },
    ApartMinusProbability { beyond: f64 },
    /// On a line: correlation coefficient of x₁ and x₂ over the runs at
    /// the end, and that minus the field's.
    Correlation,
    CorrelationMinusField,
    /// Median change of a coordinate between each run and the same run
    /// of the twin: the experiment without its potentials.
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
    pub space: Space,
    pub cells: usize,
    pub dx: f64,
    pub walls: Vec<Wall>,
    pub potentials: Vec<Potential>,
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
    pub theory: Option<Theory>,
    pub expectations: Vec<Expectation>,
}

fn inside(wall: &Wall, p: Vec2) -> bool {
    match wall {
        Wall::Slits { y, thickness, centers, width } => {
            (p[1] - y).abs() < thickness / 2.0
                && !centers.iter().any(|c| (p[0] - c).abs() < width / 2.0)
        }
        Wall::Band { axis, from, to } => (*from..*to).contains(&p[*axis]),
    }
}

fn potential(v: &Potential, p: Vec2) -> f64 {
    match v {
        Potential::Ramp { axis, from, slope } => slope * (p[*axis] - from).max(0.0),
        Potential::Harmonic { axis, omega, cap } => {
            (0.5 * omega * omega * p[*axis] * p[*axis]).min(*cap)
        }
        Potential::Bond { depth, range } => {
            let r = p[0] - p[1];
            -depth * (-r * r / (2.0 * range * range)).exp()
        }
    }
}

/// The potential's largest magnitude on a grid reaching `half` from
/// the origin, which bounds the time step.
fn potential_max(v: &Potential, half: f64) -> f64 {
    match v {
        Potential::Ramp { from, slope, .. } => slope * (half - from).max(0.0),
        Potential::Harmonic { cap, .. } => *cap,
        Potential::Bond { depth, .. } => *depth,
    }
}

fn gaussian(x: f64, width: f64) -> f64 {
    (-x * x / (4.0 * width * width)).exp()
}

impl Source {
    pub fn write(&self, wave: &mut Wave, dt: f64) {
        match self {
            Source::Packet(packet) => packet.write(wave),
            Source::Pair { a, b, symmetry } => pair::write(wave, *a, *b, *symmetry),
            Source::Epr { widths: [c, r], wavenumber: k, entangled: true } => {
                pair::fill(wave, |x1, x2| {
                    let (u, v) = ((x1 + x2) / 2.0, x1 - x2);
                    let a = gaussian(u, *c) * gaussian(v, *r);
                    [a * (k * v).cos(), -a * (k * v).sin()]
                })
            }
            Source::Epr { widths: [c, r], wavenumber: k, entangled: false } => {
                let width = (c * c + r * r / 4.0).sqrt();
                let one = Orbital { center: 0.0, width, wavenumber: -k };
                let two = Orbital { wavenumber: *k, ..one };
                pair::write(wave, one, two, Symmetry::Distinguishable);
            }
            Source::Molecule { center, widths: [w, s], relax, kick } => {
                pair::fill(wave, |x1, x2| {
                    [gaussian((x1 + x2) / 2.0 - center, *w) * gaussian(x1 - x2, *s), 0.0]
                });
                pair::relax(wave, dt, (relax / dt) as usize);
                pair::kick(wave, *kick);
            }
            Source::Channel { from, length, omega } => pair::fill(wave, |x, y| {
                let along = if (*from..from + length).contains(&y) {
                    (PI * (y - from) / length).sin()
                } else {
                    0.0
                };
                [along * (-0.5 * omega * x * x).exp(), 0.0]
            }),
        }
    }
}

/// Fraunhofer intensity per unit x on a screen at distance d from
/// slits of the given width at the given centres. The cylindrical
/// wave falls as 1/r and meets the flat screen obliquely, which
/// together give a factor cos²θ.
fn fraunhofer(x: f64, d: f64, centers: &[f64], width: f64) -> f64 {
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

    pub fn dt(&self) -> f64 {
        let half = self.extent() / 2.0;
        let v = self.potentials.iter().map(|v| potential_max(v, half)).fold(0.0, f64::max);
        stable_dt(self.dx, v)
    }

    /// The field at the start: sponge at the edges, walls cut out
    /// after it so they stay whole, potentials added, the source
    /// written.
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
                wave.potential[k] = self.potentials.iter().map(|v| potential(v, p)).sum::<f64>() as f32;
            }
        }
        self.source.write(&mut wave, self.dt());
        wave
    }

    /// The same experiment without its potentials, for measures that
    /// compare runs with and without them.
    pub fn twin(&self) -> Experiment {
        Experiment { potentials: Vec::new(), expectations: Vec::new(), ..self.clone() }
    }

    pub fn needs_twin(&self) -> bool {
        let twin = |m: &Measure| matches!(m, Measure::MovedVsTwin { .. } | Measure::MarginalVsTwin { .. });
        self.expectations.iter().any(|x| twin(&x.measure))
    }

    /// The theory's value at x on the screen, zero without a theory.
    pub fn theory(&self, x: f64) -> f64 {
        let Some(Theory::Fraunhofer) = self.theory else { return 0.0 };
        let screen = self.screen.expect("a theory needs a screen");
        let [Wall::Slits { y, thickness, centers, width }] = &self.walls[..] else {
            panic!("Fraunhofer needs exactly one wall of slits")
        };
        fraunhofer(x, screen - y - thickness / 2.0, centers, *width)
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
    fn epr_state_has_the_closed_form_correlation_and_its_control_none() {
        let wave = |entangled| {
            let mut wave = Wave::new(128, 0.1);
            Source::Epr { widths: [1.0, 0.6], wavenumber: 3.0, entangled }.write(&mut wave, 0.0);
            wave
        };
        let rho = (1.0 - 0.36 / 4.0) / (1.0 + 0.36 / 4.0);
        let correlation = pair::correlation;
        assert!((correlation(&wave(true)) - rho).abs() < 1e-3);
        assert!(correlation(&wave(false)).abs() < 1e-6);
    }
}
