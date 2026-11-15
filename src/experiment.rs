//! An experiment as data: the stage, the particles, how many runs,
//! how long, the detector and what is expected, in the engine's units
//! with the conversion to SI. No I/O.
//!
//! Engine units: ħ = m = 1, so a velocity is a wavenumber. A beam's
//! wavelength is 1, so its speed is 2π and its kinetic energy 2π².
//! The particle's mass and the stage's scale, the SI length of one
//! unit, give the units of time, energy and speed.
//!
//! Every particle is an independent run of the one-particle
//! experiment: they never interact.
use std::f64::consts::{PI, TAU};

use crate::packet::{Complex, Packet};
use crate::plane::{inside, Rect, Vec2};

pub const SPEED: f64 = TAU;
pub const ENERGY: f64 = SPEED * SPEED / 2.0;
pub const HBAR: f64 = 1.054_571_817e-34;
pub const ELECTRON_MASS: f64 = 9.109_383_7e-31;
pub const ELECTRON_VOLT: f64 = 1.602_176_634e-19;

/// The potential at a position.
pub type Potential = fn(Vec2) -> f64;
/// The amplitude of the wave at the start, at a position.
pub type State = fn(Vec2) -> Complex;
/// A closed form along the screen.
pub type Curve = fn(f64) -> f64;

/// A kind of particle: its name and its mass in kg, which with the
/// stage's scale gives the units.
#[derive(Clone, Copy)]
pub struct Kind {
    pub name: &'static str,
    pub mass: f64,
}

pub const ELECTRON: Kind = Kind { name: "electron", mass: ELECTRON_MASS };

/// A particle in a run: its kind, where it starts, how fast it moves,
/// and the spread of the packet that carries it, standard deviations
/// of where it may be along x and y. The wave at the start is the
/// Gaussian packet with these moments, unless the experiment gives
/// its `state`.
#[derive(Clone, Copy)]
pub struct Particle {
    pub kind: Kind,
    pub position: Vec2,
    pub velocity: Vec2,
    pub width: Vec2,
}

/// Where the experiment happens: the grid, what one unit of length
/// is in metres, the hard walls and the potential V(x, y), all fixed
/// before the run.
#[derive(Clone)]
pub struct Stage {
    pub cells: usize,
    pub dx: f64,
    pub scale: f64,
    pub walls: Vec<Rect>,
    pub potential: Potential,
}

/// What a run is asked about: only what a detector could record,
/// the particles' positions, counts and times.
#[derive(Clone)]
pub enum Measure {
    /// Fraction of the particles that crossed the screen.
    Arrived,
    /// Mean time of the first crossings.
    MeanArrivalTime,
    /// Rms gap between the arrivals along the screen and the theory,
    /// both scaled to the same total, as a fraction of the arrivals'
    /// peak.
    ArrivalsVsTheory,
    /// The arrivals in the bin at x, as a fraction of the peak bin.
    ArrivalsAt { x: f64 },
    /// Where the arrivals peak within the given distance of x.
    ArrivalsPeakNear { x: f64, within: f64 },
    /// Mean and standard deviation of the particles' coordinate along
    /// the axis at the end.
    Mean { axis: usize },
    Spread { axis: usize },
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
    /// The laboratory experiment or closed form this one is checked
    /// against.
    pub reference: &'static str,
    pub stage: Stage,
    /// The particle of a run.
    pub particle: Particle,
    /// The wave at the start, when it is not the particle's Gaussian
    /// packet.
    pub state: Option<State>,
    /// How many times the run is repeated, each from its own start;
    /// the runs never interact.
    pub runs: usize,
    pub duration: f64,
    /// y of the line where arrivals are counted.
    pub screen: Option<f64>,
    /// The closed form the arrivals along the screen are compared with.
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
/// between `from` and `to`.
pub fn band(axis: usize, from: f64, to: f64) -> Rect {
    let mut rect = [[f64::NEG_INFINITY, f64::INFINITY]; 2];
    rect[axis] = [from, to];
    rect
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

impl Stage {
    /// Side length of the plane, centred on the origin.
    pub fn extent(&self) -> f64 {
        self.cells as f64 * self.dx
    }

    /// Centre of cell k, counted along x and then up y.
    pub fn position(&self, k: usize) -> Vec2 {
        [k % self.cells, k / self.cells].map(|i| (i as f64 + 0.5) * self.dx - self.extent() / 2.0)
    }

    pub fn walled(&self, p: Vec2) -> bool {
        self.walls.iter().any(|w| inside(w, p))
    }
}

impl Experiment {
    /// SI value of one engine unit of time, energy and speed, from
    /// the particle's mass and the stage's scale.
    pub fn time(&self) -> f64 {
        let (m, l) = (self.particle.kind.mass, self.stage.scale);
        m * l * l / HBAR
    }

    pub fn energy(&self) -> f64 {
        let (m, l) = (self.particle.kind.mass, self.stage.scale);
        HBAR * HBAR / (m * l * l)
    }

    pub fn speed(&self) -> f64 {
        HBAR / (self.particle.kind.mass * self.stage.scale)
    }

    /// The potential's largest magnitude on the grid, which bounds
    /// the step an engine can take.
    pub fn steepest(&self) -> f64 {
        let stage = &self.stage;
        let at = |k| (stage.potential)(stage.position(k)).abs();
        (0..stage.cells * stage.cells).map(at).fold(0.0, f64::max)
    }

    /// The wave at the start: the given state, or the particle's
    /// packet.
    pub fn amplitude(&self, p: Vec2) -> Complex {
        let Particle { position, velocity, width, .. } = self.particle;
        let packet = Packet { center: position, width, wavevector: velocity };
        self.state.map_or_else(|| packet.amplitude(p), |state| state(p))
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
    fn an_electron_at_fifty_picometres_has_six_hundred_electron_volts() {
        let e = crate::experiments::single_slit::single_slit();
        let ev = ENERGY * e.energy() / ELECTRON_VOLT;
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
