//! An experiment as data: the stage, the source, how long it runs and
//! what is expected, in the engine's units with the conversion to a
//! physical realisation. No I/O.
//!
//! Engine units: ħ = m = 1 and the packet's wavelength is 1, so its
//! wavenumber and speed are both 2π and its kinetic energy 2π². A
//! realisation names the particle and its de Broglie wavelength, from
//! which the units of length, time, energy and speed follow.
use std::f64::consts::{PI, TAU};

use crate::packet::Packet;
use crate::wave::{stable_dt, Vec2, Wave};

pub const SPEED: f64 = TAU;
pub const ENERGY: f64 = SPEED * SPEED / 2.0;
pub const HBAR: f64 = 1.054_571_817e-34;
pub const ELECTRON_MASS: f64 = 9.109_383_7e-31;
pub const ELECTRON_VOLT: f64 = 1.602_176_634e-19;

/// A physical particle whose de Broglie wavelength is the length unit.
pub struct Realisation {
    pub particle: &'static str,
    /// Mass in kg and wavelength in m.
    pub mass: f64,
    pub wavelength: f64,
    pub reference: &'static str,
}

impl Realisation {
    /// SI value of one engine unit of length, time, energy and speed.
    pub fn length(&self) -> f64 {
        self.wavelength
    }

    pub fn time(&self) -> f64 {
        self.mass * self.wavelength * self.wavelength / HBAR
    }

    pub fn energy(&self) -> f64 {
        HBAR * HBAR / (self.mass * self.wavelength * self.wavelength)
    }

    pub fn speed(&self) -> f64 {
        HBAR / (self.mass * self.wavelength)
    }
}

pub enum Wall {
    /// A wall across the domain at height y, with openings of the
    /// given width at the given x.
    Slits { y: f64, thickness: f64, centers: Vec<f64>, width: f64 },
}

/// The closed form the screen's current is compared with.
pub enum Theory {
    Fraunhofer,
}

pub enum Measure {
    /// Fraction of the particles that crossed the screen.
    Arrived,
    /// That fraction minus the probability that crossed the screen.
    ArrivedMinusThrough,
    /// Rms gap between the screen's current and the theory, both
    /// scaled to the same total, as a fraction of the current's peak.
    CurrentVsTheory,
    /// Rms gap between the arrivals and the current, likewise.
    ArrivalsVsCurrent,
    /// The current at x on the screen, as a fraction of its peak.
    CurrentAt { x: f64 },
}

pub struct Expectation {
    pub what: &'static str,
    pub measure: Measure,
    pub expected: f64,
    pub tolerance: f64,
}

pub struct Experiment {
    pub name: &'static str,
    pub summary: &'static str,
    pub realisation: Realisation,
    pub cells: usize,
    pub dx: f64,
    pub walls: Vec<Wall>,
    pub source: Packet,
    pub particles: usize,
    pub duration: f64,
    /// y of the line where arrivals are counted.
    pub screen: Option<f64>,
    pub theory: Theory,
    pub expectations: Vec<Expectation>,
}

fn inside(wall: &Wall, p: Vec2) -> bool {
    match wall {
        Wall::Slits { y, thickness, centers, width } => {
            (p[1] - y).abs() < thickness / 2.0
                && !centers.iter().any(|c| (p[0] - c).abs() < width / 2.0)
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
    pub fn dt(&self) -> f64 {
        stable_dt(self.dx, 0.0)
    }

    pub fn extent(&self) -> f64 {
        self.cells as f64 * self.dx
    }

    /// The field at the start: walls cut out, sponge at the edges,
    /// the source packet written.
    pub fn build(&self) -> Wave {
        let mut wave = Wave::new(self.cells, self.dx);
        for j in 0..wave.n {
            for i in 0..wave.n {
                let p = wave.position((i, j));
                if self.walls.iter().any(|w| inside(w, p)) {
                    let k = wave.index((i, j));
                    wave.keep[k] = 0.0;
                }
            }
        }
        wave.absorb_edges(self.dt());
        self.source.write(&mut wave);
        wave
    }

    /// The theory's value at x on the screen.
    pub fn theory(&self, x: f64) -> f64 {
        let screen = self.screen.expect("a theory needs a screen");
        let Theory::Fraunhofer = self.theory;
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
        let r = Realisation { particle: "electron", mass: ELECTRON_MASS, wavelength: 50e-12, reference: "" };
        let ev = ENERGY * r.energy() / ELECTRON_VOLT;
        assert!((ev - 602.0).abs() < 2.0, "{ev}");
    }
}
