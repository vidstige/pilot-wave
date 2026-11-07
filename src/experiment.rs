//! Named experiments: a field with walls and potential, a packet, the
//! kind of particles, a detection line and how long to run. No I/O.
//!
//! The packet's wavelength is 1, so k = 2π and the packet speed is 2π.
use std::f64::consts::{PI, TAU};

use crate::grid::Vec2;
use crate::packet::Packet;
use crate::particles::Kind;
use crate::wave::{stable_dt, Wave};

pub const NAMES: [&str; 8] = [
    "free", "double-slit", "single-slit", "barrier", "disc", "harmonic", "ground", "box",
];
pub const SPEED: f64 = TAU;
pub const ENERGY: f64 = SPEED * SPEED / 2.0;
const DX: f64 = 0.1;
/// Where the packet starts and where arrivals are counted.
const START: f64 = -20.0;
const SCREEN: f64 = 30.0;
/// Walls at y = 0 with this thickness.
const WALL: f64 = 0.4;
pub const SLIT_WIDTH: f64 = 1.0;
pub const SLIT_SEPARATION: f64 = 4.0;
pub const SINGLE_SLIT_WIDTH: f64 = 2.0;
/// Barrier height relative to the packet's energy, and its width.
pub const BARRIER_HEIGHT: f64 = 1.2;
pub const BARRIER_WIDTH: f64 = 0.3;
pub const DISC_RADIUS: f64 = 2.0;
/// Trap frequency and the orbit radius of the coherent state.
pub const TRAP: f64 = 0.5;
const ORBIT: f64 = 4.0;
const BOX_HALF: f64 = 6.0;

pub struct Experiment {
    pub name: &'static str,
    pub wave: Wave,
    pub packet: Packet,
    pub kind: Kind,
    pub dt: f64,
    pub duration: f64,
    /// y of the line where arrivals are counted, if any.
    pub screen: Option<f64>,
}

fn upward(center: Vec2, width: Vec2) -> Packet {
    Packet { center, width, wavevector: [0.0, SPEED] }
}

/// Zero the field wherever the predicate holds.
fn wall(wave: &mut Wave, inside: impl Fn(Vec2) -> bool) {
    for j in 0..wave.n {
        for i in 0..wave.n {
            if inside(wave.position((i, j))) {
                let k = wave.index((i, j));
                wave.keep[k] = 0.0;
            }
        }
    }
}

fn potential(wave: &mut Wave, v: impl Fn(Vec2) -> f64) {
    for j in 0..wave.n {
        for i in 0..wave.n {
            let k = wave.index((i, j));
            wave.potential[k] = v(wave.position((i, j))) as f32;
        }
    }
}

fn in_wall(p: Vec2) -> bool {
    p[1].abs() < WALL / 2.0
}

fn slits(centers: &[f64], width: f64) -> impl Fn(Vec2) -> bool + '_ {
    move |p| in_wall(p) && !centers.iter().any(|c| (p[0] - c).abs() < width / 2.0)
}

/// How far the packet travels to the screen and back off the walls.
fn travel(wave: &Wave) -> f64 {
    (SCREEN - START + 10.0).min(wave.extent()) / SPEED
}

pub fn experiment(name: &'static str, n: usize, kind: Kind) -> Experiment {
    let mut wave = Wave::new(n, DX);
    let beam = |wave: &Wave| upward([0.0, START.max(-wave.extent() / 2.0 + 6.0)], [6.0, 2.0]);
    let (packet, potential_max, screen) = match name {
        "free" => (upward([0.0, START], [1.5, 1.5]), 0.0, None),
        "double-slit" => {
            let d = SLIT_SEPARATION / 2.0;
            wall(&mut wave, slits(&[-d, d], SLIT_WIDTH));
            (beam(&wave), 0.0, Some(SCREEN))
        }
        "single-slit" => {
            wall(&mut wave, slits(&[0.0], SINGLE_SLIT_WIDTH));
            (beam(&wave), 0.0, Some(SCREEN))
        }
        "barrier" => {
            let height = BARRIER_HEIGHT * ENERGY;
            potential(&mut wave, |p| if p[1].abs() < BARRIER_WIDTH / 2.0 { height } else { 0.0 });
            (upward([0.0, START], [4.0, 2.0]), height, Some(SCREEN))
        }
        "disc" => {
            wall(&mut wave, |p| p[0].hypot(p[1]) < DISC_RADIUS);
            (upward([0.0, START], [3.0, 2.0]), 0.0, Some(SCREEN))
        }
        "harmonic" => {
            potential(&mut wave, |p| 0.5 * TRAP * TRAP * (p[0] * p[0] + p[1] * p[1]));
            let width = 1.0 / (2.0 * TRAP).sqrt();
            let packet = Packet { center: [ORBIT, 0.0], width: [width; 2], wavevector: [0.0, TRAP * ORBIT] };
            (packet, 0.5 * TRAP * TRAP * (wave.extent() / 2.0).powi(2), None)
        }
        "ground" => {
            potential(&mut wave, |p| 0.5 * TRAP * TRAP * (p[0] * p[0] + p[1] * p[1]));
            let width = 1.0 / (2.0 * TRAP).sqrt();
            let packet = Packet { center: [0.0; 2], width: [width; 2], wavevector: [0.0; 2] };
            (packet, 0.5 * TRAP * TRAP * (wave.extent() / 2.0).powi(2), None)
        }
        "box" => {
            wall(&mut wave, |p| p[0].abs() > BOX_HALF || p[1].abs() > BOX_HALF);
            (Packet { center: [-3.0, 0.0], width: [1.0, 1.0], wavevector: [SPEED / 2.0, SPEED / 4.0] }, 0.0, None)
        }
        other => panic!("unknown experiment {other}"),
    };
    let dt = stable_dt(DX, potential_max);
    if !matches!(name, "box") {
        wave.absorb_edges(dt);
    }
    packet.write(&mut wave);
    let duration = match name {
        "harmonic" | "ground" => 2.0 * TAU / TRAP,
        "box" | "free" => 6.0,
        _ => travel(&wave),
    };
    Experiment { name, wave, packet, kind, dt, duration, screen }
}

/// Fraunhofer intensity on the screen for slits of the given width at
/// the given centres, as a function of x on the screen.
pub fn fraunhofer(x: f64, centers: &[f64], width: f64) -> f64 {
    let sin = x / x.hypot(SCREEN - WALL / 2.0);
    let beta = PI * width * sin;
    let envelope = if beta.abs() < 1e-9 { 1.0 } else { (beta.sin() / beta).powi(2) };
    let (re, im) = centers.iter().fold((0.0, 0.0), |(re, im), c| {
        let phase = TAU * c * sin;
        (re + phase.cos(), im + phase.sin())
    });
    envelope * (re * re + im * im) / (centers.len() * centers.len()) as f64
}
