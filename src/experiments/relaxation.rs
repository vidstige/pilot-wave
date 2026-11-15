//! Runs started in one corner of a box of sloshing waves forget where
//! they started: in the long run they are spread as the runs started
//! at the source are.
use std::f64::consts::{PI, TAU};

use crate::experiment::{band, free, Expectation, Experiment, Measure, Particle, Stage, ELECTRON};
use crate::packet::Complex;
use crate::plane::{Rect, Vec2};
use crate::rng::Rng;

/// The plane is |x|, |y| ≤ 6.4, one unit a nanometre; the box sits in
/// its middle, clear of the sponge.
const CELLS: usize = 128;
const DX: f64 = 0.1;
const SCALE: f64 = 1e-9;
/// The box's side and the thickness of its walls.
const SIDE: f64 = 4.0;
const WALL: f64 = 0.4;
/// Modes along each axis in the wave, sixteen in all, and the seed of
/// their phases.
const MODES: usize = 4;
const PHASES: u64 = 7;
/// Where the runs start: a square of side 1 in a corner of the box.
const CORNER: Rect = [[-1.8, -0.8], [-1.8, -0.8]];
/// About four periods of the box's slowest beat, 16/3π.
const DURATION: f64 = 40.0;
/// Side of the squares the runs are counted in.
const SQUARE: f64 = 0.5;

/// A particle at rest in the box, about 1.15 wide as a box evenly
/// filled is; the wave itself is `sloshing`.
const BOXED: Particle = Particle {
    kind: ELECTRON,
    position: [0.0; 2],
    velocity: [0.0; 2],
    width: [1.15; 2],
};

/// The first sixteen modes of the box in equal measure, each with its
/// own phase, so the density sloshes and its nodes wander.
fn sloshing([x, y]: Vec2) -> Complex {
    let mut rng = Rng::new(PHASES);
    let mode = |m: usize, x: f64| ((m + 1) as f64 * PI * (x / SIDE + 0.5)).sin();
    let inside = (x.abs().max(y.abs()) < SIDE / 2.0) as u8 as f64;
    (0..MODES * MODES).fold([0.0; 2], |[re, im], k| {
        let phase = TAU * rng.uniform();
        let a = inside * mode(k % MODES, x) * mode(k / MODES, y);
        [re + a * phase.cos(), im + a * phase.sin()]
    })
}

fn walls() -> Vec<Rect> {
    let (inner, outer) = (SIDE / 2.0, SIDE / 2.0 + WALL);
    let pair = |axis| [band(axis, inner, outer), band(axis, -outer, -inner)];
    [0, 1].into_iter().flat_map(pair).collect()
}

pub fn relaxation() -> Experiment {
    Experiment {
        name: "relaxation",
        summary: "\
A box 4 wide holds a wave made of its first sixteen modes in equal \
measure, each with its own phase, so the density sloshes about and \
the nodes of the wave wander through the box. Every other run is \
not placed where the source says but in a square of side 1 in one \
corner. Each is carried round the wandering nodes, neighbours are \
pulled apart, and the ensemble is stirred through the box until, \
counted in squares of side 0.5, the runs from the corner are spread \
as the runs from the source are. So where particles are found need \
not be put in at the start; in a wave busy enough it comes about, \
and stays.",
        reference: "Valentini and Westman 2005, Proc. R. Soc. A 461, 253: sixteen modes of a \
box, an ensemble started out of step with |ψ|² relaxing to it; no laboratory experiment",
        stage: Stage { cells: CELLS, dx: DX, scale: SCALE, walls: walls(), potential: free },
        particle: BOXED,
        state: Some(sloshing),
        start: Some(CORNER),
        runs: 20000,
        duration: DURATION,
        screen: None,
        theory: None,
        expectations: vec![Expectation {
            what: "share of the runs from the corner misplaced at the end, in squares of side 0.5",
            measure: Measure::Misplaced { side: SQUARE },
            expected: 0.0,
            tolerance: 0.1,
        }],
    }
}
