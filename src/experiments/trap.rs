//! A harmonic trap: the ground state, in which nothing moves, and the
//! coherent state, which swings like a classical particle.
use std::f64::consts::PI;

use crate::experiment::{Expectation, Experiment, Measure, Particle, Stage, ELECTRON};
use crate::wave::Vec2;

/// The plane is |x|, |y| ≤ 25.6, one unit 10 nm: an electron in a
/// quantum dot, for scale only, with ħω = 0.38 meV.
const CELLS: usize = 512;
const DX: f64 = 0.1;
const SCALE: f64 = 10e-9;
/// The trap's frequency; its ground state is a Gaussian of width
/// 1/√(2ω) = 1 and its period 2π/ω. The potential is held at CAP far
/// out, where the state is negligible, so the time step stays large.
const OMEGA: f64 = 0.5;
const WIDTH: f64 = 1.0;
const PERIOD: f64 = 2.0 * PI / OMEGA;
const CAP: f64 = 50.0;
/// How far the coherent state starts from the centre.
const SWING: f64 = 4.0;

fn trap(p: Vec2) -> f64 {
    (0.5 * OMEGA * OMEGA * (p[0] * p[0] + p[1] * p[1])).min(CAP)
}

/// The ground-state Gaussian at rest, displaced to x.
const fn at_rest(x: f64) -> Particle {
    Particle { kind: ELECTRON, position: [x, 0.0], velocity: [0.0; 2], width: [WIDTH; 2], spin: [0.0; 3] }
}

fn trapped(
    name: &'static str,
    summary: &'static str,
    particle: Particle,
    duration: f64,
    expectations: Vec<Expectation>,
) -> Experiment {
    let mut expectations = expectations;
    expectations.push(Expectation {
        what: "spread of the particles along x at the end, the ground state's width",
        measure: Measure::Spread { axis: 0 },
        expected: WIDTH,
        tolerance: 0.02,
    });
    Experiment {
        name,
        summary,
        reference: "the harmonic oscillator's ground and coherent states in closed form; no \
laboratory experiment",
        stage: Stage { cells: CELLS, dx: DX, scale: SCALE, walls: vec![], potential: trap },
        particles: vec![particle],
        state: None,
        runs: 20000,
        duration,
        screen: None,
        theory: None,
        expectations,
    }
}

pub fn ground() -> Experiment {
    trapped(
        "ground",
        "\
The ground state of a harmonic trap V = ω²r²/2 with ω = 0.5, a \
Gaussian 1 wide, for one period. A stationary state has a real \
wavefunction up to a global phase, so its phase gradient is zero \
everywhere and every particle stands still, although the state's \
momentum distribution is as wide as its position distribution. The \
density should not change either.",
        at_rest(0.0),
        PERIOD,
        vec![
            Expectation {
                what: "median distance the particles moved over a period",
                measure: Measure::Moved,
                expected: 0.0,
                tolerance: 0.001,
            },
            Expectation {
                what: "field: probability left",
                measure: Measure::Norm,
                expected: 1.0,
                tolerance: 0.0001,
            },
        ],
    )
}

pub fn coherent() -> Experiment {
    trapped(
        "coherent",
        "\
The same Gaussian displaced to x = 4 and released: a coherent state, \
which swings through the trap without changing shape, its centre \
following the classical motion x = 4 cos ωt. Its phase is linear in \
position, so every particle has the same velocity as the centre and \
moves with it: after half a period each particle has moved exactly 8, \
the mean x is −4 and the spread is still 1.",
        at_rest(SWING),
        PERIOD / 2.0,
        vec![
            Expectation {
                what: "mean x of the particles after half a period",
                measure: Measure::Mean { axis: 0 },
                expected: -SWING,
                tolerance: 0.05,
            },
            Expectation {
                what: "median distance the particles moved, twice the swing",
                measure: Measure::Moved,
                expected: 2.0 * SWING,
                tolerance: 0.05,
            },
        ],
    )
}
