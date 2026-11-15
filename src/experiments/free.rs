//! A packet travelling straight ahead: the sanity check of motion
//! and spreading.
use crate::experiment::{free, Expectation, Experiment, Measure, Particle, Stage, ELECTRON, SPEED};
use crate::plane::Rect;

/// The plane is |x|, |y| ≤ 25.6, one unit the electron's wavelength.
const CELLS: usize = 512;
const DX: f64 = 0.1;
const SCALE: f64 = 50e-12;
/// Where the packet starts, its width along both axes, the screen it
/// crosses and how long it flies.
const START: f64 = -15.0;
const WIDTH: f64 = 2.0;
const SCREEN: f64 = 0.0;
pub const DURATION: f64 = 4.0;

const PACKET: Particle = Particle {
    kind: ELECTRON,
    position: [0.0, START],
    velocity: [0.0, SPEED],
    width: [WIDTH; 2],
};

/// Width of the beam after free spreading, σ(t) = σ₀ √(1 + t²/4σ₀⁴).
fn spread(t: f64) -> f64 {
    WIDTH * (1.0 + t * t / (4.0 * WIDTH.powi(4))).sqrt()
}

/// Where the beam's centre is after flying for t.
pub fn centre(t: f64) -> f64 {
    START + SPEED * t
}

/// The packet's flight past the given walls and through a screen:
/// what is expected along the motion, and across it the same
/// whatever the walls.
pub fn flight(
    name: &'static str,
    summary: &'static str,
    walls: Vec<Rect>,
    screen: f64,
    expectations: Vec<Expectation>,
) -> Experiment {
    let mut expectations = expectations;
    expectations.extend([
        Expectation {
            what: "mean x of the particles at the end",
            measure: Measure::Mean { axis: 0 },
            expected: 0.0,
            tolerance: 0.05,
        },
        Expectation {
            what: "spread of the particles across the motion, free spreading in closed form",
            measure: Measure::Spread { axis: 0 },
            expected: spread(DURATION),
            tolerance: 0.03,
        },
    ]);
    Experiment {
        name,
        summary,
        reference: "free spreading of a Gaussian packet in closed form; no laboratory \
experiment",
        stage: Stage { cells: CELLS, dx: DX, scale: SCALE, walls, potential: free },
        particle: PACKET,
        state: None,
        start: None,
        runs: 20000,
        duration: DURATION,
        screen: Some(screen),
        theory: None,
        expectations,
    }
}

pub fn free_flight() -> Experiment {
    flight(
        "free",
        "\
A packet 2 wide along both axes starts at y = −15 and heads up at \
speed 2π for a time of 4, crossing a screen at y = 0 on the way. The \
particles should move with it: their mean y ends at −15 + 4·2π, \
every one of them crosses the screen, at a mean time of 15/2π, and \
their spread grows as a beam this narrow must, \
σ(t) = σ₀ √(1 + t²/4σ₀⁴). The tolerances along the motion are 1% of \
the flight, which is what a grid of ten cells per wavelength can \
hold.",
        vec![],
        SCREEN,
        vec![
            Expectation {
                what: "mean y of the particles at the end, −15 + 4·2π",
                measure: Measure::Mean { axis: 1 },
                expected: centre(DURATION),
                tolerance: 0.25,
            },
            Expectation {
                what: "spread of the particles along the motion, free spreading in closed form",
                measure: Measure::Spread { axis: 1 },
                expected: spread(DURATION),
                tolerance: 0.08,
            },
            Expectation {
                what: "fraction of the particles through the screen",
                measure: Measure::Arrived,
                expected: 1.0,
                tolerance: 0.001,
            },
            Expectation {
                what: "mean arrival time at the screen, 15/2π",
                measure: Measure::MeanArrivalTime,
                expected: -START / SPEED,
                tolerance: 0.03,
            },
        ],
    )
}
