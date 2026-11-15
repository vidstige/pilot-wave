//! A beam flying straight ahead: motion and spreading.
use crate::experiment::{free, Expectation, Experiment, Measure, Particle, Stage, ELECTRON, SPEED};

/// The plane is |x|, |y| ≤ 25.6, one unit the electron's wavelength.
const CELLS: usize = 512;
const DX: f64 = 0.1;
const SCALE: f64 = 50e-12;
/// Where the beam starts, its width along both axes, the screen it
/// crosses and how long it flies.
const START: f64 = -15.0;
const WIDTH: f64 = 2.0;
const SCREEN: f64 = 0.0;
pub const DURATION: f64 = 4.0;

/// Width of the beam after free spreading, σ(t) = σ₀ √(1 + t²/4σ₀⁴).
fn spread(t: f64) -> f64 {
    WIDTH * (1.0 + t * t / (4.0 * WIDTH.powi(4))).sqrt()
}

/// Where the beam's centre is after flying for t.
pub fn centre(t: f64) -> f64 {
    START + SPEED * t
}

/// What is expected across the motion, whatever lies ahead.
pub fn across() -> [Expectation; 2] {
    [
        Expectation {
            what: "mean x of the particles at the end",
            measure: Measure::Mean { axis: 0 },
            expected: 0.0,
            tolerance: 0.05,
        },
        Expectation {
            what: "spread of the particles across the motion, σ₀ √(1 + t²/4σ₀⁴)",
            measure: Measure::Spread { axis: 0 },
            expected: spread(DURATION),
            tolerance: 0.03,
        },
    ]
}

/// What is expected along the motion: 1% of the flight is allowed,
/// which is what a grid of ten cells per wavelength can hold.
fn along() -> [Expectation; 4] {
    [
        Expectation {
            what: "mean y of the particles at the end, −15 + 4·2π",
            measure: Measure::Mean { axis: 1 },
            expected: centre(DURATION),
            tolerance: 0.25,
        },
        Expectation {
            what: "spread of the particles along the motion, σ₀ √(1 + t²/4σ₀⁴)",
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
    ]
}

pub fn free_flight() -> Experiment {
    Experiment {
        name: "free",
        summary: "\
A beam 2 wide starts at y = −15 and flies up at speed 2π for a time \
of 4, through a screen at y = 0. Every particle should cross the \
screen, the beam should arrive where and when its speed says, and \
it should widen as a beam this narrow does.",
        reference: "free spreading of a Gaussian beam in closed form; no laboratory experiment",
        stage: Stage { cells: CELLS, dx: DX, scale: SCALE, walls: vec![], potential: free },
        particle: Particle {
            kind: ELECTRON,
            position: [0.0, START],
            velocity: [0.0, SPEED],
            width: [WIDTH; 2],
        },
        state: None,
        start: None,
        runs: 20000,
        duration: DURATION,
        screen: Some(SCREEN),
        theory: None,
        expectations: along().into_iter().chain(across()).collect(),
    }
}
