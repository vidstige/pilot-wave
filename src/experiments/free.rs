//! A packet travelling straight ahead, and the same packet meeting a
//! wall: the sanity checks of motion, spreading and reflection.
use crate::experiment::{
    band, free, Expectation, Experiment, Measure, Realisation, ELECTRON_MASS, SPEED,
};
use crate::packet::Packet;
use crate::wave::{effective_mass, group_velocity, Wave};

/// The plane is |x|, |y| ≤ 25.6.
const CELLS: usize = 512;
const DX: f64 = 0.1;
/// Where the packet starts, its width along both axes, the screen it
/// crosses and how long it flies.
const START: f64 = -15.0;
const WIDTH: f64 = 2.0;
const SCREEN: f64 = 0.0;
const DURATION: f64 = 4.0;
/// The wall the mirrored packet meets, and the screen behind it.
const WALL: f64 = 5.0;
const BEHIND: f64 = 8.0;

pub const ELECTRON: Realisation = Realisation {
    particle: "electron",
    mass: ELECTRON_MASS,
    length: 50e-12,
    reference: "a 600 eV electron, as in the slit experiments; no particular experiment",
};

fn packet(wave: &mut Wave) {
    Packet { center: [0.0, START], width: [WIDTH; 2], wavevector: [0.0, SPEED] }.write(wave)
}

/// Width of |ψ|² after free spreading, σ(t) = σ₀ √(1 + t²/4m²σ₀⁴),
/// with the mass the grid gives a packet along its motion.
fn spread(t: f64, mass: f64) -> f64 {
    WIDTH * (1.0 + t * t / (4.0 * mass * mass * WIDTH.powi(4))).sqrt()
}

/// Where the packet's centre is after flying for t at the grid's
/// speed for its wavenumber.
fn centre(t: f64) -> f64 {
    START + group_velocity(SPEED, DX) * t
}

fn flight(
    name: &'static str,
    summary: &'static str,
    walls: Vec<[[f64; 2]; 2]>,
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
            expected: spread(DURATION, 1.0),
            tolerance: 0.03,
        },
        Expectation {
            what: "spread of the field across the motion, likewise",
            measure: Measure::FieldSpread { axis: 0 },
            expected: spread(DURATION, 1.0),
            tolerance: 0.01,
        },
    ]);
    Experiment {
        name,
        summary,
        realisation: ELECTRON,
        cells: CELLS,
        dx: DX,
        walls,
        potential: free,
        source: packet,
        spin: [0.0; 3],
        particles: 20000,
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
their spread grows as the free spreading of |ψ|², \
σ(t) = σ₀ √(1 + t²/4σ₀⁴), each particle keeping its place in the \
packet as the trails fan out in straight lines. The grid has its own \
dispersion: at ten cells per wavelength a wave moves 0.5% slower \
than k and spreads along its motion as if 2.5% heavier, and the \
expectations along y allow for that. The particles' spread along the \
motion runs 2.5% above the field's at this resolution, ten cells per \
wavelength, and matches it at twenty. Nothing should be lost: the \
probability stays 1.",
        vec![],
        SCREEN,
        vec![
            Expectation {
                what: "mean y of the particles at the end, −15 + 4·2π at the grid's speed",
                measure: Measure::Mean { axis: 1 },
                expected: centre(DURATION),
                tolerance: 0.05,
            },
            Expectation {
                what: "spread of the particles along the motion, closed form with the grid's mass",
                measure: Measure::Spread { axis: 1 },
                expected: spread(DURATION, effective_mass(SPEED, DX)),
                tolerance: 0.08,
            },
            Expectation {
                what: "spread of the field along the motion, likewise",
                measure: Measure::FieldSpread { axis: 1 },
                expected: spread(DURATION, effective_mass(SPEED, DX)),
                tolerance: 0.03,
            },
            Expectation {
                what: "probability left in the field",
                measure: Measure::Norm,
                expected: 1.0,
                tolerance: 0.001,
            },
            Expectation {
                what: "fraction of the particles through the screen",
                measure: Measure::Arrived,
                expected: 1.0,
                tolerance: 0.001,
            },
            Expectation {
                what: "mean arrival time at the screen, 15 over the grid's speed",
                measure: Measure::MeanArrivalTime,
                expected: -START / group_velocity(SPEED, DX),
                tolerance: 0.02,
            },
            Expectation {
                what: "rms gap between the arrivals and the current, of the peak",
                measure: Measure::ArrivalsVsCurrent,
                expected: 0.0,
                tolerance: 0.03,
            },
        ],
    )
}

pub fn mirror() -> Experiment {
    flight(
        "mirror",
        "\
The same packet meets a hard wall across the plane at y = 5 and comes \
back: after a time of 4 its mean y should be the mirror image of the \
free flight's in the wall, with the spread across the motion \
as in free flight, nothing lost through the wall and no particle \
beyond it. Where the incoming and the reflected wave overlap they \
make a standing wave, and the particles slow and turn in its \
near-nodes rather than crossing them.",
        vec![band(1, WALL, WALL + 0.4)],
        BEHIND,
        vec![
            Expectation {
                what: "mean y of the particles at the end, mirrored in the wall",
                measure: Measure::Mean { axis: 1 },
                expected: 2.0 * WALL - centre(DURATION),
                tolerance: 0.1,
            },
            Expectation {
                what: "probability left in the field",
                measure: Measure::Norm,
                expected: 1.0,
                tolerance: 0.005,
            },
            Expectation {
                what: "fraction of the particles beyond the wall",
                measure: Measure::Arrived,
                expected: 0.0,
                tolerance: 0.0001,
            },
        ],
    )
}
