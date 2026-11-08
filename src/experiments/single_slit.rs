//! Electrons through one slit: the sinc² diffraction pattern.
use crate::experiment::{
    Expectation, Experiment, Measure, Realisation, Source, Space, Theory, Wall, ELECTRON_MASS,
    SPEED,
};
use crate::packet::Packet;

const WIDTH: f64 = 2.0;
const SCREEN: f64 = 30.0;
const START: f64 = -20.0;
/// The run lasts until the packet's centre would be this far past
/// the screen, so the slow, wide-angle tail crosses it too.
const TAIL: f64 = 25.0;

pub fn single_slit() -> Experiment {
    Experiment {
        name: "single-slit",
        summary: "\
A packet 6 wide across and 2 long heads up at speed 2π towards a wall \
0.4 thick at y = 0 with one slit 2 wavelengths wide; a screen at \
y = 30 counts what comes through. The current through the screen \
should follow the Fraunhofer pattern sinc²(πa sinθ) cos²θ, with its \
first zeros at sinθ = ±1/2, x = ±17.2 on the screen, and the \
particles should land where the current says. Every particle is an \
independent run of the one-particle experiment, since they share the \
field and never interact, so one run with many particles is many \
runs with one; the arrivals are recorded in order so the pattern can \
be watched building up one particle at a time.",
        realisation: Realisation {
            particle: "electron",
            mass: ELECTRON_MASS,
            length: 50e-12,
            reference: "Bach et al. 2013, 600 eV electrons; their slit was 62 nm, \
over a thousand wavelengths, so their pattern sits at angles a thousand times smaller",
        },
        space: Space::Plane,
        cells: 1024,
        dx: 0.1,
        walls: vec![Wall::Slits { y: 0.0, thickness: 0.4, centers: vec![0.0], width: WIDTH }],
        potentials: vec![],
        source: Source::Packet(Packet { center: [0.0, START], width: [6.0, 2.0], wavevector: [0.0, SPEED] }),
        spin: [0.0; 3],
        particles: 40000,
        duration: (SCREEN - START + TAIL) / SPEED,
        screen: Some(SCREEN),
        theory: Some(Theory::Fraunhofer),
        expectations: vec![
            Expectation {
                what: "fraction of particles through the screen, the beam's share across the slit, erf(1/6√2)",
                measure: Measure::Arrived,
                expected: 0.133,
                tolerance: 0.01,
            },
            Expectation {
                what: "arrived fraction minus the probability through the screen",
                measure: Measure::ArrivedMinusThrough,
                expected: 0.0,
                tolerance: 0.005,
            },
            Expectation {
                what: "rms gap between the current and Fraunhofer, of the peak",
                measure: Measure::CurrentVsTheory,
                expected: 0.0,
                tolerance: 0.03,
            },
            Expectation {
                what: "rms gap between the arrivals and the current, of the peak",
                measure: Measure::ArrivalsVsCurrent,
                expected: 0.0,
                tolerance: 0.03,
            },
            Expectation {
                what: "current at the first zero of sinc², x = 17.2, of the peak",
                measure: Measure::CurrentAt { x: 17.2 },
                expected: 0.0,
                tolerance: 0.02,
            },
        ],
    }
}
