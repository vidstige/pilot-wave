//! Electrons through two slits: cos² fringes under the sinc² envelope.
use crate::experiment::{
    Expectation, Experiment, Measure, Realisation, Source, Space, Theory, Wall, ELECTRON_MASS,
    SPEED,
};
use crate::packet::Packet;

const WIDTH: f64 = 1.0;
const SEPARATION: f64 = 4.0;
const SCREEN: f64 = 30.0;
const START: f64 = -20.0;
const TAIL: f64 = 25.0;
/// On the screen, the first-order maximum, d sinθ = 1, and the first
/// minimum, d sinθ = 1/2, of slits 4 apart 29.8 before the screen.
const FIRST_MAXIMUM: f64 = 7.7;
const FIRST_MINIMUM: f64 = 3.8;

pub fn double_slit() -> Experiment {
    Experiment {
        name: "double-slit",
        summary: "\
The same beam as the single slit meets a wall with two slits 1 wide and \
4 apart, centre to centre. The current through the screen should show \
cos²(πd sinθ) fringes under the sinc²(πa sinθ) envelope, orders at \
sinθ = m/4, so maxima at x = 0, ±7.7, ±17 and minima between them, \
and the particles should land on the bright fringes and avoid the \
dark ones, each one going through one slit or the other and never \
crossing the axis of symmetry. The arrivals are recorded in order, so \
the fringes can be watched building up one particle at a time.",
        realisation: Realisation {
            particle: "electron",
            mass: ELECTRON_MASS,
            length: 50e-12,
            reference: "Bach et al. 2013, 600 eV electrons through slits 62 nm wide and 272 nm \
apart: d/a = 4.4 against our 4, so about as many fringes fit under the envelope, but theirs \
are 0.18 mrad apart against our 0.25 rad",
        },
        space: Space::Plane,
        cells: 1024,
        dx: 0.1,
        walls: vec![Wall::Slits {
            y: 0.0,
            thickness: 0.4,
            centers: vec![-SEPARATION / 2.0, SEPARATION / 2.0],
            width: WIDTH,
        }],
        potentials: vec![],
        source: Source::Packet(Packet { center: [0.0, START], width: [6.0, 2.0], wavevector: [0.0, SPEED] }),
        spin: [0.0; 3],
        particles: 60000,
        duration: (SCREEN - START + TAIL) / SPEED,
        screen: Some(SCREEN),
        theory: Some(Theory::Fraunhofer),
        expectations: vec![
            Expectation {
                what: "fraction of particles through the screen, the beam's share across both slits",
                measure: Measure::Arrived,
                expected: 0.126,
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
                tolerance: 0.035,
            },
            Expectation {
                what: "rms gap between the arrivals and the current, of the peak",
                measure: Measure::ArrivalsVsCurrent,
                expected: 0.0,
                tolerance: 0.03,
            },
            Expectation {
                what: "x of the first-order maximum, d sinθ = 1",
                measure: Measure::MaximumNear { x: FIRST_MAXIMUM },
                expected: FIRST_MAXIMUM,
                tolerance: 0.8,
            },
            Expectation {
                what: "current at the first minimum, d sinθ = 1/2, of the peak",
                measure: Measure::CurrentAt { x: FIRST_MINIMUM },
                expected: 0.0,
                tolerance: 0.1,
            },
        ],
    }
}
