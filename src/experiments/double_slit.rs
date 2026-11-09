//! Electrons through two slits: cos² fringes under the sinc² envelope.
use crate::experiment::{fraunhofer, slits, Expectation, Experiment, Measure, SPEED};
use crate::experiments::single_slit::{stage, BEAM};

const WIDTH: f64 = 1.0;
const SEPARATION: f64 = 4.0;
const THICKNESS: f64 = 0.4;
const SCREEN: f64 = 30.0;
const START: f64 = -20.0;
const TAIL: f64 = 25.0;
/// On the screen, the first-order maximum, d sinθ = 1, and the first
/// minimum, d sinθ = 1/2, of slits 4 apart 29.8 before the screen.
const FIRST_MAXIMUM: f64 = 7.7;
const FIRST_MINIMUM: f64 = 3.8;
/// Half an order on the screen: the first maximum is sought within
/// this of where it should be.
const HALF_ORDER: f64 = 3.5;

fn theory(x: f64) -> f64 {
    fraunhofer(x, SCREEN - THICKNESS / 2.0, &[-SEPARATION / 2.0, SEPARATION / 2.0], WIDTH)
}

pub fn double_slit() -> Experiment {
    Experiment {
        name: "double-slit",
        summary: "\
The same beam as the single slit meets a wall with two slits 1 wide and \
4 apart, centre to centre. The arrivals should show cos²(πd sinθ) \
fringes under the sinc²(πa sinθ) envelope, orders at sinθ = m/4, so \
maxima at x = 0, ±7.7, ±17 and minima between them. In the trails each \
particle goes through one slit or the other and never crosses the \
axis of symmetry, which no detector could see without destroying the \
fringes, so it is not expected, only shown. The arrivals are recorded \
in order, so the fringes can be watched building up one particle at \
a time.",
        reference: "Bach et al. 2013, 600 eV electrons through slits 62 nm wide and 272 nm apart: \
d/a = 4.4 against our 4, so about as many fringes fit under the envelope, but theirs are \
0.18 mrad apart against our 0.25 rad",
        stage: stage(slits(0.0, THICKNESS, &[-SEPARATION / 2.0, SEPARATION / 2.0], WIDTH)),
        particles: vec![BEAM],
        state: None,
        runs: 60000,
        duration: (SCREEN - START + TAIL) / SPEED,
        screen: Some(SCREEN),
        theory: Some(theory),
        expectations: vec![
            Expectation {
                what: "fraction of runs through the screen, the beam's share across both slits",
                measure: Measure::Arrived,
                expected: 0.126,
                tolerance: 0.01,
            },
            Expectation {
                what: "rms gap between the arrivals and Fraunhofer, of the peak",
                measure: Measure::ArrivalsVsTheory,
                expected: 0.0,
                tolerance: 0.05,
            },
            Expectation {
                what: "x of the first-order maximum of the arrivals, d sinθ = 1",
                measure: Measure::ArrivalsPeakNear { x: FIRST_MAXIMUM, within: HALF_ORDER },
                expected: FIRST_MAXIMUM,
                tolerance: 0.8,
            },
            Expectation {
                what: "arrivals at the first minimum, d sinθ = 1/2, of the peak",
                measure: Measure::ArrivalsAt { x: FIRST_MINIMUM },
                expected: 0.0,
                tolerance: 0.1,
            },
        ],
    }
}
