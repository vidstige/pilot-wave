//! Electrons through one slit: the sinc² diffraction pattern.
use crate::experiment::{
    fraunhofer, free, slits, Expectation, Experiment, Measure, Realisation, ELECTRON_MASS, SPEED,
};
use crate::packet::Packet;
use crate::wave::Wave;

const WIDTH: f64 = 2.0;
/// The wall at y = 0, this thick, and the screen.
const THICKNESS: f64 = 0.4;
const SCREEN: f64 = 30.0;
/// Where the beam starts, and how far past the screen its centre
/// would be at the end, so the slow, wide-angle tail crosses it too.
const START: f64 = -20.0;
const TAIL: f64 = 25.0;
/// Where the first zero of sinc²(πa sinθ), sinθ = 1/2, meets the screen.
const FIRST_ZERO: f64 = 17.2;

/// The beam of the slit experiments: a packet 6 wide across and 2
/// long heading up at speed 2π.
pub fn beam(wave: &mut Wave) {
    Packet { center: [0.0, START], width: [6.0, 2.0], wavevector: [0.0, SPEED] }.write(wave)
}

fn theory(x: f64) -> f64 {
    fraunhofer(x, SCREEN - THICKNESS / 2.0, &[0.0], WIDTH)
}

pub fn single_slit() -> Experiment {
    Experiment {
        name: "single-slit",
        summary: "\
A packet 6 wide across and 2 long heads up at speed 2π towards a wall \
0.4 thick at y = 0 with one slit 2 wavelengths wide; a screen at \
y = 30 counts what comes through. The arrivals along the screen, the \
only thing a detector records, should follow the Fraunhofer pattern \
sinc²(πa sinθ) cos²θ, with its first zeros at sinθ = ±1/2, x = ±17.2 \
on the screen. Every particle is an independent run of the \
one-particle experiment, since they share the field and never \
interact, so one run with many particles is many runs with one; the \
arrivals are recorded in order so the pattern can be watched \
building up one particle at a time.",
        realisation: Realisation {
            particle: "electron",
            mass: ELECTRON_MASS,
            length: 50e-12,
            reference: "Bach et al. 2013, 600 eV electrons; their slit was 62 nm, \
over a thousand wavelengths, so their pattern sits at angles a thousand times smaller",
        },
        cells: 1024,
        dx: 0.1,
        walls: slits(0.0, THICKNESS, &[0.0], WIDTH),
        potential: free,
        source: beam,
        spin: [0.0; 3],
        particles: 40000,
        duration: (SCREEN - START + TAIL) / SPEED,
        screen: Some(SCREEN),
        theory: Some(theory),
        expectations: vec![
            Expectation {
                what: "fraction of particles through the screen, the beam's share across the slit, erf(1/6√2)",
                measure: Measure::Arrived,
                expected: 0.133,
                tolerance: 0.01,
            },
            Expectation {
                what: "rms gap between the arrivals and Fraunhofer, of the peak",
                measure: Measure::ArrivalsVsTheory,
                expected: 0.0,
                tolerance: 0.04,
            },
            Expectation {
                what: "arrivals at the first zero of sinc², x = 17.2, of the peak",
                measure: Measure::ArrivalsAt { x: FIRST_ZERO },
                expected: 0.0,
                tolerance: 0.03,
            },
        ],
    }
}
