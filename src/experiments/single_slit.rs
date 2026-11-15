//! Electrons through one slit: the sinc² diffraction pattern.
use std::f64::consts::{PI, TAU};

use crate::experiment::{free, Expectation, Experiment, Measure, Particle, Stage, ELECTRON, SPEED};
use crate::plane::Rect;

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

/// The electron of the slit experiments: a packet 6 wide across and
/// 2 long heading up at speed 2π, its wavelength 1.
pub const BEAM: Particle = Particle {
    kind: ELECTRON,
    position: [0.0, START],
    velocity: [0.0, SPEED],
    width: [6.0, 2.0],
};

/// The plane |x|, |y| ≤ 51.2, one unit the electron's wavelength of
/// 50 pm, so the beam carries 600 eV.
pub fn stage(walls: Vec<Rect>) -> Stage {
    Stage { cells: 1024, dx: 0.1, scale: 50e-12, walls, potential: free }
}

/// A wall across the plane at height y with openings of the given
/// width at the given x, as the rectangles between the openings.
pub fn slits(y: f64, thickness: f64, centers: &[f64], width: f64) -> Vec<Rect> {
    let edges: Vec<f64> = centers.iter().flat_map(|c| [c - width / 2.0, c + width / 2.0]).collect();
    let ends = [f64::NEG_INFINITY].into_iter().chain(edges).chain([f64::INFINITY]);
    let pieces: Vec<f64> = ends.collect();
    let piece = |x: &[f64]| [[x[0], x[1]], [y - thickness / 2.0, y + thickness / 2.0]];
    pieces.chunks(2).map(piece).collect()
}

/// Fraunhofer intensity per unit x on a screen at distance d from
/// slits of the given width at the given centres. The cylindrical
/// wave falls as 1/r and meets the flat screen obliquely, which
/// together give a factor cos²θ.
pub fn fraunhofer(x: f64, d: f64, centers: &[f64], width: f64) -> f64 {
    let r = x.hypot(d);
    let (sin, cos) = (x / r, d / r);
    let beta = PI * width * sin;
    let envelope = if beta.abs() < 1e-9 { 1.0 } else { (beta.sin() / beta).powi(2) };
    let (re, im) = centers.iter().fold((0.0, 0.0), |(re, im), c| {
        let phase = TAU * c * sin;
        (re + phase.cos(), im + phase.sin())
    });
    envelope * (re * re + im * im) * cos * cos
}

fn theory(x: f64) -> f64 {
    fraunhofer(x, SCREEN - THICKNESS / 2.0, &[0.0], WIDTH)
}

pub fn single_slit() -> Experiment {
    Experiment {
        name: "single-slit",
        summary: "\
A beam 6 wide and 2 long flies up at speed 2π at a wall with one slit \
2 wavelengths wide; a screen 30 beyond it counts what comes through. \
The arrivals along the screen should follow the Fraunhofer pattern \
sinc²(πa sinθ) cos²θ, its first zeros at x = ±17.2.",
        reference: "Bach et al. 2013, 600 eV electrons; their slit was 62 nm, over a thousand \
wavelengths, so their pattern sits at angles a thousand times smaller",
        stage: stage(slits(0.0, THICKNESS, &[0.0], WIDTH)),
        particle: BEAM,
        state: None,
        start: None,
        runs: 40000,
        duration: (SCREEN - START + TAIL) / SPEED,
        screen: Some(SCREEN),
        theory: Some(theory),
        expectations: vec![
            Expectation {
                what: "share of the runs through the screen, erf(1/6√2) of the beam",
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plane::inside;

    #[test]
    fn fraunhofer_has_its_first_zero_at_the_slit_angle() {
        let zero = (0.5f64).asin().tan() * 10.0;
        assert!(fraunhofer(zero, 10.0, &[0.0], 2.0) < 1e-12);
        assert!(fraunhofer(0.0, 10.0, &[0.0], 2.0) > 0.99);
    }

    #[test]
    fn slits_leave_only_their_openings_in_the_wall() {
        let walls = slits(0.0, 0.4, &[-2.0, 2.0], 1.0);
        assert_eq!(walls.len(), 3);
        let open = |x: f64| !walls.iter().any(|w| inside(w, [x, 0.1]));
        assert!(open(-2.0) && open(2.4) && !open(0.0) && !open(-9.0) && !open(9.0));
        assert!(!walls.iter().any(|w| inside(w, [0.0, 0.3])));
    }
}
