//! What a run prints and writes: the description, the units, the
//! expectations, then frames, trails, the screen and the verdicts.
use std::fs;

use crate::experiment::{Experiment, Expectation, ELECTRON_VOLT, ENERGY, SPEED};
use crate::particles::Particles;
use crate::run::{measure, Outcome, BINS};
use crate::wave::{Vec2, Wave};
use crate::{plot, png, view};

pub const PIXELS: usize = 800;
pub const FRAME_EVERY: f64 = 0.5;
const TRAILS: usize = 300;
/// Particles drawn in a frame, so the field stays visible under them.
const DOTS: usize = 2000;
const PLOT: (usize, usize) = (800, 300);

pub fn describe(e: &Experiment) {
    let r = &e.realisation;
    println!("{}\n\n{}\n", e.name, e.summary);
    println!("Realisation: {} with de Broglie wavelength {:.3e} m ({}).", r.particle, r.wavelength, r.reference);
    println!("One unit of length is {:.3e} m, of time {:.3e} s, of speed {:.3e} m/s, of energy {:.3e} J.",
        r.length(), r.time(), r.speed(), r.energy());
    println!("The packet moves at {:.3e} m/s with kinetic energy {:.1} eV; the domain is {:.3e} m across; the run lasts {:.3e} s.\n",
        SPEED * r.speed(), ENERGY * r.energy() / ELECTRON_VOLT, e.extent() * r.length(), e.duration * r.time());
    println!("Grid {} x {} cells of {}, dt {:.2e}, {} particles, duration {:.2}.", e.cells, e.cells, e.dx, e.dt(), e.particles, e.duration);
    println!("\nExpectations:");
    for x in &e.expectations {
        println!("  {:>8.4} ± {:<7.4} {}", x.expected, x.tolerance, x.what);
    }
    println!();
}

pub fn frame(e: &Experiment, wave: &Wave, particles: &Particles, path: &str) {
    let peak = e.source.peak();
    let mut image = view::field(wave, PIXELS, view::Style::Phase, peak);
    let heads: Vec<Vec<Vec2>> = particles.positions.iter().take(DOTS).map(|p| vec![*p]).collect();
    view::particles(&mut image, wave, &heads);
    png::write(path, &image).expect("write frame");
}

fn write_trails(e: &Experiment, o: &Outcome, trails: &[Vec<Vec2>], out: &str) {
    let mut glow = Wave::new(o.wave.n, o.wave.dx);
    glow.re = o.glow.iter().map(|g| g.sqrt()).collect();
    glow.keep = o.wave.keep.clone();
    let mut image = view::field(&glow, PIXELS, view::Style::Density, glow.peak());
    view::particles(&mut image, &o.wave, trails);
    png::write(format!("{out}/{}_trails.png", e.name), &image).expect("write trails");
}

fn write_screen(e: &Experiment, o: &Outcome, out: &str) {
    let Some(screen) = &o.screen else { return };
    let total: f64 = screen.flux.iter().sum();
    let scale = |v: Vec<f64>| {
        let sum: f64 = v.iter().sum::<f64>().max(1e-300);
        v.iter().map(|x| x * total / sum).collect::<Vec<f64>>()
    };
    let hits = scale(screen.hits.clone());
    let theory = scale((0..BINS).map(|b| e.theory(screen.x(b))).collect());
    let series: [(&[f64], _); 3] =
        [(&hits, [255, 255, 255]), (&screen.flux, [80, 200, 255]), (&theory, [255, 160, 30])];
    png::write(format!("{out}/{}_screen.png", e.name), &plot::curves(PLOT.0, PLOT.1, &series))
        .expect("write screen plot");
    let rows: String = (0..BINS)
        .map(|b| format!("{},{},{},{}\n", screen.x(b), hits[b], screen.flux[b], theory[b]))
        .collect();
    fs::write(format!("{out}/{}_screen.csv", e.name), format!("x,hits,flux,theory\n{rows}"))
        .expect("write screen table");
    let arrivals: String = screen.arrivals.iter().map(|(t, x)| format!("{t},{x}\n")).collect();
    fs::write(format!("{out}/{}_arrivals.csv", e.name), format!("t,x\n{arrivals}"))
        .expect("write arrivals");
}

/// Print every expectation against its measurement.
pub fn verdicts(e: &Experiment, o: &Outcome) -> bool {
    println!("Results:");
    let mut all = true;
    for Expectation { what, measure: m, expected, tolerance } in &e.expectations {
        let value = measure(e, o, m);
        let ok = (value - expected).abs() <= *tolerance;
        all &= ok;
        println!("  {} {:>8.4} expected {:>8.4} ± {:<7.4} {}", if ok { "pass" } else { "FAIL" }, value, expected, tolerance, what);
    }
    all
}

pub fn write(e: &Experiment, o: &Outcome, trails: &[Vec<Vec2>], out: &str) {
    write_trails(e, o, &trails[..trails.len().min(TRAILS)], out);
    write_screen(e, o, out);
}
