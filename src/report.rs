//! What a run prints and writes: the description, the units, the
//! expectations, then frames, trails, the screen, the runs and the
//! verdicts.
use std::fs;

use crate::experiment::{Experiment, Expectation};
use crate::particles::Particles;
use crate::run::{measure, Outcome, BINS};
use crate::wave::{Vec2, Wave};
use crate::{plot, png, view};

pub const PIXELS: usize = 800;
pub const FRAME_EVERY: f64 = 0.5;
const TRAILS: usize = 300;
/// Trails written as a table, and the steps between their rows.
const TRAIL_ROWS: usize = 40;
const TRAIL_EVERY: usize = 20;
/// Particles drawn in a frame, so the field stays visible under them.
const DOTS: usize = 2000;
const PLOT: (usize, usize) = (800, 300);

pub fn describe(e: &Experiment) {
    let r = &e.realisation;
    println!("{}\n\n{}\n", e.name, e.summary);
    println!("Realisation: {} ({}).", r.particle, r.reference);
    println!("One unit of length is {:.3e} m, of time {:.3e} s, of speed {:.3e} m/s, of energy {:.3e} J.",
        r.length, r.time(), r.speed(), r.energy());
    println!("The domain is {:.3e} m across; the run lasts {:.3e} s.\n", e.extent() * r.length, e.duration * r.time());
    println!("Grid {} x {} cells of {}, dt {:.2e}, {} particles, duration {:.2}.", e.cells, e.cells, e.dx, e.dt(), e.particles, e.duration);
    println!("\nExpectations:");
    for x in &e.expectations {
        println!("  {:>8.4} ± {:<7.4} {}", x.expected, x.tolerance, x.what);
    }
    println!();
}

/// The field with the particles on it, densities relative to `peak`.
pub fn frame(wave: &Wave, particles: &Particles, peak: f32, path: &str) {
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

/// The first trails as rows of trail, time and position.
fn write_trail_table(e: &Experiment, trails: &[Vec<Vec2>], out: &str) {
    let dt = e.dt();
    let row = |(k, trail): (usize, &Vec<Vec2>)| -> String {
        let points = trail.iter().enumerate().step_by(TRAIL_EVERY);
        points.map(|(s, p)| format!("{k},{},{},{}\n", s as f64 * dt, p[0], p[1])).collect()
    };
    let rows: String = trails.iter().take(TRAIL_ROWS).enumerate().map(row).collect();
    fs::write(format!("{out}/{}_trails.csv", e.name), format!("trail,t,x,y\n{rows}")).expect("write trails");
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

/// Where each particle started and ended, and ended in the twin if
/// there is one.
fn write_runs(e: &Experiment, o: &Outcome, out: &str) {
    let twin = o.twin.as_ref().map(|t| &t.particles.positions);
    let row = |k: usize| {
        let (s, p) = (o.start[k], o.particles.positions[k]);
        let pair = twin.map_or(String::new(), |t| format!(",{},{}", t[k][0], t[k][1]));
        format!("{},{},{},{}{pair}\n", s[0], s[1], p[0], p[1])
    };
    let rows: String = (0..o.start.len()).map(row).collect();
    let head = if twin.is_some() { "x_start,y_start,x,y,x_twin,y_twin" } else { "x_start,y_start,x,y" };
    fs::write(format!("{out}/{}_runs.csv", e.name), format!("{head}\n{rows}")).expect("write runs");
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
    write_trail_table(e, trails, out);
    write_screen(e, o, out);
    write_runs(e, o, out);
}
