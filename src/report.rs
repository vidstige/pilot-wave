//! What a run prints and writes: the description, the units, the
//! expectations, then frames, trails, the screen, the runs and the
//! verdicts.
use std::fs;

use crate::engine::Engine;
use crate::experiment::{Experiment, Expectation};
use crate::measure::measure;
use crate::plane::Vec2;
use crate::run::{Outcome, BINS};
use crate::{plot, png, view};

pub const PIXELS: usize = 800;
pub const FRAME_EVERY: f64 = 0.5;
const TRAILS: usize = 300;
/// Time between the points of a trail: a step of the wave, so a
/// medium with a finer step keeps trails no longer.
const TRAIL_STEP: f64 = 0.002;
/// Trails written as a table, and the points between their rows.
const TRAIL_ROWS: usize = 40;
const TRAIL_EVERY: usize = 20;
/// Particles drawn in a frame, so the field stays visible under them.
const DOTS: usize = 2000;
const PLOT: (usize, usize) = (800, 300);

pub fn describe(e: &Experiment, engine: &str, dt: f64) {
    let stage = &e.stage;
    println!("{}\n\n{}\n\nReference: {}.\n", e.name, e.summary, e.reference);
    let p = &e.particle;
    println!("Particle: {} of {:.3e} kg at ({}, {}), moving at ({}, {}), its packet {} by {} wide.",
        p.kind.name, p.kind.mass, p.position[0], p.position[1], p.velocity[0], p.velocity[1],
        p.width[0], p.width[1]);
    println!("One unit of length is {:.3e} m, of time {:.3e} s, of speed {:.3e} m/s, of energy {:.3e} J.",
        stage.scale, e.time(), e.speed(), e.energy());
    println!("The stage is {:.3e} m across; a run lasts {:.3e} s.\n", stage.extent() * stage.scale, e.duration * e.time());
    println!("Grid {} x {} cells of {}, dt {:.2e}, {} runs, duration {:.2}.", stage.cells, stage.cells, stage.dx, dt, e.runs, e.duration);
    println!("Engine: {engine}.");
    println!("\nExpectations:");
    for x in &e.expectations {
        println!("  {:>8.4} ± {:<7.4} {}", x.expected, x.tolerance, x.what);
    }
    println!();
}

/// The brightest the engine's field is.
pub fn peak(engine: &dyn Engine) -> f32 {
    let [re, im] = engine.field().unwrap_or([&[], &[]]);
    re.iter().zip(im).map(|(re, im)| re * re + im * im).fold(f32::MIN_POSITIVE, f32::max)
}

/// The engine's field, dark if it has none, with the particles on
/// it, densities relative to `peak`.
pub fn frame(e: &Experiment, engine: &dyn Engine, peak: f32, path: &str) {
    let dark = vec![0.0; e.stage.cells * e.stage.cells];
    let field = engine.field().unwrap_or([&dark, &dark]);
    let mut image = view::field(&e.stage, field, PIXELS, view::Style::Phase, peak);
    let heads: Vec<Vec<Vec2>> = engine.positions().iter().take(DOTS).map(|p| vec![*p]).collect();
    view::particles(&mut image, e.stage.extent(), &heads);
    png::write(path, &image).expect("write frame");
}

fn write_trails(e: &Experiment, o: &Outcome, trails: &[Vec<Vec2>], out: &str) {
    let glow: Vec<f32> = o.glow.iter().map(|g| g.sqrt()).collect();
    let peak = o.glow.iter().cloned().fold(f32::MIN_POSITIVE, f32::max);
    let field = [&glow[..], &vec![0.0; glow.len()]];
    let mut image = view::field(&e.stage, field, PIXELS, view::Style::Density, peak);
    view::particles(&mut image, e.stage.extent(), trails);
    png::write(format!("{out}/{}_trails.png", e.name), &image).expect("write trails");
}

/// Steps of dt between the points of a trail.
pub fn trail_stride(dt: f64) -> usize {
    (TRAIL_STEP / dt).round().max(1.0) as usize
}

/// The first trails as rows of trail, time and position.
fn write_trail_table(e: &Experiment, dt: f64, trails: &[Vec<Vec2>], out: &str) {
    let dt = dt * trail_stride(dt) as f64;
    let row = |(k, trail): (usize, &Vec<Vec2>)| -> String {
        let points = trail.iter().enumerate().step_by(TRAIL_EVERY);
        points.map(|(s, p)| format!("{k},{},{},{}\n", s as f64 * dt, p[0], p[1])).collect()
    };
    let rows: String = trails.iter().take(TRAIL_ROWS).enumerate().map(row).collect();
    fs::write(format!("{out}/{}_trails.csv", e.name), format!("trail,t,x,y\n{rows}")).expect("write trails");
}

fn write_screen(e: &Experiment, o: &Outcome, out: &str) {
    let Some(screen) = &o.screen else { return };
    let runs = o.start.len() as f64;
    let hits: Vec<f64> = screen.hits.iter().map(|h| h / runs).collect();
    let theory: Vec<f64> = (0..BINS).map(|b| e.theory(screen.x(b))).collect();
    let share = hits.iter().sum::<f64>() / theory.iter().sum::<f64>().max(1e-300);
    let theory: Vec<f64> = theory.iter().map(|t| t * share).collect();
    let series: [(&[f64], _); 2] = [(&hits, [255, 255, 255]), (&theory, [255, 160, 30])];
    png::write(format!("{out}/{}_screen.png", e.name), &plot::curves(PLOT.0, PLOT.1, &series))
        .expect("write screen plot");
    let rows: String = (0..BINS)
        .map(|b| format!("{},{},{}\n", screen.x(b), hits[b], theory[b]))
        .collect();
    fs::write(format!("{out}/{}_screen.csv", e.name), format!("x,hits,theory\n{rows}"))
        .expect("write screen table");
    let arrivals: String = screen.arrivals.iter().map(|(t, x)| format!("{t},{x}\n")).collect();
    fs::write(format!("{out}/{}_arrivals.csv", e.name), format!("t,x\n{arrivals}"))
        .expect("write arrivals");
}

/// Where each particle started and ended.
fn write_runs(e: &Experiment, o: &Outcome, out: &str) {
    let row = |(s, p): (&Vec2, &Vec2)| format!("{},{},{},{}\n", s[0], s[1], p[0], p[1]);
    let rows: String = o.start.iter().zip(&o.end).map(row).collect();
    fs::write(format!("{out}/{}_runs.csv", e.name), format!("x_start,y_start,x,y\n{rows}")).expect("write runs");
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
    write_trail_table(e, o.dt, trails, out);
    write_screen(e, o, out);
    write_runs(e, o, out);
}
