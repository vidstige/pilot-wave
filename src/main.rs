use std::env;

use pilot_wave::experiment::{experiment, fraunhofer, SINGLE_SLIT_WIDTH, SLIT_SEPARATION, SLIT_WIDTH};
use pilot_wave::grid::Vec2;
use pilot_wave::particles::{Kind, Particles};
use pilot_wave::rng::Rng;
use pilot_wave::screen::Screen;
use pilot_wave::{plot, png, view};

const PIXELS: usize = 800;
const GRID: usize = 1024;
/// Time between frames and the number of particles whose trails are
/// drawn.
const FRAME_EVERY: f64 = 0.5;
const TRAILS: usize = 300;
const BINS: usize = 120;
const PLOT: (usize, usize) = (800, 300);

fn main() {
    let args: Vec<String> = env::args().collect();
    let name: &'static str = pilot_wave::experiment::NAMES
        .into_iter()
        .find(|n| Some(*n) == args.get(1).map(String::as_str))
        .unwrap_or("double-slit");
    let count: usize = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(2000);
    let out = args.get(3).map(String::as_str).unwrap_or("out");
    let kind = if args.get(4).map(String::as_str) == Some("classical") { Kind::Classical } else { Kind::Quantum };
    std::fs::create_dir_all(out).expect("create output directory");

    let mut e = experiment(name, GRID, kind);
    let mut rng = Rng::new(1);
    let positions: Vec<Vec2> = (0..count).map(|_| e.packet.sample(&mut rng)).collect();
    let mut particles = Particles::new(kind, positions, e.packet.wavevector);
    let mut trails: Vec<Vec<Vec2>> = particles.positions[..count.min(TRAILS)].iter().map(|p| vec![*p]).collect();
    let half = e.wave.extent() / 2.0;
    let mut screen = e.screen.map(|y| Screen::new(y, half, BINS));
    let steps = (e.duration / e.dt).ceil() as usize;
    let per_frame = (FRAME_EVERY / e.dt) as usize;
    for step in 1..=steps {
        e.wave.step(e.dt);
        let before = particles.positions.clone();
        particles.step(&e.wave, e.dt);
        if let Some(screen) = &mut screen {
            screen.record(&before, &particles.positions);
            screen.integrate(&e.wave, e.dt);
        }
        for (trail, p) in trails.iter_mut().zip(&particles.positions) {
            trail.push(*p);
        }
        if step % per_frame == 0 || step == steps {
            let mut image = view::field(&e.wave, PIXELS, view::Style::Phase, e.packet.peak());
            let heads: Vec<Vec<Vec2>> = particles.positions.iter().map(|p| vec![*p]).collect();
            view::particles(&mut image, &e.wave, &heads);
            png::write(format!("{out}/{name}_{:04}.png", step / per_frame), &image).expect("write frame");
        }
    }
    let mut image = view::field(&e.wave, PIXELS, view::Style::Density, e.packet.peak());
    view::particles(&mut image, &e.wave, &trails);
    png::write(format!("{out}/{name}_trajectories.png"), &image).expect("write trajectories");
    if let Some(screen) = &screen {
        let theory: Vec<f64> = (0..BINS)
            .map(|b| match name {
                "double-slit" => fraunhofer(screen.x(b), screen.y, &[-SLIT_SEPARATION / 2.0, SLIT_SEPARATION / 2.0], SLIT_WIDTH),
                "single-slit" => fraunhofer(screen.x(b), screen.y, &[0.0], SINGLE_SLIT_WIDTH),
                _ => 0.0,
            })
            .collect();
        let series: [(&[f64], _); 3] = [
            (&screen.hits, [255, 255, 255]),
            (&screen.flux, [80, 200, 255]),
            (&theory, [255, 160, 30]),
        ];
        let image = plot::curves(PLOT.0, PLOT.1, &series);
        png::write(format!("{out}/{name}_screen.png"), &image).expect("write screen plot");
        let arrived: f64 = screen.hits.iter().sum();
        println!("{name}: {arrived} of {count} particles reached the screen");
    }
}
