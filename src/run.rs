//! Running an experiment and measuring its outcome. No I/O: frames
//! go to a callback.
use crate::experiment::{Experiment, Measure};
use crate::particles::{Kind, Particles};
use crate::rng::Rng;
use crate::sample::Sampler;
use crate::screen::Screen;
use crate::wave::{Vec2, Wave};

pub const BINS: usize = 128;
/// Steps between additions to the time-integrated density.
const GLOW_EVERY: usize = 10;

pub struct Outcome {
    pub wave: Wave,
    pub particles: Particles,
    pub start: Vec<Vec2>,
    pub screen: Option<Screen>,
    /// Time-integrated density, for drawing trails over.
    pub glow: Vec<f32>,
}

/// Run to the experiment's duration, calling `frame` with the field,
/// the particles and the step count at every step.
pub fn run(
    e: &Experiment,
    particles: usize,
    seed: u64,
    mut frame: impl FnMut(&Wave, &Particles, usize),
) -> Outcome {
    let mut wave = e.build();
    let mut rng = Rng::new(seed);
    let sampler = Sampler::new(&wave);
    let start: Vec<Vec2> = (0..particles).map(|_| sampler.draw(&mut rng)).collect();
    let mut particles = Particles::new(Kind::Quantum, start.clone(), &wave);
    let mut screen = e.screen.map(|y| Screen::new(y, e.extent() / 2.0, BINS));
    let mut glow = vec![0.0; wave.n * wave.n];
    let dt = e.dt();
    let steps = (e.duration / dt).ceil() as usize;
    for step in 1..=steps {
        wave.step(dt);
        let before = particles.positions.clone();
        particles.step(&wave, dt);
        if let Some(screen) = &mut screen {
            screen.record(wave.time, &before, &particles.positions);
            screen.integrate(&wave, dt);
        }
        if step % GLOW_EVERY == 0 {
            for (g, k) in glow.iter_mut().zip(0..) {
                *g += wave.density(k);
            }
        }
        frame(&wave, &particles, step);
    }
    Outcome { wave, particles, start, screen, glow }
}

/// Scale a curve to the total of another.
fn scaled_to(curve: &[f64], reference: &[f64]) -> Vec<f64> {
    let total: f64 = reference.iter().sum();
    let sum: f64 = curve.iter().sum::<f64>().max(1e-300);
    curve.iter().map(|v| v * total / sum).collect()
}

fn rms_gap(a: &[f64], b: &[f64]) -> f64 {
    (a.iter().zip(b).map(|(x, y)| (x - y).powi(2)).sum::<f64>() / a.len() as f64).sqrt()
}

pub fn measure(e: &Experiment, o: &Outcome, measure: &Measure) -> f64 {
    let screen = o.screen.as_ref().expect("the measures so far need a screen");
    let arrived = screen.hits.iter().sum::<f64>() / o.particles.positions.len() as f64;
    let through: f64 = screen.flux.iter().sum();
    match measure {
        Measure::Arrived => arrived,
        Measure::ArrivedMinusThrough => arrived - through,
        Measure::CurrentVsTheory => {
            let theory: Vec<f64> = (0..BINS).map(|b| e.theory(screen.x(b))).collect();
            rms_gap(&screen.flux, &scaled_to(&theory, &screen.flux)) / screen.peak()
        }
        Measure::ArrivalsVsCurrent => {
            rms_gap(&screen.flux, &scaled_to(&screen.hits, &screen.flux)) / screen.peak()
        }
        Measure::CurrentAt { x } => {
            let bin = ((x + screen.half) / (2.0 * screen.half) * BINS as f64) as usize;
            screen.flux[bin] / screen.peak()
        }
    }
}
