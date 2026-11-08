//! Running an experiment and measuring its outcome. No I/O: frames
//! go to a callback.
use crate::encounter::{self, Crossings};
use crate::experiment::{Experiment, Measure, TOGETHER};
use crate::pair;
use crate::particles::Particles;
use crate::rng::Rng;
use crate::sample::Sampler;
use crate::screen::Screen;
use crate::wave::{Vec2, Wave};

pub const BINS: usize = 128;
/// Steps between additions to the time-integrated density and
/// readings of the bunching.
const GLOW_EVERY: usize = 10;

pub struct Outcome {
    pub wave: Wave,
    pub particles: Particles,
    pub start: Vec<Vec2>,
    pub screen: Option<Screen>,
    /// Time-integrated density, for drawing trails over.
    pub glow: Vec<f32>,
    /// Which runs crossed the diagonal, and the most of the runs and
    /// of the probability ever within TOGETHER, when measured.
    pub crossings: Crossings,
    pub together_peak: f64,
    pub probability_together_peak: f64,
    /// The same runs of the experiment without its potential.
    pub twin: Option<Box<Outcome>>,
}

fn together(p: Vec2) -> bool {
    encounter::separation(p) < TOGETHER
}

fn run_alone(
    e: &Experiment,
    particles: usize,
    seed: u64,
    mut frame: impl FnMut(&Wave, &Particles, usize),
) -> Outcome {
    let mut wave = e.build();
    let mut rng = Rng::new(seed);
    let sampler = Sampler::new(&wave);
    let start: Vec<Vec2> = (0..particles).map(|_| sampler.draw(&mut rng)).collect();
    let mut particles = Particles::new(start.clone());
    let mut screen = e.screen.map(|y| Screen::new(y, e.extent() / 2.0, BINS, start.len()));
    let mut glow = vec![0.0; wave.n * wave.n];
    let mut crossings = Crossings::new(&start);
    let bunching = e.measures(|m| matches!(m, Measure::TogetherPeak | Measure::ProbabilityTogetherPeak));
    let (mut together_peak, mut probability_together_peak) = (0.0f64, 0.0f64);
    let dt = e.dt();
    let steps = (e.duration / dt).ceil() as usize;
    frame(&wave, &particles, 0);
    for step in 1..=steps {
        wave.step(dt);
        let before = particles.positions.clone();
        particles.step(&wave, dt);
        if let Some(screen) = &mut screen {
            screen.record(wave.time, &before, &particles.positions);
            screen.integrate(&wave, dt);
        }
        crossings.record(&particles.positions);
        if step % GLOW_EVERY == 0 {
            for (g, k) in glow.iter_mut().zip(0..) {
                *g += wave.density(k);
            }
            if bunching {
                let runs = encounter::fraction_within(&particles.positions, TOGETHER);
                together_peak = together_peak.max(runs);
                probability_together_peak = probability_together_peak.max(pair::probability(&wave, together));
            }
        }
        frame(&wave, &particles, step);
    }
    Outcome {
        wave,
        particles,
        start,
        screen,
        glow,
        crossings,
        together_peak,
        probability_together_peak,
        twin: None,
    }
}

/// Run to the experiment's duration, calling `frame` with the field,
/// the particles and the step count at the start and after every
/// step; then the twin on the same starting points, if a measure
/// needs it.
pub fn run(
    e: &Experiment,
    particles: usize,
    seed: u64,
    frame: impl FnMut(&Wave, &Particles, usize),
) -> Outcome {
    let mut outcome = run_alone(e, particles, seed, frame);
    let twin = |m: &Measure| matches!(m, Measure::MovedVsTwin { .. } | Measure::MarginalVsTwin { .. });
    if e.measures(twin) {
        outcome.twin = Some(Box::new(run_alone(&e.twin(), particles, seed, |_, _, _| {})));
    }
    outcome
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

/// Largest gap between two curves as a fraction of the first's peak.
fn largest_gap(a: &[f64], b: &[f64]) -> f64 {
    let peak = a.iter().cloned().fold(0.0, f64::max);
    a.iter().zip(b).map(|(p, q)| (p - q).abs()).fold(0.0, f64::max) / peak
}

fn on_screen(e: &Experiment, o: &Outcome, measure: &Measure) -> f64 {
    let screen = o.screen.as_ref().expect("this measure needs a screen");
    let runs = o.particles.positions.len() as f64;
    match measure {
        Measure::Arrived => screen.hits.iter().sum::<f64>() / runs,
        Measure::Through => screen.flux.iter().sum(),
        Measure::ArrivedAfter { time } => {
            screen.arrivals.iter().filter(|(t, _)| t > time).count() as f64 / runs
        }
        Measure::MeanArrivalTime => {
            screen.arrivals.iter().map(|(t, _)| t).sum::<f64>() / screen.arrivals.len() as f64
        }
        Measure::CurrentVsTheory => {
            let theory: Vec<f64> = (0..BINS).map(|b| e.theory(screen.x(b))).collect();
            rms_gap(&screen.flux, &scaled_to(&theory, &screen.flux)) / screen.peak()
        }
        Measure::ArrivalsVsCurrent => {
            rms_gap(&screen.flux, &scaled_to(&screen.hits, &screen.flux)) / screen.peak()
        }
        Measure::CurrentAt { x } => screen.flux[screen.bin_of(*x)] / screen.peak(),
        Measure::MaximumNear { x } => {
            let f = &screen.flux;
            let peaks = (1..BINS - 1).filter(|&b| f[b] > f[b - 1] && f[b] >= f[b + 1]);
            let nearest = peaks.min_by(|&a, &b| {
                (screen.x(a) - x).abs().total_cmp(&(screen.x(b) - x).abs())
            });
            screen.x(nearest.expect("a current with a maximum"))
        }
        _ => unreachable!(),
    }
}

pub fn measure(e: &Experiment, o: &Outcome, measure: &Measure) -> f64 {
    let runs = &o.particles.positions;
    let twin = || o.twin.as_ref().expect("this measure needs the twin");
    match measure {
        Measure::Crossed => o.crossings.fraction(),
        Measure::TogetherPeak => o.together_peak,
        Measure::ProbabilityTogetherPeak => o.probability_together_peak,
        Measure::Apart { beyond } => 1.0 - encounter::fraction_within(runs, *beyond),
        Measure::ProbabilityApart { beyond } => {
            pair::probability(&o.wave, |p| encounter::separation(p) > *beyond) + 1.0 - o.wave.norm()
        }
        Measure::Correlation => encounter::correlation(runs),
        Measure::FieldCorrelation => pair::correlation(&o.wave),
        Measure::MovedVsTwin { axis } => encounter::moved(runs, &twin().particles.positions, *axis),
        Measure::MarginalVsTwin { axis } => {
            largest_gap(&pair::marginals(&o.wave)[*axis], &pair::marginals(&twin().wave)[*axis])
        }
        _ => on_screen(e, o, measure),
    }
}
