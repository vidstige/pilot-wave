//! Running an experiment on an engine. No I/O: frames go to a
//! callback.
use crate::engine::Engine;
use crate::experiment::Experiment;
use crate::plane::Vec2;
use crate::rng::Rng;
use crate::sample::Sampler;
use crate::screen::Screen;

pub const BINS: usize = 128;
/// Steps between additions to the time-integrated field.
const GLOW_EVERY: usize = 10;

pub struct Outcome {
    pub start: Vec<Vec2>,
    pub end: Vec<Vec2>,
    pub screen: Option<Screen>,
    /// The engine's field squared and summed over time, for drawing
    /// trails over.
    pub glow: Vec<f32>,
    /// The step the engine took.
    pub dt: f64,
}

/// Where the runs start: where the wave at the start is strong, as
/// its square, or for every other run evenly in the experiment's
/// start rectangle.
pub fn starts(e: &Experiment, runs: usize, seed: u64) -> Vec<Vec2> {
    let sampler = Sampler::new(&e.stage, |p| e.amplitude(p).iter().map(|a| a * a).sum());
    let mut rng = Rng::new(seed);
    let place = |k: usize| match e.start {
        Some(rect) if k % 2 == 0 => rect.map(|[from, to]| from + rng.uniform() * (to - from)),
        _ => sampler.draw(&mut rng),
    };
    (0..runs).map(place).collect()
}

/// Run to the experiment's duration, calling `frame` with the engine
/// and the step count at the start and after every step.
pub fn run(
    e: &Experiment,
    engine: &mut dyn Engine,
    mut frame: impl FnMut(&dyn Engine, usize),
) -> Outcome {
    let start = engine.positions().to_vec();
    let mut screen = e.screen.map(|y| Screen::new(y, e.stage.extent() / 2.0, BINS, start.len()));
    let mut glow = vec![0.0; e.stage.cells * e.stage.cells];
    let dt = engine.dt();
    frame(engine, 0);
    for step in 1..=(e.duration / dt).ceil() as usize {
        let before = engine.positions().to_vec();
        engine.step();
        if let Some(screen) = &mut screen {
            screen.record(step as f64 * dt, &before, engine.positions());
        }
        if let (0, Some([re, im])) = (step % GLOW_EVERY, engine.field()) {
            for (g, (re, im)) in glow.iter_mut().zip(re.iter().zip(im)) {
                *g += re * re + im * im;
            }
        }
        frame(engine, step);
    }
    Outcome { start, end: engine.positions().to_vec(), screen, glow, dt }
}
