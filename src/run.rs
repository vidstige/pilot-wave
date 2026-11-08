//! Running an experiment and measuring its outcome. No I/O: frames
//! go to a callback.
use crate::experiment::{Experiment, Measure};
use crate::particles::Particles;
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
    /// Which particles ever crossed the axis x = 0.
    pub crossed: Vec<bool>,
}

/// Run to the experiment's duration, calling `frame` with the field,
/// the particles and the step count at the start and after every
/// step.
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
    let mut particles = Particles::new(start.clone());
    let mut screen = e.screen.map(|y| Screen::new(y, e.extent() / 2.0, BINS, start.len()));
    let mut glow = vec![0.0; wave.n * wave.n];
    let mut crossed = vec![false; start.len()];
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
        for ((c, s), p) in crossed.iter_mut().zip(&start).zip(&particles.positions) {
            *c |= (p[0] < 0.0) != (s[0] < 0.0);
        }
        if step % GLOW_EVERY == 0 {
            for (g, k) in glow.iter_mut().zip(0..) {
                *g += wave.density(k);
            }
        }
        frame(&wave, &particles, step);
    }
    Outcome { wave, particles, start, screen, glow, crossed }
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

fn fraction(flags: &[bool]) -> f64 {
    flags.iter().filter(|f| **f).count() as f64 / flags.len() as f64
}

/// The share of the arrived particles that started among the leading
/// starts along y, as many of those as there are arrived particles.
fn arrived_from_front(start: &[Vec2], arrived: &[bool]) -> f64 {
    let count = arrived.iter().filter(|a| **a).count();
    let mut order: Vec<usize> = (0..start.len()).collect();
    order.sort_by(|&a, &b| start[b][1].total_cmp(&start[a][1]));
    order[..count].iter().filter(|&&k| arrived[k]).count() as f64 / count as f64
}

/// Standard deviation of the field's density along an axis.
fn field_spread(wave: &Wave, axis: usize) -> f64 {
    let n = wave.n;
    let moment = |power: i32| {
        let at = |k: usize| wave.density(k) as f64 * wave.position((k % n, k / n))[axis].powi(power);
        (0..n * n).map(at).sum::<f64>()
    };
    let (total, mean) = (moment(0), moment(1) / moment(0));
    (moment(2) / total - mean * mean).sqrt()
}

fn on_screen(e: &Experiment, o: &Outcome, measure: &Measure) -> f64 {
    let screen = o.screen.as_ref().expect("this measure needs a screen");
    let runs = o.particles.positions.len() as f64;
    match measure {
        Measure::Arrived => fraction(&screen.arrived),
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
        Measure::ArrivedFromFront => arrived_from_front(&o.start, &screen.arrived),
        _ => unreachable!(),
    }
}

pub fn measure(e: &Experiment, o: &Outcome, measure: &Measure) -> f64 {
    let runs = &o.particles.positions;
    let n = runs.len() as f64;
    let mean = |axis: usize| runs.iter().map(|p| p[axis]).sum::<f64>() / n;
    match measure {
        Measure::CrossedAxis => fraction(&o.crossed),
        Measure::Mean { axis } => mean(*axis),
        Measure::Spread { axis } => {
            let m = mean(*axis);
            (runs.iter().map(|p| (p[*axis] - m).powi(2)).sum::<f64>() / n).sqrt()
        }
        Measure::Moved => {
            let mut d: Vec<f64> =
                o.start.iter().zip(runs).map(|(s, p)| (p[0] - s[0]).hypot(p[1] - s[1])).collect();
            d.sort_by(f64::total_cmp);
            d[d.len() / 2]
        }
        Measure::FieldSpread { axis } => field_spread(&o.wave, *axis),
        Measure::Norm => o.wave.norm(),
        _ => on_screen(e, o, measure),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arrivals_from_the_front_count_the_leading_starts() {
        let start = [[0.0, 1.0], [0.0, 3.0], [0.0, 2.0], [0.0, 0.0]];
        assert_eq!(arrived_from_front(&start, &[false, true, true, false]), 1.0);
        assert_eq!(arrived_from_front(&start, &[true, true, false, false]), 0.5);
    }
}
