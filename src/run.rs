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
}

/// Run to the experiment's duration, calling `frame` with the field,
/// the particles and the step count at the start and after every
/// step.
pub fn run(
    e: &Experiment,
    runs: usize,
    seed: u64,
    mut frame: impl FnMut(&Wave, &Particles, usize),
) -> Outcome {
    let mut wave = e.build();
    let mut rng = Rng::new(seed);
    let sampler = Sampler::new(&wave);
    let start: Vec<Vec2> = (0..runs).map(|_| sampler.draw(&mut rng)).collect();
    let mut particles = Particles::new(start.clone());
    let mut screen = e.screen.map(|y| Screen::new(y, e.extent() / 2.0, BINS, start.len()));
    let mut glow = vec![0.0; wave.n * wave.n];
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
    let hits = &screen.hits;
    let peak = hits.iter().cloned().fold(0.0, f64::max);
    match measure {
        Measure::Arrived => fraction(&screen.arrived),
        Measure::ArrivedAfter { time } => {
            screen.arrivals.iter().filter(|(t, _)| t > time).count() as f64 / runs
        }
        Measure::MeanArrivalTime => {
            screen.arrivals.iter().map(|(t, _)| t).sum::<f64>() / screen.arrivals.len() as f64
        }
        Measure::ArrivalsVsTheory => {
            let theory: Vec<f64> = (0..BINS).map(|b| e.theory(screen.x(b))).collect();
            rms_gap(hits, &scaled_to(&theory, hits)) / peak
        }
        Measure::ArrivalsAt { x } => hits[screen.bin_of(*x)] / peak,
        Measure::ArrivalsPeakNear { x, within } => {
            let near = (0..BINS).filter(|&b| (screen.x(b) - x).abs() < *within);
            screen.x(near.max_by(|&a, &b| hits[a].total_cmp(&hits[b])).expect("bins near x"))
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
