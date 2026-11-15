//! What a run's outcome says to each of an experiment's questions.
//! No I/O.
use crate::experiment::{Experiment, Measure};
use crate::run::{Outcome, BINS};

/// Scale a curve to the total of another.
fn scaled_to(curve: &[f64], reference: &[f64]) -> Vec<f64> {
    let total: f64 = reference.iter().sum();
    let sum: f64 = curve.iter().sum::<f64>().max(1e-300);
    curve.iter().map(|v| v * total / sum).collect()
}

fn rms_gap(a: &[f64], b: &[f64]) -> f64 {
    (a.iter().zip(b).map(|(x, y)| (x - y).powi(2)).sum::<f64>() / a.len() as f64).sqrt()
}

fn on_screen(e: &Experiment, o: &Outcome, measure: &Measure) -> f64 {
    let screen = o.screen.as_ref().expect("this measure needs a screen");
    let runs = o.end.len() as f64;
    let hits = &screen.hits;
    let peak = hits.iter().cloned().fold(0.0, f64::max);
    match measure {
        Measure::Arrived => screen.arrivals.len() as f64 / runs,
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
        _ => unreachable!(),
    }
}

pub fn measure(e: &Experiment, o: &Outcome, measure: &Measure) -> f64 {
    let runs = &o.end;
    let n = runs.len() as f64;
    let mean = |axis: usize| runs.iter().map(|p| p[axis]).sum::<f64>() / n;
    match measure {
        Measure::Mean { axis } => mean(*axis),
        Measure::Spread { axis } => {
            let m = mean(*axis);
            (runs.iter().map(|p| (p[*axis] - m).powi(2)).sum::<f64>() / n).sqrt()
        }
        _ => on_screen(e, o, measure),
    }
}
