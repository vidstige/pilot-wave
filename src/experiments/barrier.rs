//! Tunnelling through a rectangular barrier: the textbook transmission,
//! and which particles get through.
use crate::experiment::{Expectation, Experiment, Measure, ENERGY, SPEED};
use crate::experiments::free::ELECTRON;
use crate::packet::Packet;
use crate::wave::{Vec2, Wave};

/// The plane is |x|, |y| ≤ 25.6.
const CELLS: usize = 512;
const DX: f64 = 0.1;
/// The barrier across the plane from y = 0, four cells wide and 1.2
/// times the packet's energy high.
const WIDTH: f64 = 0.4;
const HEIGHT: f64 = 1.2 * ENERGY;
/// The packet: wide across so the barrier is effectively one
/// dimensional, 2 long so its wavenumber spreads 1/(2·2) = 0.25.
const START: f64 = -12.0;
const WIDTHS: [f64; 2] = [6.0, 2.0];
const SCREEN: f64 = 8.0;
const TAIL: f64 = 8.0;
/// Samples on each side of the centre when averaging over the spread.
const QUADRATURE: i32 = 4000;

fn barrier(p: Vec2) -> f64 {
    if (0.0..WIDTH).contains(&p[1]) { HEIGHT } else { 0.0 }
}

fn packet(wave: &mut Wave) {
    Packet { center: [0.0, START], width: WIDTHS, wavevector: [0.0, SPEED] }.write(wave)
}

/// Transmission through a rectangular barrier for a plane wave of
/// wavenumber k: the textbook formula, with sinh below the top of the
/// barrier and sin above it.
fn transmission(k: f64, height: f64, width: f64) -> f64 {
    let e = k * k / 2.0;
    let q = (2.0 * (height - e)).abs().sqrt() * width;
    let s = if e < height { q.sinh() } else { q.sin() };
    1.0 / (1.0 + height * height * s * s / (4.0 * e * (height - e).abs()))
}

/// The transmission averaged over a Gaussian spread of wavenumbers
/// with standard deviation `spread`, 1/(2σ) for a packet whose density
/// has standard deviation σ: the faster part tunnels more easily.
fn mean_transmission(k: f64, spread: f64, height: f64, width: f64) -> f64 {
    let (sum, total) = (-QUADRATURE..=QUADRATURE)
        .map(|i| 6.0 * i as f64 / QUADRATURE as f64)
        .map(|u| ((-u * u / 2.0).exp(), transmission(k + spread * u, height, width)))
        .fold((0.0, 0.0), |(s, t), (w, tr)| (s + w * tr, t + w));
    sum / total
}

pub fn tunnelling() -> Experiment {
    let expected = mean_transmission(SPEED, 1.0 / (2.0 * WIDTHS[1]), HEIGHT, WIDTH);
    Experiment {
        name: "barrier",
        summary: "\
A packet 6 wide across and 2 long heads up at speed 2π at a barrier \
0.4 wide and 1.2 times its energy high, across the whole plane, so \
the problem is the one-dimensional one of the textbooks: the \
transmission of a plane wave is T = 1/(1 + V² sinh²(κa)/4E(V − E)) \
with κ = √(2(V − E)), 0.227 at k = 2π, and averaged over the packet's \
spread of wavenumbers 0.230. A screen beyond the barrier counts what \
tunnels through; on the grid it comes out a few percent less, the \
grid's dispersion acting like a heavier mass. Which particles tunnel is \
decided by where they start: trajectories in one dimension cannot \
cross, so the ones that get through are exactly the ones that \
started at the front of the packet, the leading share equal to the \
transmission.",
        realisation: ELECTRON,
        cells: CELLS,
        dx: DX,
        walls: vec![],
        potential: barrier,
        source: packet,
        spin: [0.0; 3],
        particles: 20000,
        duration: (SCREEN - START + TAIL) / SPEED,
        screen: Some(SCREEN),
        theory: None,
        expectations: vec![
            Expectation {
                what: "fraction of the particles through the barrier, the averaged textbook transmission",
                measure: Measure::Arrived,
                expected,
                tolerance: 0.02,
            },
            Expectation {
                what: "fraction of the transmitted particles that started in the front of the packet",
                measure: Measure::ArrivedFromFront,
                expected: 1.0,
                tolerance: 0.01,
            },
        ],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn transmission_is_whole_without_a_barrier() {
        assert!((transmission(SPEED, 0.0, WIDTH) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn transmission_is_continuous_over_the_barrier_top() {
        let top = (2.0 * HEIGHT).sqrt();
        let limit = 1.0 / (1.0 + HEIGHT * WIDTH * WIDTH / 2.0);
        for k in [top - 1e-6, top + 1e-6] {
            assert!((transmission(k, HEIGHT, WIDTH) - limit).abs() < 1e-4);
        }
    }

    #[test]
    fn a_narrow_spread_averages_to_the_plane_wave() {
        let plane = transmission(SPEED, HEIGHT, WIDTH);
        assert!((mean_transmission(SPEED, 1e-4, HEIGHT, WIDTH) - plane).abs() < 1e-6);
        assert!(mean_transmission(SPEED, 0.25, HEIGHT, WIDTH) > plane);
    }
}
