//! Positions drawn from |ψ|² on the grid, for fields with no closed
//! form to sample. No I/O.
use crate::grid::Vec2;
use crate::rng::Rng;
use crate::wave::Wave;

/// Positions drawn cell by cell from the density, each placed
/// uniformly within its cell.
pub fn positions(wave: &Wave, count: usize, rng: &mut Rng) -> Vec<Vec2> {
    let n = wave.n;
    let cumulative: Vec<f64> = (0..n * n)
        .scan(0.0, |sum, k| {
            *sum += wave.density(k) as f64;
            Some(*sum)
        })
        .collect();
    let total = cumulative[n * n - 1];
    let mut draw = || {
        let u = rng.uniform() * total;
        let k = cumulative.partition_point(|&c| c < u);
        let p = wave.position((k % n, k / n));
        let jitter = |rng: &mut Rng| (rng.uniform() - 0.5) * wave.dx;
        [p[0] + jitter(rng), p[1] + jitter(rng)]
    };
    (0..count).map(|_| draw()).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::packet::Packet;

    #[test]
    fn samples_spread_like_the_density() {
        let mut wave = Wave::new(128, 0.1);
        let packet = Packet { center: [1.0, -2.0], width: [0.5, 2.0], wavevector: [3.0, 0.0] };
        packet.write(&mut wave);
        let count = 8000;
        let samples = positions(&wave, count, &mut Rng::new(5));
        for axis in 0..2 {
            let values = samples.iter().map(|p| p[axis]);
            let mean = values.clone().sum::<f64>() / count as f64;
            let var = values.map(|v| (v - mean).powi(2)).sum::<f64>() / count as f64;
            let width = packet.width[axis];
            assert!((mean - packet.center[axis]).abs() < 0.05 * width);
            assert!((var.sqrt() - width).abs() < 0.05 * width);
        }
    }
}
