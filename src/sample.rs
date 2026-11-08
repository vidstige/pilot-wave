//! Draw positions from |ψ|² on the grid: a cell by its probability,
//! then anywhere within it. No I/O.
use crate::rng::Rng;
use crate::wave::{Vec2, Wave};

pub struct Sampler {
    cumulative: Vec<f64>,
    n: usize,
    dx: f64,
    half: f64,
}

impl Sampler {
    pub fn new(wave: &Wave) -> Sampler {
        let mut total = 0.0;
        let cumulative = (0..wave.n * wave.n)
            .map(|k| {
                total += wave.density(k) as f64;
                total
            })
            .collect();
        Sampler { cumulative, n: wave.n, dx: wave.dx, half: wave.extent() / 2.0 }
    }

    pub fn draw(&self, rng: &mut Rng) -> Vec2 {
        let total = *self.cumulative.last().unwrap();
        let u = rng.uniform() * total;
        let k = self.cumulative.partition_point(|&c| c < u).min(self.n * self.n - 1);
        let (i, j) = (k % self.n, k / self.n);
        [
            (i as f64 + rng.uniform()) * self.dx - self.half,
            (j as f64 + rng.uniform()) * self.dx - self.half,
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::packet::Packet;

    #[test]
    fn draws_spread_like_the_density() {
        let mut wave = Wave::new(128, 0.1);
        let packet = Packet { center: [1.0, -2.0], width: [0.5, 1.5], wavevector: [0.0; 2] };
        packet.write(&mut wave);
        let sampler = Sampler::new(&wave);
        let mut rng = Rng::new(3);
        let n = 4000;
        let draws: Vec<Vec2> = (0..n).map(|_| sampler.draw(&mut rng)).collect();
        for axis in 0..2 {
            let mean = draws.iter().map(|p| p[axis]).sum::<f64>() / n as f64;
            let var = draws.iter().map(|p| (p[axis] - mean).powi(2)).sum::<f64>() / n as f64;
            assert!((mean - packet.center[axis]).abs() < 0.05 * packet.width[axis]);
            assert!((var.sqrt() - packet.width[axis]).abs() < 0.05 * packet.width[axis]);
        }
    }
}
