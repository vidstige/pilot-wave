//! Draw positions from a density over the stage: a cell by its share,
//! then anywhere within it. No I/O.
use crate::experiment::Stage;
use crate::plane::Vec2;
use crate::rng::Rng;

pub struct Sampler {
    cumulative: Vec<f64>,
    n: usize,
    dx: f64,
    half: f64,
}

impl Sampler {
    pub fn new(stage: &Stage, density: impl Fn(Vec2) -> f64) -> Sampler {
        let mut total = 0.0;
        let cumulative = (0..stage.cells * stage.cells)
            .map(|k| {
                total += density(stage.position(k));
                total
            })
            .collect();
        Sampler { cumulative, n: stage.cells, dx: stage.dx, half: stage.extent() / 2.0 }
    }

    pub fn draw(&self, rng: &mut Rng) -> Vec2 {
        let total = *self.cumulative.last().unwrap();
        let u = rng.uniform() * total;
        let k = self.cumulative.partition_point(|&c| c < u).min(self.n * self.n - 1);
        [k % self.n, k / self.n].map(|i| (i as f64 + rng.uniform()) * self.dx - self.half)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::experiment::free;
    use crate::packet::Packet;

    #[test]
    fn draws_spread_like_the_density() {
        let stage = Stage { cells: 128, dx: 0.1, scale: 1.0, walls: vec![], potential: free };
        let packet = Packet { center: [1.0, -2.0], width: [0.5, 1.5], wavevector: [0.0; 2] };
        let sampler = Sampler::new(&stage, |p| packet.amplitude(p)[0].powi(2));
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
