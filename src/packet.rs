//! Gaussian wave packets: the field's starting shape and the matching
//! draw of particle positions from |ψ|². No I/O.
use std::f64::consts::PI;

use crate::grid::Vec2;
use crate::rng::Rng;
use crate::wave::Wave;

#[derive(Clone, Copy, Debug)]
pub struct Packet {
    pub center: Vec2,
    /// Standard deviations of |ψ|² along x and y.
    pub width: Vec2,
    pub wavevector: Vec2,
}

impl Packet {
    /// Write ψ = exp(-(x-c)²/4σ²) exp(i k·x) over the grid, normalised.
    pub fn write(&self, wave: &mut Wave) {
        let amplitude = 1.0 / (2.0 * PI * self.width[0] * self.width[1]).sqrt();
        for j in 0..wave.n {
            for i in 0..wave.n {
                let p = wave.position((i, j));
                let d = [p[0] - self.center[0], p[1] - self.center[1]];
                let envelope = (-d[0] * d[0] / (4.0 * self.width[0] * self.width[0])
                    - d[1] * d[1] / (4.0 * self.width[1] * self.width[1]))
                    .exp();
                let phase = self.wavevector[0] * p[0] + self.wavevector[1] * p[1];
                let k = wave.index((i, j));
                wave.re[k] = (amplitude * envelope * phase.cos()) as f32;
                wave.im[k] = (amplitude * envelope * phase.sin()) as f32;
            }
        }
    }

    /// A position drawn from |ψ|², which is Gaussian with these widths.
    pub fn sample(&self, rng: &mut Rng) -> Vec2 {
        [
            self.center[0] + self.width[0] * rng.normal(),
            self.center[1] + self.width[1] * rng.normal(),
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn samples_spread_like_the_packet() {
        let packet = Packet { center: [1.0, -2.0], width: [0.5, 2.0], wavevector: [0.0; 2] };
        let mut rng = Rng::new(3);
        let n = 4000;
        let samples: Vec<Vec2> = (0..n).map(|_| packet.sample(&mut rng)).collect();
        for axis in 0..2 {
            let mean = samples.iter().map(|p| p[axis]).sum::<f64>() / n as f64;
            let var = samples.iter().map(|p| (p[axis] - mean).powi(2)).sum::<f64>() / n as f64;
            assert!((mean - packet.center[axis]).abs() < 0.05 * packet.width[axis]);
            assert!((var.sqrt() - packet.width[axis]).abs() < 0.05 * packet.width[axis]);
        }
    }
}
