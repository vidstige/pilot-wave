//! Gaussian wave packets written onto the grid. No I/O.
use std::f64::consts::PI;

use crate::wave::{Vec2, Wave};

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


    /// The density at the packet's centre.
    pub fn peak(&self) -> f32 {
        (1.0 / (2.0 * PI * self.width[0] * self.width[1])) as f32
    }
}
