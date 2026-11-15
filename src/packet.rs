//! A Gaussian wave packet. No I/O.
use crate::plane::Vec2;

/// A wave's two quadratures at a point, as a complex amplitude.
pub type Complex = [f64; 2];

#[derive(Clone, Copy, Debug)]
pub struct Packet {
    pub center: Vec2,
    /// Standard deviations of the squared amplitude along x and y.
    pub width: Vec2,
    pub wavevector: Vec2,
}

impl Packet {
    /// exp(-(x-c)²/4σ²) exp(i k·x) at p.
    pub fn amplitude(&self, p: Vec2) -> Complex {
        let along = |a: usize| ((p[a] - self.center[a]) / (2.0 * self.width[a])).powi(2);
        let envelope = (-along(0) - along(1)).exp();
        let phase = self.wavevector[0] * p[0] + self.wavevector[1] * p[1];
        [envelope * phase.cos(), envelope * phase.sin()]
    }
}
