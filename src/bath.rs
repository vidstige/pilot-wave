//! What a bath must offer the droplets: a wave height, its slope,
//! and a way to add this bounce's impacts.
pub type Vec2 = [f64; 2];

pub trait Bath {
    fn height(&self, p: Vec2) -> f64;
    fn slope(&self, p: Vec2) -> Vec2;
    /// Age every wave by one bounce and add this bounce's impacts.
    fn remember(&mut self, impacts: &[Vec2]);
}

pub fn distance(a: Vec2, b: Vec2) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sqrt()
}
