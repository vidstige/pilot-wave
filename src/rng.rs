//! Tiny xorshift generator with normal deviates; enough for noise.
use std::f64::consts::PI;

pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Rng {
        Rng(seed.wrapping_mul(0x9E3779B97F4A7C15) | 1)
    }

    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    /// Uniform in (0, 1].
    pub fn uniform(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64 + f64::EPSILON
    }

    /// Standard normal deviate by Box-Muller.
    pub fn normal(&mut self) -> f64 {
        let r = (-2.0 * self.uniform().ln()).sqrt();
        r * (2.0 * PI * self.uniform()).cos()
    }
}
