//! An unbounded bath: the sum of every droplet's wake.
use crate::bath::{Bath, Vec2};
use crate::wake::Wake;

pub struct Open {
    wakes: Vec<Wake>,
}

impl Open {
    pub fn new(memory: f64, droplets: usize) -> Open {
        Open { wakes: (0..droplets).map(|_| Wake::new(memory)).collect() }
    }
}

impl Bath for Open {
    fn height(&self, p: Vec2) -> f64 {
        self.wakes.iter().map(|w| w.height(p)).sum()
    }

    fn slope(&self, p: Vec2) -> Vec2 {
        self.wakes.iter().map(|w| w.slope(p)).fold([0.0; 2], |g, s| [g[0] + s[0], g[1] + s[1]])
    }

    fn remember(&mut self, impacts: &[Vec2]) {
        for (wake, &p) in self.wakes.iter_mut().zip(impacts) {
            wake.remember(p);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bath_is_the_sum_of_its_wakes() {
        let mut open = Open::new(10.0, 2);
        let mut alone = [Wake::new(10.0), Wake::new(10.0)];
        for i in 0..20 {
            let a = [i as f64 * 0.1, 0.0];
            let b = [0.0, -(i as f64) * 0.1];
            open.remember(&[a, b]);
            alone[0].remember(a);
            alone[1].remember(b);
        }
        let p = [0.3, 0.4];
        let sum = alone[0].height(p) + alone[1].height(p);
        assert!((open.height(p) - sum).abs() < 1e-9);
        let s = open.slope(p);
        let (a, b) = (alone[0].slope(p), alone[1].slope(p));
        assert!((s[0] - a[0] - b[0]).abs() < 1e-9 && (s[1] - a[1] - b[1]).abs() < 1e-9);
    }
}
