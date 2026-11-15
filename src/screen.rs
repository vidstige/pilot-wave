//! A detection line: where and when particles first cross it, binned
//! along x. No I/O.
use crate::plane::Vec2;

pub struct Screen {
    pub y: f64,
    pub half: f64,
    pub bins: usize,
    pub hits: Vec<f64>,
    /// Each particle's first crossing as (time, x), in the order they
    /// happened.
    pub arrivals: Vec<(f64, f64)>,
    pub arrived: Vec<bool>,
}

impl Screen {
    pub fn new(y: f64, half: f64, bins: usize, particles: usize) -> Screen {
        Screen {
            y,
            half,
            bins,
            hits: vec![0.0; bins],
            arrivals: Vec::new(),
            arrived: vec![false; particles],
        }
    }

    fn bin(&self, x: f64) -> Option<usize> {
        let t = (x + self.half) / (2.0 * self.half);
        (0.0..1.0).contains(&t).then(|| (t * self.bins as f64) as usize)
    }

    /// The bin holding x, which must lie on the screen.
    pub fn bin_of(&self, x: f64) -> usize {
        self.bin(x).expect("x on the screen")
    }

    /// Count particles that crossed the line upward this step for
    /// the first time.
    pub fn record(&mut self, time: f64, before: &[Vec2], after: &[Vec2]) {
        for (k, (a, b)) in before.iter().zip(after).enumerate() {
            if a[1] < self.y && b[1] >= self.y && !self.arrived[k] {
                self.arrived[k] = true;
                self.arrivals.push((time, b[0]));
                if let Some(bin) = self.bin(b[0]) {
                    self.hits[bin] += 1.0;
                }
            }
        }
    }

    pub fn x(&self, bin: usize) -> f64 {
        -self.half + (bin as f64 + 0.5) * 2.0 * self.half / self.bins as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn only_first_upward_crossings_count_and_land_in_their_bin() {
        let mut screen = Screen::new(1.0, 2.0, 4, 3);
        screen.record(1.0, &[[0.5, 0.9], [-1.5, 1.2], [1.5, 0.0]], &[[0.5, 1.1], [-1.5, 0.8], [1.5, 1.5]]);
        screen.record(2.0, &[[0.5, 0.9], [-1.5, 0.8]], &[[0.5, 1.1], [-1.5, 1.1]]);
        assert_eq!(screen.arrivals, vec![(1.0, 0.5), (1.0, 1.5), (2.0, -1.5)]);
        assert_eq!(screen.hits, vec![1.0, 0.0, 1.0, 1.0]);
    }
}
