//! A detection line: where particles cross it, and the probability
//! current the field carries through it, binned along x. No I/O.
use crate::grid::Vec2;
use crate::wave::Wave;

pub struct Screen {
    pub y: f64,
    pub half: f64,
    pub bins: usize,
    pub hits: Vec<f64>,
    /// Time-integrated probability current through the line.
    pub flux: Vec<f64>,
}

impl Screen {
    pub fn new(y: f64, half: f64, bins: usize) -> Screen {
        Screen { y, half, bins, hits: vec![0.0; bins], flux: vec![0.0; bins] }
    }

    fn bin(&self, x: f64) -> Option<usize> {
        let t = (x + self.half) / (2.0 * self.half);
        (0.0..1.0).contains(&t).then(|| (t * self.bins as f64) as usize)
    }

    /// Count particles that crossed the line upward this step.
    pub fn record(&mut self, before: &[Vec2], after: &[Vec2]) {
        for (a, b) in before.iter().zip(after) {
            if a[1] < self.y && b[1] >= self.y {
                if let Some(bin) = self.bin(b[0]) {
                    self.hits[bin] += 1.0;
                }
            }
        }
    }

    /// Add this step's probability current through the line.
    pub fn integrate(&mut self, wave: &Wave, dt: f64) {
        let (_, row) = wave.cell([0.0, self.y]);
        let j = row.round() as usize;
        for i in 1..wave.n - 1 {
            let k = wave.index((i, j));
            let x = wave.position((i, j))[0];
            if let Some(bin) = self.bin(x) {
                let v = wave.velocity([x, self.y])[1];
                self.flux[bin] += wave.density(k) as f64 * v * dt * wave.dx;
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
    fn only_upward_crossings_count_and_land_in_their_bin() {
        let mut screen = Screen::new(1.0, 2.0, 4);
        screen.record(&[[0.5, 0.9], [-1.5, 1.2], [1.5, 0.0]], &[[0.5, 1.1], [-1.5, 0.8], [1.5, 1.5]]);
        assert_eq!(screen.hits, vec![0.0, 0.0, 1.0, 1.0]);
    }
}
