//! A scalar field sampled on a square grid of cells around a center.
use std::sync::Mutex;
use std::thread;

use crate::bath::{Bath, Vec2};

pub struct Grid {
    pub n: usize,
    pub center: Vec2,
    pub half: f64,
    pub values: Vec<f64>,
}

impl Grid {
    pub fn zeros(center: Vec2, half: f64, n: usize) -> Grid {
        Grid { n, center, half, values: vec![0.0; n * n] }
    }

    /// Bath heights at every cell, rows handed out to all cores.
    pub fn sample(bath: &(impl Bath + ?Sized), center: Vec2, half: f64, n: usize) -> Grid {
        let mut grid = Grid { n, center, half, values: Vec::new() };
        let mut values = vec![0.0; n * n];
        let rows = Mutex::new(values.chunks_mut(n).enumerate());
        let threads = thread::available_parallelism().map_or(1, |t| t.get());
        thread::scope(|s| {
            for _ in 0..threads {
                s.spawn(|| loop {
                    // A while let would hold the lock through the row.
                    let Some((j, row)) = rows.lock().unwrap().next() else { break };
                    for (i, v) in row.iter_mut().enumerate() {
                        *v = bath.height(grid.world(i, j));
                    }
                });
            }
        });
        grid.values = values;
        grid
    }

    /// World position of the center of cell (i, j); j grows upward.
    pub fn world(&self, i: usize, j: usize) -> Vec2 {
        let step = 2.0 * self.half / self.n as f64;
        [
            self.center[0] - self.half + (i as f64 + 0.5) * step,
            self.center[1] - self.half + (j as f64 + 0.5) * step,
        ]
    }

    /// Fractional cell coordinates of a world position.
    pub fn cell(&self, p: Vec2) -> (f64, f64) {
        let scale = self.n as f64 / (2.0 * self.half);
        (
            (p[0] - self.center[0] + self.half) * scale - 0.5,
            (p[1] - self.center[1] + self.half) * scale - 0.5,
        )
    }

    pub fn at(&self, i: usize, j: usize) -> f64 {
        self.values[j * self.n + i]
    }

    pub fn peak(&self) -> f64 {
        self.values.iter().fold(0.0, |m, v| m.max(v.abs()))
    }

    /// Count a position into its cell, ignoring positions outside.
    pub fn count(&mut self, p: Vec2) {
        let (x, y) = self.cell(p);
        let (i, j) = (x.round(), y.round());
        if i >= 0.0 && j >= 0.0 && i < self.n as f64 && j < self.n as f64 {
            self.values[j as usize * self.n + i as usize] += 1.0;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::open::Open;

    #[test]
    fn cell_inverts_world() {
        let g = Grid::zeros([1.0, -2.0], 3.0, 50);
        let (x, y) = g.cell(g.world(7, 31));
        assert!((x - 7.0).abs() < 1e-9 && (y - 31.0).abs() < 1e-9);
    }

    #[test]
    fn sampled_cells_hold_the_bath_height_at_their_centers() {
        let mut bath = Open::new(5.0, 1);
        bath.remember(&[[0.2, -0.1]]);
        let g = Grid::sample(&bath, [0.3, 0.1], 1.5, 37);
        for (j, i) in [(0, 0), (5, 36), (36, 20), (18, 18)] {
            assert_eq!(g.at(i, j), bath.height(g.world(i, j)));
        }
    }

    #[test]
    fn counted_position_lands_in_its_cell() {
        let mut g = Grid::zeros([0.0, 0.0], 1.0, 4);
        g.count(g.world(2, 3));
        assert_eq!(g.at(2, 3), 1.0);
        assert_eq!(g.values.iter().sum::<f64>(), 1.0);
    }
}
