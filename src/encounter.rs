//! Statistics of two-particle runs on a line: whether each run's
//! point crossed the diagonal x₁ = x₂, how far apart the particles
//! are, how correlated their positions are, and how much a coordinate
//! changed against paired runs of a twin. No I/O.
use crate::wave::Vec2;

pub fn separation(p: Vec2) -> f64 {
    (p[0] - p[1]).abs()
}

/// Fraction of the runs whose particles are within d of each other.
pub fn fraction_within(points: &[Vec2], d: f64) -> f64 {
    let close = points.iter().filter(|p| separation(**p) < d).count();
    close as f64 / points.len() as f64
}

/// Correlation coefficient of x₁ and x₂ over the runs.
pub fn correlation(points: &[Vec2]) -> f64 {
    let mean = |f: fn(&Vec2) -> f64| points.iter().map(f).sum::<f64>() / points.len() as f64;
    let (x, y) = (mean(|p| p[0]), mean(|p| p[1]));
    let xy = mean(|p| p[0] * p[1]) - x * y;
    let xx = mean(|p| p[0] * p[0]) - x * x;
    let yy = mean(|p| p[1] * p[1]) - y * y;
    xy / (xx * yy).sqrt()
}

/// Median change of a coordinate between paired runs.
pub fn moved(a: &[Vec2], b: &[Vec2], axis: usize) -> f64 {
    let mut d: Vec<f64> = a.iter().zip(b).map(|(p, q)| (p[axis] - q[axis]).abs()).collect();
    d.sort_by(f64::total_cmp);
    d[d.len() / 2]
}

/// Which runs have ever crossed the diagonal.
pub struct Crossings {
    /// Whether x₁ < x₂ at the start.
    side: Vec<bool>,
    crossed: Vec<bool>,
}

impl Crossings {
    pub fn new(points: &[Vec2]) -> Crossings {
        Crossings { side: points.iter().map(|p| p[0] < p[1]).collect(), crossed: vec![false; points.len()] }
    }

    pub fn record(&mut self, points: &[Vec2]) {
        for (k, p) in points.iter().enumerate() {
            self.crossed[k] |= (p[0] < p[1]) != self.side[k];
        }
    }

    /// Fraction of the runs that ever crossed.
    pub fn fraction(&self) -> f64 {
        self.crossed.iter().filter(|c| **c).count() as f64 / self.crossed.len() as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crossing_is_remembered_after_coming_back() {
        let start = [[-1.0, 1.0], [2.0, 0.0]];
        let mut crossings = Crossings::new(&start);
        crossings.record(&[[0.5, 0.2], [1.5, 0.5]]);
        crossings.record(&start);
        assert_eq!(crossings.fraction(), 0.5);
    }

    #[test]
    fn points_on_a_line_are_fully_correlated_either_way() {
        let up: Vec<Vec2> = (0..5).map(|i| [i as f64, 2.0 * i as f64 + 1.0]).collect();
        let down: Vec<Vec2> = up.iter().map(|p| [p[0], -p[1]]).collect();
        assert!((correlation(&up) - 1.0).abs() < 1e-12);
        assert!((correlation(&down) + 1.0).abs() < 1e-12);
    }

    #[test]
    fn moved_is_the_median_change_and_within_counts_close_pairs() {
        let a = [[0.0, 0.0], [1.0, 5.0], [2.0, -1.0]];
        let b = [[0.5, 9.0], [1.0, 5.0], [-1.0, 3.0]];
        assert_eq!(moved(&a, &b, 0), 0.5);
        assert_eq!(fraction_within(&a, 1.5), 1.0 / 3.0);
    }
}
