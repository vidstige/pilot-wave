//! Runs of a two-particle experiment on a line followed through time:
//! whether each run's point crossed the diagonal x₁ = x₂, how close
//! its two particles came, and how far apart they are. No I/O.
use crate::grid::Vec2;

pub fn separation(p: Vec2) -> f64 {
    (p[0] - p[1]).abs()
}

/// Fraction of the runs whose particles are farther apart than d.
pub fn fraction_apart(points: &[Vec2], d: f64) -> f64 {
    let apart = points.iter().filter(|p| separation(**p) > d).count();
    apart as f64 / points.len() as f64
}

/// Mean and standard deviation of the separations.
pub fn spread(points: &[Vec2]) -> (f64, f64) {
    let n = points.len() as f64;
    let mean = points.iter().map(|p| separation(*p)).sum::<f64>() / n;
    let var = points.iter().map(|p| (separation(*p) - mean).powi(2)).sum::<f64>() / n;
    (mean, var.sqrt())
}

pub struct Encounter {
    /// Whether x₁ < x₂ at the start.
    side: Vec<bool>,
    crossed: Vec<bool>,
    /// Each run's least separation so far.
    pub closest: Vec<f64>,
}

impl Encounter {
    pub fn new(points: &[Vec2]) -> Encounter {
        Encounter {
            side: points.iter().map(|p| p[0] < p[1]).collect(),
            crossed: vec![false; points.len()],
            closest: points.iter().map(|p| separation(*p)).collect(),
        }
    }

    pub fn record(&mut self, points: &[Vec2]) {
        for (k, p) in points.iter().enumerate() {
            self.crossed[k] |= (p[0] < p[1]) != self.side[k];
            self.closest[k] = self.closest[k].min(separation(*p));
        }
    }

    /// Fraction of the runs that ever crossed the diagonal.
    pub fn crossed(&self) -> f64 {
        self.crossed.iter().filter(|c| **c).count() as f64 / self.crossed.len() as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crossing_is_remembered_after_coming_back() {
        let mut encounter = Encounter::new(&[[-1.0, 1.0], [2.0, 0.0]]);
        encounter.record(&[[0.5, 0.2], [1.5, 0.5]]);
        encounter.record(&[[-1.0, 1.0], [2.0, 0.0]]);
        assert_eq!(encounter.crossed(), 0.5);
        assert!((encounter.closest[0] - 0.3).abs() < 1e-12);
        assert_eq!(encounter.closest[1], 1.0);
    }

    #[test]
    fn spread_is_unchanged_by_exchange_and_translation() {
        let points = [[0.0, 1.0], [3.0, 0.0], [-2.0, -4.0]];
        let moved = points.map(|p| [p[1] + 5.0, p[0] + 5.0]);
        assert_eq!(spread(&points), spread(&moved));
        assert_eq!(spread(&points), (2.0, (2.0f64 / 3.0).sqrt()));
        assert_eq!(fraction_apart(&points, 1.5), 2.0 / 3.0);
    }
}
