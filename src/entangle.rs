//! Entanglement of two particles on a line, over runs or over the
//! field: the correlation of x₁ and x₂, the mutual information of
//! their binned joint distribution, and the comparison of paired runs
//! with and without a wall on particle 2. No I/O.
use crate::grid::Vec2;
use crate::wave::Wave;

/// Width of the bins along each axis, a whole number of cells of 0.1
/// or 0.05, so the line holds a whole number of them.
pub const BIN: f64 = 0.4;
/// Fewest runs a bin must expect to count in a χ² against a density.
const MIN_EXPECTED: f64 = 5.0;

/// A configuration and its weight: one run, or one cell of the field
/// weighted by its density.
pub type Weighted = (Vec2, f64);

pub fn runs(points: &[Vec2]) -> Vec<Weighted> {
    points.iter().map(|p| (*p, 1.0)).collect()
}

pub fn field(wave: &Wave) -> Vec<Weighted> {
    let n = wave.n;
    let at = |k: usize| (wave.position((k % n, k / n)), wave.density(k) as f64);
    (0..n * n).map(at).collect()
}

/// Correlation coefficient of x₁ and x₂.
pub fn correlation(points: &[Weighted]) -> f64 {
    let total: f64 = points.iter().map(|(_, w)| w).sum();
    let mean = |f: fn(Vec2) -> f64| {
        points.iter().map(|(p, w)| w * f(*p)).sum::<f64>() / total
    };
    let (x, y) = (mean(|p| p[0]), mean(|p| p[1]));
    let xy = mean(|p| p[0] * p[1]) - x * y;
    let xx = mean(|p| p[0] * p[0]) - x * x;
    let yy = mean(|p| p[1] * p[1]) - y * y;
    xy / (xx * yy).sqrt()
}

/// Bins of BIN across a line |x| ≤ half.
pub fn bins(half: f64) -> usize {
    (2.0 * half / BIN).round() as usize
}

/// Bin of x on a line |x| ≤ half, the ends taking what lies beyond.
fn bin(x: f64, half: f64) -> usize {
    (((x + half) / BIN) as usize).min(bins(half) - 1)
}

/// Mutual information of x₁ and x₂ in nats, from their joint
/// distribution binned in squares of BIN on the line |x| ≤ half: the
/// Kullback-Leibler divergence of the binned joint from the product
/// of its marginals, zero exactly when it is that product.
pub fn mutual_information(points: &[Weighted], half: f64) -> f64 {
    let n = bins(half);
    let mut joint = vec![0.0; n * n];
    let mut margins = [vec![0.0; n], vec![0.0; n]];
    for (p, w) in points {
        let (i, j) = (bin(p[0], half), bin(p[1], half));
        joint[j * n + i] += w;
        margins[0][i] += w;
        margins[1][j] += w;
    }
    let total: f64 = margins[0].iter().sum();
    let term = |k: usize| {
        let p = joint[k] / total;
        let q = margins[0][k % n] * margins[1][k / n] / (total * total);
        if p > 0.0 { p * (p / q).ln() } else { 0.0 }
    };
    (0..n * n).map(term).sum()
}

/// Correlation over the runs and in the field, then their mutual
/// information.
pub fn measure(wave: &Wave, points: &[Vec2]) -> [f64; 4] {
    let half = wave.extent() / 2.0;
    let (runs, field) = (runs(points), field(wave));
    [
        correlation(&runs),
        correlation(&field),
        mutual_information(&runs, half),
        mutual_information(&field, half),
    ]
}

/// Runs per bin of one coordinate on the line |x| ≤ half.
pub fn histogram(points: &[Vec2], axis: usize, half: f64) -> Vec<f64> {
    let mut counts = vec![0.0; bins(half)];
    for p in points {
        counts[bin(p[axis], half)] += 1.0;
    }
    counts
}

/// The largest difference between two histograms of the same runs,
/// as a share of the runs, and its statistical error in that bin
/// were the two independent samples of one distribution.
pub fn largest_difference(a: &[f64], b: &[f64]) -> (f64, f64) {
    let total: f64 = a.iter().sum();
    let gap = |k: usize| (a[k] - b[k]).abs();
    let k = (0..a.len()).max_by(|&i, &j| gap(i).total_cmp(&gap(j))).unwrap();
    (gap(k) / total, (a[k] + b[k]).sqrt() / total)
}

/// χ² of two histograms of the same number of runs, over the bins
/// either has runs in, and the number of those bins: were the two
/// independent samples of one distribution, χ² would be about the
/// number of bins.
pub fn chi_square(a: &[f64], b: &[f64]) -> (f64, usize) {
    let full = a.iter().zip(b).filter(|(p, q)| *p + *q > 0.0);
    full.fold((0.0, 0), |(chi, n), (p, q)| (chi + (p - q).powi(2) / (p + q), n + 1))
}

/// Runs expected in each bin from a density per cell of width dx.
pub fn expected(density: &[f64], dx: f64, runs: usize) -> Vec<f64> {
    let cells = (BIN / dx).round() as usize;
    let bin = |c: &[f64]| c.iter().sum::<f64>() * dx * runs as f64;
    density.chunks(cells).map(bin).collect()
}

/// χ² of a histogram against the runs expected in each bin, over the
/// bins expecting at least MIN_EXPECTED, and the number of those bins,
/// which χ² should be about if the runs follow the density.
pub fn goodness(h: &[f64], expected: &[f64]) -> (f64, usize) {
    let full = h.iter().zip(expected).filter(|(_, e)| **e >= MIN_EXPECTED);
    full.fold((0.0, 0), |(chi, n), (h, e)| (chi + (h - e).powi(2) / e, n + 1))
}

/// Largest and median change of a coordinate between paired runs.
pub fn moved(a: &[Vec2], b: &[Vec2], axis: usize) -> (f64, f64) {
    let mut d: Vec<f64> = a.iter().zip(b).map(|(p, q)| (p[axis] - q[axis]).abs()).collect();
    d.sort_by(f64::total_cmp);
    (d[d.len() - 1], d[d.len() / 2])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pair;

    const HALF: f64 = 6.4;

    /// Weights on the centres of 9 × 9 bins around the origin.
    fn grid(f: impl Fn(f64, f64) -> f64) -> Vec<Weighted> {
        let x = |i: i32| (i as f64 - 4.0) * BIN;
        let cell = |k: i32| ([x(k % 9), x(k / 9)], f(x(k % 9), x(k / 9)));
        (0..81).map(cell).collect()
    }

    #[test]
    fn points_on_a_line_are_fully_correlated_either_way() {
        let up: Vec<Vec2> = (0..5).map(|i| [i as f64, 2.0 * i as f64 + 1.0]).collect();
        let down: Vec<Vec2> = up.iter().map(|p| [p[0], -p[1]]).collect();
        assert!((correlation(&runs(&up)) - 1.0).abs() < 1e-12);
        assert!((correlation(&runs(&down)) + 1.0).abs() < 1e-12);
    }

    #[test]
    fn product_weights_share_no_information() {
        let points = grid(|x, y| (-x * x).exp() * (1.0 + y * y));
        assert!(mutual_information(&points, HALF).abs() < 1e-12);
        assert!(correlation(&points).abs() < 1e-12);
    }

    #[test]
    fn a_diagonal_over_k_bins_shares_ln_k() {
        let points = grid(|x, y| if x == y { 1.0 } else { 0.0 });
        let information = mutual_information(&points, HALF);
        assert!((information - 9f64.ln()).abs() < 1e-12);
    }

    #[test]
    fn gaussian_pair_matches_its_closed_form() {
        let mut wave = Wave::new(128, 0.1);
        let (c, r) = (1.0f64, 0.6f64);
        pair::fill(&mut wave, |x1, x2| {
            let (u, v) = ((x1 + x2) / 2.0, x1 - x2);
            [(-u * u / (4.0 * c * c) - v * v / (4.0 * r * r)).exp(), 0.0]
        });
        let rho = (c * c - r * r / 4.0) / (c * c + r * r / 4.0);
        let points = field(&wave);
        assert!((correlation(&points) - rho).abs() < 1e-3);
        let gaussian = -0.5 * (1.0 - rho * rho).ln();
        let information = mutual_information(&points, wave.extent() / 2.0);
        assert!((information - gaussian).abs() < 0.1 * gaussian, "{information}");
    }

    #[test]
    fn moved_runs_report_the_largest_and_median_change() {
        let a = [[0.0, 0.0], [1.0, 5.0], [2.0, -1.0]];
        let b = [[0.5, 9.0], [1.0, 5.0], [-1.0, 3.0]];
        assert_eq!(moved(&a, &b, 0), (3.0, 0.5));
        let h = histogram(&a, 0, HALF);
        assert_eq!(largest_difference(&h, &h).0, 0.0);
        let g = histogram(&b, 0, HALF);
        let (gap, error) = largest_difference(&h, &g);
        assert_eq!((gap, error), (1.0 / 3.0, 1.0 / 3.0));
    }
}
