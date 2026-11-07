//! Bessel J0 and J1, rational approximations from Abramowitz & Stegun
//! 9.4, accurate to about 1e-8.
use std::f64::consts::PI;

const SWITCH: f64 = 8.0;

/// Horner evaluation, coefficients from the highest power down.
fn poly(c: &[f64], x: f64) -> f64 {
    c.iter().fold(0.0, |acc, &a| acc * x + a)
}

fn asymptotic(x: f64, phase: f64, p: &[f64], q: &[f64]) -> f64 {
    let z = SWITCH / x;
    let y = z * z;
    let xx = x - phase;
    let amplitude = (2.0 / (PI * x)).sqrt();
    amplitude * (xx.cos() * poly(p, y) - z * xx.sin() * poly(q, y))
}

pub fn j0(x: f64) -> f64 {
    let ax = x.abs();
    if ax < SWITCH {
        let y = ax * ax;
        let p = poly(
            &[-184.9052456, 77392.33017, -11214424.18, 651619640.7,
              -13362590354.0, 57568490574.0], y);
        let q = poly(
            &[1.0, 267.8532712, 59272.64853, 9494680.718,
              1029532985.0, 57568490411.0], y);
        p / q
    } else {
        asymptotic(
            ax, 0.785398164,
            &[0.2093887211e-6, -0.2073370639e-5, 0.2734510407e-4,
              -0.1098628627e-2, 1.0],
            &[-0.934935152e-7, 0.7621095161e-6, -0.6911147651e-5,
              0.1430488765e-3, -0.1562499995e-1])
    }
}

pub fn j1(x: f64) -> f64 {
    let ax = x.abs();
    let value = if ax < SWITCH {
        let y = ax * ax;
        let p = ax * poly(
            &[-30.16036606, 15704.48260, -2972611.439, 242396853.1,
              -7895059235.0, 72362614232.0], y);
        let q = poly(
            &[1.0, 376.9991397, 99447.43394, 18583304.74,
              2300535178.0, 144725228442.0], y);
        p / q
    } else {
        asymptotic(
            ax, 2.356194491,
            &[-0.240337019e-6, 0.2457520174e-5, -0.3516396496e-4,
              0.183105e-2, 1.0],
            &[0.105787412e-6, -0.88228987e-6, 0.8449199096e-5,
              -0.2002690873e-3, 0.04687499995])
    };
    value * x.signum()
}

/// J1(x)/x, continuous through x = 0 where it equals 1/2.
pub fn j1_over_x(x: f64) -> f64 {
    if x.abs() < 1e-6 { 0.5 } else { j1(x) / x }
}

#[cfg(test)]
mod tests {
    use super::*;

    const H: f64 = 1e-5;

    #[test]
    fn j0_derivative_is_minus_j1() {
        let worst = (0..2001)
            .map(|i| -30.0 + 0.03 * i as f64)
            .map(|x| ((j0(x + H) - j0(x - H)) / (2.0 * H) + j1(x)).abs())
            .fold(0.0, f64::max);
        assert!(worst < 1e-6, "{worst}");
    }

    #[test]
    fn j0_vanishes_at_its_first_zero() {
        assert!(j0(2.404825557695773).abs() < 1e-7);
    }

    #[test]
    fn j1_over_x_is_continuous_at_zero() {
        assert!((j1_over_x(1e-4) - 0.5).abs() < 1e-6);
    }
}
