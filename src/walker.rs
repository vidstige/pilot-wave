//! Droplets stepping on the field. No I/O.
//!
//! A bouncer at rest sits on top of its own wave hill. Once
//! WAVENUMBER² · FORCE · memory² > 2 · DRAG that hill is steep enough
//! that the slightest noise grows into steady walking.
use crate::bath::{Bath, Vec2};
use crate::rng::Rng;

pub const FORCE: f64 = 6e-4;
pub const DRAG: f64 = 0.3;
const NOISE: f64 = 1e-4;

pub struct Walkers {
    pub positions: Vec<Vec2>,
    pub velocities: Vec<Vec2>,
}

impl Walkers {
    pub fn new(positions: Vec<Vec2>, velocities: Vec<Vec2>) -> Walkers {
        Walkers { positions, velocities }
    }
}

/// Advance every droplet by one bounce under an external force.
pub fn step(
    bath: &mut (impl Bath + ?Sized),
    walkers: &mut Walkers,
    external: impl Fn(Vec2) -> Vec2,
    rng: &mut Rng,
) {
    for (p, v) in walkers.positions.iter_mut().zip(&mut walkers.velocities) {
        let g = bath.slope(*p);
        let e = external(*p);
        for k in 0..2 {
            v[k] = (1.0 - DRAG) * v[k] - FORCE * g[k] + e[k]
                + NOISE * rng.normal();
            p[k] += v[k];
        }
    }
    bath.remember(&walkers.positions);
}

/// Bounce droplets that left the disc back off its rim.
pub fn confine(walkers: &mut Walkers, radius: f64) {
    for (p, v) in walkers.positions.iter_mut().zip(&mut walkers.velocities) {
        let r = (p[0] * p[0] + p[1] * p[1]).sqrt();
        if r > radius {
            let n = [p[0] / r, p[1] / r];
            let outward = v[0] * n[0] + v[1] * n[1];
            v[0] -= 2.0 * outward * n[0];
            v[1] -= 2.0 * outward * n[1];
            *p = [n[0] * radius, n[1] * radius];
        }
    }
}

pub fn spring(stiffness: f64) -> impl Fn(Vec2) -> Vec2 {
    move |p| [-stiffness * p[0], -stiffness * p[1]]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bath::distance;
    use crate::open::Open;
    use std::f64::consts::TAU;

    fn run(memory: f64, w: &mut Walkers, bounces: usize) -> Vec<Vec<Vec2>> {
        let mut f = Open::new(memory, w.positions.len());
        let mut rng = Rng::new(1);
        (0..bounces)
            .map(|_| {
                step(&mut f, w, spring(0.0), &mut rng);
                w.positions.clone()
            })
            .collect()
    }

    fn speeds(t: &[Vec<Vec2>], i: usize) -> Vec<f64> {
        t.windows(2).map(|w| distance(w[0][i], w[1][i])).collect()
    }

    #[test]
    fn short_memory_bouncer_stays_put() {
        let mut w = Walkers::new(vec![[0.0; 2]], vec![[0.0; 2]]);
        let t = run(3.0, &mut w, 1500);
        assert!(distance(t[1499][0], [0.0; 2]) < 0.05);
    }

    #[test]
    fn long_memory_bouncer_walks_at_steady_speed() {
        let mut w = Walkers::new(vec![[0.0; 2]], vec![[0.0; 2]]);
        let t = run(10.0, &mut w, 1500);
        let v = &speeds(&t, 0)[1200..];
        let mean = v.iter().sum::<f64>() / v.len() as f64;
        let var = v.iter().map(|x| (x - mean).powi(2)).sum::<f64>()
            / v.len() as f64;
        assert!(mean > 0.02, "{mean}");
        assert!(var.sqrt() < 0.1 * mean);
    }

    fn wrap(angle: f64) -> f64 {
        angle - (angle / TAU).round() * TAU
    }

    #[test]
    fn nudged_close_pair_locks_into_orbit() {
        let mut w = Walkers::new(
            vec![[-0.35, 0.0], [0.35, 0.0]],
            vec![[0.0, 0.02], [0.0, -0.02]],
        );
        let t = run(8.0, &mut w, 4000);
        let sep: Vec<f64> =
            t[3000..].iter().map(|p| distance(p[0], p[1])).collect();
        let mean = sep.iter().sum::<f64>() / sep.len() as f64;
        assert!(sep.iter().all(|s| (s - mean).abs() < 0.02 * mean));
        let angle = |p: &Vec<Vec2>| {
            (p[1][1] - p[0][1]).atan2(p[1][0] - p[0][0])
        };
        let turned: f64 =
            t.windows(2).map(|w| wrap(angle(&w[1]) - angle(&w[0]))).sum();
        assert!(turned.abs() > 2.0 * TAU, "{turned}");
    }
}
