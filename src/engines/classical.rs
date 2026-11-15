//! The classical engine: Newton's particles and no wave. A run flies
//! through the potential, keeping its energy, and bounces off the
//! walls. No I/O.
//!
//! A beam's runs differ in velocity as well as in position: each
//! starts with the particle's velocity and a share of the beam's
//! divergence, a Gaussian ½ over the packet's width along each axis,
//! in units where ħ = m = 1. That is the spread a beam this narrow is
//! seen to have, taken here as a property of the source.
use crate::engine::Engine;
use crate::experiment::{Experiment, Particle, Stage};
use crate::plane::Vec2;
use crate::rng::Rng;

/// The time a step takes: a beam's particle moves an eighth of a cell.
const STEP: f64 = 0.002;
/// Seed of the runs' velocities.
const SEED: u64 = 2;

pub struct Classical {
    stage: Stage,
    positions: Vec<Vec2>,
    velocities: Vec<Vec2>,
}

/// One particle a step on, an axis at a time. Where the potential
/// rises the particle pays for it in speed along that axis, and it
/// turns back where it cannot pay, as it does at a wall: its energy
/// is kept exactly, over a step in the potential as on a slope. Only
/// a particle exactly at rest along an axis never starts along it.
fn advance(stage: &Stage, p: &mut Vec2, v: &mut Vec2) {
    for axis in 0..2 {
        let mut moved = *p;
        moved[axis] += STEP * v[axis];
        let climb = (stage.potential)(moved) - (stage.potential)(*p);
        let left = v[axis] * v[axis] - 2.0 * climb;
        if stage.walled(moved) || left < 0.0 {
            v[axis] = -v[axis];
        } else {
            v[axis] = left.sqrt().copysign(v[axis]);
            *p = moved;
        }
    }
}

impl Classical {
    pub fn new(e: &Experiment, positions: Vec<Vec2>) -> Classical {
        let Particle { velocity, width, .. } = e.particle;
        let mut rng = Rng::new(SEED);
        let mut draw = || [0, 1].map(|a| velocity[a] + rng.normal() / (2.0 * width[a]));
        let velocities = positions.iter().map(|_| draw()).collect();
        Classical { stage: e.stage.clone(), positions, velocities }
    }
}

impl Engine for Classical {
    fn dt(&self) -> f64 {
        STEP
    }

    fn step(&mut self) {
        for (p, v) in self.positions.iter_mut().zip(&mut self.velocities) {
            advance(&self.stage, p, v);
        }
    }

    fn positions(&self) -> &[Vec2] {
        &self.positions
    }

    fn field(&self) -> Option<[&[f32]; 2]> {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::experiment::band;

    /// A barrier 10 high from y = 0 to 0.4, and a wall from y = 3.
    fn barrier(p: Vec2) -> f64 {
        if (0.0..0.4).contains(&p[1]) { 10.0 } else { 0.0 }
    }

    fn stage() -> Stage {
        let walls = vec![band(1, 3.0, 3.4)];
        Stage { cells: 128, dx: 0.1, scale: 1.0, walls, potential: barrier }
    }

    /// Where a particle from y = −1 is and how it moves after a time.
    fn flown(speed: f64, time: f64) -> (Vec2, Vec2) {
        let (mut p, mut v) = ([0.0, -1.0], [0.0, speed]);
        for _ in 0..(time / STEP) as usize {
            advance(&stage(), &mut p, &mut v);
        }
        (p, v)
    }

    #[test]
    fn a_particle_with_the_energy_crosses_the_barrier_and_keeps_its_speed() {
        let (p, v) = flown(5.0, 0.5);
        assert!(p[1] > 0.4 && (v[1] - 5.0).abs() < 1e-9, "{p:?} {v:?}");
    }

    #[test]
    fn a_particle_short_of_the_barriers_height_turns_back() {
        let (p, v) = flown(4.0, 0.5);
        assert!(p[1] < 0.0 && (v[1] + 4.0).abs() < 1e-9, "{p:?} {v:?}");
    }

    #[test]
    fn a_particle_bounces_off_a_wall() {
        let (p, v) = flown(5.0, 1.2);
        assert!((0.4..3.0).contains(&p[1]) && (v[1] + 5.0).abs() < 1e-9, "{p:?} {v:?}");
    }
}
