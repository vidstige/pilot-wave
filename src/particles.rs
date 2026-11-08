//! Particles with definite positions riding the field's guiding
//! velocity. No I/O.
use std::thread;

use crate::wave::{Vec2, Wave};

/// Longest move in one substep, in cells, and the most substeps a
/// step may take. Next to a node of ψ the guiding velocity is huge,
/// so the particle creeps up to it and turns back as it should,
/// instead of jumping across.
const SUBSTEP: f64 = 0.25;
const MAX_SUBSTEPS: usize = 64;
/// Particles per thread: fewer is not worth a thread's start-up.
const PARTICLES_PER_THREAD: usize = 2000;

pub struct Particles {
    pub positions: Vec<Vec2>,
}

impl Particles {
    pub fn new(positions: Vec<Vec2>) -> Particles {
        Particles { positions }
    }

    /// Every particle takes its guided step, in chunks on all cores.
    pub fn step(&mut self, wave: &Wave, dt: f64) {
        let cores = thread::available_parallelism().map_or(1, |t| t.get());
        let chunk = (self.positions.len() / cores).max(PARTICLES_PER_THREAD);
        thread::scope(|s| {
            for positions in self.positions.chunks_mut(chunk) {
                s.spawn(move || {
                    for p in positions {
                        guide(wave, p, dt);
                    }
                });
            }
        });
    }
}

/// Midpoint steps along the guiding velocity field, shortened
/// wherever the field would move the particle more than a fraction
/// of a cell.
fn guide(wave: &Wave, p: &mut Vec2, dt: f64) {
    let reach = SUBSTEP * wave.dx;
    let mut left = dt;
    for substep in 0..MAX_SUBSTEPS {
        let v = wave.velocity(*p);
        let speed = v[0].hypot(v[1]).max(f64::MIN_POSITIVE);
        let h = if substep + 1 == MAX_SUBSTEPS { left } else { left.min(reach / speed) };
        let mid = [p[0] + 0.5 * h * v[0], p[1] + 0.5 * h * v[1]];
        let v = wave.velocity(mid);
        p[0] += h * v[0];
        p[1] += h * v[1];
        left -= h;
        if left <= 0.0 {
            break;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::packet::Packet;

    #[test]
    fn quantum_particles_in_a_moving_packet_move_with_it() {
        let mut wave = Wave::new(128, 0.1);
        let packet = Packet { center: [-2.0, 0.0], width: [1.0, 1.0], wavevector: [4.0, 0.0] };
        packet.write(&mut wave);
        let mut particles = Particles::new(vec![[-2.0, 0.0], [-2.5, 0.4]]);
        let dt = 0.002;
        for _ in 0..200 {
            wave.step(dt);
            particles.step(&wave, dt);
        }
        for p in &particles.positions {
            assert!((p[0] - (-2.0 + 4.0 * dt * 200.0)).abs() < 0.6, "{p:?}");
        }
    }
}
