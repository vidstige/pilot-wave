//! Particles carried along a velocity field. No I/O.
use std::thread;

use crate::plane::Vec2;

/// Longest move in one substep, in cells, and the most substeps a
/// step may take. Next to a node of the wave the velocity is huge, so
/// the particle creeps up to it and turns back as it should, instead
/// of jumping across.
const SUBSTEP: f64 = 0.25;
const MAX_SUBSTEPS: usize = 64;
/// Particles per thread: fewer is not worth a thread's start-up.
const PARTICLES_PER_THREAD: usize = 2000;

/// Midpoint steps along the velocity field, shortened wherever the
/// field would move the particle further than `reach`.
fn guide(velocity: impl Fn(Vec2) -> Vec2, reach: f64, p: &mut Vec2, dt: f64) {
    let mut left = dt;
    for substep in 0..MAX_SUBSTEPS {
        let v = velocity(*p);
        let speed = v[0].hypot(v[1]).max(f64::MIN_POSITIVE);
        let h = if substep + 1 == MAX_SUBSTEPS { left } else { left.min(reach / speed) };
        let mid = [p[0] + 0.5 * h * v[0], p[1] + 0.5 * h * v[1]];
        let v = velocity(mid);
        p[0] += h * v[0];
        p[1] += h * v[1];
        left -= h;
        if left <= 0.0 {
            break;
        }
    }
}

/// Every particle takes its step along the velocity field of a grid
/// of cells dx wide, in chunks on all cores.
pub fn carry(positions: &mut [Vec2], dx: f64, dt: f64, velocity: impl Fn(Vec2) -> Vec2 + Sync) {
    let cores = thread::available_parallelism().map_or(1, |t| t.get());
    let chunk = (positions.len() / cores).max(PARTICLES_PER_THREAD);
    let velocity = &velocity;
    thread::scope(|s| {
        for positions in positions.chunks_mut(chunk) {
            s.spawn(move || {
                for p in positions {
                    guide(velocity, SUBSTEP * dx, p, dt);
                }
            });
        }
    });
}

#[cfg(test)]
mod tests {
    use super::super::wave::Wave;
    use super::*;
    use crate::packet::Packet;

    #[test]
    fn particles_in_a_moving_packet_move_with_it() {
        let mut wave = Wave::new(128, 0.1);
        let packet = Packet { center: [-2.0, 0.0], width: [1.0, 1.0], wavevector: [4.0, 0.0] };
        wave.fill(|p| packet.amplitude(p));
        let mut positions = [[-2.0, 0.0], [-2.5, 0.4]];
        let dt = 0.002;
        for _ in 0..200 {
            wave.step(dt);
            carry(&mut positions, wave.dx, dt, |p| wave.velocity(p));
        }
        for p in positions {
            assert!((p[0] - (-2.0 + 4.0 * dt * 200.0)).abs() < 0.6, "{p:?}");
        }
    }
}
