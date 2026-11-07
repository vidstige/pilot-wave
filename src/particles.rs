//! Particles with definite positions. Quantum ones ride the field's
//! phase gradient; classical ones obey Newton in the potential. No I/O.
use crate::grid::Vec2;
use crate::wave::Wave;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Kind {
    Quantum,
    Classical,
}

pub struct Particles {
    pub kind: Kind,
    pub positions: Vec<Vec2>,
    /// Used by classical particles only.
    pub velocities: Vec<Vec2>,
}

impl Particles {
    pub fn new(kind: Kind, positions: Vec<Vec2>, velocity: Vec2) -> Particles {
        let velocities = vec![velocity; positions.len()];
        Particles { kind, positions, velocities }
    }

    /// Midpoint step along the guiding velocity field.
    fn guide(&mut self, wave: &Wave, dt: f64) {
        for p in &mut self.positions {
            let v = wave.velocity(*p);
            let mid = [p[0] + 0.5 * dt * v[0], p[1] + 0.5 * dt * v[1]];
            let v = wave.velocity(mid);
            p[0] += dt * v[0];
            p[1] += dt * v[1];
        }
    }

    /// Newton's law in the potential, with walls as hard reflectors.
    fn push(&mut self, wave: &Wave, dt: f64) {
        for (p, v) in self.positions.iter_mut().zip(&mut self.velocities) {
            let f = force(wave, *p);
            v[0] += dt * f[0];
            v[1] += dt * f[1];
            let next = [p[0] + dt * v[0], p[1] + dt * v[1]];
            if wall(wave, next) {
                let (x, y) = (wall(wave, [next[0], p[1]]), wall(wave, [p[0], next[1]]));
                if x || !y { v[0] = -v[0]; }
                if y || !x { v[1] = -v[1]; }
            } else {
                *p = next;
            }
        }
    }

    pub fn step(&mut self, wave: &Wave, dt: f64) {
        match self.kind {
            Kind::Quantum => self.guide(wave, dt),
            Kind::Classical => self.push(wave, dt),
        }
    }
}

fn nearest(wave: &Wave, p: Vec2) -> usize {
    let (x, y) = wave.cell(p);
    let clamp = |c: f64| c.round().clamp(0.0, (wave.n - 1) as f64) as usize;
    wave.index((clamp(x), clamp(y)))
}

fn wall(wave: &Wave, p: Vec2) -> bool {
    wave.keep[nearest(wave, p)] == 0.0
}

/// Minus the potential gradient from neighbouring cells.
fn force(wave: &Wave, p: Vec2) -> Vec2 {
    let k = nearest(wave, p);
    let (n, v) = (wave.n, &wave.potential);
    let scale = -1.0 / (2.0 * wave.dx);
    [
        scale * (v[k + 1] - v[k - 1]) as f64,
        scale * (v[k + n] - v[k - n]) as f64,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::packet::Packet;

    fn walled_wave() -> Wave {
        let mut wave = Wave::new(64, 0.1);
        for j in 0..64 {
            for i in 40..64 {
                let k = wave.index((i, j));
                wave.keep[k] = 0.0;
            }
        }
        wave
    }

    #[test]
    fn classical_particle_reflects_off_a_wall() {
        let wave = walled_wave();
        let mut particles = Particles::new(Kind::Classical, vec![[0.0, 0.0]], [5.0, 0.0]);
        for _ in 0..200 {
            particles.step(&wave, 0.01);
        }
        let p = particles.positions[0];
        assert!(p[0] < 1.0 && particles.velocities[0][0] < 0.0, "{p:?}");
    }

    #[test]
    fn quantum_particles_in_a_moving_packet_move_with_it() {
        let mut wave = Wave::new(128, 0.1);
        let packet = Packet { center: [-2.0, 0.0], width: [1.0, 1.0], wavevector: [4.0, 0.0] };
        packet.write(&mut wave);
        let mut particles = Particles::new(Kind::Quantum, vec![[-2.0, 0.0], [-2.5, 0.4]], [0.0; 2]);
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
