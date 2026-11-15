//! The quantum engine: a wavefunction stepped by the Schrödinger
//! equation, and particles carried by its probability current. No I/O.
use super::guided::carry;
use super::wave::{self, Wave};
use crate::engine::Engine;
use crate::experiment::Experiment;
use crate::plane::Vec2;

pub struct Quantum {
    wave: Wave,
    positions: Vec<Vec2>,
    dt: f64,
}

/// The wave at the start, to be stepped by dt: sponge at the edges,
/// walls cut out after it so they stay whole, the potential, then the
/// experiment's amplitude.
fn build(e: &Experiment, dt: f64) -> Wave {
    let stage = &e.stage;
    let mut wave = Wave::new(stage.cells, stage.dx);
    wave.absorb_edges(dt);
    for k in 0..stage.cells * stage.cells {
        let p = stage.position(k);
        if stage.walled(p) {
            wave.keep[k] = 0.0;
        }
        wave.potential[k] = (stage.potential)(p) as f32;
    }
    wave.fill(|p| e.amplitude(p));
    wave
}

impl Quantum {
    pub fn new(e: &Experiment, positions: Vec<Vec2>) -> Quantum {
        let dt = wave::stable_dt(e.stage.dx, e.steepest());
        Quantum { wave: build(e, dt), positions, dt }
    }
}

impl Engine for Quantum {
    fn dt(&self) -> f64 {
        self.dt
    }

    fn step(&mut self) {
        self.wave.step(self.dt);
        let wave = &self.wave;
        carry(&mut self.positions, wave.dx, self.dt, |p| wave.velocity(p));
    }

    fn positions(&self) -> &[Vec2] {
        &self.positions
    }

    fn field(&self) -> Option<[&[f32]; 2]> {
        Some([&self.wave.re, &self.wave.im])
    }
}
