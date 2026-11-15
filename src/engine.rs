//! What a run needs of the physics: particles that move a step at a
//! time, and a field to show if there is one. No I/O.
use crate::engines::classical::Classical;
use crate::engines::quantum::Quantum;
use crate::experiment::Experiment;
use crate::plane::Vec2;

pub trait Engine {
    /// The time one step takes.
    fn dt(&self) -> f64;
    fn step(&mut self);
    /// Where every run's particle is.
    fn positions(&self) -> &[Vec2];
    /// The engine's field for the pictures, two quadratures to each
    /// cell of the stage: its square is drawn as brightness, its
    /// phase as hue.
    fn field(&self) -> Option<[&[f32]; 2]>;
}

pub const NAMES: [&str; 2] = ["quantum", "classical"];

/// The named engine set up for an experiment, its runs starting at
/// the given positions.
pub fn named(name: &str, e: &Experiment, starts: Vec<Vec2>) -> Option<Box<dyn Engine>> {
    match name {
        "quantum" => Some(Box::new(Quantum::new(e, starts))),
        "classical" => Some(Box::new(Classical::new(e, starts))),
        _ => None,
    }
}
