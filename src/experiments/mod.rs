//! The experiments, one concept per file, found by name.
use crate::experiment::Experiment;

pub mod collision;
pub mod double_slit;
pub mod molecule;
pub mod single_slit;

pub fn all() -> Vec<Experiment> {
    vec![
        single_slit::single_slit(),
        double_slit::double_slit(),
        collision::distinguishable(),
        collision::bosons(),
        collision::fermions(),
        molecule::molecule(),
        molecule::molecule_break(),
    ]
}

pub fn named(name: &str) -> Option<Experiment> {
    all().into_iter().find(|e| e.name == name)
}
