//! The experiments, one per file, found by name.
use crate::experiment::Experiment;

pub mod single_slit;

pub fn all() -> Vec<Experiment> {
    vec![single_slit::single_slit()]
}

pub fn named(name: &str) -> Option<Experiment> {
    all().into_iter().find(|e| e.name == name)
}
