//! The experiments, one concept per file, found by name.
use crate::experiment::Experiment;

pub mod barrier;
pub mod double_slit;
pub mod free;
pub mod mirror;
pub mod single_slit;

pub fn all() -> Vec<Experiment> {
    vec![
        free::free_flight(),
        mirror::mirror(),
        barrier::tunnelling(),
        single_slit::single_slit(),
        double_slit::double_slit(),
    ]
}

pub fn named(name: &str) -> Option<Experiment> {
    all().into_iter().find(|e| e.name == name)
}
