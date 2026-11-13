//! The experiments, one concept per file, found by name.
use crate::experiment::Experiment;

pub mod barrier;
pub mod double_slit;
pub mod free;
pub mod single_slit;
pub mod spin;

pub fn all() -> Vec<Experiment> {
    vec![
        free::free_flight(),
        free::mirror(),
        barrier::tunnelling(),
        single_slit::single_slit(),
        double_slit::double_slit(),
        spin::spin_in_plane(),
        spin::spin_normal(),
    ]
}

pub fn named(name: &str) -> Option<Experiment> {
    all().into_iter().find(|e| e.name == name)
}
