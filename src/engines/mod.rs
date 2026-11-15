//! The engines, one concept per file. Nothing outside this folder
//! knows what is in it but `engine`, which hands them out by name.
pub mod classical;
pub mod guided;
pub mod quantum;
pub mod wave;
