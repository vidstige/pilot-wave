//! The live viewer's running simulation and its frames, free of UI.
use crate::frame;
use crate::image::Image;
use crate::scene::{scene, Simulation};

pub const SCENES: [&str; 4] = ["single", "pair", "trap", "corral"];
/// Cells sampled for the top view, drawn a pixel each.
const TOP_SAMPLES: usize = 180;
/// The oblique view interpolates the surface, so it samples coarser
/// than it draws.
const OBLIQUE_SAMPLES: usize = 120;
const OBLIQUE_PIXELS: usize = 360;
/// One ray per pixel keeps the live oblique view fluid.
const OBLIQUE_SAMPLES_PER_PIXEL: usize = 1;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum View {
    Top,
    Oblique,
}

pub struct Live {
    pub name: &'static str,
    pub memory: f64,
    pub bounces: usize,
    simulation: Simulation,
}

impl Live {
    pub fn new(name: &'static str, memory: f64) -> Live {
        let simulation = Simulation::new(scene(name, memory));
        Live { name, memory, bounces: 0, simulation }
    }

    /// Start over if the scene or memory differ from the running ones.
    pub fn select(&mut self, name: &'static str, memory: f64) {
        if (name, memory) != (self.name, self.memory) {
            *self = Live::new(name, memory);
        }
    }

    pub fn advance(&mut self, bounces: usize) {
        for _ in 0..bounces {
            self.simulation.bounce();
        }
        self.bounces += bounces;
    }

    pub fn frame(&self, view: View) -> Image {
        let sim = &self.simulation;
        match view {
            View::Top => frame::top(sim, TOP_SAMPLES),
            View::Oblique => {
                frame::oblique(
                    sim, OBLIQUE_SAMPLES, OBLIQUE_PIXELS, OBLIQUE_SAMPLES_PER_PIXEL)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn advanced(name: &'static str, bounces: usize) -> Live {
        let mut live = Live::new(name, 20.0);
        live.advance(bounces);
        live
    }

    #[test]
    fn only_a_new_scene_or_memory_starts_over() {
        let mut live = advanced("pair", 30);
        live.select("pair", 20.0);
        assert_eq!(live.bounces, 30);
        live.select("pair", 21.0);
        assert_eq!(live.bounces, 0);
        live.advance(5);
        live.select("trap", 21.0);
        assert_eq!((live.name, live.bounces), ("trap", 0));
    }

    #[test]
    fn frames_change_as_the_simulation_advances() {
        let mut live = advanced("corral", 50);
        let before = live.frame(View::Top);
        live.advance(10);
        assert_ne!(before.pixels, live.frame(View::Top).pixels);
    }
}
