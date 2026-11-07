//! A simulation's current state as an image, in either view.
use crate::colormap::COOLWARM;
use crate::grid::Grid;
use crate::image::Image;
use crate::scene::Simulation;
use crate::topdown;

fn grid(sim: &Simulation, n: usize) -> Grid {
    let s = &sim.scene;
    Grid::sample(&*s.bath, sim.center(), s.half, n)
}

/// Top view of the bath sampled on n by n cells, a pixel per cell.
pub fn top(sim: &Simulation, n: usize) -> Image {
    let grid = grid(sim, n);
    let peak = grid.peak();
    let mut image = topdown::paint(&grid, &COOLWARM, -peak, peak, 1);
    topdown::droplets(&mut image, &grid, &sim.trails);
    image
}

/// Oblique view of the bath sampled on n by n cells.
pub fn oblique(sim: &Simulation, n: usize, pixels: usize) -> Image {
    let droplets = &sim.scene.walkers.positions;
    crate::oblique::render(&grid(sim, n), droplets, pixels, pixels)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scene::scene;
    use crate::topdown::DROPLET_COLORS;

    fn walked(name: &str, bounces: usize) -> Simulation {
        let mut sim = Simulation::new(scene(name, 20.0));
        for _ in 0..bounces {
            sim.bounce();
        }
        sim
    }

    #[test]
    fn top_view_keeps_a_followed_walker_in_the_middle() {
        let image = top(&walked("single", 300), 60);
        let middle = image.pixels[30 * 60 + 30];
        assert_eq!(middle, DROPLET_COLORS[0]);
    }

    #[test]
    fn oblique_size_follows_pixels_not_samples() {
        let image = oblique(&walked("pair", 50), 20, 48);
        assert_eq!((image.width, image.height), (48, 48));
    }
}
