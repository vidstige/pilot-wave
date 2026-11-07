use std::env;

use pilot_wave::grid::Grid;
use pilot_wave::scene::{scene, Simulation};
use pilot_wave::{colormap, frame, png, topdown};

const PIXELS: usize = 600;
/// Rays per pixel in the oblique view.
const STILL_SAMPLES: usize = 4;
const HISTOGRAM_BINS: usize = 100;
const HISTOGRAM_SKIP: usize = 5000;

fn write_frame(sim: &Simulation, view: &str, path: &str) {
    let image = match view {
        "top" => frame::top(sim, PIXELS),
        "oblique" => frame::oblique(sim, PIXELS, PIXELS, STILL_SAMPLES),
        other => panic!("unknown view {other}"),
    };
    png::write(path, &image).expect("write frame");
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let name = args.get(1).map(String::as_str).unwrap_or("single");
    let memory: f64 = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(20.0);
    let bounces: usize = args.get(3).map(|s| s.parse().unwrap()).unwrap_or(1200);
    let every: usize = args.get(4).map(|s| s.parse().unwrap()).unwrap_or(bounces);
    let out = args.get(5).map(String::as_str).unwrap_or("out");
    let view = args.get(6).map(String::as_str).unwrap_or("top");
    std::fs::create_dir_all(out).expect("create output directory");

    let mut sim = Simulation::new(scene(name, memory));
    let mut histogram = Grid::zeros([0.0; 2], sim.scene.half, HISTOGRAM_BINS);
    for bounce in 1..=bounces {
        sim.bounce();
        if bounce > HISTOGRAM_SKIP {
            histogram.count(sim.scene.walkers.positions[0]);
        }
        if bounce % every == 0 {
            write_frame(&sim, view, &format!("{out}/{name}_{:05}.png", bounce / every));
        }
    }
    let scale = PIXELS / HISTOGRAM_BINS;
    let image = topdown::paint(&histogram, &colormap::VIRIDIS, 0.0, histogram.peak(), scale);
    png::write(format!("{out}/{name}_histogram.png"), &image).expect("write histogram");
}
