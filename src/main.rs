mod bessel;
mod colormap;
mod field;
mod grid;
mod image;
mod oblique;
mod png;
mod rng;
mod topdown;
mod walker;

use std::env;

use field::{Field, Vec2};
use grid::Grid;
use rng::Rng;
use walker::{free, spring, step, Walkers};

const PIXELS: usize = 600;
const TRAIL: usize = 400;
const SPRING: f64 = 0.002;
const HISTOGRAM_SKIP: usize = 5000;

struct Scene {
    walkers: Walkers,
    half: f64,
    follow: bool,
    trap: bool,
}

fn scene(name: &str) -> Scene {
    match name {
        "single" => Scene {
            walkers: Walkers::new(vec![[0.0; 2]], vec![[0.0; 2]]),
            half: 3.0,
            follow: true,
            trap: false,
        },
        "pair" => Scene {
            walkers: Walkers::new(
                vec![[-0.35, 0.0], [0.35, 0.0]],
                vec![[0.0, 0.02], [0.0, -0.02]],
            ),
            half: 2.5,
            follow: true,
            trap: false,
        },
        "trap" => Scene {
            walkers: Walkers::new(vec![[0.3, 0.0]], vec![[0.0; 2]]),
            half: 3.0,
            follow: false,
            trap: true,
        },
        other => panic!("unknown scene {other}"),
    }
}

fn centroid(points: &[Vec2]) -> Vec2 {
    let n = points.len() as f64;
    let sum = points.iter().fold([0.0; 2], |s, p| [s[0] + p[0], s[1] + p[1]]);
    [sum[0] / n, sum[1] / n]
}

fn frame(field: &Field, scene: &Scene, trails: &[Vec<Vec2>], view: &str, path: &str) {
    let center = if scene.follow { centroid(&scene.walkers.positions) } else { [0.0; 2] };
    let grid = Grid::sample(field, center, scene.half, PIXELS);
    let image = match view {
        "top" => {
            let peak = grid.peak();
            let mut image = topdown::paint(&grid, &colormap::COOLWARM, -peak, peak);
            topdown::droplets(&mut image, &grid, trails);
            image
        }
        "oblique" => oblique::render(&grid, &scene.walkers.positions, PIXELS, PIXELS),
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

    let mut scene = scene(name);
    let mut field = Field::new(memory, scene.walkers.positions.len());
    let mut rng = Rng::new(1);
    let mut trails: Vec<Vec<Vec2>> = scene.walkers.positions.iter().map(|p| vec![*p]).collect();
    let mut histogram = Grid::zeros([0.0; 2], scene.half, PIXELS / 4);
    for bounce in 1..=bounces {
        if scene.trap {
            step(&mut field, &mut scene.walkers, spring(SPRING), &mut rng);
        } else {
            step(&mut field, &mut scene.walkers, free, &mut rng);
        }
        for (trail, p) in trails.iter_mut().zip(&scene.walkers.positions) {
            trail.push(*p);
            if trail.len() > TRAIL {
                trail.remove(0);
            }
        }
        if bounce > HISTOGRAM_SKIP {
            histogram.count(scene.walkers.positions[0]);
        }
        if bounce % every == 0 {
            frame(&field, &scene, &trails, view, &format!("{out}/{name}_{:05}.png", bounce / every));
        }
    }
    if scene.trap {
        let image = topdown::paint(&histogram, &colormap::VIRIDIS, 0.0, histogram.peak());
        png::write(format!("{out}/{name}_histogram.png"), &image).expect("write histogram");
    }
}
