use std::env;

use pilot_wave::bath::{Bath, Vec2};
use pilot_wave::corral::Corral;
use pilot_wave::field::Field;
use pilot_wave::grid::Grid;
use pilot_wave::rng::Rng;
use pilot_wave::walker::{confine, spring, step, Walkers};
use pilot_wave::{colormap, oblique, png, topdown};

const PIXELS: usize = 600;
const TRAIL: usize = 400;
const SPRING: f64 = 0.002;
const CORRAL_RADIUS: f64 = 2.0;
const HISTOGRAM_BINS: usize = 100;
const HISTOGRAM_SKIP: usize = 5000;

struct Scene {
    bath: Box<dyn Bath>,
    walkers: Walkers,
    half: f64,
    follow: bool,
    spring: f64,
    rim: f64,
}

fn scene(name: &str, memory: f64) -> Scene {
    let free = |walkers: Walkers, half: f64| Scene {
        bath: Box::new(Field::new(memory, walkers.positions.len())),
        walkers,
        half,
        follow: true,
        spring: 0.0,
        rim: f64::INFINITY,
    };
    match name {
        "single" => free(Walkers::new(vec![[0.0; 2]], vec![[0.0; 2]]), 3.0),
        "pair" => free(
            Walkers::new(
                vec![[-0.35, 0.0], [0.35, 0.0]],
                vec![[0.0, 0.02], [0.0, -0.02]],
            ),
            2.5,
        ),
        "trap" => Scene {
            follow: false,
            spring: SPRING,
            ..free(Walkers::new(vec![[0.3, 0.0]], vec![[0.0; 2]]), 3.0)
        },
        "corral" => Scene {
            bath: Box::new(Corral::new(CORRAL_RADIUS, memory)),
            walkers: Walkers::new(vec![[0.3, 0.1]], vec![[0.0; 2]]),
            half: CORRAL_RADIUS,
            follow: false,
            spring: 0.0,
            rim: CORRAL_RADIUS,
        },
        other => panic!("unknown scene {other}"),
    }
}

fn centroid(points: &[Vec2]) -> Vec2 {
    let n = points.len() as f64;
    let sum = points.iter().fold([0.0; 2], |s, p| [s[0] + p[0], s[1] + p[1]]);
    [sum[0] / n, sum[1] / n]
}

fn frame(scene: &Scene, trails: &[Vec<Vec2>], view: &str, path: &str) {
    let center = if scene.follow { centroid(&scene.walkers.positions) } else { [0.0; 2] };
    let grid = Grid::sample(&*scene.bath, center, scene.half, PIXELS);
    let image = match view {
        "top" => {
            let peak = grid.peak();
            let mut image = topdown::paint(&grid, &colormap::COOLWARM, -peak, peak, 1);
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

    let mut scene = scene(name, memory);
    let mut rng = Rng::new(1);
    let mut trails: Vec<Vec<Vec2>> = scene.walkers.positions.iter().map(|p| vec![*p]).collect();
    let mut histogram = Grid::zeros([0.0; 2], scene.half, HISTOGRAM_BINS);
    for bounce in 1..=bounces {
        step(&mut *scene.bath, &mut scene.walkers, spring(scene.spring), &mut rng);
        confine(&mut scene.walkers, scene.rim);
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
            frame(&scene, &trails, view, &format!("{out}/{name}_{:05}.png", bounce / every));
        }
    }
    let scale = PIXELS / HISTOGRAM_BINS;
    let image = topdown::paint(&histogram, &colormap::VIRIDIS, 0.0, histogram.peak(), scale);
    png::write(format!("{out}/{name}_histogram.png"), &image).expect("write histogram");
}
