//! Named starting setups and their bounce by bounce run. No I/O.
use crate::bath::{Bath, Vec2};
use crate::corral::Corral;
use crate::open::Open;
use crate::rng::Rng;
use crate::walker::{confine, spring, step, Walkers};

/// Bounces of history kept in each droplet's trail.
const TRAIL: usize = 400;
const SPRING: f64 = 0.002;
const CORRAL_RADIUS: f64 = 2.0;

pub struct Scene {
    pub bath: Box<dyn Bath>,
    pub walkers: Walkers,
    /// Half width of the square the views show.
    pub half: f64,
    /// Whether views follow the droplets' centroid.
    pub follow: bool,
    pub spring: f64,
    pub rim: f64,
}

pub fn scene(name: &str, memory: f64) -> Scene {
    let free = |walkers: Walkers, half: f64| Scene {
        bath: Box::new(Open::new(memory, walkers.positions.len())),
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

/// A scene in motion, with the recent trail of every droplet.
pub struct Simulation {
    pub scene: Scene,
    pub trails: Vec<Vec<Vec2>>,
    rng: Rng,
}

impl Simulation {
    pub fn new(scene: Scene) -> Simulation {
        let positions = &scene.walkers.positions;
        let trails = positions.iter().map(|p| vec![*p]).collect();
        Simulation { scene, trails, rng: Rng::new(1) }
    }

    pub fn bounce(&mut self) {
        let s = &mut self.scene;
        let force = spring(s.spring);
        step(&mut *s.bath, &mut s.walkers, force, &mut self.rng);
        confine(&mut s.walkers, s.rim);
        for (trail, p) in self.trails.iter_mut().zip(&s.walkers.positions) {
            trail.push(*p);
            if trail.len() > TRAIL {
                trail.remove(0);
            }
        }
    }

    /// Where the views look.
    pub fn center(&self) -> Vec2 {
        let s = &self.scene;
        if s.follow { centroid(&s.walkers.positions) } else { [0.0; 2] }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(name: &str, bounces: usize) -> Simulation {
        let mut sim = Simulation::new(scene(name, 20.0));
        for _ in 0..bounces {
            sim.bounce();
        }
        sim
    }

    #[test]
    fn trails_end_at_droplets_and_keep_only_recent_bounces() {
        let sim = run("pair", TRAIL + 50);
        let positions = &sim.scene.walkers.positions;
        for (trail, p) in sim.trails.iter().zip(positions) {
            assert_eq!(trail.len(), TRAIL);
            assert_eq!(trail.last(), Some(p));
        }
    }

    #[test]
    fn corral_droplet_never_leaves_the_rim() {
        let sim = run("corral", 2000);
        let r = |p: &Vec2| (p[0] * p[0] + p[1] * p[1]).sqrt();
        assert!(sim.trails[0].iter().all(|p| r(p) <= CORRAL_RADIUS));
    }
}
