//! The live viewer's running experiment and its pictures, free of UI.
use crate::experiment::{
    experiment, fraunhofer, Experiment, SINGLE_SLIT_WIDTH, SLIT_SEPARATION, SLIT_WIDTH,
};
use crate::grid::Vec2;
use crate::image::{Color, Image};
use crate::particles::{Kind, Particles};
use crate::rng::Rng;
use crate::screen::Screen;
use crate::wave::Wave;
use crate::{oblique, plot, view};

/// The field is drawn a pixel per cell up to this size.
const FIELD_PIXELS: usize = 512;
/// The oblique view samples the density coarser than it draws, with
/// one ray per pixel to stay fluid.
const OBLIQUE_CELLS: usize = 160;
const OBLIQUE_PIXELS: usize = 360;
const OBLIQUE_SAMPLES: usize = 1;
/// Least distance from the screen to the domain edge, clear of the
/// sponge, for grids too small to reach the experiment's screen.
const SCREEN_MARGIN: f64 = 5.0;
/// Particles closer than this to the domain's edge are dropped.
const EDGE: f64 = 1.0;
const BINS: usize = 120;
const SEED: u64 = 1;
/// The plot's curves: arrivals, integrated current and Fraunhofer.
const LEGEND: [(&str, Color); 3] = [
    ("arrivals", [255, 255, 255]),
    ("current", [80, 200, 255]),
    ("Fraunhofer", [255, 160, 30]),
];

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum View {
    Phase,
    Density,
    Oblique,
}

/// What a run starts from; changing any of it starts over.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Setup {
    pub name: &'static str,
    pub count: usize,
    pub kind: Kind,
    /// Grid cells along a side.
    pub n: usize,
}

pub struct Live {
    pub setup: Setup,
    pub steps: usize,
    /// The fraction of a step played but not yet run.
    owed: f64,
    experiment: Experiment,
    particles: Particles,
    screen: Option<Screen>,
}

/// Slit centres and width for the experiments that have slits.
fn slits(name: &str) -> Option<(Vec<f64>, f64)> {
    let d = SLIT_SEPARATION / 2.0;
    match name {
        "double-slit" => Some((vec![-d, d], SLIT_WIDTH)),
        "single-slit" => Some((vec![0.0], SINGLE_SLIT_WIDTH)),
        _ => None,
    }
}

/// Drop the particles that strayed within EDGE of the domain's edge.
/// The run goes on indefinitely, so particles do leave, and the
/// field has nothing to say about them out there.
fn crop(particles: &mut Particles, wave: &Wave) {
    let half = wave.extent() / 2.0 - EDGE;
    let inside = |p: &Vec2| p[0].abs() < half && p[1].abs() < half;
    let velocities = particles.velocities.iter().copied();
    let pairs = particles.positions.iter().copied().zip(velocities);
    (particles.positions, particles.velocities) = pairs.filter(|(p, _)| inside(p)).unzip();
}

impl Live {
    pub fn new(setup: Setup) -> Live {
        let e = experiment(setup.name, setup.n, setup.kind);
        let mut rng = Rng::new(SEED);
        let positions = (0..setup.count).map(|_| e.packet.sample(&mut rng)).collect();
        let mut particles = Particles::new(setup.kind, positions, e.packet.wavevector);
        crop(&mut particles, &e.wave);
        let half = e.wave.extent() / 2.0;
        let screen = e.screen.map(|y| Screen::new(y.min(half - SCREEN_MARGIN), half, BINS));
        Live { setup, steps: 0, owed: 0.0, experiment: e, particles, screen }
    }

    /// Start over if the setup differs from the running one.
    pub fn select(&mut self, setup: Setup) {
        if setup != self.setup {
            *self = Live::new(setup);
        }
    }

    pub fn advance(&mut self, steps: usize) {
        let e = &mut self.experiment;
        for _ in 0..steps {
            e.wave.step(e.dt);
            let before = self.particles.positions.clone();
            self.particles.step(&e.wave, e.dt);
            if let Some(screen) = &mut self.screen {
                screen.record(&before, &self.particles.positions);
                screen.integrate(&e.wave, e.dt);
            }
            crop(&mut self.particles, &e.wave);
        }
        self.steps += steps;
    }

    /// Play a frame of dt seconds at speed time units per second, at
    /// most limit steps, and return the steps run. A fraction of a
    /// step carries to the next frame; beyond the limit the replay
    /// falls behind rather than catching up.
    pub fn play(&mut self, speed: f64, dt: f64, limit: usize) -> usize {
        self.owed += speed * dt / self.experiment.dt;
        let steps = (self.owed as usize).min(limit);
        self.owed = (self.owed - steps as f64).fract();
        self.advance(steps);
        steps
    }

    pub fn time(&self) -> f64 {
        self.experiment.wave.time
    }

    /// Particles that crossed the screen, if there is one.
    pub fn arrivals(&self) -> Option<f64> {
        self.screen.as_ref().map(|s| s.hits.iter().sum())
    }

    pub fn frame(&self, view: View) -> Image {
        let wave = &self.experiment.wave;
        let style = match view {
            View::Phase => view::Style::Phase,
            View::Density => view::Style::Density,
            View::Oblique => {
                let grid = view::relief(wave, OBLIQUE_CELLS);
                let (w, h) = (OBLIQUE_PIXELS, OBLIQUE_PIXELS);
                return oblique::render(&grid, &[], w, h, OBLIQUE_SAMPLES);
            }
        };
        let peak = self.experiment.packet.peak();
        let mut image = view::field(wave, wave.n.min(FIELD_PIXELS), style, peak);
        let heads: Vec<Vec<Vec2>> =
            self.particles.positions.iter().map(|p| vec![*p]).collect();
        view::particles(&mut image, wave, &heads);
        image
    }

    /// Fraunhofer intensity at the screen's bins. The formula is for
    /// the experiment's screen, so x scales by the ratio of distances
    /// when the screen sits closer.
    fn theory(&self, screen: &Screen) -> Option<Vec<f64>> {
        let (centers, width) = slits(self.setup.name)?;
        let ratio = self.experiment.screen? / screen.y;
        let at = |b| fraunhofer(screen.x(b) * ratio, &centers, width);
        Some((0..screen.bins).map(at).collect())
    }

    /// Arrivals, integrated current and, for slits, the Fraunhofer
    /// curve, if there is a screen.
    pub fn plot(&self, width: usize, height: usize) -> Option<Image> {
        let screen = self.screen.as_ref()?;
        let theory = self.theory(screen);
        let color = |k: usize| LEGEND[k].1;
        let mut series: Vec<(&[f64], Color)> =
            vec![(&screen.hits, color(0)), (&screen.flux, color(1))];
        series.extend(theory.as_deref().map(|t| (t, color(2))));
        Some(plot::curves(width, height, &series))
    }

    /// Names and colours of the plotted curves.
    pub fn legend(&self) -> &'static [(&'static str, Color)] {
        match (&self.screen, slits(self.setup.name)) {
            (None, _) => &[],
            (Some(_), None) => &LEGEND[..2],
            (Some(_), Some(_)) => &LEGEND,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn setup(name: &'static str) -> Setup {
        Setup { name, count: 50, kind: Kind::Quantum, n: 128 }
    }

    fn advanced(name: &'static str, steps: usize) -> Live {
        let mut live = Live::new(setup(name));
        live.advance(steps);
        live
    }

    #[test]
    fn only_a_new_setup_starts_over() {
        let mut live = advanced("free", 10);
        live.select(setup("free"));
        assert_eq!(live.steps, 10);
        live.select(Setup { count: 51, ..setup("free") });
        assert_eq!(live.steps, 0);
        live.advance(5);
        live.select(Setup { kind: Kind::Classical, ..setup("free") });
        assert_eq!(live.steps, 0);
    }

    #[test]
    fn played_time_matches_speed_until_the_limit() {
        let mut live = Live::new(setup("free"));
        let dt = live.experiment.dt;
        let steps: usize = (0..64).map(|_| live.play(0.5, 1.0 / 64.0, 1000)).sum();
        assert!((steps as f64 * dt - 0.5).abs() <= dt, "{steps}");
        assert_eq!(live.play(100.0, 1.0, 7), 7);
        assert!(live.owed < 1.0);
    }

    #[test]
    fn slow_speeds_run_a_step_every_few_frames() {
        let mut live = Live::new(setup("free"));
        let speed = live.experiment.dt * 8.0;
        let played: Vec<usize> =
            (0..64).map(|_| live.play(speed, 1.0 / 64.0, 100)).collect();
        assert_eq!(live.steps, 8);
        assert!(played.iter().all(|&n| n <= 1));
    }

    #[test]
    fn classical_particles_bouncing_off_the_slits_leave_the_domain() {
        let classical = Setup { kind: Kind::Classical, n: 256, ..setup("double-slit") };
        let mut live = Live::new(classical);
        live.advance(2000);
        let half = live.experiment.wave.extent() / 2.0 - EDGE;
        let positions = &live.particles.positions;
        assert!(positions.len() < classical.count);
        assert!(positions.iter().all(|p| p[0].abs() < half && p[1].abs() < half));
    }

    #[test]
    fn small_grids_keep_the_screen_inside_the_domain() {
        let live = Live::new(setup("double-slit"));
        let half = live.experiment.wave.extent() / 2.0;
        let screen = live.screen.as_ref().unwrap();
        assert!(screen.y <= half - SCREEN_MARGIN);
        assert!(live.plot(64, 32).is_some());
        assert_eq!(live.legend().len(), 3);
        let free = Live::new(setup("free"));
        assert!(free.plot(64, 32).is_none() && free.legend().is_empty());
    }

    #[test]
    fn theory_peaks_mid_screen_and_is_symmetric() {
        let live = Live::new(setup("double-slit"));
        let screen = live.screen.as_ref().unwrap();
        let theory = live.theory(screen).unwrap();
        let peak = theory.iter().cloned().fold(0.0, f64::max);
        assert!(theory[BINS / 2] > 0.9 * peak);
        for b in 0..BINS {
            assert!((theory[b] - theory[BINS - 1 - b]).abs() < 1e-9);
        }
    }

    #[test]
    fn frames_change_as_the_field_evolves() {
        let mut live = advanced("free", 1);
        let before = live.frame(View::Phase);
        live.advance(20);
        assert_ne!(before.pixels, live.frame(View::Phase).pixels);
        assert_eq!(live.frame(View::Oblique).width, OBLIQUE_PIXELS);
    }
}
