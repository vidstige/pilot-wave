//! Named experiments: a field with walls and potential, a packet, the
//! kind of particles, a detection line and how long to run. No I/O.
//!
//! The packet's wavelength is 1, so k = 2π and the packet speed is 2π.
use std::f64::consts::{PI, TAU};

use crate::grid::Vec2;
use crate::packet::Packet;
use crate::particles::Kind;
use crate::wave::{stable_dt, Wave};

pub const NAMES: [&str; 8] = [
    "free", "double-slit", "single-slit", "barrier", "disc", "harmonic", "ground", "box",
];
pub const SPEED: f64 = TAU;
pub const ENERGY: f64 = SPEED * SPEED / 2.0;
const DX: f64 = 0.1;
/// Samples on each side of the centre in the quadratures below.
const QUADRATURE: i32 = 4000;
/// Where the packet starts and where arrivals are counted, pulled in
/// on small grids to stay this far from the edge.
const START: f64 = -20.0;
const SCREEN: f64 = 30.0;
const MARGIN: f64 = 6.0;
/// Extra distance run after the packet centre reaches the screen, so
/// that the slow and wide-angle tail of the wave crosses it too.
const TAIL: f64 = 25.0;
/// Walls at y = 0 with this thickness.
const WALL: f64 = 0.4;
pub const SLIT_WIDTH: f64 = 1.0;
pub const SLIT_SEPARATION: f64 = 4.0;
pub const SINGLE_SLIT_WIDTH: f64 = 2.0;
/// Barrier height relative to the packet's energy, and its width.
/// The width is a whole number of cells, so the grid holds exactly it.
pub const BARRIER_HEIGHT: f64 = 1.2;
pub const BARRIER_WIDTH: f64 = 0.4;
pub const DISC_RADIUS: f64 = 4.0;
/// Width of the beam around the disc.
pub const DISC_BEAM: f64 = 8.0;
/// Trap frequency and the orbit radius of the coherent state.
pub const TRAP: f64 = 0.5;
const ORBIT: f64 = 4.0;
const BOX_HALF: f64 = 6.0;

pub struct Experiment {
    pub name: &'static str,
    pub wave: Wave,
    pub packet: Packet,
    pub kind: Kind,
    pub dt: f64,
    pub duration: f64,
    /// y of the line where arrivals are counted, if any.
    pub screen: Option<f64>,
}

fn upward(center: Vec2, width: Vec2) -> Packet {
    Packet { center, width, wavevector: [0.0, SPEED] }
}

/// Zero the field wherever the predicate holds.
fn wall(wave: &mut Wave, inside: impl Fn(Vec2) -> bool) {
    for j in 0..wave.n {
        for i in 0..wave.n {
            if inside(wave.position((i, j))) {
                let k = wave.index((i, j));
                wave.keep[k] = 0.0;
            }
        }
    }
}

fn potential(wave: &mut Wave, v: impl Fn(Vec2) -> f64) {
    for j in 0..wave.n {
        for i in 0..wave.n {
            let k = wave.index((i, j));
            wave.potential[k] = v(wave.position((i, j))) as f32;
        }
    }
}

fn in_wall(p: Vec2) -> bool {
    p[1].abs() < WALL / 2.0
}

fn slits(centers: &[f64], width: f64) -> impl Fn(Vec2) -> bool + '_ {
    move |p| in_wall(p) && !centers.iter().any(|c| (p[0] - c).abs() < width / 2.0)
}

pub fn experiment(name: &'static str, n: usize, kind: Kind) -> Experiment {
    let mut wave = Wave::new(n, DX);
    let half = wave.extent() / 2.0;
    let start = START.max(-half + MARGIN);
    let screen_y = SCREEN.min(half - MARGIN);
    let (packet, potential_max, screen) = match name {
        "free" => (upward([0.0, start], [1.5, 1.5]), 0.0, None),
        "double-slit" => {
            let d = SLIT_SEPARATION / 2.0;
            wall(&mut wave, slits(&[-d, d], SLIT_WIDTH));
            (upward([0.0, start], [6.0, 2.0]), 0.0, Some(screen_y))
        }
        "single-slit" => {
            wall(&mut wave, slits(&[0.0], SINGLE_SLIT_WIDTH));
            (upward([0.0, start], [6.0, 2.0]), 0.0, Some(screen_y))
        }
        "barrier" => {
            let height = BARRIER_HEIGHT * ENERGY;
            potential(&mut wave, |p| if p[1].abs() < BARRIER_WIDTH / 2.0 { height } else { 0.0 });
            (upward([0.0, start], [4.0, 2.0]), height, Some(screen_y))
        }
        "disc" => {
            wall(&mut wave, |p| p[0].hypot(p[1]) < DISC_RADIUS);
            (upward([0.0, start], [DISC_BEAM, 2.0]), 0.0, Some(screen_y))
        }
        "harmonic" => {
            potential(&mut wave, |p| 0.5 * TRAP * TRAP * (p[0] * p[0] + p[1] * p[1]));
            let width = 1.0 / (2.0 * TRAP).sqrt();
            let packet = Packet { center: [ORBIT, 0.0], width: [width; 2], wavevector: [0.0, TRAP * ORBIT] };
            (packet, 0.5 * TRAP * TRAP * (wave.extent() / 2.0).powi(2), None)
        }
        "ground" => {
            potential(&mut wave, |p| 0.5 * TRAP * TRAP * (p[0] * p[0] + p[1] * p[1]));
            let width = 1.0 / (2.0 * TRAP).sqrt();
            let packet = Packet { center: [0.0; 2], width: [width; 2], wavevector: [0.0; 2] };
            (packet, 0.5 * TRAP * TRAP * (wave.extent() / 2.0).powi(2), None)
        }
        "box" => {
            wall(&mut wave, |p| p[0].abs() > BOX_HALF || p[1].abs() > BOX_HALF);
            (Packet { center: [-3.0, 0.0], width: [1.0, 1.0], wavevector: [SPEED / 2.0, SPEED / 4.0] }, 0.0, None)
        }
        other => panic!("unknown experiment {other}"),
    };
    let dt = stable_dt(DX, potential_max);
    if !matches!(name, "box") {
        wave.absorb_edges(dt);
    }
    packet.write(&mut wave);
    let duration = match name {
        "harmonic" | "ground" => 2.0 * TAU / TRAP,
        "box" | "free" => 6.0,
        _ => (screen_y - start + TAIL) / SPEED,
    };
    Experiment { name, wave, packet, kind, dt, duration, screen }
}

/// Fraunhofer intensity for slits of the given width at the given
/// centres, per unit of x on a screen at the given y. The angular
/// intensity of the cylindrical wave falls as 1/r and meets the flat
/// screen obliquely, which together give a factor cos²θ.
pub fn fraunhofer(x: f64, y: f64, centers: &[f64], width: f64) -> f64 {
    let r = x.hypot(y - WALL / 2.0);
    let (sin, cos) = (x / r, (y - WALL / 2.0) / r);
    let beta = PI * width * sin;
    let envelope = if beta.abs() < 1e-9 { 1.0 } else { (beta.sin() / beta).powi(2) };
    let (re, im) = centers.iter().fold((0.0, 0.0), |(re, im), c| {
        let phase = TAU * c * sin;
        (re + phase.cos(), im + phase.sin())
    });
    let array = (re * re + im * im) / (centers.len() * centers.len()) as f64;
    cos * cos * envelope * array
}

/// Transmission through a rectangular barrier for a plane wave of
/// wavenumber k: the textbook formula, with sinh below the top of the
/// barrier and sin above it.
pub fn transmission(k: f64, height: f64, width: f64) -> f64 {
    let e = k * k / 2.0;
    let q = (2.0 * (height - e)).abs().sqrt() * width;
    let s = if e < height { q.sinh() } else { q.sin() };
    1.0 / (1.0 + height * height * s * s / (4.0 * e * (height - e).abs()))
}

/// The barrier's transmission averaged over a Gaussian spread of
/// wavenumbers with standard deviation `spread`, which is 1/(2σ) for
/// a packet whose density has standard deviation σ.
pub fn mean_transmission(k: f64, spread: f64, height: f64, width: f64) -> f64 {
    let (sum, total) = (-QUADRATURE..=QUADRATURE)
        .map(|i| 6.0 * i as f64 / QUADRATURE as f64)
        .map(|u| ((-u * u / 2.0).exp(), transmission(k + spread * u, height, width)))
        .fold((0.0, 0.0), |(s, t), (w, tr)| (s + w * tr, t + w));
    sum / total
}

/// Fresnel diffraction of the beam around the disc, as intensity
/// against x on the screen: the beam's Gaussian profile, blocked
/// across the disc's diameter at y = 0, carried paraxially to the
/// screen. In two dimensions the disc is a strip in cross-section, so
/// its bright spot is a line on the axis.
pub fn shadow(x: f64, y: f64) -> f64 {
    let ds = 5.0 * DISC_BEAM / QUADRATURE as f64;
    let (re, im) = (-QUADRATURE..=QUADRATURE)
        .map(|i| i as f64 * ds)
        .filter(|s| s.abs() > DISC_RADIUS)
        .fold((0.0, 0.0), |(re, im), s| {
            let a = (-s * s / (4.0 * DISC_BEAM * DISC_BEAM)).exp();
            let phase = PI * (x - s) * (x - s) / y;
            (re + a * phase.cos(), im + a * phase.sin())
        });
    re * re + im * im
}

/// The expected shape of the arrivals at x on a screen at y, or zero
/// where there is no closed form.
pub fn pattern(name: &str, x: f64, y: f64) -> f64 {
    let d = SLIT_SEPARATION / 2.0;
    match name {
        "double-slit" => fraunhofer(x, y, &[-d, d], SLIT_WIDTH),
        "single-slit" => fraunhofer(x, y, &[0.0], SINGLE_SLIT_WIDTH),
        "disc" => shadow(x, y),
        _ => 0.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const HEIGHT: f64 = BARRIER_HEIGHT * ENERGY;

    #[test]
    fn transmission_is_whole_without_a_barrier() {
        assert!((transmission(SPEED, 0.0, BARRIER_WIDTH) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn transmission_is_continuous_over_the_barrier_top() {
        let top = (2.0 * HEIGHT).sqrt();
        let limit = 1.0 / (1.0 + HEIGHT * BARRIER_WIDTH * BARRIER_WIDTH / 2.0);
        for k in [top - 1e-6, top + 1e-6] {
            assert!((transmission(k, HEIGHT, BARRIER_WIDTH) - limit).abs() < 1e-4);
        }
    }

    #[test]
    fn narrow_spread_averages_to_the_plane_wave() {
        let plane = transmission(SPEED, HEIGHT, BARRIER_WIDTH);
        let mean = mean_transmission(SPEED, 1e-4, HEIGHT, BARRIER_WIDTH);
        assert!((mean - plane).abs() < 1e-6);
    }

    #[test]
    fn disc_shadow_is_symmetric_with_a_bright_axis() {
        assert!((shadow(1.3, SCREEN) - shadow(-1.3, SCREEN)).abs() < 1e-9 * shadow(1.3, SCREEN));
        let beside = (1..=10).map(|i| shadow(0.2 * i as f64, SCREEN)).fold(f64::MAX, f64::min);
        assert!(shadow(0.0, SCREEN) > 2.0 * beside);
    }
}
