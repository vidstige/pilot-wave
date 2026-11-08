//! Named experiments: a field with walls and potential, a packet, the
//! kind of particles, a detection line and how long to run. No I/O.
//!
//! The packet's wavelength is 1, so k = 2π and the packet speed is 2π.
//! Two particles on a line are one point (x₁, x₂) in the plane.
use std::f64::consts::{PI, TAU};

use crate::grid::Vec2;
use crate::packet::Packet;
use crate::pair::{self, Orbital, Symmetry};
use crate::particles::Kind;
use crate::wave::{stable_dt, Wave};

pub const NAMES: [&str; 22] = [
    "free", "double-slit", "single-slit", "barrier", "disc", "harmonic", "ground", "vortex",
    "antivortex", "orbital", "box", "flux-trap", "flux-slit", "fermions", "bosons",
    "distinguishable", "molecule", "molecule-break", "entangled", "unentangled",
    "entangled-wall", "unentangled-wall",
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
/// Two particles on a line |x| ≤ 12.8, whatever the grid asked for.
const LINE_CELLS: usize = 256;
/// The colliding packets: where they start, their width and speed.
const PAIR_START: f64 = 4.0;
const PAIR_WIDTH: f64 = 1.0;
const PAIR_SPEED: f64 = 2.0;
const COLLISION_TIME: f64 = 5.0;
/// The molecule's Gaussian well U(x₁ − x₂), which binds exactly one
/// state of the relative motion, 2.8 deep; a second one would need
/// a well half again as deep.
pub const BOND_DEPTH: f64 = 6.0;
pub const BOND_RANGE: f64 = 0.4;
/// The hard wall on the line, which only particle 1 feels: a wall
/// both feel is a right-angled corner in configuration space, which
/// reflects the pair like a corner reflector, mirroring the relative
/// motion so that the bond survives any kick. Then where the molecule
/// starts, the widths of its centre and of its separation before
/// relaxing, and how long it relaxes in imaginary time.
pub const LINE_WALL: f64 = 5.0;
const MOLECULE_START: f64 = -2.0;
const MOLECULE_WIDTH: [f64; 2] = [1.0, 0.5];
const RELAX_TIME: f64 = 4.0;
/// Total wavenumbers of the gentle and the hard kick, and how long
/// each runs. Particle 1 bouncing alone gives the relative motion
/// the pair's kinetic energy k²/4, which breaks the bond once it
/// passes the binding energy.
pub const KICK: f64 = 2.0;
pub const BREAK_KICK: f64 = 5.0;
const MOLECULE_TIME: f64 = 15.0;
const BREAK_TIME: f64 = 6.0;
/// The entangled pair on a line twice as long, |x| ≤ 25.6, so that
/// it flies well apart, in cells half as wide: the widths of |ψ|² in
/// the pair's centre and in its separation, the wavenumber carrying
/// each particle away from the other, where the wall only particle 2
/// feels begins, the force with which it pushes back, and how long it
/// runs.
///
/// The wall is a steep ramp, which turns particle 2 back within half
/// a unit. Cells held at zero let a few runs step into them, as the
/// interpolated velocity stays finite up to such cells and is zero
/// inside, freezing particle 1 with them; a step a hundred high made
/// ψ fall off within a cell. Where particle 2 turns, the wave coming
/// back interferes with the wave still coming in in fringes 0.5 to 1
/// apart, and on cells of 0.1 the runs that crossed them bunched in x₁
/// where the field did not (χ² of 112 over 56 bins against the
/// field's marginal at the end); on cells of 0.05 they do not (55
/// over 57) and each run's x₁ changes as before.
const EPR_CELLS: usize = 1024;
const EPR_DX: f64 = 0.05;
const EPR_WIDTHS: [f64; 2] = [1.0, 0.6];
const EPR_KICK: f64 = 3.0;
const EPR_WALL: f64 = 5.0;
const EPR_FORCE: f64 = 20.0;
const EPR_TIME: f64 = 3.5;

/// Where the experiment's particles live: one particle in the plane,
/// or two on a line, the grid then being their configuration space.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Space {
    Plane,
    Line,
}
/// Radius of the hard core shielding the flux line in the trap, the
/// imaginary time its state relaxes for, and the grid it relaxes on:
/// the state is below 1e-11 beyond that grid, and relaxing on a 1024
/// grid took 34 s against 1.2 s.
const FLUX_CORE: f64 = 0.2;
const RELAX: f64 = 10.0;
const RELAX_GRID: usize = 256;

pub struct Experiment {
    pub name: &'static str,
    pub wave: Wave,
    pub packet: Packet,
    pub kind: Kind,
    pub dt: f64,
    pub duration: f64,
    /// y of the line where arrivals are counted, if any.
    pub screen: Option<f64>,
    /// The starting density's peak, the reference for brightness.
    pub peak: f32,
    /// For a line, `packet` holds the two particles' packets as one
    /// in configuration space.
    pub space: Space,
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

/// Write ψ from a function of position, then normalise.
fn write(wave: &mut Wave, psi: impl Fn(Vec2) -> (f64, f64)) {
    for j in 0..wave.n {
        for i in 0..wave.n {
            let (re, im) = psi(wave.position((i, j)));
            let k = wave.index((i, j));
            wave.re[k] = re as f32;
            wave.im[k] = im as f32;
        }
    }
    let scale = (1.0 / wave.norm()).sqrt() as f32;
    for k in 0..wave.n * wave.n {
        wave.re[k] *= scale;
        wave.im[k] *= scale;
    }
}

/// The trap's n = 1 states: angular momentum m = ±1 circulates, and
/// their real combination is the p orbital along x, which stands still.
fn trap_state(wave: &mut Wave, m: i32) {
    write(wave, |p| {
        let envelope = (-0.5 * TRAP * (p[0] * p[0] + p[1] * p[1])).exp();
        match m {
            0 => (p[0] * envelope, 0.0),
            _ => (p[0] * envelope, m as f64 * p[1] * envelope),
        }
    });
}

/// ψ of a smaller grid of the same cells copied into the middle of
/// a larger one, zero around it.
fn embed(small: &Wave, wave: &mut Wave) {
    let offset = (wave.n - small.n) / 2;
    wave.re.fill(0.0);
    wave.im.fill(0.0);
    for j in 0..small.n {
        for i in 0..small.n {
            let (k, s) = (wave.index((i + offset, j + offset)), small.index((i, j)));
            (wave.re[k], wave.im[k]) = (small.re[s], small.im[s]);
        }
    }
}

/// Half a flux quantum through a shielded core at the trap's centre,
/// and the lowest state of angular momentum 1/2 there: the ground
/// Gaussian times e^{iθ/2}, θ counterclockwise from the cut, relaxed.
fn flux_trap_state(wave: &mut Wave, dt: f64) {
    wall(wave, |p| p[0].hypot(p[1]) < FLUX_CORE);
    wave.flux_line([0.0; 2]);
    if wave.n > RELAX_GRID {
        embed(&experiment("flux-trap", RELAX_GRID, Kind::Quantum).wave, wave);
    } else {
        write(wave, |p| {
            let envelope = (-0.5 * TRAP * (p[0] * p[0] + p[1] * p[1])).exp();
            let half = (-p[0]).atan2(p[1]).rem_euclid(TAU) / 2.0;
            (envelope * half.cos(), envelope * half.sin())
        });
        for _ in 0..(RELAX / dt) as usize {
            wave.relax(dt);
        }
    }
}

fn in_wall(p: Vec2) -> bool {
    p[1].abs() < WALL / 2.0
}

fn slits(centers: &[f64], width: f64) -> impl Fn(Vec2) -> bool + '_ {
    move |p| in_wall(p) && !centers.iter().any(|c| (p[0] - c).abs() < width / 2.0)
}

/// Two packets heading for each other.
fn collision(wave: &mut Wave, symmetry: Symmetry) -> Packet {
    let (x, k) = (PAIR_START, PAIR_SPEED);
    let a = Orbital { center: -x, width: PAIR_WIDTH, wavenumber: k };
    let b = Orbital { center: x, width: PAIR_WIDTH, wavenumber: -k };
    pair::write(wave, a, b, symmetry);
    Packet { center: [-x, x], width: [PAIR_WIDTH; 2], wavevector: [k, -k] }
}

pub fn bond(r: f64) -> f64 {
    -BOND_DEPTH * (-r * r / (2.0 * BOND_RANGE * BOND_RANGE)).exp()
}

/// Standard deviation of x₁.
fn marginal_width(wave: &Wave) -> f64 {
    let [m, _] = pair::marginals(wave);
    let x = |i: usize| wave.position((i, 0))[0];
    let moment = |p| wave.dx * (0..wave.n).map(|i| m[i] * x(i).powi(p)).sum::<f64>();
    (moment(2) - moment(1).powi(2)).sqrt()
}

/// The bound pair in the ground state of its relative motion, a
/// Gaussian in its centre, kicked towards particle 1's wall with
/// total wavenumber k. The packet is either particle's, centred on the pair.
fn molecule(wave: &mut Wave, dt: f64, k: f64) -> Packet {
    pair::interaction(wave, bond);
    wall(wave, |p| p[0] > LINE_WALL);
    let [w, s] = MOLECULE_WIDTH;
    pair::fill(wave, |x1, x2| {
        let (c, r) = ((x1 + x2) / 2.0 - MOLECULE_START, x1 - x2);
        [(-c * c / (4.0 * w * w) - r * r / (4.0 * s * s)).exp(), 0.0]
    });
    pair::relax(wave, dt, (RELAX_TIME / dt) as usize);
    pair::kick(wave, k);
    let width = marginal_width(wave);
    let center = [MOLECULE_START; 2];
    Packet { center, width: [width; 2], wavevector: [k / 2.0; 2] }
}

/// Two particles on a line, any wall set after the sponge so the
/// sponge leaves it whole.
fn line(name: &'static str, kind: Kind) -> Experiment {
    let mut wave = Wave::new(LINE_CELLS, DX);
    let dt = stable_dt(DX, BOND_DEPTH);
    wave.absorb_edges(dt);
    let (packet, duration) = match name {
        "fermions" => (collision(&mut wave, Symmetry::Fermions), COLLISION_TIME),
        "bosons" => (collision(&mut wave, Symmetry::Bosons), COLLISION_TIME),
        "distinguishable" => {
            (collision(&mut wave, Symmetry::Distinguishable), COLLISION_TIME)
        }
        "molecule" => (molecule(&mut wave, dt, KICK), MOLECULE_TIME),
        _ => (molecule(&mut wave, dt, BREAK_KICK), BREAK_TIME),
    };
    let space = Space::Line;
    let peak = wave.peak();
    Experiment { name, wave, packet, kind, dt, duration, screen: None, peak, space }
}

/// The entangled pair: together, far wider in its centre
/// c = (x₁ + x₂)/2 than in its separation r = x₁ − x₂, and flying
/// apart, particle 1 to the left: ψ ∝ e^{−c²/4σ_c² − r²/4σ_r² − ikr}.
/// Unentangled, it is the product of the same two marginal packets.
fn epr(wave: &mut Wave, entangled: bool) {
    let ([c, r], k) = (EPR_WIDTHS, EPR_KICK);
    if entangled {
        pair::fill(wave, |x1, x2| {
            let (u, v) = ((x1 + x2) / 2.0, x1 - x2);
            let a = (-u * u / (4.0 * c * c) - v * v / (4.0 * r * r)).exp();
            [a * (k * v).cos(), -a * (k * v).sin()]
        });
    } else {
        let width = (c * c + r * r / 4.0).sqrt();
        let one = Orbital { center: 0.0, width, wavenumber: -k };
        let two = Orbital { wavenumber: k, ..one };
        pair::write(wave, one, two, Symmetry::Distinguishable);
    }
}

/// The pair flying apart, entangled or not, with or without a wall
/// far out that only particle 2 feels.
fn flight(name: &'static str, kind: Kind) -> Experiment {
    let mut wave = Wave::new(EPR_CELLS, EPR_DX);
    // The same step with the wall and without, for runs side by side.
    let dt = stable_dt(EPR_DX, EPR_FORCE * (wave.extent() / 2.0 - EPR_WALL));
    wave.absorb_edges(dt);
    if name.ends_with("-wall") {
        potential(&mut wave, |p| EPR_FORCE * (p[1] - EPR_WALL).max(0.0));
    }
    epr(&mut wave, name.starts_with("entangled"));
    let width = marginal_width(&wave);
    let wavevector = [-EPR_KICK, EPR_KICK];
    let packet = Packet { center: [0.0; 2], width: [width; 2], wavevector };
    let (duration, peak, space) = (EPR_TIME, wave.peak(), Space::Line);
    Experiment { name, wave, packet, kind, dt, duration, screen: None, peak, space }
}

/// The same experiment without the wall on particle 2, if it has one.
pub fn twin(name: &'static str) -> Option<&'static str> {
    name.strip_suffix("-wall")
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
        "ground" | "vortex" | "antivortex" | "orbital" | "flux-trap" => {
            potential(&mut wave, |p| 0.5 * TRAP * TRAP * (p[0] * p[0] + p[1] * p[1]));
            let width = 1.0 / (2.0 * TRAP).sqrt();
            let packet = Packet { center: [0.0; 2], width: [width; 2], wavevector: [0.0; 2] };
            (packet, 0.5 * TRAP * TRAP * (wave.extent() / 2.0).powi(2), None)
        }
        "box" => {
            wall(&mut wave, |p| p[0].abs() > BOX_HALF || p[1].abs() > BOX_HALF);
            (Packet { center: [-3.0, 0.0], width: [1.0, 1.0], wavevector: [SPEED / 2.0, SPEED / 4.0] }, 0.0, None)
        }
        "fermions" | "bosons" | "distinguishable" | "molecule" | "molecule-break" => {
            return line(name, kind)
        }
        "flux-slit" => {
            let d = SLIT_SEPARATION / 2.0;
            wall(&mut wave, slits(&[-d, d], SLIT_WIDTH));
            wave.flux_line([0.0; 2]);
            (upward([0.0, start], [6.0, 2.0]), 0.0, Some(screen_y))
        }
        "entangled" | "unentangled" | "entangled-wall" | "unentangled-wall" => {
            return flight(name, kind)
        }
        other => panic!("unknown experiment {other}"),
    };
    let dt = stable_dt(DX, potential_max);
    if !matches!(name, "box") {
        wave.absorb_edges(dt);
    }
    packet.write(&mut wave);
    match name {
        "vortex" => trap_state(&mut wave, 1),
        "antivortex" => trap_state(&mut wave, -1),
        "orbital" => trap_state(&mut wave, 0),
        "flux-trap" => flux_trap_state(&mut wave, dt),
        _ => {}
    }
    let duration = match name {
        "harmonic" | "ground" | "vortex" | "antivortex" | "orbital" => 2.0 * TAU / TRAP,
        "flux-trap" => 2.0 * TAU / TRAP,
        "box" | "free" => 6.0,
        _ => (screen_y - start + TAIL) / SPEED,
    };
    let peak = wave.peak();
    Experiment { name, wave, packet, kind, dt, duration, screen, peak, space: Space::Plane }
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
        // Half a flux quantum between the slits makes their waves meet
        // with opposite signs, and |a - b|² = 2|a|² + 2|b|² - |a + b|²:
        // the double slit's fringes turned over in the same envelope.
        "flux-slit" => {
            fraunhofer(x, y, &[0.0], SLIT_WIDTH) - fraunhofer(x, y, &[-d, d], SLIT_WIDTH)
        }
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

    #[test]
    fn molecule_starts_bound_with_the_kick_on_top() {
        let e = experiment("molecule", 0, Kind::Quantum);
        // The relative motion's one bound state is 2.82 deep.
        let relative = pair::energy(&e.wave) - KICK * KICK / 4.0;
        assert!((-2.9..-2.6).contains(&relative), "{relative}");
    }

    #[test]
    fn flux_slit_fringes_sit_half_an_order_from_the_double_slits() {
        let x = |order: f64| (SCREEN - WALL / 2.0) * (order / SLIT_SEPARATION).asin().tan();
        let flux = |order| pattern("flux-slit", x(order), SCREEN);
        let plain = |order| pattern("double-slit", x(order), SCREEN);
        assert!(flux(0.0) < 1e-12 && plain(0.5) < 1e-12);
        assert!(flux(0.5) > 0.9 * plain(0.0));
    }

    fn largest_gap(a: &[f64], b: &[f64]) -> f64 {
        let peak = a.iter().cloned().fold(0.0, f64::max);
        a.iter().zip(b).map(|(p, q)| (p - q).abs()).fold(0.0, f64::max) / peak
    }

    #[test]
    fn entangled_pair_is_correlated_and_its_control_has_its_marginals() {
        let [e, u] = ["entangled", "unentangled"].map(|n| experiment(n, 0, Kind::Quantum));
        let [c, r] = EPR_WIDTHS;
        let rho = (c * c - r * r / 4.0) / (c * c + r * r / 4.0);
        let correlation = |w: &Wave| crate::entangle::correlation(&crate::entangle::field(w));
        assert!((correlation(&e.wave) - rho).abs() < 1e-3);
        assert!(correlation(&u.wave).abs() < 1e-6);
        let [m, n] = [&e.wave, &u.wave].map(pair::marginals);
        assert!(largest_gap(&m[0], &n[0]) < 1e-3 && largest_gap(&m[1], &n[1]) < 1e-3);
    }

    /// Densities of x₁ and x₂ from the scheme's conserved density
    /// R² + I₋I₊, the imaginary part half a step behind R before the
    /// step and half a step ahead after it. R² + I² is off by about
    /// dt H ρ, and where ψ sits on particle 2's wall that error reaches
    /// x₁'s density too.
    fn conserved_marginals(wave: &Wave, behind: &[f32]) -> [Vec<f64>; 2] {
        let n = wave.n;
        let mut m = [vec![0.0; n], vec![0.0; n]];
        for k in 0..n * n {
            let d = (wave.re[k] * wave.re[k] + behind[k] * wave.im[k]) as f64 * wave.dx;
            m[0][k % n] += d;
            m[1][k / n] += d;
        }
        m
    }

    /// The pair on a small grid flown for a while, with particle 2's
    /// wall close by or without it: the marginals, and runs started on
    /// a row of points across the pair.
    fn flown(entangled: bool, walled: bool) -> ([Vec<f64>; 2], Vec<Vec2>) {
        let mut wave = Wave::new(128, 0.1);
        if walled {
            potential(&mut wave, |p| EPR_FORCE * (p[1] - 1.0).max(0.0));
        }
        epr(&mut wave, entangled);
        let start = (0..20).map(|i| [0.1 * i as f64 - 1.0, 0.3 - 0.05 * i as f64]).collect();
        let mut particles = crate::particles::Particles::new(Kind::Quantum, start, &wave);
        let dt = stable_dt(0.1, EPR_FORCE * 5.4);
        let mut behind = wave.im.clone();
        for _ in 0..(1.0 / dt) as usize {
            behind.clone_from(&wave.im);
            wave.step(dt);
            particles.step(&wave, dt);
        }
        (conserved_marginals(&wave, &behind), particles.positions)
    }

    #[test]
    fn a_wall_on_particle_two_moves_particle_one_only_if_entangled() {
        for entangled in [true, false] {
            let [(m, p), (n, q)] = [true, false].map(|w| flown(entangled, w));
            let gap = largest_gap(&m[0], &n[0]);
            let other = largest_gap(&m[1], &n[1]);
            let moved = p.iter().zip(&q).map(|(p, q)| (p[0] - q[0]).abs()).fold(0.0, f64::max);
            assert!(gap < 1e-4 && other > 0.5, "{entangled} {gap} {other}");
            assert_eq!(moved > 0.1, entangled, "{moved}");
        }
    }
}
