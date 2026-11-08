use std::env;
use std::f64::consts::TAU;
use std::fs;

use pilot_wave::encounter::{self, Encounter};
use pilot_wave::experiment::{
    experiment, mean_transmission, pattern, transmission, Experiment, Space, BARRIER_HEIGHT,
    BARRIER_WIDTH, ENERGY, SLIT_SEPARATION, SPEED,
};
use pilot_wave::grid::Vec2;
use pilot_wave::image::{self, Image};
use pilot_wave::particles::{Kind, Particles};
use pilot_wave::rng::Rng;
use pilot_wave::sample::Sampler;
use pilot_wave::screen::Screen;
use pilot_wave::wave::Wave;
use pilot_wave::{pair, plot, png, view};

const PIXELS: usize = 800;
const GRID: usize = 1024;
/// The traps and the box fit a smaller domain.
const SMALL_GRID: usize = 256;
/// Time between frames and the number of particles whose trails are
/// drawn.
const FRAME_EVERY: f64 = 0.5;
const TRAILS: usize = 300;
/// Time between samples of the density added into the
/// time-integrated picture behind the trails.
const GLOW_EVERY: f64 = 0.05;
/// Bins of a whole number of cells, eight on the 1024 grid, so each
/// gathers the same share of the current.
const BINS: usize = 128;
/// Classical particles stop this far inside the domain's edge, within
/// the sponge that absorbs the wave there.
const EDGE: f64 = 1.0;
const PLOT: (usize, usize) = (800, 300);
/// Two particles on a line: the separation beyond which a pair counts
/// as apart and within which as close, the line view's height and
/// the time between reports of the separation.
const APART: f64 = 3.0;
const NEAR: f64 = 0.5;
const LINE_HEIGHT: usize = 200;
const REPORT_EVERY: f64 = 1.0;
/// Local maxima on the screen below this fraction of the largest are
/// taken as noise rather than fringes.
const FRINGE_FLOOR: f64 = 0.05;

type Trail = Vec<Vec2>;

fn distance(a: Vec2, b: Vec2) -> f64 {
    (a[0] - b[0]).hypot(a[1] - b[1])
}

fn densities(wave: &Wave) -> Vec<f32> {
    (0..wave.n * wave.n).map(|k| wave.density(k)).collect()
}

/// The time-integrated density as a field to draw, walls included.
fn glow_field(wave: &Wave, glow: &[f32]) -> Wave {
    let mut field = Wave::new(wave.n, wave.dx);
    field.re = glow.iter().map(|g| g.sqrt()).collect();
    field.keep = wave.keep.clone();
    field
}

/// Stop classical particles that reach the sponge, as the wave stops
/// there, before they leave the grid their force is read from.
fn stop_at_edges(particles: &mut Particles, half: f64) {
    for (p, v) in particles.positions.iter().zip(&mut particles.velocities) {
        if p[0].abs().max(p[1].abs()) > half - EDGE {
            *v = [0.0; 2];
        }
    }
}

/// Standard deviations of the density along x and y.
fn spread(wave: &Wave) -> Vec2 {
    let mut sums = [0.0; 5];
    for k in 0..wave.n * wave.n {
        let p = wave.position((k % wave.n, k / wave.n));
        let w = wave.density(k) as f64;
        let terms = [w, w * p[0], w * p[1], w * p[0] * p[0], w * p[1] * p[1]];
        for (s, t) in sums.iter_mut().zip(terms) {
            *s += t;
        }
    }
    let [_, x, y, xx, yy] = sums.map(|s| s / sums[0]);
    [(xx - x * x).sqrt(), (yy - y * y).sqrt()]
}

fn largest_displacement(start: &[Vec2], end: &[Vec2]) -> f64 {
    start.iter().zip(end).map(|(a, b)| distance(*a, *b)).fold(0.0, f64::max)
}

/// Largest departure of a trail from the circle that the packet
/// centre, starting at `center`, traces around the trap's middle,
/// shifted by the particle's starting offset from that centre.
fn orbit_error(trails: &[Trail], center: Vec2) -> f64 {
    let radius = center[0].hypot(center[1]);
    let error = |t: &Trail| {
        let middle = [t[0][0] - center[0], t[0][1] - center[1]];
        t.iter().map(|p| (distance(*p, middle) - radius).abs()).fold(0.0, f64::max)
    };
    trails.iter().map(error).fold(0.0, f64::max)
}

/// Cells in a wall or next to one, where the guidance is least
/// accurate and the density vanishes.
fn at_wall(wave: &Wave) -> Vec<bool> {
    let n = wave.n;
    let wall = |k: usize| wave.keep[k] == 0.0;
    let at = |k: usize| {
        let (i, j) = (k % n, k / n);
        let inner = i > 0 && j > 0 && i < n - 1 && j < n - 1;
        wall(k) || inner && [k - 1, k + 1, k - n, k + n].into_iter().any(wall)
    };
    (0..n * n).map(at).collect()
}

/// Particles in or beside walls against the probability there, which
/// equivariance says should agree.
fn print_walls(name: &str, wave: &Wave, positions: &[Vec2]) {
    let at = at_wall(wave);
    let cell = |p: &Vec2| {
        let (x, y) = wave.cell(*p);
        let clamp = |c: f64| c.round().clamp(0.0, (wave.n - 1) as f64) as usize;
        wave.index((clamp(x), clamp(y)))
    };
    let stuck = positions.iter().filter(|p| at[cell(p)]).count();
    let fraction = stuck as f64 / positions.len() as f64;
    let probability: f64 = (0..wave.n * wave.n)
        .filter(|&k| at[k])
        .map(|k| wave.density(k) as f64 * wave.dx * wave.dx)
        .sum();
    println!("{name}: {fraction:.4} of the particles at walls, probability there {probability:.4}");
}

fn print_screen(name: &str, screen: &Screen, count: usize) {
    let arrived: f64 = screen.hits.iter().sum();
    let through: f64 = screen.flux.iter().sum();
    let fraction = arrived / count as f64;
    println!(
        "{name}: {arrived} of {count} particles reached the screen ({fraction:.4}), \
         probability through it {through:.4}"
    );
}

/// The barrier's transmission against the analytic one, and where in
/// the packet the transmitted particles started.
fn print_barrier(e: &Experiment, start: &[Vec2], end: &[Vec2]) {
    let height = BARRIER_HEIGHT * ENERGY;
    let spread = 1.0 / (2.0 * e.packet.width[1]);
    let mean = mean_transmission(SPEED, spread, height, BARRIER_WIDTH);
    let plane = transmission(SPEED, height, BARRIER_WIDTH);
    println!("barrier: analytic transmission {mean:.4} over the packet, {plane:.4} at its centre");
    let started = |through: bool| {
        let pairs = start.iter().zip(end);
        pairs.filter(move |(_, b)| (b[1] > 0.0) == through).map(|(a, _)| a[1])
    };
    let rear = started(true).fold(f64::MAX, f64::min);
    let front = started(false).fold(f64::MIN, f64::max);
    println!("barrier: transmitted started at y >= {rear:.3}, reflected at y <= {front:.3}");
}

/// Mean r × v of the particles, the orbital angular momentum each one
/// carries about the trap's centre.
fn angular_momentum(wave: &Wave, positions: &[Vec2]) -> f64 {
    let l = |p: &Vec2| {
        let v = wave.velocity(*p);
        p[0] * v[1] - p[1] * v[0]
    };
    positions.iter().map(l).sum::<f64>() / positions.len() as f64
}

fn print_summary(e: &Experiment, start: &[Vec2], end: &[Vec2], trails: &[Trail], initial: &[f32]) {
    let s = spread(&e.wave);
    match e.name {
        "harmonic" => println!(
            "harmonic: spread {:.4} x {:.4}, orbit error {:.4}, largest displacement {:.4}",
            s[0],
            s[1],
            orbit_error(trails, e.packet.center),
            largest_displacement(start, end)
        ),
        "ground" | "vortex" | "antivortex" | "orbital" | "flux-trap" => {
            let peak = initial.iter().cloned().fold(0.0, f32::max);
            let now = e.wave.centred_densities(e.dt);
            let change = now.iter().zip(initial).map(|(a, b)| (a - b).abs()).fold(0.0, f32::max);
            println!(
                "{}: energy {:.4}, angular momentum {:.4} per particle, largest displacement \
                 {:.2e}, largest density change {:.2e} of the peak",
                e.name,
                e.wave.energy(),
                angular_momentum(&e.wave, end),
                largest_displacement(start, end),
                change / peak
            );
        }
        "free" => {
            let w = e.packet.width;
            let t = e.wave.time;
            let expected = w.map(|w| w * (1.0 + (t / (2.0 * w * w)).powi(2)).sqrt());
            println!(
                "free: spread {:.4} x {:.4}, expected {:.4} x {:.4}",
                s[0], s[1], expected[0], expected[1]
            );
        }
        "box" => {
            println!("box: norm {:.4}", e.wave.norm());
            print_walls(e.name, &e.wave, end);
        }
        "barrier" => print_barrier(e, start, end),
        "double-slit" | "single-slit" | "disc" | "flux-slit" => print_walls(e.name, &e.wave, end),
        _ => {}
    }
}

/// Local maxima of a curve over the screen above FRINGE_FLOOR of its
/// largest, each refined by a parabola through its bin and the two
/// beside it, as fringe orders d sinθ.
fn fringe_orders(screen: &Screen, v: &[f64]) -> Vec<f64> {
    let top = v.iter().cloned().fold(0.0, f64::max);
    let width = screen.x(1) - screen.x(0);
    let order = |x: f64| SLIT_SEPARATION * x / x.hypot(screen.y);
    let refined = |b: usize| {
        let shift = 0.5 * (v[b - 1] - v[b + 1]) / (v[b - 1] - 2.0 * v[b] + v[b + 1]);
        order(screen.x(b) + shift * width)
    };
    (1..v.len() - 1)
        .filter(|&b| v[b] > v[b - 1] && v[b] >= v[b + 1] && v[b] > FRINGE_FLOOR * top)
        .map(refined)
        .collect()
}

/// The fringes of the current and of the theory in orders d sinθ,
/// and the two bins at the centre against the peak for the current,
/// the arrivals and the theory: half a flux quantum between the
/// slits moves every maximum by half an order and darkens the centre.
fn print_fringes(name: &str, screen: &Screen) {
    let theory: Vec<f64> = (0..BINS).map(|b| pattern(name, screen.x(b), screen.y)).collect();
    let orders = |v: &[f64]| -> String {
        fringe_orders(screen, v).iter().map(|o| format!(" {o:.2}")).collect()
    };
    let centre = |v: &[f64]| {
        let top = v.iter().cloned().fold(0.0, f64::max);
        (v[BINS / 2 - 1] + v[BINS / 2]) / (2.0 * top)
    };
    println!("{name}: current maxima at orders{}", orders(&screen.flux));
    println!("{name}: theory maxima at orders{}", orders(&theory));
    println!(
        "{name}: centre {:.4} of the peak for the current, {:.4} for the arrivals, \
         {:.4} for the theory",
        centre(&screen.flux),
        centre(&screen.hits),
        centre(&theory)
    );
}

/// Write the screen as a plot and a table, arrivals and theory scaled
/// to the same total as the probability current.
fn write_screen(out: &str, name: &str, screen: &Screen) {
    let total: f64 = screen.flux.iter().sum();
    let scale = |v: Vec<f64>| {
        let sum: f64 = v.iter().sum();
        v.iter().map(|x| x * total / sum.max(1e-300)).collect::<Vec<f64>>()
    };
    let hits = scale(screen.hits.clone());
    let theory = scale((0..BINS).map(|b| pattern(name, screen.x(b), screen.y)).collect());
    let series: [(&[f64], _); 3] = [
        (&hits, [255, 255, 255]),
        (&screen.flux, [80, 200, 255]),
        (&theory, [255, 160, 30]),
    ];
    let image = plot::curves(PLOT.0, PLOT.1, &series);
    png::write(format!("{out}/{name}_screen.png"), &image).expect("write screen plot");
    let rows: String = (0..BINS)
        .map(|b| format!("{},{},{},{}\n", screen.x(b), hits[b], screen.flux[b], theory[b]))
        .collect();
    fs::write(format!("{out}/{name}_screen.csv"), format!("x,hits,flux,theory\n{rows}"))
        .expect("write screen table");
}

/// The configuration-space field with the diagonal and the first
/// run, over the line view.
fn line_frame(e: &Experiment, positions: &[Vec2], run: &[Vec2]) -> Image {
    let mut image = view::field(&e.wave, PIXELS, view::Style::Phase, e.peak);
    let heads: Vec<Trail> = positions.iter().map(|p| vec![*p]).collect();
    view::particles(&mut image, &e.wave, &heads);
    view::diagonal(&mut image);
    view::run(&mut image, &e.wave, run);
    let peak = 1.0 / (TAU.sqrt() * e.packet.width[0]);
    let line = view::line(&e.wave, positions, PIXELS, LINE_HEIGHT, peak);
    image::stack(&image, &line)
}

/// Crossings and closest approach, the most probability and runs
/// ever found close together, and which side of the diagonal the
/// runs and the field end on.
fn print_encounter(e: &Experiment, encounter: &Encounter, near: [f64; 2], end: &[Vec2]) {
    let name = e.name;
    let mut closest = encounter.closest.clone();
    closest.sort_by(f64::total_cmp);
    let (least, median) = (closest[0], closest[closest.len() / 2]);
    let crossed = encounter.crossed();
    println!(
        "{name}: {crossed:.4} of the runs crossed the diagonal, \
         closest approach {least:.4}, median {median:.3}"
    );
    println!(
        "{name}: within {NEAR} of each other at most {:.4} of the probability, {:.4} of the runs",
        near[0], near[1]
    );
    let past = end.iter().filter(|p| p[0] > p[1]).count() as f64 / end.len() as f64;
    let field = pair::probability(&e.wave, |p| p[0] > p[1]);
    let left = e.wave.norm();
    println!(
        "{name}: x1 > x2 at the end in {past:.4} of the runs, \
         probability {field:.4} of the {left:.4} left"
    );
}

/// The pair's separation over time and how many runs end apart,
/// against the probability apart and that absorbed, which only
/// fragments flying off reach.
fn print_separations(e: &Experiment, separations: &[(f64, f64, f64)], end: &[Vec2]) {
    let name = e.name;
    for (t, mean, spread) in separations {
        println!("{name}: t = {t:4.1}, separation {mean:.3} ± {spread:.3}");
    }
    let runs = encounter::fraction_apart(end, APART);
    let field = pair::probability(&e.wave, |p| encounter::separation(p) > APART);
    let absorbed = 1.0 - e.wave.norm();
    println!(
        "{name}: apart beyond {APART} at the end {runs:.4} of the runs, \
         probability {field:.4}, absorbed {absorbed:.4}"
    );
}

/// Two particles on a line: frames of configuration space over the
/// line, trails, and the checks.
fn run_line(mut e: Experiment, count: usize, out: &str) {
    let name = e.name;
    let sampler = Sampler::new(&e.wave);
    let mut rng = Rng::new(1);
    let start: Vec<Vec2> = (0..count).map(|_| sampler.draw(&mut rng)).collect();
    let mut particles = Particles::new(e.kind, start.clone(), &e.wave);
    let mut trails: Vec<Trail> = start[..count.min(TRAILS)].iter().map(|p| vec![*p]).collect();
    let mut encounter = Encounter::new(&start);
    let mut near = [0.0f64; 2];
    let mut separations = vec![];
    let mut glow = vec![0.0; e.wave.n * e.wave.n];
    let steps = (e.duration / e.dt).ceil() as usize;
    let [per_frame, per_glow, per_report] =
        [FRAME_EVERY, GLOW_EVERY, REPORT_EVERY].map(|t| (t / e.dt) as usize);
    for step in 0..=steps {
        if step > 0 {
            e.wave.step(e.dt);
            particles.step(&e.wave, e.dt);
        }
        let positions = &particles.positions;
        encounter.record(positions);
        for (trail, p) in trails.iter_mut().zip(positions) {
            trail.push(*p);
        }
        if step % per_glow == 0 {
            for (g, d) in glow.iter_mut().zip(densities(&e.wave)) {
                *g += d;
            }
            let field = pair::probability(&e.wave, |p| encounter::separation(p) < NEAR);
            let runs = 1.0 - encounter::fraction_apart(positions, NEAR);
            near = [near[0].max(field), near[1].max(runs)];
        }
        if step % per_report == 0 {
            let (mean, spread) = encounter::spread(positions);
            separations.push((e.wave.time, mean, spread));
        }
        if step % per_frame == 0 || step == steps {
            let image = line_frame(&e, positions, &trails[0]);
            png::write(format!("{out}/{name}_{:04}.png", step / per_frame), &image)
                .expect("write frame");
        }
    }
    let glow = glow_field(&e.wave, &glow);
    let mut image = view::field(&glow, PIXELS, view::Style::Density, glow.peak());
    view::particles(&mut image, &e.wave, &trails);
    view::diagonal(&mut image);
    view::run(&mut image, &e.wave, &trails[0]);
    png::write(format!("{out}/{name}_trajectories.png"), &image).expect("write trajectories");
    let end = &particles.positions;
    print_encounter(&e, &encounter, near, end);
    if name.starts_with("molecule") {
        print_separations(&e, &separations, end);
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let name: &'static str = pilot_wave::experiment::NAMES
        .into_iter()
        .find(|n| Some(*n) == args.get(1).map(String::as_str))
        .unwrap_or("double-slit");
    let count: usize = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(2000);
    let out = args.get(3).map(String::as_str).unwrap_or("out");
    let kind = if args.get(4).map(String::as_str) == Some("classical") { Kind::Classical } else { Kind::Quantum };
    fs::create_dir_all(out).expect("create output directory");

    let small = matches!(name, "harmonic" | "ground" | "box" | "flux-trap");
    let n = if small { SMALL_GRID } else { GRID };
    let mut e = experiment(name, n, kind);
    if e.space == Space::Line {
        return run_line(e, count, out);
    }
    let mut rng = Rng::new(1);
    let sampler = Sampler::new(&e.wave);
    let start: Vec<Vec2> = (0..count).map(|_| sampler.draw(&mut rng)).collect();
    let mut particles = Particles::new(kind, start.clone(), &e.wave);
    let mut trails: Vec<Trail> = start[..count.min(TRAILS)].iter().map(|p| vec![*p]).collect();
    let half = e.wave.extent() / 2.0;
    let mut screen = e.screen.map(|y| Screen::new(y, half, BINS));
    let initial = e.wave.centred_densities(e.dt);
    let mut glow = vec![0.0; n * n];
    let steps = (e.duration / e.dt).ceil() as usize;
    let per_frame = (FRAME_EVERY / e.dt) as usize;
    let per_glow = (GLOW_EVERY / e.dt) as usize;
    for step in 1..=steps {
        e.wave.step(e.dt);
        let before = particles.positions.clone();
        particles.step(&e.wave, e.dt);
        stop_at_edges(&mut particles, half);
        if let Some(screen) = &mut screen {
            screen.record(&before, &particles.positions);
            screen.integrate(&e.wave, e.dt);
        }
        for (trail, p) in trails.iter_mut().zip(&particles.positions) {
            trail.push(*p);
        }
        if step % per_glow == 0 {
            for (g, d) in glow.iter_mut().zip(densities(&e.wave)) {
                *g += d;
            }
        }
        if step % per_frame == 0 || step == steps {
            let mut image = view::field(&e.wave, PIXELS, view::Style::Phase, e.peak);
            let heads: Vec<Trail> = particles.positions.iter().map(|p| vec![*p]).collect();
            view::particles(&mut image, &e.wave, &heads);
            png::write(format!("{out}/{name}_{:04}.png", step / per_frame), &image).expect("write frame");
        }
    }
    let glow = glow_field(&e.wave, &glow);
    let mut image = view::field(&glow, PIXELS, view::Style::Density, glow.peak());
    view::particles(&mut image, &e.wave, &trails);
    png::write(format!("{out}/{name}_trajectories.png"), &image).expect("write trajectories");
    if let Some(screen) = &screen {
        write_screen(out, name, screen);
        print_screen(name, screen, count);
        if matches!(name, "double-slit" | "flux-slit") {
            print_fringes(name, screen);
        }
    }
    print_summary(&e, &start, &particles.positions, &trails, &initial);
}
