//! Two particles bound by a Gaussian well, sent at a wall only one
//! of them feels: the pair bounces as a pair, or breaks up.
use crate::experiment::{band, Expectation, Experiment, Measure, Realisation, ELECTRON_MASS};
use crate::pair;
use crate::wave::{stable_dt, Vec2, Wave};

/// The line is |x| ≤ 12.8.
const CELLS: usize = 256;
const DX: f64 = 0.1;
/// The well U(x₁ − x₂), which binds exactly one state of the relative
/// motion, 2.8 deep; a second would need a well half again as deep.
const DEPTH: f64 = 6.0;
const RANGE: f64 = 0.4;
/// The hard wall particle 1 feels, where the molecule starts, the
/// widths of its centre and separation before relaxing, and how long
/// it relaxes in imaginary time.
const WALL: f64 = 5.0;
const START: f64 = -2.0;
const WIDTHS: [f64; 2] = [1.0, 0.5];
const RELAX: f64 = 4.0;
/// The gentle and the hard kick, total wavenumbers of the pair.
const KICK: f64 = 2.0;
const HARD_KICK: f64 = 5.0;
/// Particles count as apart beyond this separation.
const APART: f64 = 3.0;

const REALISATION: Realisation = Realisation {
    particle: "two electrons on a line",
    mass: ELECTRON_MASS,
    length: 10e-9,
    reference: "for scale only, a nanowire; no particular experiment",
};

fn bond(p: Vec2) -> f64 {
    let r = p[0] - p[1];
    -DEPTH * (-r * r / (2.0 * RANGE * RANGE)).exp()
}

fn gaussian(x: f64, width: f64) -> f64 {
    (-x * x / (4.0 * width * width)).exp()
}

/// The pair in the ground state of its bond, Gaussian in its centre.
fn bound(wave: &mut Wave) {
    let [w, s] = WIDTHS;
    pair::fill(wave, |x1, x2| [gaussian((x1 + x2) / 2.0 - START, w) * gaussian(x1 - x2, s), 0.0]);
    let dt = stable_dt(DX, DEPTH);
    pair::relax(wave, dt, (RELAX / dt) as usize);
}

fn kicked(wave: &mut Wave) {
    bound(wave);
    pair::kick(wave, KICK);
}

fn kicked_hard(wave: &mut Wave) {
    bound(wave);
    pair::kick(wave, HARD_KICK);
}

fn molecule_at(
    name: &'static str,
    summary: &'static str,
    source: fn(&mut Wave),
    duration: f64,
    apart: f64,
    tolerance: f64,
) -> Experiment {
    Experiment {
        name,
        summary,
        realisation: REALISATION,
        cells: CELLS,
        dx: DX,
        walls: vec![band(0, WALL, f64::INFINITY)],
        potential: bond,
        source,
        spin: [0.0; 3],
        particles: 20000,
        duration,
        screen: None,
        theory: None,
        expectations: vec![
            Expectation {
                what: "fraction of the runs more than 3 apart at the end",
                measure: Measure::Apart { beyond: APART },
                expected: apart,
                tolerance,
            },
            Expectation {
                what: "probability apart at the end, counting what left the line",
                measure: Measure::ProbabilityApart { beyond: APART },
                expected: apart,
                tolerance: tolerance.max(0.02),
            },
        ],
    }
}

pub fn molecule() -> Experiment {
    molecule_at(
        "molecule",
        "\
Two particles on a line bound by a Gaussian well U = −6 exp(−r²/0.32) \
in their separation r = x₁ − x₂, which binds one state of the relative \
motion, 2.8 deep. The pair starts as a Gaussian relaxed in imaginary \
time into that bound state, then gets total wavenumber 2, kinetic \
energy 1, a third of the binding energy, towards a hard wall at x = 5 \
that only particle 1 feels: in configuration space the half plane \
x₁ > 5. A wall both felt would be a corner, which reflects the pair \
like a corner reflector and never breaks it. Particle 1 bounces and \
drags particle 2 round with it; the pair stays bound.",
        kicked,
        15.0,
        0.0,
        0.002,
    )
}

pub fn molecule_break() -> Experiment {
    molecule_at(
        "molecule-break",
        "\
The same molecule kicked with total wavenumber 5, kinetic energy 6.25, \
2.2 times the binding energy. Particle 1 turns at the wall, particle \
2 runs on through where the wall is not, and the relative motion \
takes up the pair's kinetic energy: most pairs break, flying off as a \
ring of fragments around the bound part that still bounces back along \
the diagonal, and a fraction stays bound although the kick carried \
more than twice the binding energy.",
        kicked_hard,
        6.0,
        0.85,
        0.03,
    )
}
