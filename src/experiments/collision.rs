//! Two particles on a line heading for each other: fermions, bosons
//! and distinguishable particles, the exchange symmetry alone telling
//! them apart; and distinguishable particles scattering off each
//! other into an entangled pair.
use crate::experiment::{free, Expectation, Experiment, Measure, Realisation, ELECTRON_MASS};
use crate::pair::{self, Orbital, Symmetry};
use crate::wave::{Vec2, Wave};

/// The line is |x| ≤ 12.8.
const CELLS: usize = 256;
const DX: f64 = 0.1;
/// Where the packets start, their width and their wavenumber.
const START: f64 = 4.0;
const WIDTH: f64 = 1.0;
const WAVENUMBER: f64 = 2.0;
const DURATION: f64 = 5.0;
/// The Gaussian bump between the scattering particles, which at their
/// relative energy of 4 lets about half the probability through, and
/// that share, measured from the field.
const BUMP: f64 = 4.8;
const BUMP_RANGE: f64 = 0.4;
const PASSED: f64 = 0.5;
/// The scattered pair runs on a line twice as long, |x| ≤ 25.6, for
/// longer, so that the slow tails of both branches clear the origin.
const SCATTER_CELLS: usize = 512;
const SCATTER_DURATION: f64 = 8.0;

const REALISATION: Realisation = Realisation {
    particle: "two electrons on a line",
    mass: ELECTRON_MASS,
    length: 10e-9,
    reference: "for scale only, a nanowire; no particular experiment",
};

fn orbitals() -> [Orbital; 2] {
    [
        Orbital { center: -START, width: WIDTH, wavenumber: WAVENUMBER },
        Orbital { center: START, width: WIDTH, wavenumber: -WAVENUMBER },
    ]
}

fn distinguishable_pair(wave: &mut Wave) {
    let [a, b] = orbitals();
    pair::write(wave, a, b, Symmetry::Distinguishable)
}

fn boson_pair(wave: &mut Wave) {
    let [a, b] = orbitals();
    pair::write(wave, a, b, Symmetry::Bosons)
}

fn fermion_pair(wave: &mut Wave) {
    let [a, b] = orbitals();
    pair::write(wave, a, b, Symmetry::Fermions)
}

fn bump(p: Vec2) -> f64 {
    let r = p[0] - p[1];
    BUMP * (-r * r / (2.0 * BUMP_RANGE * BUMP_RANGE)).exp()
}

fn meeting(
    name: &'static str,
    summary: &'static str,
    source: fn(&mut Wave),
    expectations: Vec<Expectation>,
) -> Experiment {
    Experiment {
        name,
        summary,
        realisation: REALISATION,
        cells: CELLS,
        dx: DX,
        walls: vec![],
        potential: free,
        source,
        spin: [0.0; 3],
        particles: 20000,
        duration: DURATION,
        screen: None,
        theory: None,
        expectations,
    }
}

/// Crossing, bunching and the bunching's probability, the peaks from
/// free propagation of the relative motion in closed form.
fn free_expectations(crossed: f64, together: f64) -> Vec<Expectation> {
    vec![
        Expectation {
            what: "fraction of the runs in which the particles pass each other",
            measure: Measure::Crossed,
            expected: crossed,
            tolerance: 0.005,
        },
        Expectation {
            what: "most of the runs ever within 0.5 of each other, free relative motion in closed form",
            measure: Measure::TogetherPeak,
            expected: together,
            tolerance: 0.01,
        },
        Expectation {
            what: "most probability ever within 0.5, likewise",
            measure: Measure::ProbabilityTogetherPeak,
            expected: together,
            tolerance: 0.005,
        },
    ]
}

pub fn distinguishable() -> Experiment {
    let mut expectations = free_expectations(1.0, 0.197);
    expectations.push(Expectation {
        what: "mutual information of x₁ and x₂ in the field at the end: a product state",
        measure: Measure::MutualInformation,
        expected: 0.0,
        tolerance: 0.001,
    });
    meeting(
        "distinguishable",
        "\
Two packets 1 wide start at x = ∓4 on a line and head for each other \
at wavenumber ±2, as the product a(x₁)b(x₂): in configuration space \
one packet moving towards the diagonal x₁ = x₂. Nothing acts between \
the particles, so they pass straight through each other and stay a \
product; the probability of finding them within 0.5 of each other \
peaks at t ≈ 2 at 0.197, the free spreading of the relative motion in \
closed form.",
        distinguishable_pair,
        expectations,
    )
}

pub fn bosons() -> Experiment {
    meeting(
        "bosons",
        "\
The same two packets as identical bosons, a(x₁)b(x₂) + b(x₁)a(x₂): \
two mirrored packets in configuration space. The symmetric ψ guides \
mirrored points in mirrored ways, so the diagonal is a trajectory of \
its own and no run crosses it, yet the density is largest there and \
the particles bunch: the probability of finding them within 0.5 of \
each other peaks at 0.287, 1.46 times the distinguishable value.",
        boson_pair,
        free_expectations(0.0, 0.287),
    )
}

pub fn fermions() -> Experiment {
    meeting(
        "fermions",
        "\
The same two packets as identical fermions, a(x₁)b(x₂) − b(x₁)a(x₂). \
The antisymmetric ψ vanishes on the whole diagonal, and that node is \
all there is to Pauli exclusion here: no force acts, yet the guiding \
velocity turns every run back before the particles meet, and the \
probability of finding them within 0.5 of each other peaks at only \
0.109, 0.56 times the distinguishable value.",
        fermion_pair,
        free_expectations(0.0, 0.109),
    )
}

pub fn scattered() -> Experiment {
    let entropy = -PASSED * PASSED.ln() - (1.0 - PASSED) * (1.0 - PASSED).ln();
    let mut e = meeting(
        "scattered",
        "\
The two distinguishable packets again, now with a Gaussian bump of \
height 4 and range 0.4 between the particles, which at their relative \
energy of 4 lets about half the wave through and reflects the rest. \
Entangled by nothing but having met, the pair leaves in two branches: \
passed, particle 1 on the right and particle 2 on the left, or \
bounced, each back where it came from. Either way the particles end \
on opposite sides of the origin, so a measurement of one, which side \
it is found on, tells which side the other is on, however far apart \
they are, and the mutual information of the joint density stays at \
the entropy of the two branches, about ln 2, while they separate. In \
this theory the measurement reads a position the particle already \
has; the branch the pair is not in goes on as an empty wave. No \
laboratory can do this: it would need two particles in one dimension \
that interact only with each other.",
        distinguishable_pair,
        vec![
            Expectation {
                what: "fraction of the runs in which the particles pass each other, the field's transmission",
                measure: Measure::Crossed,
                expected: PASSED,
                tolerance: 0.03,
            },
            Expectation {
                what: "fraction of the runs ending with the particles on opposite sides of the origin",
                measure: Measure::OppositeSides,
                expected: 1.0,
                tolerance: 0.005,
            },
            Expectation {
                what: "mutual information of x₁ and x₂ in the field at the end, the entropy of the two branches",
                measure: Measure::MutualInformation,
                expected: entropy,
                tolerance: 0.05,
            },
        ],
    );
    e.potential = bump;
    e.cells = SCATTER_CELLS;
    e.duration = SCATTER_DURATION;
    e
}
