//! Two particles on a line heading for each other: fermions, bosons
//! and distinguishable particles, the exchange symmetry alone telling
//! them apart.
use crate::experiment::{free, Expectation, Experiment, Measure, Realisation, ELECTRON_MASS};
use crate::pair::{self, Orbital, Symmetry};
use crate::wave::Wave;

/// The line is |x| ≤ 12.8.
const CELLS: usize = 256;
const DX: f64 = 0.1;
/// Where the packets start, their width and their wavenumber.
const START: f64 = 4.0;
const WIDTH: f64 = 1.0;
const WAVENUMBER: f64 = 2.0;
const DURATION: f64 = 5.0;
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
        free_expectations(1.0, 0.197),
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

