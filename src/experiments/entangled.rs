//! An entangled pair flying apart: the correlation that survives
//! separation, and a wall only one particle feels moving the other
//! without any observable sign of it.
use crate::experiment::{
    Expectation, Experiment, Measure, Potential, Realisation, Source, Space, ELECTRON_MASS,
};

/// The line is |x| ≤ 25.6 in cells of 0.05: on cells of 0.1 the runs
/// that crossed the fringes of the wave coming back from the wall
/// bunched where the field did not.
const CELLS: usize = 1024;
const DX: f64 = 0.05;
/// Widths of |ψ|² in the pair's centre and separation, the wavenumber
/// carrying each particle away from the other, and how long it flies.
const WIDTHS: [f64; 2] = [1.0, 0.6];
const WAVENUMBER: f64 = 3.0;
const DURATION: f64 = 3.5;
/// The wall only particle 2 feels: a steep ramp from x₂ = 5 on, which
/// turns it back within half a unit. Cells held at zero let a few runs
/// step into them and froze particle 1 with them.
const WALL: f64 = 5.0;
const SLOPE: f64 = 20.0;

const REALISATION: Realisation = Realisation {
    particle: "two electrons on a line",
    mass: ELECTRON_MASS,
    length: 10e-9,
    reference: "for scale only, a nanowire; no particular experiment",
};

fn flight(
    name: &'static str,
    summary: &'static str,
    entangled: bool,
    walled: bool,
    expectations: Vec<Expectation>,
) -> Experiment {
    let ramp = Potential::Ramp { axis: 1, from: WALL, slope: SLOPE };
    Experiment {
        name,
        summary,
        realisation: REALISATION,
        space: Space::Line,
        cells: CELLS,
        dx: DX,
        walls: vec![],
        potentials: if walled { vec![ramp] } else { vec![] },
        source: Source::Epr { widths: WIDTHS, wavenumber: WAVENUMBER, entangled },
        spin: [0.0; 3],
        particles: 20000,
        duration: DURATION,
        screen: None,
        theory: None,
        expectations,
    }
}

/// Correlation of x₁ and x₂ in free flight, from the centre's and the
/// separation's spreading: (Var c − Var r/4)/(Var c + Var r/4).
fn free_correlation(t: f64) -> f64 {
    let [c, r] = WIDTHS;
    let (centre, relative) = (c * c + t * t / (16.0 * c * c), r * r + t * t / (r * r));
    (centre - relative / 4.0) / (centre + relative / 4.0)
}

pub fn entangled() -> Experiment {
    flight(
        "entangled",
        "\
Two particles start together and fly apart, particle 1 to the left: \
with the centre c = (x₁ + x₂)/2 and the separation r = x₁ − x₂, \
ψ ∝ exp(−c²/4 − r²/1.44) e^{−3ir}, positions correlated and momenta \
anticorrelated. Each particle's velocity depends on where the other \
is, however far away. The correlation of x₁ and x₂ starts at 0.835, \
passes through zero at t = 1.2, when |ψ|² is momentarily a product \
and the entanglement sits wholly in the phase, and settles negative, \
−0.659 at t = 3.5 with the particles 21 apart: a particle 1 that went \
far left is matched by a particle 2 that went far right. The runs \
follow the field's correlation, and the field the closed form of free \
spreading.",
        true,
        false,
        vec![
            Expectation {
                what: "correlation of x₁ and x₂ over the runs at the end, free spreading in closed form",
                measure: Measure::Correlation,
                expected: free_correlation(DURATION),
                tolerance: 0.02,
            },
            Expectation {
                what: "that correlation minus the field's",
                measure: Measure::CorrelationMinusField,
                expected: 0.0,
                tolerance: 0.01,
            },
        ],
    )
}

fn wall_expectations(moved: f64, tolerance: f64) -> Vec<Expectation> {
    vec![
        Expectation {
            what: "median change of x₁ against the same run without the wall",
            measure: Measure::MovedVsTwin { axis: 0 },
            expected: moved,
            tolerance,
        },
        Expectation {
            what: "largest change of x₁'s density, of its peak: no signalling",
            measure: Measure::MarginalVsTwin { axis: 0 },
            expected: 0.0,
            tolerance: 0.001,
        },
        Expectation {
            what: "largest change of x₂'s density, of its peak",
            measure: Measure::MarginalVsTwin { axis: 1 },
            expected: 0.95,
            tolerance: 0.1,
        },
    ]
}

pub fn entangled_wall() -> Experiment {
    flight(
        "entangled-wall",
        "\
The same entangled pair with a wall only particle 2 feels, a steep \
ramp from x₂ = 5 on, which it meets at t ≈ 1.7 with particle 1 about \
10 away on the other side and nothing in its way. The same starting \
points are run again without the wall. Each run's particle 1 ends \
somewhere else once particle 2 has bounced, since its velocity \
depends on where particle 2 is; yet the density of x₁ is the same \
with the wall and without, so no observer of particle 1 alone can \
tell: nonlocal guidance, no signalling.",
        true,
        true,
        wall_expectations(0.75, 0.15),
    )
}

pub fn unentangled_wall() -> Experiment {
    flight(
        "unentangled-wall",
        "\
The control: the product of the entangled state's two marginal \
packets, each 1.04 wide with wavenumber ∓3, so the same x₁ and x₂ \
distributions at the start and no entanglement, with the same wall \
on particle 2. Now particle 1 goes the same way with the wall and \
without, up to the time stepping, which keeps a product only nearly \
a product.",
        false,
        true,
        wall_expectations(0.0, 0.01),
    )
}
