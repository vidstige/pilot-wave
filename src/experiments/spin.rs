//! Spin without a magnetic field: a spin-1/2 particle released from a
//! box into a channel, its arrival times at a line down the channel
//! depending on which way its spin points (Das and Dürr 2019).
use crate::experiment::{
    Expectation, Experiment, Measure, Potential, Realisation, Source, Space, Wall, ELECTRON_MASS,
};

/// The plane is |x|, |y| ≤ 51.2; the channel runs up it, with the
/// sponge 49 past the arrival line. Closer, the slow tail's long
/// wavelengths reflect off the sponge, and arrivals drop by a third.
const CELLS: usize = 1024;
const DX: f64 = 0.1;
/// The closed end of the channel, the box's length and the harmonic
/// confinement across the channel, ½ω²x², held at CAP far out.
const END: f64 = -40.0;
const BOX: f64 = 4.0;
const OMEGA: f64 = 4.0;
const CAP: f64 = 50.0;
/// The arrival line, ten box lengths down the channel.
const SCREEN: f64 = END + 10.0 * BOX;
const DURATION: f64 = 120.0;
/// Spin normal to the plane leaves no arrivals after this time.
const LAST: f64 = 90.0;

const REALISATION: Realisation = Realisation {
    particle: "electron in a waveguide, one unit a quarter of the box",
    mass: ELECTRON_MASS,
    length: 12.5e-6,
    reference: "Das and Dürr 2019, Sci. Rep. 9, 2242: released from a box 50 μm long into a \
cylindrical waveguide; this is the plane version, their spin along the axis being our spin in \
the plane and their spin across it our spin normal to the plane",
};

fn channel(
    name: &'static str,
    summary: &'static str,
    spin: [f64; 3],
    arrived: f64,
    after: f64,
    mean: f64,
) -> Experiment {
    Experiment {
        name,
        summary,
        realisation: REALISATION,
        space: Space::Plane,
        cells: CELLS,
        dx: DX,
        walls: vec![Wall::Band { axis: 1, from: END - 0.4, to: END }],
        potentials: vec![Potential::Harmonic { axis: 0, omega: OMEGA, cap: CAP }],
        source: Source::Channel { from: END, length: BOX, omega: OMEGA },
        spin,
        particles: 20000,
        duration: DURATION,
        screen: Some(SCREEN),
        theory: None,
        expectations: vec![
            Expectation {
                what: "fraction of the particles arriving by t = 120",
                measure: Measure::Arrived,
                expected: arrived,
                tolerance: 0.01,
            },
            Expectation {
                what: "fraction arriving after t = 90",
                measure: Measure::ArrivedAfter { time: LAST },
                expected: after,
                tolerance: 0.01,
            },
            Expectation {
                what: "mean arrival time",
                measure: Measure::MeanArrivalTime,
                expected: mean,
                tolerance: 1.0,
            },
        ],
    }
}

pub fn spin_in_plane() -> Experiment {
    channel(
        "spin-in-plane",
        "\
A particle sits in the ground state of a box 4 long at the closed end \
of a channel along y, confined across it by ½ω²x² with ω = 4, so \
|ψ|² is 0.35 wide. At t = 0 the box's far wall is removed and the \
wave spreads down the channel; a line 40 further on, ten box \
lengths, records when each particle first crosses it; behind it the \
sponge takes the wave away. In the de \
Broglie-Bohm theory of the Pauli equation the velocity is the \
current over the density plus Holland's term ∇×(ρ ŝ)/2ρ, where ŝ is \
the spin direction. With the spin in the plane that term has no \
component in the plane, the trajectories are those of a spinless \
particle, and the arrival times are the quantum flux through the \
line: a main lobe and a long tail of slow arrivals falling as t⁻⁴, \
a tenth of them still on their way at t = 120.",
        [0.0, 1.0, 0.0],
        0.90,
        0.10,
        60.5,
    )
}

pub fn spin_normal() -> Experiment {
    channel(
        "spin-normal",
        "\
The same release with the spin normal to the plane. Nothing in the \
Hamiltonian knows about spin and |ψ|² evolves exactly as before, since \
∇×(ρ ŝ) has no divergence; but Holland's term now lies in the plane, \
(∂ᵧρ, −∂ₓρ)/2ρ, and the trajectories change. Across the channel the \
density is Gaussian, so the term speeds particles on one side of the \
axis down the channel and slows those on the other, by ωx; along the \
channel it drifts the particles sideways where the density changes \
behind the front, so the slow tail is carried to the fast side. The \
arrivals pinch off: all have arrived by t = 90, where the spinless \
tail goes on past 120, and the arrival-time distribution, the one \
thing measured here, depends on the spin although no field acts on \
it.",
        [0.0, 0.0, 1.0],
        1.0,
        0.0,
        55.0,
    )
}
