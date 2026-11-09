//! Spin without a magnetic field: a spin-1/2 particle released from a
//! box into a channel, its arrival times at a line down the channel
//! depending on which way its spin points (Das and Dürr 2019).
use std::f64::consts::PI;

use crate::experiment::{
    band, Expectation, Experiment, Measure, Particle, Stage, ELECTRON,
};
use crate::wave::{Vec2, Wave};

/// The plane is |x|, |y| ≤ 51.2; the channel runs up it, with the
/// sponge 49 past the arrival line. Closer, the slow tail's long
/// wavelengths reflect off the sponge, and arrivals drop by a third.
/// One unit is a quarter of the paper's 50 μm box.
const CELLS: usize = 1024;
const DX: f64 = 0.1;
const SCALE: f64 = 12.5e-6;
/// The closed end of the channel, the box's length and the harmonic
/// confinement across the channel, ½ω²x², held at CAP far out where
/// the state is negligible so that the time step stays large.
const END: f64 = -40.0;
const BOX: f64 = 4.0;
const OMEGA: f64 = 4.0;
const CAP: f64 = 50.0;
/// The arrival line, ten box lengths down the channel.
const SCREEN: f64 = END + 10.0 * BOX;
const DURATION: f64 = 120.0;
/// Spin normal to the plane leaves no arrivals after this time.
const LAST: f64 = 90.0;

/// Width of |ψ|² in the box's ground state along it and across it,
/// L √(1/12 − 1/2π²) and 1/√(2ω); the state itself is written by
/// `boxed`.
const WIDTHS: [f64; 2] = [0.354, 0.723];

/// The electron at rest in the box, with the given spin.
const fn electron(spin: [f64; 3]) -> Particle {
    Particle { kind: ELECTRON, position: [0.0, END + BOX / 2.0], velocity: [0.0; 2], width: WIDTHS, spin }
}

fn confinement(p: Vec2) -> f64 {
    (0.5 * OMEGA * OMEGA * p[0] * p[0]).min(CAP)
}

/// The ground state of the box at the closed end of the channel,
/// Gaussian across it; the box's far wall is removed at t = 0.
fn boxed(wave: &mut Wave) {
    wave.fill(|[x, y]| {
        let along = if (END..END + BOX).contains(&y) { (PI * (y - END) / BOX).sin() } else { 0.0 };
        [along * (-0.5 * OMEGA * x * x).exp(), 0.0]
    })
}

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
        reference: "Das and Dürr 2019, Sci. Rep. 9, 2242: an electron released from a box 50 μm long into a cylindrical waveguide; this is the plane version, their spin along the axis being our spin in the plane and their spin across it our spin normal to the plane",
        stage: Stage {
            cells: CELLS,
            dx: DX,
            scale: SCALE,
            walls: vec![band(1, END - 0.4, END)],
            potential: confinement,
        },
        particles: vec![electron(spin)],
        state: Some(boxed),
        runs: 20000,
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
sponge takes the wave away. In the de Broglie-Bohm theory of the \
Pauli equation the velocity is the current over the density plus \
Holland's term ∇×(ρ ŝ)/2ρ, where ŝ is the spin direction. With the \
spin in the plane that term has no component in the plane, the \
trajectories are those of a spinless particle, and the arrival times \
are the quantum flux through the line: a main lobe and a long tail of \
slow arrivals falling as t⁻⁴, a tenth of them still on their way at \
t = 120.",
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
(∂ᵧρ, −∂ₓρ)/2ρ, a circulation along the contours of the density, and \
the particles loop as they go: up the channel on one side of the axis \
and down on the other, by ωx, and across it where the density changes \
along the channel, so the slow tail behind the front is carried to \
the fast side. The arrivals pinch off: all have arrived by t = 90, \
where the spinless tail goes on past 120, and the arrival-time \
distribution, the one thing measured here, depends on the spin \
although no field acts on it.",
        [0.0, 0.0, 1.0],
        1.0,
        0.0,
        55.0,
    )
}
