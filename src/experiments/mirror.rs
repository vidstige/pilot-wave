//! The free packet meeting a wall: the sanity check of reflection.
use crate::experiment::{band, Expectation, Experiment, Measure};
use crate::experiments::free::{centre, flight, DURATION};

/// The wall the packet meets, its thickness, and the screen behind it.
const WALL: f64 = 5.0;
const THICKNESS: f64 = 0.4;
const BEHIND: f64 = 8.0;

pub fn mirror() -> Experiment {
    flight(
        "mirror",
        "\
The packet of the free flight meets a hard wall across the plane at \
y = 5 and comes back. After a time of 4 its mean y should be the \
mirror image of the free flight's in the wall, with the spread across \
the motion as in free flight, and no particle beyond it.",
        vec![band(1, WALL, WALL + THICKNESS)],
        BEHIND,
        vec![
            Expectation {
                what: "mean y of the particles at the end, mirrored in the wall",
                measure: Measure::Mean { axis: 1 },
                expected: 2.0 * WALL - centre(DURATION),
                tolerance: 0.25,
            },
            Expectation {
                what: "fraction of the particles beyond the wall",
                measure: Measure::Arrived,
                expected: 0.0,
                tolerance: 0.0001,
            },
        ],
    )
}
