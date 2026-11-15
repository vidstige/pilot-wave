//! The free beam meeting a wall: reflection.
use crate::experiment::{Expectation, Experiment, Measure, Stage};
use crate::experiments::free::{across, centre, free_flight, DURATION};
use crate::plane::band;

/// The wall the beam meets, its thickness, and the screen behind it.
const WALL: f64 = 5.0;
const THICKNESS: f64 = 0.4;
const BEHIND: f64 = 8.0;

/// What is expected along the motion.
fn along() -> [Expectation; 2] {
    [
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
    ]
}

/// The free flight with a wall in the way and the screen behind it.
pub fn mirror() -> Experiment {
    let free = free_flight();
    Experiment {
        name: "mirror",
        summary: "\
The beam of the free flight meets a hard wall across the plane at \
y = 5. It should come back as if mirrored in the wall, as wide \
across as in free flight, with no particle beyond the wall.",
        reference: "reflection off a hard wall; no laboratory experiment",
        stage: Stage { walls: vec![band(1, WALL, WALL + THICKNESS)], ..free.stage },
        screen: Some(BEHIND),
        expectations: along().into_iter().chain(across()).collect(),
        ..free
    }
}
