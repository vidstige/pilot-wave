//! Rates of events, such as frames, smoothed over recent time. No I/O.

/// Seconds over which older frames fade by a factor e.
const WINDOW: f64 = 1.0;

/// Events and seconds both decay with time, so their ratio is a
/// rate weighted towards recent frames, exact from the first frame
/// when steady.
#[derive(Default)]
pub struct Rate {
    events: f64,
    seconds: f64,
}

impl Rate {
    /// Count n events in a frame that took dt seconds.
    pub fn tick(&mut self, n: f64, dt: f64) {
        let keep = (-dt / WINDOW).exp();
        self.events = keep * self.events + n;
        self.seconds = keep * self.seconds + dt;
    }

    /// Events per second, zero before the first frame.
    pub fn rate(&self) -> f64 {
        if self.seconds > 0.0 {
            self.events / self.seconds
        } else {
            0.0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn frames(times: &[f64]) -> Rate {
        let mut rate = Rate::default();
        for &dt in times {
            rate.tick(1.0, dt);
        }
        rate
    }

    #[test]
    fn steady_frames_give_their_rate_from_the_first_frame() {
        assert_eq!(frames(&[0.02]).rate(), 50.0);
        assert!((frames(&[0.02; 10]).rate() - 50.0).abs() < 1e-9);
    }

    #[test]
    fn one_slow_frame_dents_the_rate_and_heals() {
        let mut times = vec![0.02; 10];
        times.push(0.2);
        let dented = frames(&times).rate();
        times.extend([0.02; 500]);
        let healed = frames(&times).rate();
        assert!(10.0 < dented && dented < 50.0, "{dented}");
        assert!((healed - 50.0).abs() < 0.1, "{healed}");
    }

    #[test]
    fn sparse_events_average_over_frames() {
        let mut rate = Rate::default();
        for i in 0..256 {
            rate.tick(if i % 8 == 0 { 1.0 } else { 0.0 }, 1.0 / 64.0);
        }
        assert!((rate.rate() - 8.0).abs() < 1.0, "{}", rate.rate());
    }
}
