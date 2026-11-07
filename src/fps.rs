//! Frame rate smoothed over recent frames. No I/O.

/// Weight of the newest frame time in the running average.
const SMOOTHING: f64 = 0.1;

#[derive(Default)]
pub struct Fps {
    seconds: Option<f64>,
}

impl Fps {
    /// Count a frame that took dt seconds.
    pub fn tick(&mut self, dt: f64) {
        let s = self.seconds.unwrap_or(dt);
        self.seconds = Some(s + SMOOTHING * (dt - s));
    }

    /// Frames per second, zero before the first frame.
    pub fn rate(&self) -> f64 {
        self.seconds.map_or(0.0, |s| 1.0 / s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ticked(times: &[f64]) -> Fps {
        let mut fps = Fps::default();
        for &dt in times {
            fps.tick(dt);
        }
        fps
    }

    #[test]
    fn steady_frames_give_their_rate_from_the_first_frame() {
        assert_eq!(ticked(&[0.02]).rate(), 50.0);
        assert!((ticked(&[0.02; 10]).rate() - 50.0).abs() < 1e-9);
    }

    #[test]
    fn one_slow_frame_dents_the_rate_and_heals() {
        let mut times = vec![0.02; 10];
        times.push(0.2);
        let dented = ticked(&times).rate();
        times.extend([0.02; 100]);
        let healed = ticked(&times).rate();
        assert!(10.0 < dented && dented < 50.0, "{dented}");
        assert!((healed - 50.0).abs() < 0.1, "{healed}");
    }
}
