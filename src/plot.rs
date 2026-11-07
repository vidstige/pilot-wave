//! Curves on a dark background, each scaled to its own peak.
use crate::image::{Color, Image};

const BACKGROUND: Color = [12, 12, 16];
const MARGIN: usize = 8;

/// Draw each series as a polyline, scaled to its own peak so shapes
/// can be compared.
pub fn curves(width: usize, height: usize, series: &[(&[f64], Color)]) -> Image {
    let mut image = Image::new(width, height, BACKGROUND);
    let inner = (width - 2 * MARGIN, height - 2 * MARGIN);
    for (values, color) in series {
        let peak = values.iter().cloned().fold(0.0, f64::max).max(1e-300);
        let point = |i: usize| {
            let x = MARGIN as f64 + i as f64 * inner.0 as f64 / (values.len() - 1).max(1) as f64;
            let y = MARGIN as f64 + inner.1 as f64 * (1.0 - values[i] / peak);
            (x, y)
        };
        for i in 1..values.len() {
            line(&mut image, point(i - 1), point(i), *color);
        }
    }
    image
}

fn line(image: &mut Image, a: (f64, f64), b: (f64, f64), color: Color) {
    let steps = (b.0 - a.0).abs().max((b.1 - a.1).abs()).ceil().max(1.0) as usize;
    for s in 0..=steps {
        let t = s as f64 / steps as f64;
        let (x, y) = (a.0 + t * (b.0 - a.0), a.1 + t * (b.1 - a.1));
        image.set(x.round() as i64, y.round() as i64, color);
    }
}
