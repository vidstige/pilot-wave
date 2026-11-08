//! The field as an image: hue from the phase, brightness from the
//! density; walls grey; particles and their trails on top.
use std::f64::consts::TAU;

use crate::colormap::{sample, VIRIDIS};
use crate::image::{mix, Color, Image};
use crate::wave::{Vec2, Wave};

const WALL: Color = [70, 70, 76];
const PARTICLE: Color = [255, 255, 255];
const TRAIL: Color = [255, 220, 160];
/// Brightness spans this many decades of density below the reference
/// peak, so faint fringes still show and frames stay comparable.
const DECADES: f32 = 3.0;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Style {
    /// Hue from the phase, brightness from the density.
    Phase,
    /// Density alone through a sequential colormap.
    Density,
}

/// A fully saturated colour around the hue wheel, t in 0..1.
fn hue(t: f64) -> Color {
    let channel = |offset: f64| {
        let x = ((t + offset).rem_euclid(1.0) * 6.0 - 3.0).abs() - 1.0;
        (x.clamp(0.0, 1.0) * 255.0).round() as u8
    };
    [channel(0.0), channel(2.0 / 3.0), channel(1.0 / 3.0)]
}

/// Brightness of a density relative to the reference peak.
fn brightness(density: f32, peak: f32) -> f64 {
    (1.0 + (density / peak).max(1e-30).log10() / DECADES).clamp(0.0, 1.0) as f64
}

/// The field at every pixel, densities relative to the given peak.
pub fn field(wave: &Wave, pixels: usize, style: Style, peak: f32) -> Image {
    let mut image = Image::new(pixels, pixels, [0; 3]);
    for y in 0..pixels {
        for x in 0..pixels {
            let i = x * wave.n / pixels;
            let j = (pixels - 1 - y) * wave.n / pixels;
            let k = wave.index((i, j));
            let brightness = brightness(wave.density(k), peak);
            let color = if wave.keep[k] == 0.0 {
                WALL
            } else {
                match style {
                    Style::Phase => {
                        let phase = (wave.im[k].atan2(wave.re[k]) as f64) / TAU + 0.5;
                        mix([0; 3], hue(phase), brightness)
                    }
                    Style::Density => sample(&VIRIDIS, brightness),
                }
            };
            image.set(x as i64, y as i64, color);
        }
    }
    image
}

/// Pixel of a position in an image of the whole domain.
pub fn pixel(wave: &Wave, pixels: usize, p: Vec2) -> (i64, i64) {
    let half = wave.extent() / 2.0;
    let s = pixels as f64 / wave.extent();
    (((p[0] + half) * s) as i64, ((half - p[1]) * s) as i64)
}

/// The diagonal x₁ = x₂ of a configuration space, dotted.
pub fn diagonal(image: &mut Image) {
    for t in (0..image.width).step_by(4) {
        image.set(t as i64, (image.width - 1 - t) as i64, WALL);
    }
}

pub fn particles(image: &mut Image, wave: &Wave, trails: &[Vec<Vec2>]) {
    for trail in trails {
        for p in trail {
            let (x, y) = pixel(wave, image.width, *p);
            image.set(x, y, TRAIL);
        }
        let (x, y) = pixel(wave, image.width, *trail.last().unwrap());
        image.disc(x as f64, y as f64, 1.5, PARTICLE);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn brightness_spans_the_decades_below_the_peak() {
        assert_eq!(brightness(1.0, 1.0), 1.0);
        assert_eq!(brightness(2.0, 1.0), 1.0);
        assert!((brightness(1e-3, 1.0) - 0.0).abs() < 1e-6);
        assert!((brightness(1e-1, 1.0) - 2.0 / 3.0).abs() < 1e-6);
    }

    #[test]
    fn hue_wheel_passes_through_the_primaries_and_wraps() {
        assert_eq!(hue(0.0), [255, 0, 0]);
        assert_eq!(hue(1.0 / 3.0), [0, 255, 0]);
        assert_eq!(hue(2.0 / 3.0), [0, 0, 255]);
        assert_eq!(hue(1.0), hue(0.0));
    }
}
