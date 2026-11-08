//! The field as an image: hue from the phase, brightness from the
//! density; walls grey; particles and their trails on top. For two
//! particles on a line, the line itself with both particles on it.
use std::f64::consts::TAU;

use crate::colormap::{sample, VIRIDIS};
use crate::grid::{Grid, Vec2};
use crate::image::{mix, Color, Image};
use crate::pair;
use crate::wave::Wave;

const WALL: Color = [70, 70, 76];
const PARTICLE: Color = [255, 255, 255];
const TRAIL: Color = [255, 220, 160];
/// Brightness spans this many decades of density below the reference
/// peak, so faint fringes still show and frames stay comparable.
const DECADES: f32 = 3.0;
/// Particle 1 and 2 on the line, the one run followed, the diagonal
/// x₁ = x₂ and the line view's background.
pub const ONE: Color = [255, 140, 50];
pub const TWO: Color = [70, 170, 255];
const RUN: Color = [255, 60, 200];
const DIAGONAL: Color = [150, 150, 150];
const BACKGROUND: Color = [12, 12, 16];
/// Share of the line view's height the marginals take, the rest
/// going to the two rows of dots, and how strongly a marginal tints
/// what is under it.
const CURVES: f64 = 0.6;
const TINT: f64 = 0.45;
/// The top of the line view's bands, over the reference peak.
const HEADROOM: f64 = 3.0;

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

/// The diagonal x₁ = x₂ of a configuration-space image, dotted.
pub fn diagonal(image: &mut Image) {
    for x in (0..image.width as i64).step_by(4) {
        image.set(x, image.height as i64 - 1 - x, DIAGONAL);
    }
}

/// The one run to follow in configuration space: its trail and a
/// ring where it is now.
pub fn run(image: &mut Image, wave: &Wave, trail: &[Vec2]) {
    for p in trail {
        let (x, y) = pixel(wave, image.width, *p);
        image.disc(x as f64, y as f64, 1.0, RUN);
    }
    let (x, y) = pixel(wave, image.width, *trail.last().unwrap());
    image.disc(x as f64, y as f64, 4.5, RUN);
    image.disc(x as f64, y as f64, 2.5, [0; 3]);
}

/// Cells of the line where the particle along the axis is in a wall
/// wherever the other one is.
fn walls(wave: &Wave, axis: usize) -> Vec<bool> {
    let n = wave.n;
    let cell = |i, j| if axis == 0 { (i, j) } else { (j, i) };
    let walled = |i| (0..n).all(|j| wave.keep[wave.index(cell(i, j))] == 0.0);
    (0..n).map(walled).collect()
}

/// Each particle's walls grey in its row and above the baseline, then
/// bands up to each marginal density, scaled to the given peak.
fn marginals(image: &mut Image, wave: &Wave, baseline: usize, peak: f64) {
    let n = wave.n;
    let row = (image.height - baseline) / 2;
    let rows = [0..baseline, baseline..baseline + row, baseline + row..image.height];
    let walls = [walls(wave, 0), walls(wave, 1)];
    let marginals = pair::marginals(wave);
    for x in 0..image.width {
        let i = x * n / image.width;
        for axis in (0..2).filter(|&a| walls[a][i]) {
            for y in rows[0].clone().chain(rows[axis + 1].clone()) {
                image.set(x as i64, y as i64, WALL);
            }
        }
        for (m, color) in marginals.iter().zip([ONE, TWO]) {
            let h = ((m[i] / (HEADROOM * peak)).min(1.0) * (baseline - 1) as f64) as usize;
            for y in baseline - h..baseline {
                let k = y * image.width + x;
                image.pixels[k] = mix(image.pixels[k], color, TINT);
            }
            image.set(x as i64, (baseline - 1 - h) as i64, color);
        }
    }
}

/// Two particles on a line: their marginal densities as bands, and
/// below them x₁ and x₂ of every point as a row of dots each, spread
/// over the row's height to show how they crowd, the first point's
/// large: the one run to follow.
pub fn line(wave: &Wave, points: &[Vec2], width: usize, height: usize, peak: f64) -> Image {
    let mut image = Image::new(width, height, BACKGROUND);
    let baseline = (height as f64 * CURVES) as usize;
    marginals(&mut image, wave, baseline, peak);
    let row = (height - baseline) / 2;
    let x = |v: f64| (v / wave.extent() + 0.5) * width as f64;
    let y = |axis: usize, k: usize| (baseline + axis * row + k * 7919 % row) as f64;
    for (k, p) in points.iter().enumerate() {
        for (axis, color) in [ONE, TWO].into_iter().enumerate() {
            image.set(x(p[axis]) as i64, y(axis, k) as i64, color);
        }
    }
    if let Some(p) = points.first() {
        for axis in 0..2 {
            let middle = (baseline + axis * row + row / 2) as f64;
            image.disc(x(p[axis]), middle, 3.5, RUN);
        }
    }
    image
}

/// The density as a heightfield for the oblique view, with a square
/// root so the packet's tail shows.
pub fn relief(wave: &Wave, cells: usize) -> Grid {
    let mut grid = Grid::zeros([0.0; 2], wave.extent() / 2.0, cells);
    for j in 0..cells {
        for i in 0..cells {
            let k = wave.index((i * wave.n / cells, j * wave.n / cells));
            grid.values[j * cells + i] = (wave.density(k) as f64).sqrt();
        }
    }
    grid
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

    #[test]
    fn a_wall_only_particle_one_feels_is_only_on_its_line() {
        let mut wave = Wave::new(32, 0.1);
        for k in 0..32 * 32 {
            if wave.position((k % 32, k / 32))[0] > 1.0 {
                wave.keep[k] = 0.0;
            }
        }
        let count = |walls: Vec<bool>| walls.iter().filter(|w| **w).count();
        assert_eq!(count(walls(&wave, 0)), 6);
        assert_eq!(count(walls(&wave, 1)), 0);
    }
}
