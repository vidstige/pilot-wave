//! Top-down view: grid values through a colormap, droplets on top.
use crate::colormap::{sample, Colormap};
use crate::grid::{Grid, Vec2};
use crate::image::{Color, Image};

pub const DROPLET_RADIUS: f64 = 4.0;
pub const DROPLET_COLORS: [Color; 3] =
    [[255, 255, 255], [255, 160, 30], [110, 200, 255]];

/// Paint the grid with values from lo to hi spanning the colormap,
/// each cell as a square of scale by scale pixels.
pub fn paint(grid: &Grid, map: &Colormap, lo: f64, hi: f64, scale: usize) -> Image {
    let side = grid.n * scale;
    let mut image = Image::new(side, side, [0; 3]);
    for y in 0..side {
        for x in 0..side {
            let (i, j) = (x / scale, grid.n - 1 - y / scale);
            let t = (grid.at(i, j) - lo) / (hi - lo);
            image.set(x as i64, y as i64, sample(map, t));
        }
    }
    image
}

pub fn pixel(grid: &Grid, p: Vec2) -> (f64, f64) {
    let (x, y) = grid.cell(p);
    (x, grid.n as f64 - 1.0 - y)
}

pub fn droplets(image: &mut Image, grid: &Grid, trails: &[Vec<Vec2>]) {
    for (trail, color) in trails.iter().zip(DROPLET_COLORS.iter().cycle()) {
        for &p in trail {
            let (x, y) = pixel(grid, p);
            image.set(x.round() as i64, y.round() as i64, *color);
        }
        let (x, y) = pixel(grid, *trail.last().unwrap());
        image.disc(x, y, DROPLET_RADIUS + 1.0, [0; 3]);
        image.disc(x, y, DROPLET_RADIUS, *color);
    }
}
