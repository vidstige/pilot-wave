//! Oblique view: the grid as a lit liquid surface, droplets as spheres.
use std::thread;

use crate::bath::Vec2;
use crate::grid::Grid;
use crate::image::{Color, Image};
use crate::topdown::DROPLET_COLORS;

type Vec3 = [f64; 3];

/// Peak wave height, in grid half widths.
const RELIEF: f64 = 0.07;
/// Waves fade out between these radii, in grid half widths.
const FADE: (f64, f64) = (0.75, 1.0);
/// Camera angle from vertical, in radians.
const TILT: f64 = 1.05;
/// Camera distance in grid half widths, and tan of half the view angle.
const DISTANCE: f64 = 4.0;
const LENS: f64 = 0.17;
/// In wavelengths, like the grid.
const DROPLET_RADIUS: f64 = 0.09;
const KEY: Vec3 = [-0.7, 0.5, 0.5];
const KEY_COLOR: Vec3 = [12.0, 11.0, 9.5];
const FILL: f64 = 0.35;
const LIQUID: Vec3 = [0.02, 0.07, 0.12];
const ZENITH: Vec3 = [0.02, 0.03, 0.07];
const HORIZON: Vec3 = [1.3, 1.05, 0.8];
const SHININESS: f64 = 60.0;
const EXPOSURE: f64 = 1.5;
const GOLDEN: f64 = 0.618_033_988_749_895;

fn add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale(a: Vec3, s: f64) -> Vec3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn normalize(a: Vec3) -> Vec3 {
    scale(a, 1.0 / dot(a, a).sqrt())
}

fn smoothstep(a: f64, b: f64, x: f64) -> f64 {
    let t = ((x - a) / (b - a)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

fn bilinear(grid: &Grid, p: Vec2) -> f64 {
    let m = (grid.n - 2) as f64;
    let (x, y) = grid.cell(p);
    let (x, y) = (x.clamp(0.0, m), y.clamp(0.0, m));
    let (i, j) = (x as usize, y as usize);
    let (fx, fy) = (x - i as f64, y - j as f64);
    let row = |j| grid.at(i, j) + (grid.at(i + 1, j) - grid.at(i, j)) * fx;
    row(j) + (row(j + 1) - row(j)) * fy
}

/// Surface height at p, scaled by k and faded out towards the edge.
fn elevation(grid: &Grid, k: f64, p: Vec2) -> f64 {
    let d = [p[0] - grid.center[0], p[1] - grid.center[1]];
    let r = (d[0] * d[0] + d[1] * d[1]).sqrt() / grid.half;
    k * bilinear(grid, p) * (1.0 - smoothstep(FADE.0, FADE.1, r))
}

fn normal(grid: &Grid, k: f64, p: Vec2) -> Vec3 {
    let e = grid.half / grid.n as f64;
    let h = |dx, dy| elevation(grid, k, [p[0] + dx, p[1] + dy]);
    normalize([h(-e, 0.0) - h(e, 0.0), h(0.0, -e) - h(0.0, e), 2.0 * e])
}

/// Distance along the ray to the surface, marching one cell at a time.
fn march(grid: &Grid, k: f64, o: Vec3, d: Vec3) -> f64 {
    let top = RELIEF * grid.half;
    let dt = 2.0 * grid.half / grid.n as f64 / d[0].hypot(d[1]);
    let gap = |t: f64| {
        o[2] + t * d[2] - elevation(grid, k, [o[0] + t * d[0], o[1] + t * d[1]])
    };
    let mut t = (top - o[2]) / d[2];
    let mut g = gap(t);
    loop {
        let next = gap(t + dt);
        if next <= 0.0 {
            return t + dt * g / (g - next);
        }
        t += dt;
        g = next;
    }
}

/// Distance along a ray to a sphere, if it hits.
fn sphere(o: Vec3, d: Vec3, c: Vec3, r: f64) -> Option<f64> {
    let oc = sub(o, c);
    let b = dot(oc, d);
    let disc = b * b - dot(oc, oc) + r * r;
    (disc > 0.0).then(|| -b - disc.sqrt())
}

fn sky(r: Vec3) -> Vec3 {
    let t = r[2].max(0.0).sqrt();
    add(scale(HORIZON, 1.0 - t), scale(ZENITH, t))
}

/// Lambert and Blinn-Phong from the key light, sky reflection by Fresnel.
fn shade(n: Vec3, v: Vec3, base: Vec3) -> Vec3 {
    let l = normalize(KEY);
    let h = normalize(add(l, v));
    let diffuse = dot(n, l).max(0.0) + FILL * (0.5 + 0.5 * n[2]);
    let specular = dot(n, h).max(0.0).powf(SHININESS);
    let fresnel = 0.04 + 0.96 * (1.0 - dot(n, v).max(0.0)).powi(5);
    let r = sub(scale(n, 2.0 * dot(n, v)), v);
    let mirror = add(sky(r), scale(KEY_COLOR, specular));
    add(scale(base, diffuse * (1.0 - fresnel)), scale(mirror, fresnel))
}

fn linear(c: Color) -> Vec3 {
    c.map(|x| (x as f64 / 255.0).powf(2.2))
}

/// Tone map and gamma encode.
fn srgb(c: Vec3) -> Color {
    let tone = |x: f64| 1.0 - (-EXPOSURE * x).exp();
    c.map(|x| (tone(x).powf(1.0 / 2.2) * 255.0).round() as u8)
}

/// Sample k of n within a pixel: stratified in x, golden ratio in y.
fn subpixel(k: usize, n: usize) -> (f64, f64) {
    let k = k as f64;
    ((k + 0.5) / n as f64, (0.5 + k * GOLDEN).fract())
}

fn trace(grid: &Grid, k: f64, droplets: &[Vec3], o: Vec3, d: Vec3) -> Vec3 {
    let v = scale(d, -1.0);
    let r = DROPLET_RADIUS;
    let hit = droplets.iter().zip(DROPLET_COLORS.iter().cycle())
        .filter_map(|(&c, &color)| Some((sphere(o, d, c, r)?, c, color)))
        .min_by(|a, b| a.0.total_cmp(&b.0));
    let t = march(grid, k, o, d);
    match hit {
        Some((s, c, color)) if s < t => {
            let p = add(o, scale(d, s));
            let n = scale(sub(p, c), 1.0 / r);
            shade(n, v, linear(color))
        }
        _ => {
            let p = [o[0] + t * d[0], o[1] + t * d[1]];
            shade(normal(grid, k, p), v, LIQUID)
        }
    }
}

/// Oblique view of the grid with droplets at the given positions,
/// averaging samples rays per pixel.
pub fn render(
    grid: &Grid,
    droplets: &[Vec2],
    width: usize,
    height: usize,
    samples: usize,
) -> Image {
    let k = RELIEF * grid.half / grid.peak().max(f64::MIN_POSITIVE);
    let spheres: Vec<Vec3> = droplets.iter()
        .map(|&p| [p[0], p[1], elevation(grid, k, p) + DROPLET_RADIUS])
        .collect();
    let forward = [0.0, TILT.sin(), -TILT.cos()];
    let up = [0.0, TILT.cos(), TILT.sin()];
    let target = [grid.center[0], grid.center[1], 0.0];
    let eye = add(target, scale(forward, -DISTANCE * grid.half));
    let s = 2.0 * LENS / height as f64;
    let ray = |x: f64, y: f64| {
        let sx = (x - 0.5 * width as f64) * s;
        let sy = (0.5 * height as f64 - y) * s;
        normalize(add(forward, add([sx, 0.0, 0.0], scale(up, sy))))
    };
    let pixel = |i: usize| {
        let (x, y) = ((i % width) as f64, (i / width) as f64);
        let c = (0..samples)
            .map(|i| subpixel(i, samples))
            .map(|(dx, dy)| trace(grid, k, &spheres, eye, ray(x + dx, y + dy)))
            .fold([0.0; 3], add);
        srgb(scale(c, 1.0 / samples as f64))
    };
    let mut pixels = vec![[0; 3]; width * height];
    let threads = thread::available_parallelism().map_or(1, |n| n.get());
    let band = height.div_ceil(threads) * width;
    thread::scope(|s| {
        for (b, chunk) in pixels.chunks_mut(band).enumerate() {
            let pixel = &pixel;
            s.spawn(move || {
                for (i, c) in chunk.iter_mut().enumerate() {
                    *c = pixel(b * band + i);
                }
            });
        }
    });
    Image { width, height, pixels }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bump(amplitude: f64) -> Grid {
        let mut grid = Grid::zeros([0.0; 2], 1.0, 40);
        for j in 0..grid.n {
            for i in 0..grid.n {
                let [x, y] = grid.world(i, j);
                let h = amplitude * (-(x * x + y * y) / 0.05).exp();
                grid.values[j * grid.n + i] = h;
            }
        }
        grid
    }

    #[test]
    fn bump_is_brighter_on_the_side_facing_the_key_light() {
        let grid = bump(1.0);
        let k = RELIEF / grid.peak();
        let up = [0.0, 0.0, 1.0];
        let lit = |x| shade(normal(&grid, k, [x, 0.0]), up, LIQUID);
        assert!(lit(-0.2).iter().sum::<f64>() > lit(0.2).iter().sum());
    }

    #[test]
    fn image_is_independent_of_wave_amplitude() {
        // A power of two scales the heights exactly.
        let a = render(&bump(1.0), &[], 32, 32, 4);
        let b = render(&bump(4.0), &[], 32, 32, 4);
        assert_eq!(a.pixels, b.pixels);
    }

    #[test]
    fn droplet_appears_where_it_sits() {
        let grid = Grid::zeros([0.0; 2], 1.0, 40);
        let calm = render(&grid, &[], 32, 32, 4);
        let drop = render(&grid, &[[0.0; 2]], 32, 32, 4);
        let at = |x: usize, y: usize| y * 32 + x;
        assert_ne!(calm.pixels[at(16, 16)], drop.pixels[at(16, 16)]);
        assert_eq!(calm.pixels[at(0, 0)], drop.pixels[at(0, 0)]);
    }
}
