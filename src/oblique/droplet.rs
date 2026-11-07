//! Droplets: dark glossy beads resting in dimples of the bath.
use super::{add, dot, fresnel, normalize, reflect, scale, sky, sub};
use super::{Vec3, KEY, KEY_COLOR};
use crate::bath::Vec2;

/// In wavelengths, like the grid.
pub const RADIUS: f64 = 0.09;
/// How deep a droplet sinks below the surface under it, in radii.
pub const SINK: f64 = 0.05;
/// Depth and Gaussian width of the dimple under a droplet, in radii.
const DIMPLE: (f64, f64) = (0.3, 0.8);
/// Brightness of the oil's identity color under the gloss.
const BODY: f64 = 0.04;
const GLOSS: f64 = 120.0;

pub struct Droplet {
    pub center: Vec3,
    /// Identity color, linear.
    pub tint: Vec3,
}

/// Depth at p of the dimple under a droplet at c.
pub fn dimple(c: Vec2, p: Vec2) -> f64 {
    let r2 = (p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2);
    let w = DIMPLE.1 * RADIUS;
    DIMPLE.0 * RADIUS * (-0.5 * r2 / (w * w)).exp()
}

/// Distance along a ray to a droplet, if it hits.
fn sphere(o: Vec3, d: Vec3, c: Vec3) -> Option<f64> {
    let oc = sub(o, c);
    let b = dot(oc, d);
    let disc = b * b - dot(oc, oc) + RADIUS * RADIUS;
    (disc > 0.0).then(|| -b - disc.sqrt()).filter(|&t| t > 0.0)
}

/// The nearest droplet along a ray, and the distance to it.
fn hit(droplets: &[Droplet], o: Vec3, d: Vec3) -> Option<(f64, &Droplet)> {
    droplets.iter()
        .filter_map(|x| Some((sphere(o, d, x.center)?, x)))
        .min_by(|a, b| a.0.total_cmp(&b.0))
}

/// Fraction of the hemisphere above p, with normal n, the droplet
/// covers: the solid angle of the sphere weighted by cosine.
fn cover(x: &Droplet, p: Vec3, n: Vec3) -> f64 {
    let d = sub(x.center, p);
    let d2 = dot(d, d);
    dot(n, d).max(0.0) * RADIUS * RADIUS / (d2 * d2.sqrt())
}

/// Ambient occlusion at p with normal n, from 0 open to 1 covered.
pub fn occlusion(droplets: &[Droplet], p: Vec3, n: Vec3) -> f64 {
    droplets.iter().map(|x| cover(x, p, n)).sum::<f64>().min(1.0)
}

/// What a droplet mirrors in direction r: the sky with the key in it,
/// and below the horizon the sky again, mirrored by a flat bath.
fn environment(r: Vec3) -> Vec3 {
    let f = if r[2] < 0.0 { fresnel(-r[2]) } else { 1.0 };
    let glint = dot(r, normalize(KEY)).max(0.0).powf(GLOSS);
    let mirror = scale(sky([r[0], r[1], r[2].abs()]), f);
    add(mirror, scale(KEY_COLOR, glint))
}

/// A dark glossy droplet at normal n seen from v: its oil's color,
/// glowing towards the bottom where the oil focuses the sky, under a
/// Fresnel mirror.
fn shade(x: &Droplet, n: Vec3, v: Vec3) -> Vec3 {
    let f = fresnel(dot(n, v));
    let body = scale(x.tint, BODY * (1.0 - n[2]) * (1.0 - f));
    add(body, scale(environment(reflect(n, v)), f))
}

/// Distance to and color of the nearest droplet along a ray.
pub fn trace(droplets: &[Droplet], o: Vec3, d: Vec3) -> Option<(f64, Vec3)> {
    let (t, x) = hit(droplets, o, d)?;
    let n = scale(sub(add(o, scale(d, t)), x.center), 1.0 / RADIUS);
    Some((t, shade(x, n, scale(d, -1.0))))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn resting() -> Droplet {
        Droplet { center: [0.0, 0.0, RADIUS], tint: [1.0; 3] }
    }

    #[test]
    fn occlusion_fades_with_distance_from_the_contact() {
        let x = resting();
        let up = [0.0, 0.0, 1.0];
        let at = |r| cover(&x, [r, 0.0, 0.0], up);
        assert!((at(0.0) - 1.0).abs() < 1e-9);
        assert!(at(RADIUS) > at(2.0 * RADIUS));
        assert!(at(2.0 * RADIUS) > at(4.0 * RADIUS));
    }
}
