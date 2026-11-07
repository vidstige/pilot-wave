//! Droplets as spheres resting on the bath.
use super::{add, dot, fresnel, normalize, reflect, scale, sky, sub, Vec3};
use super::{FILL, KEY, KEY_COLOR};
use crate::bath::Vec2;

/// In wavelengths, like the grid.
pub const RADIUS: f64 = 0.09;
/// How deep a droplet sinks below the surface under it, in radii.
pub const SINK: f64 = 0.05;
/// Depth and Gaussian width of the dimple under a droplet, in radii.
const DIMPLE: (f64, f64) = (0.3, 0.8);
/// Albedo of the dark body, and how much identity color it shows.
const BODY: f64 = 0.03;
const TINT: f64 = 0.5;
/// Key light focused through the droplet onto its far side.
const GLOW: f64 = 0.3;
const FOCUS: f64 = 8.0;
const GLOSS: f64 = 300.0;

pub struct Droplet {
    pub center: Vec3,
    /// Identity color, linear.
    pub tint: Vec3,
}

/// Depth of the droplet's dimple at p.
pub fn dimple(x: &Droplet, p: Vec2) -> f64 {
    let r2 = (p[0] - x.center[0]).powi(2) + (p[1] - x.center[1]).powi(2);
    let w = DIMPLE.1 * RADIUS;
    DIMPLE.0 * RADIUS * (-0.5 * r2 / (w * w)).exp()
}

/// Distance along a ray to a droplet, if it hits.
fn sphere(o: Vec3, d: Vec3, c: Vec3) -> Option<f64> {
    let oc = sub(o, c);
    let b = dot(oc, d);
    let disc = b * b - dot(oc, oc) + RADIUS * RADIUS;
    (disc > 0.0).then(|| -b - disc.sqrt())
}

/// The nearest droplet along a ray, and the distance to it.
pub fn hit(droplets: &[Droplet], o: Vec3, d: Vec3) -> Option<(f64, &Droplet)> {
    droplets.iter()
        .filter_map(|x| Some((sphere(o, d, x.center)?, x)))
        .min_by(|a, b| a.0.total_cmp(&b.0))
}

/// What a droplet mirrors: the sky, and below the horizon the bath
/// as a flat mirror of the sky.
fn environment(r: Vec3) -> Vec3 {
    let f = if r[2] < 0.0 { fresnel(-r[2]) } else { 1.0 };
    scale(sky([r[0], r[1], r[2].abs()]), f)
}

/// A dark glassy droplet at normal n seen from v: a tinted body lit
/// and back lit by the key, under a Fresnel mirror.
pub fn shade(x: &Droplet, n: Vec3, v: Vec3) -> Vec3 {
    let l = normalize(KEY);
    let h = normalize(add(l, v));
    let diffuse = dot(n, l).max(0.0) + FILL * (0.5 + 0.5 * n[2]);
    let glow = (-dot(n, l)).max(0.0).powf(FOCUS);
    let tint = add(scale([1.0; 3], 1.0 - TINT), scale(x.tint, TINT));
    let body = scale(tint, BODY * diffuse + GLOW * glow);
    let specular = dot(n, h).max(0.0).powf(GLOSS);
    let mirror = add(environment(reflect(n, v)), scale(KEY_COLOR, specular));
    let f = fresnel(dot(n, v));
    add(scale(body, 1.0 - f), scale(mirror, f))
}
