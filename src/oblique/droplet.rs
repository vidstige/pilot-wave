//! Droplets as spheres resting on the bath.
use super::{dot, sub, Vec3};
use crate::bath::Vec2;

/// In wavelengths, like the grid.
pub const RADIUS: f64 = 0.09;
/// How deep a droplet sinks below the surface under it, in radii.
pub const SINK: f64 = 0.05;
/// Depth and Gaussian width of the dimple under a droplet, in radii.
const DIMPLE: (f64, f64) = (0.3, 0.8);

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
