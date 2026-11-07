//! Droplets as spheres resting on the bath.
use super::{dot, sub, Vec3};

/// In wavelengths, like the grid.
pub const RADIUS: f64 = 0.09;

pub struct Droplet {
    pub center: Vec3,
    /// Identity color, linear.
    pub tint: Vec3,
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
