//! Droplets as spheres resting on the bath.
use super::{add, dot, fresnel, mul, normalize, reflect, scale};
use super::{sky, smoothstep, sub, Vec3};
use super::{FILL, KEY, KEY_COLOR, LIQUID};
use crate::bath::Vec2;

/// In wavelengths, like the grid.
pub const RADIUS: f64 = 0.09;
/// How deep a droplet sinks below the surface under it, in radii.
pub const SINK: f64 = 0.05;
/// Depth and Gaussian width of the dimple under a droplet, in radii.
const DIMPLE: (f64, f64) = (0.3, 0.8);
/// Refractive index of silicone oil.
const IOR: f64 = 1.4;
/// Fraction of light passing through, and how much identity color
/// it picks up.
const CLARITY: f64 = 0.15;
const TINT: f64 = 0.75;
const GLOSS: f64 = 120.0;
/// Angular radius of the key light, as a tangent.
const PENUMBRA: f64 = 0.25;

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

/// Fraction of the key light reaching p past the droplet, softened
/// by how far behind it p lies.
fn shadow(x: &Droplet, p: Vec3) -> f64 {
    let d = sub(x.center, p);
    let along = dot(d, normalize(KEY));
    let miss = (dot(d, d) - along * along).sqrt();
    let w = PENUMBRA * along;
    if along > 0.0 { smoothstep(RADIUS - w, RADIUS + w, miss) } else { 1.0 }
}

/// Fraction of the key light reaching p past all droplets.
pub fn sunlight(droplets: &[Droplet], p: Vec3) -> f64 {
    droplets.iter().map(|x| shadow(x, p)).product()
}

/// What a droplet mirrors: the sky, and below the horizon a flat lit
/// bath mirroring the sky.
fn environment(r: Vec3) -> Vec3 {
    let f = if r[2] < 0.0 { fresnel(-r[2]) } else { 1.0 };
    let lit = normalize(KEY)[2] + FILL;
    let mirror = scale(sky([r[0], r[1], r[2].abs()]), f);
    add(mirror, scale(LIQUID, lit * (1.0 - f)))
}

/// Direction d bent into a surface with normal n facing it, for a
/// ratio eta of refractive indices.
fn refract(d: Vec3, n: Vec3, eta: f64) -> Vec3 {
    let c = -dot(n, d);
    let k = 1.0 - eta * eta * (1.0 - c * c);
    add(scale(d, eta), scale(n, eta * c - k.sqrt()))
}

/// Direction a view ray from v leaves the droplet after entering at
/// normal n. It exits where the chord meets the sphere again.
fn through(n: Vec3, v: Vec3) -> Vec3 {
    let d = refract(scale(v, -1.0), n, 1.0 / IOR);
    let exit = sub(n, scale(d, 2.0 * dot(n, d)));
    refract(d, scale(exit, -1.0), IOR)
}

/// The environment and the key as seen in direction r.
fn light(r: Vec3, specular: f64) -> Vec3 {
    let glint = dot(r, normalize(KEY)).max(0.0).powf(specular);
    add(environment(r), scale(KEY_COLOR, glint))
}

/// A dark glassy droplet at normal n seen from v: the scene behind it
/// refracted, upside down, through tinted oil, under a Fresnel mirror.
fn shade(x: &Droplet, n: Vec3, v: Vec3) -> Vec3 {
    let tint = add(scale([1.0; 3], 1.0 - TINT), scale(x.tint, TINT));
    let body = scale(mul(tint, light(through(n, v), GLOSS)), CLARITY);
    let mirror = light(reflect(n, v), GLOSS);
    let f = fresnel(dot(n, v));
    add(scale(body, 1.0 - f), scale(mirror, f))
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

    /// Ground point at distance r from the contact, towards the key.
    fn towards_key(r: f64) -> Vec3 {
        let l = normalize([KEY[0], KEY[1], 0.0]);
        [l[0] * r, l[1] * r, 0.0]
    }

    #[test]
    fn shadow_falls_away_from_the_key_light() {
        let x = resting();
        let r = 1.7 * RADIUS;
        assert!(shadow(&x, towards_key(-r)) < 0.1);
        assert_eq!(shadow(&x, towards_key(r)), 1.0);
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
