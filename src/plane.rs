//! The plane: positions and rectangles in it. No I/O.

pub type Vec2 = [f64; 2];
/// A rectangle [[x₀, x₁], [y₀, y₁]].
pub type Rect = [[f64; 2]; 2];

pub fn inside(rect: &Rect, p: Vec2) -> bool {
    (rect[0][0]..rect[0][1]).contains(&p[0]) && (rect[1][0]..rect[1][1]).contains(&p[1])
}

/// The rectangle of everything with the coordinate along the axis
/// between `from` and `to`.
pub fn band(axis: usize, from: f64, to: f64) -> Rect {
    let mut rect = [[f64::NEG_INFINITY, f64::INFINITY]; 2];
    rect[axis] = [from, to];
    rect
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_band_spans_the_other_axis() {
        let wall = band(0, 5.0, f64::INFINITY);
        assert!(inside(&wall, [6.0, -1e9]) && !inside(&wall, [4.9, 0.0]));
    }
}
