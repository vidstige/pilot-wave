//! The plane: positions and rectangles in it. No I/O.

pub type Vec2 = [f64; 2];
/// A rectangle [[x₀, x₁], [y₀, y₁]].
pub type Rect = [[f64; 2]; 2];

pub fn inside(rect: &Rect, p: Vec2) -> bool {
    (rect[0][0]..rect[0][1]).contains(&p[0]) && (rect[1][0]..rect[1][1]).contains(&p[1])
}
