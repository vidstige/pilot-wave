//! An RGB raster with the few drawing primitives the views need.
pub type Color = [u8; 3];

pub struct Image {
    pub width: usize,
    pub height: usize,
    pub pixels: Vec<Color>,
}

impl Image {
    pub fn new(width: usize, height: usize, fill: Color) -> Image {
        Image { width, height, pixels: vec![fill; width * height] }
    }

    pub fn set(&mut self, x: i64, y: i64, color: Color) {
        if x >= 0 && y >= 0 && (x as usize) < self.width && (y as usize) < self.height {
            self.pixels[y as usize * self.width + x as usize] = color;
        }
    }

    pub fn disc(&mut self, cx: f64, cy: f64, radius: f64, color: Color) {
        let r = radius.ceil() as i64;
        for dy in -r..=r {
            for dx in -r..=r {
                if ((dx * dx + dy * dy) as f64) <= radius * radius {
                    self.set(cx.round() as i64 + dx, cy.round() as i64 + dy, color);
                }
            }
        }
    }
}

/// Linear blend of two colors, t in 0..1.
pub fn mix(a: Color, b: Color, t: f64) -> Color {
    let lerp = |x: u8, y: u8| (x as f64 + (y as f64 - x as f64) * t).round() as u8;
    [lerp(a[0], b[0]), lerp(a[1], b[1]), lerp(a[2], b[2])]
}
