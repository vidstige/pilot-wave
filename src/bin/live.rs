use std::time::Instant;

use eframe::egui;
use egui::{ColorImage, TextureHandle, TextureOptions};
use pilot_wave::fps::Fps;
use pilot_wave::image::Image;
use pilot_wave::live::{Live, View, SCENES};

const SCENE: &str = "pair";
const MEMORY: f64 = 20.0;
const MEMORY_RANGE: (f64, f64) = (1.0, 400.0);
const BOUNCES: usize = 10;
const MAX_BOUNCES: usize = 200;

struct App {
    live: Live,
    bounces: usize,
    view: View,
    paused: bool,
    fps: Fps,
    last: Instant,
    texture: TextureHandle,
}

fn color_image(image: &Image) -> ColorImage {
    let size = [image.width, image.height];
    ColorImage::from_rgb(size, image.pixels.as_flattened())
}

impl App {
    fn new(ctx: &egui::Context) -> App {
        let live = Live::new(SCENE, MEMORY);
        let image = color_image(&live.frame(View::Top));
        App {
            live,
            bounces: BOUNCES,
            view: View::Top,
            paused: false,
            fps: Fps::default(),
            last: Instant::now(),
            texture: ctx.load_texture("frame", image, TextureOptions::LINEAR),
        }
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        let (mut name, mut memory) = (self.live.name, self.live.memory);
        egui::ComboBox::from_label("scene")
            .selected_text(name)
            .show_ui(ui, |ui| {
                for scene in SCENES {
                    ui.selectable_value(&mut name, scene, scene);
                }
            });
        let range = MEMORY_RANGE.0..=MEMORY_RANGE.1;
        let slider = egui::Slider::new(&mut memory, range);
        ui.add(slider.logarithmic(true).text("memory"));
        self.live.select(name, memory);
        let slider = egui::Slider::new(&mut self.bounces, 1..=MAX_BOUNCES);
        ui.add(slider.text("bounces per frame"));
        ui.horizontal(|ui| {
            ui.selectable_value(&mut self.view, View::Top, "top");
            ui.selectable_value(&mut self.view, View::Oblique, "oblique");
        });
        ui.horizontal(|ui| {
            ui.toggle_value(&mut self.paused, "pause");
            if ui.button("reset").clicked() {
                self.live = Live::new(name, memory);
            }
        });
        ui.separator();
        ui.label(format!("bounce {}", self.live.bounces));
        ui.label(format!("{:.1} fps", self.fps.rate()));
    }

    fn picture(&self, ui: &mut egui::Ui) {
        let side = ui.available_size().min_elem();
        let size = egui::vec2(side, side);
        ui.centered_and_justified(|ui| ui.image((self.texture.id(), size)));
    }
}

impl eframe::App for App {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        let now = Instant::now();
        self.fps.tick((now - self.last).as_secs_f64());
        self.last = now;
        egui::Panel::left("controls").show(ui, |ui| self.controls(ui));
        if !self.paused {
            self.live.advance(self.bounces);
        }
        let image = color_image(&self.live.frame(self.view));
        self.texture.set(image, TextureOptions::LINEAR);
        egui::CentralPanel::default_margins().show(ui, |ui| self.picture(ui));
        ui.ctx().request_repaint();
    }
}

fn main() -> eframe::Result {
    eframe::run_native(
        "pilot-wave",
        eframe::NativeOptions::default(),
        Box::new(|cc| Ok(Box::new(App::new(&cc.egui_ctx)))),
    )
}
