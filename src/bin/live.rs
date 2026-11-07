use eframe::egui;
use egui::{ColorImage, TextureHandle, TextureOptions};
use pilot_wave::rate::Rate;
use pilot_wave::image::Image;
use pilot_wave::live::{Live, View};
use pilot_wave::scene::{default_memory, SCENES};

const SCENE: &str = "pair";
const MEMORY_RANGE: (f64, f64) = (1.0, 400.0);
/// Replay speed in bounces per second.
const SPEED: f64 = 200.0;
const SPEED_RANGE: (f64, f64) = (1.0, 3000.0);

struct App {
    live: Live,
    speed: f64,
    view: View,
    paused: bool,
    fps: Rate,
    bounce_rate: Rate,
    texture: TextureHandle,
}

fn color_image(image: &Image) -> ColorImage {
    let size = [image.width, image.height];
    ColorImage::from_rgb(size, image.pixels.as_flattened())
}

impl App {
    fn new(ctx: &egui::Context) -> App {
        let live = Live::new(SCENE, default_memory(SCENE));
        let image = color_image(&live.frame(View::Top));
        App {
            live,
            speed: SPEED,
            view: View::Top,
            paused: false,
            fps: Rate::default(),
            bounce_rate: Rate::default(),
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
        if name != self.live.name {
            memory = default_memory(name);
        }
        let range = MEMORY_RANGE.0..=MEMORY_RANGE.1;
        let slider = egui::Slider::new(&mut memory, range);
        ui.add(slider.logarithmic(true).text("memory"));
        self.live.select(name, memory);
        let range = SPEED_RANGE.0..=SPEED_RANGE.1;
        let slider = egui::Slider::new(&mut self.speed, range);
        ui.add(slider.logarithmic(true).integer().text("bounces per second"));
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
        ui.label(format!(
            "{:.1} fps, {:.0} bounces per second",
            self.fps.rate(),
            self.bounce_rate.rate()
        ));
    }

    fn picture(&self, ui: &mut egui::Ui) {
        let side = ui.available_size().min_elem();
        let size = egui::vec2(side, side);
        ui.centered_and_justified(|ui| ui.image((self.texture.id(), size)));
    }
}

impl eframe::App for App {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        // We repaint every frame, so egui's stable dt is the wall clock.
        let dt = ui.ctx().input(|i| i.stable_dt) as f64;
        egui::Panel::left("controls").show(ui, |ui| self.controls(ui));
        let bounces = if self.paused {
            0
        } else {
            self.live.play(self.speed, dt)
        };
        self.fps.tick(1.0, dt);
        self.bounce_rate.tick(bounces as f64, dt);
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
