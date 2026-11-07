use std::env;
use std::time::Instant;

use eframe::egui;
use egui::{Color32, ColorImage, RichText, TextureHandle, TextureOptions, ViewportCommand};
use pilot_wave::experiment::NAMES;
use pilot_wave::image::{Color, Image};
use pilot_wave::live::{Live, Setup, View};
use pilot_wave::particles::Kind;
use pilot_wave::png;
use pilot_wave::rate::Rate;

const SETUP: Setup = Setup { name: "double-slit", count: 500, kind: Kind::Quantum, n: 512 };
const COUNTS: (usize, usize) = (1, 20000);
const GRIDS: [usize; 3] = [256, 512, 1024];
/// Replay speed in simulation time units per second.
const SPEED: f64 = 0.5;
const SPEEDS: (f64, f64) = (0.01, 10.0);
/// Seconds of stepping a frame may spend before the replay falls
/// behind.
const BUDGET: f64 = 0.1;
const PLOT: (usize, usize) = (720, 160);
const WINDOW: [f32; 2] = [1000.0, 860.0];
/// Seconds to run before `--screenshot <path>` captures the window.
const SCREENSHOT_AFTER: f64 = 2.0;

struct App {
    screenshot: Option<String>,
    requested: bool,
    live: Live,
    speed: f64,
    view: View,
    paused: bool,
    fps: Rate,
    step_rate: Rate,
    /// Steps per second of stepping alone, which sets the frame's limit.
    throughput: Rate,
    field: TextureHandle,
    plot: Option<TextureHandle>,
}

fn color_image(image: &Image) -> ColorImage {
    let size = [image.width, image.height];
    ColorImage::from_rgb(size, image.pixels.as_flattened())
}

fn color32(c: Color) -> Color32 {
    Color32::from_rgb(c[0], c[1], c[2])
}

impl App {
    fn new(ctx: &egui::Context, screenshot: Option<String>) -> App {
        let live = Live::new(SETUP);
        let image = color_image(&live.frame(View::Phase));
        App {
            screenshot,
            requested: false,
            live,
            speed: SPEED,
            view: View::Phase,
            paused: false,
            fps: Rate::default(),
            step_rate: Rate::default(),
            throughput: Rate::default(),
            field: ctx.load_texture("field", image, TextureOptions::LINEAR),
            plot: None,
        }
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        let mut setup = self.live.setup;
        egui::ComboBox::from_label("experiment")
            .selected_text(setup.name)
            .show_ui(ui, |ui| {
                for name in NAMES {
                    ui.selectable_value(&mut setup.name, name, name);
                }
            });
        let slider = egui::Slider::new(&mut setup.count, COUNTS.0..=COUNTS.1);
        ui.add(slider.logarithmic(true).text("particles"));
        ui.horizontal(|ui| {
            ui.selectable_value(&mut setup.kind, Kind::Quantum, "quantum");
            ui.selectable_value(&mut setup.kind, Kind::Classical, "classical");
        });
        egui::ComboBox::from_label("grid")
            .selected_text(setup.n.to_string())
            .show_ui(ui, |ui| {
                for n in GRIDS {
                    ui.selectable_value(&mut setup.n, n, n.to_string());
                }
            });
        self.live.select(setup);
        let slider = egui::Slider::new(&mut self.speed, SPEEDS.0..=SPEEDS.1);
        ui.add(slider.logarithmic(true).text("time per second"));
        ui.horizontal(|ui| {
            ui.selectable_value(&mut self.view, View::Phase, "phase");
            ui.selectable_value(&mut self.view, View::Density, "density");
            ui.selectable_value(&mut self.view, View::Oblique, "oblique");
        });
        ui.horizontal(|ui| {
            ui.toggle_value(&mut self.paused, "pause");
            if ui.button("reset").clicked() {
                self.live = Live::new(setup);
            }
        });
        ui.separator();
        ui.label(format!("time {:.2}", self.live.time()));
        ui.label(format!("{:.1} fps", self.fps.rate()));
        ui.label(format!("{:.0} steps per second", self.step_rate.rate()));
        if let Some(arrivals) = self.live.arrivals() {
            ui.label(format!("{arrivals} of {} arrived", setup.count));
            ui.separator();
        }
        for &(text, color) in self.live.legend() {
            ui.label(RichText::new(text).color(color32(color)));
        }
    }

    /// Run this frame's steps, at most as many as fit in the budget
    /// at the measured throughput.
    fn step(&mut self, dt: f64) -> usize {
        let limit = (self.throughput.rate() * BUDGET).max(1.0) as usize;
        let start = Instant::now();
        let steps = if self.paused { 0 } else { self.live.play(self.speed, dt, limit) };
        if steps > 0 {
            self.throughput.tick(steps as f64, start.elapsed().as_secs_f64());
        }
        steps
    }

    /// Ask for the window's pixels once it has run a while, then save
    /// them and quit once they arrive.
    fn screenshot(&mut self, ctx: &egui::Context) {
        let Some(path) = &self.screenshot else { return };
        if !self.requested && ctx.input(|i| i.time) > SCREENSHOT_AFTER {
            ctx.send_viewport_cmd(ViewportCommand::Screenshot(Default::default()));
            self.requested = true;
        }
        let shot = ctx.input(|i| {
            i.events.iter().find_map(|e| match e {
                egui::Event::Screenshot { image, .. } => Some(image.clone()),
                _ => None,
            })
        });
        if let Some(shot) = shot {
            let pixels = shot.pixels.iter().map(|c| [c.r(), c.g(), c.b()]).collect();
            let image = Image { width: shot.size[0], height: shot.size[1], pixels };
            png::write(path, &image).expect("write screenshot");
            ctx.send_viewport_cmd(ViewportCommand::Close);
        }
    }

    /// The field as large a square as fits, the plot as wide below it.
    fn pictures(&self, ui: &mut egui::Ui) {
        let available = ui.available_size();
        let aspect = PLOT.1 as f32 / PLOT.0 as f32;
        let below = self.plot.as_ref().map_or(0.0, |_| aspect);
        let spacing = ui.spacing().item_spacing.y;
        let side = available.x.min((available.y - spacing) / (1.0 + below));
        ui.vertical_centered(|ui| {
            ui.image((self.field.id(), egui::vec2(side, side)));
            if let Some(texture) = &self.plot {
                ui.image((texture.id(), egui::vec2(side, side * aspect)));
            }
        });
    }
}

impl eframe::App for App {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        // We repaint every frame, so egui's stable dt is the wall clock.
        let dt = ui.ctx().input(|i| i.stable_dt) as f64;
        egui::Panel::left("controls").show(ui, |ui| self.controls(ui));
        let steps = self.step(dt);
        self.fps.tick(1.0, dt);
        self.step_rate.tick(steps as f64, dt);
        let image = color_image(&self.live.frame(self.view));
        self.field.set(image, TextureOptions::LINEAR);
        let plot = self.live.plot(PLOT.0, PLOT.1).map(|p| color_image(&p));
        self.plot = plot.map(|p| match self.plot.take() {
            Some(mut texture) => {
                texture.set(p, TextureOptions::LINEAR);
                texture
            }
            None => ui.ctx().load_texture("plot", p, TextureOptions::LINEAR),
        });
        egui::CentralPanel::default_margins().show(ui, |ui| self.pictures(ui));
        self.screenshot(ui.ctx());
        ui.ctx().request_repaint();
    }
}

fn main() -> eframe::Result {
    let args: Vec<String> = env::args().collect();
    let screenshot = (args.get(1).map(String::as_str) == Some("--screenshot"))
        .then(|| args[2].clone());
    let viewport = egui::ViewportBuilder::default().with_inner_size(WINDOW);
    eframe::run_native(
        "pilot-wave",
        eframe::NativeOptions { viewport, ..Default::default() },
        Box::new(|cc| Ok(Box::new(App::new(&cc.egui_ctx, screenshot)))),
    )
}
