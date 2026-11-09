use std::env;
use std::process::ExitCode;

use pilot_wave::experiments;
use pilot_wave::report::{self, FRAME_EVERY};
use pilot_wave::run;
use pilot_wave::wave::Vec2;

const USAGE: &str = "pilot-wave list | describe <experiment> | run <experiment> [runs] [out]";

fn main() -> ExitCode {
    let args: Vec<String> = env::args().collect();
    let name = args.get(2).map(String::as_str).unwrap_or("");
    match args.get(1).map(String::as_str) {
        Some("list") => {
            for e in experiments::all() {
                println!("{:<14} {}", e.name, e.summary.split(". ").next().unwrap_or(""));
            }
            ExitCode::SUCCESS
        }
        Some("describe") => match experiments::named(name) {
            Some(e) => {
                report::describe(&e);
                ExitCode::SUCCESS
            }
            None => {
                eprintln!("unknown experiment {name}");
                ExitCode::FAILURE
            }
        },
        Some("run") => match experiments::named(name) {
            Some(e) => {
                let runs = args.get(3).map(|s| s.parse().unwrap()).unwrap_or(e.runs);
                let out = args.get(4).map(String::as_str).unwrap_or("out");
                if run_and_report(&e, runs, out) { ExitCode::SUCCESS } else { ExitCode::FAILURE }
            }
            None => {
                eprintln!("unknown experiment {name}");
                ExitCode::FAILURE
            }
        },
        _ => {
            eprintln!("{USAGE}");
            ExitCode::FAILURE
        }
    }
}

fn run_and_report(e: &pilot_wave::experiment::Experiment, runs: usize, out: &str) -> bool {
    std::fs::create_dir_all(out).expect("create output directory");
    report::describe(e);
    let per_frame = (FRAME_EVERY / e.dt()) as usize;
    let mut trails: Vec<Vec<Vec2>> = Vec::new();
    let mut peak = 0.0;
    let outcome = run::run(e, runs, 1, |wave, particles, step| {
        if step == 0 {
            peak = wave.peak();
            trails = particles.positions.iter().take(300).map(|p| vec![*p]).collect();
        } else {
            for (trail, p) in trails.iter_mut().zip(&particles.positions) {
                trail.push(*p);
            }
        }
        if step % per_frame == 0 {
            report::frame(wave, particles, peak, &format!("{out}/{}_{:04}.png", e.name, step / per_frame));
        }
    });
    report::write(e, &outcome, &trails, out);
    report::verdicts(e, &outcome)
}
