use std::env;
use std::process::ExitCode;

use pilot_wave::engine::{self, Engine};
use pilot_wave::experiment::Experiment;
use pilot_wave::experiments;
use pilot_wave::plane::Vec2;
use pilot_wave::report::{self, FRAME_EVERY};
use pilot_wave::run;

const USAGE: &str = "pilot-wave list | describe <experiment> | run <experiment> [runs] [out]
  --engine=<name>  the physics that moves the particles, quantum unless given";
const ENGINE: &str = "quantum";
const SEED: u64 = 1;

fn main() -> ExitCode {
    let (flags, args): (Vec<String>, Vec<String>) = env::args().partition(|a| a.starts_with("--"));
    let engine = flags.iter().find_map(|f| f.strip_prefix("--engine=")).unwrap_or(ENGINE);
    let word = |k: usize| args.get(k).map_or("", String::as_str);
    match (word(1), experiments::named(word(2))) {
        ("describe", Some(e)) => {
            report::describe(&e, engine, build(engine, &e, 0).dt());
            ExitCode::SUCCESS
        }
        ("run", Some(e)) => {
            let runs = args.get(3).map_or(e.runs, |s| s.parse().unwrap());
            let out = args.get(4).map_or("out", String::as_str);
            let passed = run_and_report(&e, engine, runs, out);
            if passed { ExitCode::SUCCESS } else { ExitCode::FAILURE }
        }
        ("list", _) => {
            list();
            ExitCode::SUCCESS
        }
        _ => {
            eprintln!("{USAGE}");
            ExitCode::FAILURE
        }
    }
}

fn list() {
    println!("Engines: {}\n\nExperiments:", engine::NAMES.join(", "));
    for e in experiments::all() {
        println!("  {:<14} {}", e.name, e.summary.split(". ").next().unwrap_or(""));
    }
}

fn build(name: &str, e: &Experiment, runs: usize) -> Box<dyn Engine> {
    let starts = run::starts(e, runs, SEED);
    engine::named(name, e, starts).unwrap_or_else(|| panic!("no engine named {name}"))
}

fn run_and_report(e: &Experiment, name: &str, runs: usize, out: &str) -> bool {
    std::fs::create_dir_all(out).expect("create output directory");
    let mut engine = build(name, e, runs);
    let dt = engine.dt();
    report::describe(e, name, dt);
    let (per_frame, per_point) = ((FRAME_EVERY / dt) as usize, report::trail_stride(dt));
    let mut trails: Vec<Vec<Vec2>> = Vec::new();
    let mut peak = 0.0;
    let outcome = run::run(e, &mut *engine, |engine, step| {
        if step == 0 {
            peak = report::peak(engine);
            trails = engine.positions().iter().take(300).map(|p| vec![*p]).collect();
        } else if step % per_point == 0 {
            for (trail, p) in trails.iter_mut().zip(engine.positions()) {
                trail.push(*p);
            }
        }
        if step % per_frame == 0 {
            let path = format!("{out}/{}_{:04}.png", e.name, step / per_frame);
            report::frame(e, engine, peak, &path);
        }
    });
    report::write(e, &outcome, &trails, out);
    report::verdicts(e, &outcome)
}
