# pilot-wave

A deterministic pilot-wave (de Broglie-Bohm) simulation in two
dimensions, written from scratch in Rust with no dependencies. A
wavefunction on a square grid is stepped by the Schrödinger equation,
and particles with definite positions ride its phase gradient.
Positions drawn from |ψ|² at the start stay distributed as |ψ|², so an
ensemble of deterministic trajectories reproduces the quantum
statistics. Experiments are plain data: a stage, a source, a duration
and a list of expectations, each a named measure with an expected
value and a tolerance, checked after every run.

![Single slit: the field, the trails, the build-up and the screen](images/single-slit.png)

## Running it

```
cargo run --release -- list
cargo run --release -- describe single-slit
cargo run --release -- run single-slit [particles] [out]
cargo test --release
```

`describe` prints the experiment's summary, its physical realisation
with the unit conversions, the grid, and the expectations. `run`
prints the same, runs, and ends with a verdict per expectation; it
exits non-zero if any fails. It also writes frames
`<out>/<experiment>_NNNN.png` every half time unit (hue is the phase,
brightness the density), the trails of 300 particles over the
time-integrated density, the screen as a plot and a table, and the
arrivals in order. The figures come from

```
cargo run --release -- run single-slit 40000 out/single-slit
.venv/bin/python tools/figures.py single-slit
```

with a Python venv holding numpy and matplotlib at `.venv`.

## Units

The engine uses ħ = m = 1 and a packet wavelength of 1, so a packet's
wavenumber and speed are both 2π and its kinetic energy 2π². Each
experiment names a physical realisation, a particle and its de Broglie
wavelength λ, and the units follow: length λ, time mλ²/ħ, speed
ħ/(mλ), energy ħ²/(mλ²). For a 600 eV electron (λ = 50 pm) one unit
of length is 50 pm, of time 2.16·10⁻¹⁷ s, of speed 2.3·10⁶ m/s and of
energy 30.5 eV; `describe` prints these for every experiment.

## The model

Each step applies the Hamiltonian, a fourth order Laplacian plus the
potential, to the imaginary part to advance the real part, then the
reverse, in one fused sweep per band of rows on all cores
(`wave.rs`). Hard walls are cells kept at zero and an absorbing sponge
at the edges swallows what reaches it. The guiding velocity is the
probability current over the density, each interpolated between the
surrounding cells before dividing, with fourth order gradients, and
particles move in substeps of at most a quarter cell, so they turn
back from the near-nodes of standing waves instead of crossing them
(`particles.rs`). Positions are drawn from |ψ|² on the grid
(`sample.rs`). A screen bins the particles' crossings and the
time-integrated current through the line (`screen.rs`).

An experiment (`experiment.rs`) is a `Realisation`, a grid, walls, a
source packet, a particle count, a duration, an optional screen with
the closed form to compare it with, and expectations. The runner
(`run.rs`) builds the field, samples the particles, steps both, and
measures; `report.rs` prints and writes.

## Experiments

### Single slit

A packet 6 wide across and 2 long heads up at speed 2π towards a wall
0.4 thick with one slit 2 wavelengths wide; a screen 30 further on
counts what comes through. Every particle is an independent run of the
one-particle experiment, since all share the field and never interact,
so one run with 40000 particles is 40000 runs with one, and the
arrivals are recorded in order so the pattern can be watched building
up one particle at a time, as in the electron build-up films.

Realisation: a 600 eV electron, as in Bach et al. (2013), whose slit
was 62 nm, over a thousand wavelengths, so their pattern sits at angles
a thousand times smaller than ours; the shape, sinc² with a cos²
obliquity factor, is the same.

| expectation | expected | measured |
|---|---|---|
| fraction through the screen, the beam's share across the slit | 0.133 ± 0.010 | 0.1264 |
| arrived fraction minus the probability through the screen | 0 ± 0.005 | −0.0007 |
| rms gap between the current and Fraunhofer, of the peak | 0 ± 0.03 | 0.0177 |
| rms gap between the arrivals and the current, of the peak | 0 ± 0.03 | 0.0257 |
| current at the first zero of sinc², x = 17.2, of the peak | 0 ± 0.02 | 0.0095 |

The arrived fraction sits 5% under the geometric share because the
slit's edges diffract some of the wave back; the current's first
minimum falls to 1% of the peak at the predicted place.

### Double slit

The same beam meets a wall with two slits 1 wide and 4 apart, centre
to centre: cos²(πd sinθ) fringes under the sinc²(πa sinθ) envelope,
orders at sinθ = m/4. Each particle goes through one slit or the
other and never crosses the axis of symmetry, yet the ensemble lands
on the fringes, which build up one particle at a time as in the
electron films. The realisation is again the 600 eV electron of Bach
et al. (2013), whose slits were 62 nm wide and 272 nm apart: d/a = 4.4
against our 4, so about as many fringes fit under the envelope, but
theirs are 0.18 mrad apart against our 0.25 rad.

![Double slit: the field, the trails, the build-up and the screen](images/double-slit.png)

| expectation | expected | measured |
|---|---|---|
| fraction through the screen, the beam's share across both slits | 0.126 ± 0.010 | 0.1184 |
| arrived fraction minus the probability through the screen | 0 ± 0.005 | −0.0011 |
| rms gap between the current and Fraunhofer, of the peak | 0 ± 0.035 | 0.0282 |
| rms gap between the arrivals and the current, of the peak | 0 ± 0.03 | 0.0225 |
| x of the first-order maximum, d sinθ = 1 | 7.7 ± 0.8 | 7.6 |
| current at the first minimum, d sinθ = 1/2, of the peak | 0 ± 0.1 | 0.019 |
