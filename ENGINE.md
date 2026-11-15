# The engines

What an engine is, how each of the two works, the decisions behind
them, and which knobs an experiment turns against which are fixed in
the code.

## What an engine is

An engine is the physics of a run, and the only part of the program
that knows it. `src/engine.rs` holds the whole contract:

    trait Engine {
        fn dt(&self) -> f64;                    // the time a step takes
        fn step(&mut self);
        fn positions(&self) -> &[Vec2];         // every run's particle
        fn field(&self) -> Option<[&[f32]; 2]>; // for the pictures
    }

An engine is built from an `Experiment` and the runs' starting
positions and owns everything else it needs: its grid, its step, its
threads. The field is two numbers to each cell of the stage, drawn as
brightness (their squares) and hue (their phase), and no measure
reads it.

Around the engines, and knowing none of them: `src/experiment.rs` is
what an experiment says, `src/sample.rs` draws the starting
positions, `src/run.rs` steps the engine and keeps the detector,
`src/screen.rs` is the detection line, `src/measure.rs` answers the
expectations, `src/view.rs` and `src/report.rs` draw and write.
`src/engines/` holds the engines: `quantum.rs` on `wave.rs` and
`guided.rs`, and `classical.rs`. A new engine is a file there and a line in
`engine::named`; it passes when `run <experiment> --engine=<name>`
does.

What every engine is given is the same: the stage's walls and
potential, the particle's position, velocity and width, and the
wave at the start as a function of position, `Experiment::amplitude`,
which an engine with no wave may ignore. The runs' starts are drawn
from the square of that function.

## The quantum engine

A wavefunction ψ on an n × n grid of square cells of side dx, and a
set of points moving in the same plane, each an independent run of
the experiment. Each step advances ψ by the
Schrödinger equation with ħ = m = 1,

    i ∂ψ/∂t = −½ ∇²ψ + V ψ,

and moves every point along the guiding velocity

    v = Im(ψ* ∇ψ) / |ψ|²,

the probability current over the density. Positions drawn
from |ψ|² at the start stay distributed as |ψ|² for all time, so an
ensemble of these deterministic trajectories reproduces the quantum
statistics, and one particle's trajectory is one run of the
experiment.

The particles never interact: they share the field, and the field
does not know they are there. N interacting particles would share
one wavefunction over their joint configuration space, of dimension
N times the dimensions each moves in, and the grid holds two of those
dimensions. So the plane can also be read as the configuration space
(x₁, x₂) of two particles on a line, with a potential U(x₁ − x₂) as
their interaction and a point of the plane as the pair's two
positions; experiments of that kind were done here and set aside (see
the README). Two particles both free in the plane would need a
four-dimensional ψ, which is a different engine.

### The field

**Scheme.** Real and imaginary parts are stored separately in `f32`
and stepped alternately, the imaginary part living half a step ahead
(Visscher 1991): R ← R + dt·H·I, then I ← I − dt·H·R. The scheme is
explicit, second order in time, conserves the staggered density
R² + I₋I₊ exactly, and is stable while dt < 2 / (largest eigenvalue
of H). An implicit scheme (Crank-Nicolson) would allow larger steps
but needs a linear solve every step; the explicit scheme's step is one
sweep over the grid and runs on all cores with no synchronisation but
the band halos.

**Stencil.** The Laplacian is the fourth order five-point stencil
along each axis, [−1/12, 4/3, −5/2, 4/3, −1/12]/dx². Its dispersion
is known in closed form: at ten cells per wavelength, the usual
k = 2π on dx = 0.1, a wave moves 0.5% slower than k and a packet
spreads along its motion as if 2.5% heavier. The `free` experiment's
tolerances along the motion, 1% of the flight, allow for that. The
second order stencil was off by 3% in speed, enough to move
the double slit's fringes. The stencil's largest eigenvalue is
16/(3dx²) per axis, which sets the step.

**Time step.** `stable_dt(dx, Vmax) = 0.5 · 2 / (16/(3dx²) + Vmax)`:
half the stability limit, with the potential's largest magnitude on
the grid added to the kinetic bound. For dx = 0.1 and no potential
that is dt = 1.875·10⁻³. An experiment never chooses dt; it follows
from dx and the potential, and a potential that is large somewhere
the state never goes (the harmonic confinement of the channel far
from its axis) should be capped so it does not shrink the step.

**Walls.** A hard wall is a cell held at zero: every step multiplies
each cell by a `keep` factor, 0 in walls, 1 elsewhere. This is the
Dirichlet condition ψ = 0 on the wall's cells, so the wall is as
rough as the grid; walls are drawn as rectangles in world units and
whichever cells have their centres inside are wall. Walls hold −0.0
and the phase there is undefined; nothing reads the phase in a wall.

**Sponge.** The outermost 24 cells on every side are an absorbing
layer: `keep = exp(−40 s² dt)` with s rising from 0 to 1 over the
layer, so what reaches the edge dies away instead of reflecting. It
absorbs wavelengths of a few cells well and long wavelengths badly: a
slow tail with λ ≈ 10 partly reflects. The sponge is set
before the walls, so a wall running into the edge stays whole.

**Threads.** Rows are split into bands and the bands handed out
through a queue to `std::thread::scope` threads, two bands per thread
so a fast core takes over from a slow one, with halos of two rows
(for R) and four (for I) copied around each band. One 1024² step takes
0.45 ms on an M1; a 512² step 0.125 ms. There is no GPU path: the grid
sizes here are small enough that the CPU, 150 steps a second on
1024², keeps every run under three minutes.

### The particles

**Guiding velocity.** The current Im(ψ*∇ψ) and the density |ψ|² are
computed at the four cells around the point with the same fourth
order gradient as the field, each bilinearly interpolated to the
point, and only then divided. This matters next to a node: the
velocity there must soar, so that a particle slows and turns before
the node, and dividing after interpolating keeps that. Interpolating
a per-cell velocity, or using the phase difference between cells
(bounded by π/2dx), let particles cross the near-nodes of standing
waves: the barrier's particle transmission came out 0.39 against the
field's 0.22, and with this form 0.218. Finer substeps did not fix
it; the fix was spatial.

**Substeps.** A particle moves by midpoint steps, each shortened so
that it moves at most a quarter of a cell, up to 64 substeps per
field step; the last substep takes whatever time is left. Next to a
node the velocity is huge and the particle creeps up to it instead of
jumping across. Two known limits: a node that lies exactly along a
line of cells lets one run in twenty thousand step over it, and at
ten cells per wavelength the ensemble spreads along its motion 2.5%
more than the field does, which vanishes at twenty cells per
wavelength. That one is the step's doing rather than
the cells': the real and imaginary parts are half a step apart where
the guiding velocity reads them, and a step four times finer on the
same cells removes it.

**Sampling.** Each run's starting position is drawn from the square
of the wave at the start, on the stage's grid: a cell by its share,
then uniformly within it. This is the runner's doing and the same for
every engine. The seed is 1 for every experiment, so a run is
reproducible.

**Threads.** Particles step in chunks of at least 2000 on all cores,
each thread reading the field.

**Screen.** A horizontal line at the experiment's `screen` height.
Each particle's first upward crossing is binned along x and logged
with its time. The screen is the runner's too.

**Exact positions.** Every particle's position is known, so a run
could be asked things no laboratory can: whether the particles that
tunnelled were the ones at the front of the packet, as they are here,
or which slit each went through. The experiments do not ask, since
another engine need not agree and no detector could tell.

## The classical engine

`--engine=classical` is Newton's particles and no wave
(`src/engines/classical.rs`), here to show what the experiments do
not let through.

**Motion.** A particle moves an axis at a time. Where the potential
rises over the move it pays in speed along that axis, v² − 2ΔV, and
where it cannot pay, or the move would end in a wall, it turns back.
Its energy is kept exactly, over a step in the potential as on a
slope, which a force taken from the potential's gradient does not do
at a barrier's edge: a first version with kicks lost 4% of the speed
crossing one. The price is that a particle exactly at rest along an
axis never starts along it.

**The source.** The runs start where every engine's do. Their
velocities are the particle's plus a Gaussian share of the beam's
divergence, ½ over the packet's width along each axis: the spread a
beam that narrow is seen to have, taken as a property of the source,
since without it not even the free flight would spread.

**What it passes.** The free flight, the mirror and the relaxation,
where bouncing about the box mixes the runs as well as the wave does.
It fails the barrier (0.008 of the runs have the energy, against
0.230) and both slits (shadows of the openings, no lobes or fringes).

## Decisions

- **Engines behind one trait, experiments that name none.** The
  experiments are the fixed point and the physics is what is tried
  against them, so nothing outside `src/engines/` may depend on what
  an engine is made of: a wave is not needed to place the runs, to
  detect them or to judge them.
- **Deterministic trajectories on a grid field** for the quantum
  engine, not a path-memory model: the droplet model of the earlier
  project reproduced some quantum-like statistics qualitatively, the
  Bohmian field reproduces them exactly, double slit included.
- **Zero dependencies.** The solver, PNG writer, plotting, random
  numbers and threading are written here; `std::thread::scope` is
  enough for the parallelism, and the only tool with a dependency is
  the figure script, in Python.
- **`f32` fields, `f64` particles.** The field is memory bound and
  half the bytes are half the time; the guided particles are few and
  their arithmetic near nodes wants the precision.
- **Hard walls as zeroed cells** rather than a large potential: a
  potential step a hundred high made ψ fall off within a cell and
  shrank the step. A steep ramp potential is the alternative where a
  wall must be soft.
- **First crossings on the screen**, so arrival-time distributions are
  the Bohmian first-arrival ones.
- **Expect only what a detector records.** No field is ever observed,
  so the expectations are on the particles: where they land, how
  many, when. The closed form is drawn beside the arrivals for the
  eye. No measure reads an engine's field or where a run started.
- **Tolerances for the physics, not for a solver.** An expectation is
  the closed form's value, never one corrected for a grid; where ten
  cells per wavelength cost half a percent the tolerance says so.
- **Experiments as data** with the source and potential as named
  functions in the experiment's file rather than enums of shapes: the
  file reads top to bottom as the setup, and nothing in the framework
  has to know what a bump or a channel is.

## Knobs an experiment turns

Every field of `Experiment`, grouped as the definition is:

| field | what it sets |
|---|---|
| `name`, `summary`, `reference` | what it is, what should happen and why, and the paper or closed form it is checked against |
| `stage.cells`, `stage.dx` | the grid: n × n cells of side dx, the plane |x|, |y| ≤ n·dx/2 |
| `stage.scale` | the SI length of one unit, which with the particle's mass gives the units |
| `stage.walls` | rectangles of hard wall, built with `slits` and `band` or written out |
| `stage.potential` | V(x, y), a plain function, fixed for the run; its largest magnitude on the grid bounds dt |
| `particle` | the particle of a run: kind (name, mass), position, velocity, packet width |
| `state` | the wave at the start as a function of position, when it is not the particle's Gaussian packet |
| `start` | a rectangle every other run starts in, evenly, instead of at the source: what `Misplaced` compares |
| `runs` | how many times the one-particle run is repeated, each from its own start |
| `duration` | how long to step, in engine time |
| `screen` | height of the detection line, if any |
| `theory` | the closed form along the screen the arrivals are compared with, if any |
| `expectations` | named measures with expected values and tolerances |

`dt` is the engine's, derived, not chosen. Every experiment uses
dx = 0.1, ten cells per wavelength for the k = 2π beam; the sanity
checks run on 512 cells, the box of the relaxation on 128, the slits
and the channel on 1024. The engine is chosen when running, not by
the experiment: `--engine=<name>`.

## Knobs fixed in the code

| constant | value | where | why |
|---|---|---|---|
| `STENCIL` | fourth order, five points | `engines/wave.rs` | 0.5% speed error at ten cells per wavelength |
| `SPONGE`, `SPONGE_RATE` | 24 cells, rate 40 | `wave.rs` | absorbs short wavelengths without reflection; long ones partly reflect |
| step safety | 0.5 of the stability limit | `stable_dt` | margin for the potential and the sponge |
| `CELLS_PER_THREAD` | 65536 | `wave.rs` | below this a thread costs more to start than it saves |
| `BANDS_PER_THREAD` | 2 | `wave.rs` | lets fast cores take over bands from slow ones |
| `SUBSTEP`, `MAX_SUBSTEPS` | ¼ cell, 64 | `engines/guided.rs` | creeping up to nodes |
| `PARTICLES_PER_THREAD` | 2000 | `engines/guided.rs` | thread start-up cost |
| `classical::STEP` | 0.002 | `engines/classical.rs` | a beam's particle moves an eighth of a cell |
| seed | 1 | `main.rs` | reproducible runs |
| `BINS` | 128 across the plane | `run.rs` | screen resolution: 0.8 wide on a plane of 102.4 |
| `GLOW_EVERY` | 10 steps | `run.rs` | how often the engine's field is added to the glow under the trails |
| `DECADES` | 3 | `view.rs` | the density range the frames' brightness spans |
| `PIXELS`, `FRAME_EVERY`, `TRAILS`, `DOTS` | 800, 0.5, 300, 2000 | `report.rs` | what a run writes |
| `TRAIL_STEP` | 0.002 | `report.rs` | time between the points of a trail, whatever the engine's step |

## What is not modelled

- No spin and no magnetic field. Both were done earlier, the spin as
  Holland's term in the guiding velocity and the Aharonov-Bohm flux
  line as a sign-flipping cut in the grid, and can come back when an
  experiment done in a laboratory asks for them.
- No time-dependent potentials or walls: a barrier that is removed at
  t = 0 is written as an initial state that does not extend past it.
- No more than one particle to a run, so no interacting particles in
  the plane; two on a line fit the quantum engine's plane as
  configuration space, see above.
- No particle acting on its field: the wave does not know the
  particles are there.
- No measurement apparatus: a measurement is reading a particle's
  position, and a detector that disturbs ψ would be a wall or a
  potential in the experiment.
