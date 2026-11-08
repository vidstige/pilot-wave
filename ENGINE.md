# The engine

How the pilot-wave engine works, the decisions behind it, and which
knobs an experiment turns against which are fixed in the code. The
engine is `src/wave.rs` (the field), `src/particles.rs` (the guided
particles), `src/sample.rs` (drawing starting positions) and `src/screen.rs` (the
detection line);
`src/experiment.rs` is what an experiment may say to it, and
`src/run.rs` drives a run.

## What it computes

A wavefunction ψ on an n × n grid of square cells of side dx, and a
set of points moving in the same plane, each an independent run of
the experiment. Each step advances ψ by the
Schrödinger equation with ħ = m = 1,

    i ∂ψ/∂t = −½ ∇²ψ + V ψ,

and moves every point along the guiding velocity

    v = Im(ψ* ∇ψ) / |ψ|²  +  s (∂ᵧρ, −∂ₓρ) / 2ρ,

the probability current over the density, plus Holland's spin term
when the field carries a spin s normal to the plane. Positions drawn
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

## The field

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
is known in closed form (`group_velocity`, `effective_mass`): at ten
cells per wavelength, the usual k = 2π on dx = 0.1, a wave moves 0.5%
slower than k and a packet spreads along its motion as if 2.5%
heavier. The `free` experiment checks the field against exactly
that. The second order stencil was off by 3% in speed, enough to move
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
slow tail with λ ≈ 10 partly reflects, which is why the spin
channel's arrival line sits 49 before the sponge. The sponge is set
before the walls, so a wall running into the edge stays whole.

**Threads.** Rows are split into bands and the bands handed out
through a queue to `std::thread::scope` threads, two bands per thread
so a fast core takes over from a slow one, with halos of two rows
(for R) and four (for I) copied around each band. One 1024² step takes
0.45 ms on an M1; a 512² step 0.125 ms. There is no GPU path: the grid
sizes here are small enough that the CPU, 150 steps a second on
1024², keeps every run under three minutes.

## The particles

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
more than the field does (`free` measures both), which vanishes at
twenty cells per wavelength.

**Spin.** The field carries a spin s, the component normal to the
plane of the experiment's spin direction, and the current gets
s (∂ᵧρ, −∂ₓρ)/2, the two-dimensional form of Holland's
∇×(ρ ŝ)/2. It has no divergence, so the density evolves exactly as
without it and the position statistics are untouched, but the
trajectories differ: it circulates along the contours of the density.
A spin in the plane contributes nothing in two dimensions.

**Sampling.** Starting positions are drawn from |ψ|² on the grid: a
cell by its probability, then uniformly within it. The seed is 1 for
every run, so a run is reproducible.

**Threads.** Particles step in chunks of at least 2000 on all cores,
each thread reading the field.

**Screen.** A horizontal line at the experiment's `screen` height.
Each particle's first upward crossing is binned along x and logged
with its time; the probability current through the line is
integrated in the same bins every step. The two should agree when
particles cross once; with spin a particle may cross back and forth,
and only its first crossing counts.

**Exact positions.** Since every particle's position is known, the
measures can ask things no laboratory can: whether any particle ever
crossed the axis of symmetry (`CrossedAxis`), whether the particles
that tunnelled were the ones at the front of the packet
(`ArrivedFromFront`), how far each one moved (`Moved`), and whether
the ensemble's spread matches the field's (`Spread`, `FieldSpread`).

## Decisions

- **Deterministic trajectories on a grid field**, not a path-memory
  model: the droplet model of the earlier project reproduced some
  quantum-like statistics qualitatively, the Bohmian field reproduces
  them exactly, double slit included, because ψ is the same ψ.
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
- **Experiments as data** with the source and potential as named
  functions in the experiment's file rather than enums of shapes: the
  file reads top to bottom as the setup, and nothing in the framework
  has to know what a bump or a channel is.

## Knobs an experiment turns

Every field of `Experiment`:

| field | what it sets |
|---|---|
| `cells`, `dx` | the grid: n × n cells of side dx, the plane |x|, |y| ≤ n·dx/2 |
| `walls` | rectangles of hard wall, built with `slits` and `band` or written out |
| `potential` | V(x, y), a plain function; its largest magnitude on the grid bounds dt |
| `source` | writes ψ at t = 0: a `Packet`, or `Wave::fill` with any closed form |
| `spin` | the spin direction; only its component normal to the plane acts |
| `particles` | how many trajectories, each a run |
| `duration` | how long to step, in engine time |
| `screen` | height of the detection line, if any |
| `theory` | the closed form along the screen the current is compared with, if any |
| `realisation` | the particle and the SI length of one unit, for the printed conversions |
| `expectations` | named measures with expected values and tolerances |

`dt` is derived, not chosen. Every experiment uses dx = 0.1, ten
cells per wavelength for the k = 2π beam; the sanity checks run on
512 cells, the slits and the channel on 1024.

## Knobs fixed in the code

| constant | value | where | why |
|---|---|---|---|
| `STENCIL` | fourth order, five points | `wave.rs` | 0.5% speed error at ten cells per wavelength |
| `SPONGE`, `SPONGE_RATE` | 24 cells, rate 40 | `wave.rs` | absorbs short wavelengths without reflection; long ones partly reflect |
| step safety | 0.5 of the stability limit | `stable_dt` | margin for the potential and the sponge |
| `CELLS_PER_THREAD` | 65536 | `wave.rs` | below this a thread costs more to start than it saves |
| `BANDS_PER_THREAD` | 2 | `wave.rs` | lets fast cores take over bands from slow ones |
| `SUBSTEP`, `MAX_SUBSTEPS` | ¼ cell, 64 | `particles.rs` | creeping up to nodes |
| `PARTICLES_PER_THREAD` | 2000 | `particles.rs` | thread start-up cost |
| seed | 1 | `main.rs` | reproducible runs |
| `BINS` | 128 across the plane | `run.rs` | screen resolution: 0.8 wide on a plane of 102.4 |
| `GLOW_EVERY` | 10 steps | `run.rs` | how often the time-integrated density and the bunching are read |
| `DECADES` | 3 | `view.rs` | the density range the frames' brightness spans |
| `PIXELS`, `FRAME_EVERY`, `TRAILS`, `DOTS` | 800, 0.5, 300, 2000 | `report.rs` | what a run writes |

## What is not modelled

- No magnetic field: spin enters only through the guiding law. The
  Aharonov-Bohm flux line was done earlier as a sign-flipping cut in
  the grid and can come back.
- No time-dependent potentials or walls: a barrier that is removed at
  t = 0 is written as an initial state that does not extend past it.
- No more than one field, so no interacting particles in the plane;
  two on a line fit the plane as configuration space, see above.
- No measurement apparatus: a measurement is reading a particle's
  position, and a detector that disturbs ψ would be a wall or a
  potential in the experiment.
