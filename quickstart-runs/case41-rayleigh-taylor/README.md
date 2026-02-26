# Case 41: Rayleigh-Taylor Instability — Heavy-Over-Light Fluid

## Overview

Rayleigh-Taylor (RT) instability occurs whenever a denser fluid is supported above
a lighter fluid against gravity. The configuration is gravitationally unstable: any
small perturbation at the interface grows exponentially in time, eventually producing
the characteristic mushroom-shaped "fingers" as the heavy fluid sinks and the light
fluid rises. RT instability appears throughout nature and engineering — supernova
shell collapse, inertial confinement fusion implosions, ocean overturning, and
volcanic intrusions all share the same underlying instability mechanism.

This simulation uses MOOSE's Boussinesq finite-volume Navier-Stokes solver with
temperature acting as a density marker. Cold fluid (T = 0) occupies the upper half
of the domain and hot fluid (T = 1) occupies the lower half, separated by a thin
diffuse interface at y = 1. A small sinusoidal perturbation seeds the instability.
The Boussinesq buoyancy force causes cold (denser) fluid to sink and hot (lighter)
fluid to rise, driving the roll-up of the interface into descending cold spikes and
ascending hot bubbles.

---

## The Physics

### Gravitational Instability

Consider two immiscible fluids at rest, stacked horizontally. The heavier fluid
occupies the top half (y > H) and the lighter fluid the bottom half (y < H). The
hydrostatic pressure balances gravity locally in each layer, but the interface is
mechanically unstable: if a patch of heavy fluid is displaced slightly downward, it
finds itself surrounded by lighter fluid. The net buoyancy force is now downward
rather than restoring — the displacement grows. Conversely, displaced light fluid
experiences an upward net buoyancy and accelerates upward. This positive-feedback
mechanism is the Rayleigh-Taylor instability.

### Governing Equations

The flow obeys the incompressible Navier-Stokes equations coupled to an energy
(temperature) equation under the Boussinesq approximation:

    Continuity:   div(v) = 0

    Momentum:     rho * Dv/Dt = -grad p  +  mu * laplacian(v)  +  f_buoyancy

    Energy:       rho * cp * DT/Dt = k * laplacian(T)

where D/Dt = d/dt + (v . grad) is the material derivative.

### The Boussinesq Approximation

The full variable-density problem is replaced by the Boussinesq approximation:
density is constant everywhere (rho_0 = 1) except in the gravitational body force,
where it varies linearly with temperature:

    rho_eff(T) = rho_0 * [1 - alpha * (T - T_ref)]

With alpha = 1 and T_ref = 0.5, the buoyancy body force added to the y-momentum
equation is:

    f_y = -rho_0 * alpha * (T - T_ref) * g_y

Since g_y = -1 (gravity points downward):

    f_y = rho_0 * alpha * (T - T_ref)

For T > T_ref (hot fluid): f_y > 0, i.e., upward buoyancy (lighter than average).
For T < T_ref (cold fluid): f_y < 0, i.e., downward body force (heavier than average).

By placing cold fluid (T = 0) on top and hot fluid (T = 1) on the bottom, the
density stratification is inverted: the heavy (cold) layer sits above the light
(hot) layer — the unstable RT configuration.

### Initial Condition

The temperature field is initialized as a smooth hyperbolic-tangent profile
with a superimposed sinusoidal perturbation:

    T(x, y, t=0) = 0.5 * [1 - tanh((y - 1) / delta)]
                   + epsilon * cos(2*pi*x) * exp(-(y-1)^2 / sigma^2)

where:
- delta = 0.02 — interface half-thickness (diffuse transition layer)
- epsilon = 0.05 — perturbation amplitude (5% of the temperature jump)
- sigma^2 = 0.01 — Gaussian width of the perturbation envelope

As y increases above the interface (y > 1): T -> 0 (cold, heavy). Below the
interface (y < 1): T -> 1 (hot, light). The cosine perturbation with wavenumber
k = 2*pi corresponds to a single wavelength fitting exactly across the domain
width [0, 1].

### Linear Growth Rate

In the inviscid limit, the linear stability analysis of a sharp interface gives
the RT growth rate (Chandrasekhar, 1961):

    sigma = sqrt(A * g * k)

where:
- A = (rho_heavy - rho_light) / (rho_heavy + rho_light) — Atwood number
- g = 1 — gravitational acceleration (non-dimensional)
- k = 2*pi — perturbation wavenumber

With the Boussinesq approximation, the effective Atwood number using the
temperature jump is:

    A_eff = alpha * delta_T / 2 = 1.0 * 1.0 / 2 = 0.5

Thus:

    sigma = sqrt(0.5 * 1 * 2*pi) = sqrt(pi) ~ 1.77 s^-1

The perturbation amplitude grows as exp(sigma * t), so at t = 1 the amplitude
has grown by a factor of exp(1.77) ~ 5.9. Nonlinear saturation and the
mushroom-finger pattern emerge shortly after.

Viscosity and diffusivity damp short-wavelength modes. The viscous cutoff
wavenumber is approximately k_max ~ sqrt(g * A / nu^2)^(1/3). With nu = 0.001,
the dominant instability mode remains at k = 2*pi for the selected parameters.

### Dimensionless Parameters

| Parameter | Symbol | Value | Notes |
|-----------|--------|-------|-------|
| Density | rho | 1.0 | Non-dimensional |
| Dynamic viscosity | mu | 0.001 | Re ~ 1000 at peak velocities |
| Thermal conductivity | k | 0.001 | Pe_T ~ 1000 |
| Specific heat | cp | 1.0 | Non-dimensional |
| Thermal expansion | alpha | 1.0 | Full density contrast |
| Reference temperature | T_ref | 0.5 | Midpoint of [0, 1] range |
| Gravity | g | (0, -1, 0) | Non-dimensional |
| Interface half-width | delta | 0.02 | ~1 cell width at 50 cells |
| Perturbation amplitude | epsilon | 0.05 | 5% seed |
| Wavenumber | k | 2*pi | One wavelength in x |
| Growth rate (linear) | sigma | ~1.77 | Inviscid estimate |

---

## Input File Walkthrough

### `[Mesh]` Block

A 25x50 quadrilateral mesh on the tall domain [0,1] x [0,2]. The aspect ratio
2:1 (height:width) ensures the descending cold fingers have space to develop
without prematurely hitting the bottom wall. With 25 cells across the interface
wavelength and 50 cells in the vertical, the initial interface thickness (~0.02)
spans approximately one cell — sufficient to resolve the diffuse layer.

### `[Modules/NavierStokesFV]` Block

The NavierStokesFV action creates all FV variables (vel_x, vel_y, pressure,
T_fluid), all kernels, Rhie-Chow interpolation, and boundary conditions
automatically.

Key parameters specific to this case:

**`add_energy_equation = true`**: Activates the advection-diffusion equation for
T_fluid. Temperature here is not a true thermodynamic temperature — it serves as a
passive density marker that is advected with the flow and diffuses slowly.

**`boussinesq_approximation = true`** with **`gravity = '0 -1 0'`**,
**`ref_temperature = 0.5`**, **`thermal_expansion = 'alpha'`**: Activates the
buoyancy coupling. The body force per unit volume applied to y-momentum is:

    f_y = rho * alpha * (T - T_ref) * |g|

This single term couples the temperature field to the velocity field, driving the
instability.

**`energy_wall_types = 'heatflux heatflux heatflux heatflux'`** with
**`energy_wall_functors = '0 0 0 0'`**: All four walls are thermally insulating
(zero heat flux). No temperature is imposed at any boundary, so the total thermal
energy is conserved. The temperature field evolves only through advection and
internal diffusion.

**`momentum_wall_types = 'noslip noslip noslip noslip'`**: No-slip on all walls.
This damps the velocity field near the boundaries and prevents artificial
recirculation from the walls.

**`pin_pressure = true`**: Enforces zero mean pressure. In a closed cavity with
no-slip walls, pressure is determined only up to a constant, and this removes
the nullspace.

**`momentum_advection_interpolation = 'upwind'`** and
**`energy_advection_interpolation = 'upwind'`**: Upwind differencing adds
numerical diffusion that stabilizes the high-Re flow but slightly smooths the
interface. For a sharper interface simulation, a MUSCL or central scheme would
be preferred but requires a finer mesh and tighter tolerances.

### `[FunctorMaterials]` Block

Five constant material properties are defined via `ADGenericFunctorMaterial`.
Using `mu_val = 0.001` and `k_val = 0.001` (both equal) gives a Prandtl number
Pr = mu * cp / k = 1.0. The thermal and momentum boundary layers therefore grow
at the same rate, keeping the physics symmetric between velocity and temperature.

### `[Functions]` and `[FVICs]` Blocks

The initial temperature profile is prescribed by a `ParsedFunction`:

    0.5*(1.0 - tanh((y-1.0)/0.02)) + 0.05*cos(2*pi*x)*exp(-((y-1.0)^2)/0.01)

`FVFunctionIC` samples this function at each cell center and assigns it to
T_fluid. The velocity field starts from near-zero (1e-15) — effectively at
rest — so the only initial motion seed is the unstable buoyancy distribution.

### `[Postprocessors]` Block

Six postprocessors track the simulation evolution:

- **`max_vel_y`**, **`min_vel_y`**: Maximum upward and downward vertical velocities.
  Rising hot fingers produce positive max_vel_y; descending cold spikes produce
  negative min_vel_y. Their magnitudes grow exponentially during the linear phase
  (t < 0.5) and then saturate nonlinearly.

- **`max_vel_x`**: Maximum horizontal velocity. Grows as fingers develop sideways
  bulges and the mushroom caps form.

- **`avg_T`**: Average temperature, nominally 0.5 (conserved by insulated walls).
  Deviations from 0.5 indicate numerical diffusion across the walls.

- **`max_T`**, **`min_T`**: Extremes of the temperature field. In the inviscid
  limit, these should remain at 1 and 0 throughout (pure advection). Diffusion
  and numerical dissipation will gradually mix these toward 0.5.

### `[Executioner]` Block

**`type = Transient`** with **`end_time = 3.0`** and adaptive time-stepping from
dt = 0.05: the simulation runs 3 non-dimensional time units, capturing both the
linear exponential growth phase (t < 0.5), the nonlinear roll-up phase (0.5 < t
< 1.5), and the late-time turbulent mixing (t > 1.5).

The `IterationAdaptiveDT` time stepper adjusts dt to keep Newton iterations near
the optimal count (8). During rapid roll-up, Newton iterations increase and dt is
cut back automatically. During slower late-time mixing, dt grows (up to 1.2x per
step) to advance efficiently.

**`solve_type = NEWTON`** with **`petsc_options '-pc_type lu'`**: Direct LU
factorization is robust for this coupled nonlinear problem. At 25x50 cells with
4 unknowns per cell (vel_x, vel_y, p, T), the system has about 5000 DOF —
small enough for direct solvers. A larger 100x200 mesh would require switching
to an iterative preconditioner such as `-pc_type hypre -pc_hypre_type boomeramg`.

---

## Running the Simulation

```bash
# From the moose-next root directory:
cd quickstart-runs/case41-rayleigh-taylor

# Run with the combined application (adjust path as needed)
../../modules/combined/combined-opt -i case41_rayleigh_taylor.i

# Docker equivalent (from inside the container):
cd /home/user/projects/moose/quickstart-runs/case41-rayleigh-taylor
/home/user/projects/moose/modules/combined/combined-opt -i case41_rayleigh_taylor.i
```

The simulation writes:
- `case41_rayleigh_taylor_out.e` — Exodus mesh + fields at each time step
- `case41_rayleigh_taylor_out.csv` — postprocessor time series

Estimated wall time: 10-30 minutes at 25x50 resolution (varies with hardware).

To visualize in ParaView:
1. Open `case41_rayleigh_taylor_out.e`
2. Color by `T_fluid` to see hot (light) and cold (heavy) regions
3. Apply `Glyph` filter on (vel_x, vel_y) to see velocity arrows
4. Step through time frames to watch finger development and mushroom formation
5. Use `Threshold` filter (T > 0.6 or T < 0.4) to isolate each fluid layer

---

## Expected Results

### Phase 1: Linear Growth (t = 0 to ~0.5)

The sinusoidal perturbation at the interface grows exponentially with growth
rate sigma ~ 1.77. Velocities are small but measurable:

```
t = 0.0:  max|v_y| ~ 0      (at rest)
t = 0.3:  max|v_y| ~ 0.02   (linear phase, exp growth)
t = 0.5:  max|v_y| ~ 0.15   (linear phase ending)
```

The temperature field retains its original tanh shape, deformed only slightly
by the growing perturbation. A cross-section at y = 1 shows the interface
displaced into a cosine wave of growing amplitude.

### Phase 2: Nonlinear Roll-Up (t = 0.5 to ~1.5)

The cosine perturbation becomes asymmetric. Cold spikes sharpen into downward
plumes (falling at accelerating speed) while hot regions broaden into rounded
bubbles (rising). The classical RT "mushroom cap" structure forms as the cold
spikes develop lateral overhangs.

```
t = 1.0:  max|v_y| ~ 0.5-1.0   (nonlinear)
          Cold spike descends to y ~ 0.5
          Hot bubble ascends to y ~ 1.5
```

The interface has stretched and folded significantly. The domain cross-section
at y = 1 no longer shows a simple wave — it is a broken, multi-valued surface.

### Phase 3: Late-Time Mixing (t = 1.5 to 3.0)

Cold fluid pools near the bottom; hot fluid pools near the top (the stable
final state — they have traded places). A turbulent mixing zone occupies the
middle of the domain. avg_T remains near 0.5 (conservation check). min_T and
max_T approach each other as diffusion blends the two layers.

```
Domain cross-section at t = 3.0:
  y = 2: cold (T ~ 0)          pooled heavy fluid (stable)
  y = 1: mixed (T ~ 0.3-0.7)   turbulent mixing zone
  y = 0: hot (T ~ 1)           pooled light fluid (stable)
```

The mixing zone width grows approximately as h ~ alpha_mix * A * g * t^2
(quadratic in time) during the fully nonlinear turbulent phase, consistent
with experimental measurements (Dimonte & Schneider 2000).

### Schematic of Finger Development

```
  t=0.0          t=0.5          t=1.0          t=2.0
  -------        -------        -------        -------
  COLD           C   C          C     C        C  C  C
  (heavy)         \ /          | |   | |       |  |  |
  ~ ~ ~ ~         V           |  | | |  |     mixing
  HOT             ^            | |   | |       |  |  |
  (light)        / \           |       |       H  H  H
  -------        -------        -------        -------
  flat IC     growing wave   mushroom caps   turbulent
```

---

## Key Takeaways

| Concept | Where in Input |
|---------|----------------|
| RT instability seeded by IC perturbation | `ParsedFunction` with tanh + cosine |
| Boussinesq buoyancy drives instability | `boussinesq_approximation = true` |
| Heavy (cold) fluid on top — unstable | T=0 above interface, T=1 below |
| Insulated walls — no thermal BC forcing | `energy_wall_types = heatflux`, functors=0 |
| Adaptive time-stepping through rapid roll-up | `IterationAdaptiveDT` with cutback |
| Temperature as passive density marker | `add_energy_equation = true`, Pr=1 |
| Conservation check via avg_T ~ 0.5 | `ElementAverageValue` postprocessor |
| Linear growth rate sigma ~ sqrt(A*g*k) | Matches inviscid theory at early times |
