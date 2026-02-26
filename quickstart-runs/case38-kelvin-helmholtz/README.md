# Case 38: Kelvin-Helmholtz Instability — Shear Layer Rollup

## Overview

The Kelvin-Helmholtz instability is one of the most visually striking phenomena in
fluid mechanics: a thin shear layer between two streams moving in opposite directions
spontaneously rolls up into a row of spiraling vortices — the "KH billows." It occurs
at the boundary between any two flows moving at different speeds: ocean thermoclines,
atmospheric jet streams, the edges of volcanic eruption columns, and the corona of the
Sun. It is also the mechanism by which wind raises waves on water.

This case initialises a tanh velocity profile across a 2D channel — rightward above
the midplane, leftward below — with a small sinusoidal perturbation in the transverse
velocity. The NavierStokesFV action integrates the incompressible Navier-Stokes
equations forward in time. A passive scalar (the temperature field, repurposed as dye)
marks the upper and lower streams so that the rollup is visible in the scalar field.

The physics follows Rieutord, *Fluid Dynamics* (Springer, 2015), Chapter 6, Section
6.3.1, which gives the dispersion relation and growth-rate analysis for the KH
instability of a tanh shear layer.

---

## The Physics

### Governing Equations

The simulation solves the 2D incompressible Navier-Stokes equations together with a
passive scalar transport equation:

```
Continuity:   div(v) = 0

Momentum:     rho * (dv/dt + v . grad v) = -grad p  +  mu * lap v

Energy (dye): rho * cp * (dT/dt + v . grad T)  =  k * lap T
```

The temperature field T is used purely as a passive scalar (dye marker), not to
drive buoyancy. With cp = 1 and k = mu (Prandtl number Pr = 1), the dye diffuses
at the same rate as momentum.

### The Tanh Shear Layer

The initial velocity profile is:

```
vel_x(y) = tanh((y - 0.5) / delta)
```

where delta = 0.05 is the shear layer half-thickness. This gives:

```
y = 1.0 (top):    vel_x → +1  (upper stream, rightward)
y = 0.5 (center): vel_x  = 0  (shear layer midplane)
y = 0.0 (bottom): vel_x → -1  (lower stream, leftward)
```

The velocity jump across the shear layer is dU = 2.

The passive scalar is initialised with the same tanh shape:

```
T(y) = 0.5 * (1 + tanh((y - 0.5) / delta))
```

so T = 1 in the upper stream and T = 0 in the lower stream. As vortices develop,
the rolled-up dye pattern makes the billows visible.

### The Perturbation

A small sinusoidal transverse velocity seeds the instability:

```
vel_y(x) = A * sin(2*pi*x / lambda)
```

with A = 0.01 and lambda = 1. Two full wavelengths fit in the domain [0, 2], so
two KH billows are expected to develop simultaneously. The perturbation amplitude
A = 0.01 is 0.5% of dU = 2, ensuring the flow starts in the linear growth regime.

### Linear Stability Analysis

The KH growth rate for a tanh profile (Rayleigh's equation, Drazin & Reid 1981) is
maximum at wavenumber k*delta ~ 0.4, giving:

```
sigma_max ~ 0.2 * dU / delta = 0.2 * 2 / 0.05 = 8
```

in units of inverse time. Starting from A = 0.01, the transverse velocity during
the linear phase grows as:

```
max|vel_y| ~ 0.01 * exp(sigma_max * t)
```

| t    | Predicted max|vel_y| | Regime         |
|------|----------------------|----------------|
| 0.0  | 0.010                | Initial        |
| 0.1  | 0.022                | Linear growth  |
| 0.2  | 0.049                | Linear growth  |
| 0.3  | 0.110                | Transitional   |
| 0.5  | 0.545                | Nonlinear      |
| 1.0  | saturated            | Vortex rollup  |

Saturation occurs when max|vel_y| is O(dU) = O(1) and the perturbation is no longer
small. By t ~ 0.5-1.0, two spiraling billows should be visible in the scalar field.

### Reynolds Number

```
Re = rho * dU * L / mu = 1 * 2 * 1 / 0.001 = 2000
```

where L = 1 is the domain height. The shear-layer Reynolds number is:

```
Re_delta = rho * dU * delta / mu = 1 * 2 * 0.05 / 0.001 = 100
```

At Re_delta = 100 the KH instability is strongly unstable. The simulation does not
resolve the full Kolmogorov turbulence cascade (which would require much finer grids),
but does capture the primary rollup and the large-scale vortex structure.

---

## MOOSE Implementation

### Key Design Choices

| Feature | Choice | Reason |
|---------|--------|--------|
| NavierStokesFV action | `[Modules/NavierStokesFV]` | Handles all FV NS kernels and BCs in one block |
| Energy equation | `add_energy_equation = true` | Passive scalar transport at no extra kernel cost |
| Left boundary | `inlet` with tanh functors | Continuous injection of the shear layer profile |
| Right boundary | `outlet` fixed pressure | Free outflow; fluid exits without reflection |
| Top/bottom | `slip` walls | No tangential friction, no normal flow |
| Advection scheme | `upwind` for momentum/energy | Suppresses oscillations near the steep tanh layer |
| Initial vel_y | `FVFunctionIC` with `sin(2*pi*x)` | Seeds two KH billows simultaneously |

### Object Map

| MOOSE Object | Role |
|---|---|
| `GeneratedMeshGenerator` 80x20, [0,2]x[0,1] | 2D domain, 1600 elements |
| `NavierStokesFV` action (incompressible + energy) | Navier-Stokes + passive scalar |
| `ADGenericFunctorMaterial` (rho, mu, k, cp) | Constant fluid properties |
| `ParsedFunction` (tanh profiles, sin perturbation) | IC and inlet functors |
| `FVFunctionIC` for vel_x, vel_y, T_fluid | Overrides NavierStokesFV defaults |
| `ADElementExtremeFunctorValue` (vel_x, vel_y) | Tracks velocity extremes over time |
| `ElementAverageValue` (T_fluid) | Verifies scalar conservation |
| `Transient`, Newton, LU, dt=0.01 | Time integration |
| `IterationAdaptiveDT` | Adaptive timestep control |

### Mesh Resolution

The domain is [0, 2] x [0, 1] with 80 x 20 = 1600 elements, giving:

```
dx = 2/80 = 0.025
dy = 1/20 = 0.05
```

The shear layer half-thickness delta = 0.05 equals exactly one cell in y. This is
the minimum resolution needed to represent the initial profile. The rollup structures
that develop are larger (O(lambda) = O(1)) and are well resolved. For production runs,
nx = 160, ny = 40 (double resolution) gives significantly smoother vortex cores.

---

## How to Run

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case38-kelvin-helmholtz \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case38_kelvin_helmholtz.i 2>&1 | tail -30'
```

The `MSYS_NO_PATHCONV=1` flag prevents MSYS/Git Bash from mangling the Docker volume
mount path on Windows.

Output files produced:
- `case38_kelvin_helmholtz_out.e` — Exodus file with vel_x, vel_y, pressure, T_fluid at every timestep
- `case38_kelvin_helmholtz_out.csv` — max_vel_x, min_vel_x, max_vel_y, min_vel_y, avg_T vs. time

Typical run time: 5-15 minutes (200 timesteps at dt=0.01, adaptive growth).

---

## Expected Results

### Postprocessor Time Series

| Time | max_vel_x | max_vel_y | avg_T | Notes |
|------|-----------|-----------|-------|-------|
| 0.0  | ~+1.0     | ~+0.01    | ~0.5  | Initial condition |
| 0.2  | ~+1.0     | ~+0.05    | ~0.5  | Linear growth phase |
| 0.5  | ~+1.0     | ~+0.3-0.5 | ~0.5  | Transitional phase |
| 1.0  | ~+1.5     | ~+0.8     | ~0.5  | Vortex rollup, nonlinear |
| 2.0  | varies    | varies    | ~0.5  | Mature billows, mixing |

Key observations:
- `avg_T` stays near 0.5 throughout: the passive scalar is conserved (only redistributed)
- `max_vel_x` can briefly exceed 1 inside the vortex cores (angular momentum concentration)
- `max_vel_y` grows exponentially then saturates — the signature of convective instability
- `min_vel_y` mirrors `max_vel_y` with opposite sign (two-billow symmetry)

### Scalar Field Evolution

The passive scalar T (visualised in ParaView or via `visualize_all.py`) shows:

```
t = 0:   Flat horizontal layers (T=1 above, T=0 below, smooth tanh interface)

t ~ 0.3: Gentle waviness along the interface — sinusoidal undulation visible

t ~ 0.7: Interface curling up at the crests and down at the troughs, two
          distinct spiral arms beginning to form

t ~ 1.5: Two fully rolled-up billows — each a tight spiral of alternating
          T=1 and T=0 fluid, with thin interleaved sheets of mixed scalar

t ~ 2.0: Billows have grown and shifted; beginnings of braid thinning and
          possible secondary instabilities in the thin connecting sheets
```

### Instability Growth Rate

Plotting ln(max_vel_y) vs. t during 0 < t < 0.3 should show a linear trend with
slope approximately equal to the theoretical sigma_max = 8. Any slope between 4
and 12 is consistent with numerical damping from the coarse (80x20) mesh and
the upwind advection scheme. Refining to 160x40 gives a slope closer to the theory.

---

## ASCII Domain Sketch

```
y = 1.0  +----------------------------------------------------------+
         | vel_x → +1 (upper stream, rightward)       T = 1        |
         | vel_x → +1                                              |
  slip   |                                                          | outlet
  wall   |----- shear layer: vel_x = tanh((y-0.5)/0.05) ----------| p = 0
         |                                                          |
         | vel_x → -1 (lower stream, leftward)        T = 0        |
         | vel_x → -1                                              |
y = 0.0  +----------------------------------------------------------+
         ^                                                          ^
      inlet                                                      outlet
   (tanh profile)                                              (p = 0)
      x=0                                                          x=2

Domain: [0, 2] x [0, 1]
Mesh:   80 x 20 = 1600 quadrilateral cells
Cell:   dx = 0.025,  dy = 0.05  (shear layer delta = 0.05 = 1 cell)
```

---

## Key Takeaways

| Concept | Location in Input |
|---|---|
| FV NS action with energy for passive scalar | `[Modules/NavierStokesFV]` |
| Tanh inlet profile via ParsedFunction functor | `vel_x_inlet`, `T_inlet` |
| FVFunctionIC to set non-trivial initial conditions | `[FVICs]` block |
| Sinusoidal perturbation seeds specific wavenumber | `vel_y_init_func` |
| Slip walls for frictionless top/bottom | `momentum_wall_types = 'slip slip'` |
| Upwind advection for stability near shear layer | `momentum_advection_interpolation = 'upwind'` |
| IterationAdaptiveDT for efficient time integration | `[TimeStepper]` block |
| ADElementExtremeFunctorValue tracks perturbation growth | `max_vel_y`, `min_vel_y` |
| avg_T verifies scalar mass conservation | `avg_T` postprocessor |

---

## Experiments to Try

### Experiment 1: Vary the Perturbation Amplitude

Change the amplitude of the vel_y perturbation by editing `vel_y_init_func`:

```
# Larger perturbation — skips linear phase, rolls up faster
expression = '0.1*sin(2*pi*x)'

# Smaller perturbation — longer linear phase, cleaner growth curve
expression = '0.001*sin(2*pi*x)'
```

For A = 0.001, the linear growth phase extends to t ~ 0.7 and the log-linear growth
of max_vel_y is easier to measure. For A = 0.1, rollup begins almost immediately and
the scalar structures are more strongly asymmetric due to the nonlinear start.

### Experiment 2: Change the Shear Layer Thickness

The shear layer thickness delta controls the length scale and growth rate:

```
# Thicker layer: slower growth, larger billows
expression = 'tanh((y-0.5)/0.15)'   # delta = 0.15

# Thinner layer: faster growth, smaller billows (requires finer mesh)
expression = 'tanh((y-0.5)/0.025)'  # delta = 0.025 — needs ny >= 40
```

Apply the same change to `vel_x_inlet`, `vel_x_init_func`, `T_inlet`, and `T_init_func`.
The theoretical growth rate scales as dU/delta, so doubling delta halves sigma_max.

### Experiment 3: Single Billow vs. Two Billows

Change the perturbation wavelength to see one billow at a time:

```
# Single billow (one wavelength in [0, 2])
expression = '0.01*sin(pi*x)'

# Four billows (four wavelengths in [0, 2])
expression = '0.01*sin(4*pi*x)'
```

Not all wavenumbers are equally unstable. The tanh profile is stable for very short
wavelengths (k*delta > 1). Confirm that the single-billow (k = pi, k*delta = 0.157)
grows faster than the four-billow (k = 4*pi, k*delta = 0.628) perturbation.

### Experiment 4: Increase Reynolds Number

Lower the viscosity and conductivity to increase Re:

```
mu_val = 0.0001   # Re_full = 20000, Re_delta = 1000
k_val  = 0.0001
```

At this Re the vortex cores are thinner and the braid sheets between billows develop
secondary instabilities. The coarse 80x20 mesh will show numerical diffusion effects;
refine to 160x40 or 200x50 for cleaner results. Direct LU factorisation may run out
of memory at very fine meshes — switch to `-pc_type ilu` or `-pc_type hypre` if needed.

### Experiment 5: Unequal Stream Speeds

Modify the initial and inlet profiles so the two streams do not have equal and
opposite velocities:

```
# Upper stream faster, lower stream slower (not symmetric)
expression = '0.5 + 0.5*tanh((y-0.5)/0.05)'   # range [0, 1]
```

An asymmetric shear layer rolls up into billows that convect at the mean velocity
(0.5) rather than staying stationary. Observe how the billows drift rightward across
the domain and exit through the outlet boundary.
