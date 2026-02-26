# Case 24: Charge Drift-Diffusion Between Parallel Plates

## Overview

This case models **unipolar charge injection** from an electrode into a dielectric
gap, a central topic in Melcher's *Continuum Electromechanics* (MIT Press, 1981),
Chapter 5, sections 5.5-5.7. Positive ions are continuously injected at the left
electrode, drift across the gap under the applied electric field, and are collected
at the right electrode. As the charge front propagates, Poisson's equation is solved
simultaneously to show how the accumulating space charge distorts the potential away
from the linear Laplace solution.

The drift-diffusion equation (also called the Nernst-Planck or drift-dominated
Fokker-Planck equation) appears throughout electromechanics, semiconductor device
physics, electrochemistry, and plasma physics. The MOOSE implementation here uses
`ConservativeAdvection` for the drift term, `ADMatDiffusion` for Fickian spreading,
and `ADHeatConduction` + `ADCoupledForce` for the coupled Poisson equation.

### Why a prescribed drift velocity?

A fully self-consistent simulation would extract the electric field E = -grad(phi) at
each Newton iteration and feed it directly into the drift velocity for the charge
transport equation. This coupling creates a strongly nonlinear system and introduces
convergence challenges that obscure the fundamental physics for a first study.

The prescribed-field approach used here — v = mu_i * E_applied, with E_applied = V/d
(the uniform Laplace-solution field) — is valid when space charge is weak relative to
the applied field. Quantitatively, the approximation holds when the space charge density
rho_e is small compared to epsilon * E_applied / d. For the parameters here
(rho_e ~ 1, eps = 1, E = 1, d = 1) space charge is not negligible at steady state,
so the Poisson solve shows a measurable potential distortion, but the drift velocity
is still well-approximated by the applied-field value for the initial transient. This
is the standard starting point for space-charge-limited conduction analysis in
Melcher's treatment.

New concepts introduced in this case:

- **`ConservativeAdvection`**: MOOSE kernel for the conservative form div(v * rho_e),
  with full upwinding for stability at sharp charge fronts.
- **Coupled Poisson + transport**: solving two coupled PDEs (rho_e and phi) in a
  single nonlinear system with full Jacobian coupling via `[Preconditioning] [SMP]`.
- **`GenericConstantVectorMaterial`**: supplying a constant vector material property
  (the drift velocity vector) consumed by the advection kernel.
- **Quasi-1D geometry**: a 60x4 slab that makes the physics effectively one-dimensional
  while remaining formally two-dimensional for MOOSE.
- **Upwinding for advection stability**: `upwinding_type = full` prevents the numerical
  oscillations (Gibbs ringing) that appear at sharp fronts when central differencing
  is used for hyperbolic transport.

---

## The Physics

### Physical Problem in Plain English

Picture two parallel metal plates separated by a gap of d = 1 m filled with a weakly
conducting dielectric fluid. The left plate (anode) is connected to a 1 V source; the
right plate (cathode) is grounded. At time t = 0, positive ions begin to be injected
from the left electrode into the fluid with charge density rho_e = 1 C/m³.

The applied field E = V/d = 1 V/m pushes the positive ions rightward. The ions move
at drift velocity v = mu_i * E = 1.0 m/s. Simultaneously, random thermal motion
spreads the charge via Fickian diffusion with diffusivity D_i = 0.01 m²/s.

The result is a charge wave that advances from left to right. The front is not a
perfect step because diffusion rounds it. At the same time, the accumulated space
charge creates an additional electric field that adds to (near the injection electrode)
or subtracts from (in the bulk) the applied field — this is captured by the Poisson
equation, which shows the potential bending away from the linear Laplace solution.

After one transit time t_transit = d / v = 1/1 = 1 s, the charge front reaches the
collecting electrode and a quasi-steady state develops: a roughly linear charge
profile from rho_e = 1 at the left to rho_e = 0 at the right, with a corresponding
distorted potential profile.

### Governing Equations

**Drift-diffusion (charge transport):**

```
∂ρ_e/∂t  +  div( v · ρ_e )  =  D_i · ∇²ρ_e    in Omega = [0,1] x [0,0.067]

ρ_e = 1    at x = 0  (injection electrode)
ρ_e = 0    at x = 1  (collecting electrode)
```

where:
- `ρ_e` — free charge density [C/m³]
- `v = μ_i · E_applied · x̂ = (1.0, 0, 0)` m/s — prescribed drift velocity
- `D_i = 0.01` m²/s — ion diffusivity
- `div(v · ρ_e)` — conservative divergence form of the drift flux

**Poisson's equation (electrostatics):**

```
-div( ε · grad(φ) )  =  ρ_e    in Omega

φ = 1    at x = 0  (anode)
φ = 0    at x = 1  (cathode)
```

where:
- `φ` — electric potential [V]
- `ε = 1.0` — permittivity [C²/(N·m²)]

**Electric field (derived, not solved directly):**

```
E = -grad(φ)
```

At steady state with uniform rho_e, the Poisson equation gives a parabolic potential
profile rather than the linear Laplace solution. The space-charge correction shifts
the field: stronger near the injecting electrode, weaker near the collecting electrode.

### Dimensionless Groups

The two characteristic time scales are:

```
t_drift    = d / (μ_i · E)  = 1 / (1 · 1) = 1.0 s     (convective transit time)
t_diffusion = d² / D_i       = 1² / 0.01   = 100 s     (diffusive spread time)
```

Their ratio is the Peclet number:

```
Pe = t_diffusion / t_drift = v · d / D_i = 1 · 1 / 0.01 = 100
```

At Pe = 100, drift strongly dominates diffusion. The charge front advances as a
near-step, broadened only slightly by diffusion. Full upwinding in the advection
kernel is essential at Pe >> 1.

The space-charge parameter (ratio of space-charge field to applied field):

```
C_sc = ρ_e · d / (ε · E_applied) = 1 · 1 / (1 · 1) = 1
```

C_sc = 1 means space charge significantly distorts the field at steady state,
justifying the simultaneous Poisson solve.

### Domain Diagram

```
x=0 (anode, φ=1)                        x=1 (cathode, φ=0)
Injection: ρ_e=1                          Collection: ρ_e=0
     |                                         |
y=0.067 +=========================================+
     |    → → →  drift velocity v = 1 m/s → → →  |
     |                                             |
     |   charge front propagates rightward         |
     |   with diffuse leading edge                 |
     |                                             |
y=0  +=========================================+
     |← ————————  d = 1 m  ————————————————→|

Mesh: 60 elements (x) × 4 elements (y)
Boundary sidesets: left, right (physics), top, bottom (passive — no BC needed)

Time evolution of ρ_e at centre line (schematic):

t=0.0 s:  ρ_e = 0 everywhere (charge-free gap)
             |___________________________

t=0.2 s:  front at x ≈ 0.2, diffuse leading edge
             |▓▓▓▓░░___________________

t=0.6 s:  front at x ≈ 0.6
             |▓▓▓▓▓▓▓▓▓▓▓▓░░__________

t=1.0 s:  front reaches cathode, quasi-steady linear profile
             |▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓_|

Legend: ▓ = ρ_e near 1, ░ = diffuse front, _ = ρ_e near 0
```

---

## Input File Walkthrough

The input file is `case24_drift_diffusion.i`.

### Header Parameters

```
mu_ion = 1.0     # ion mobility [m²/(V·s)]
D_ion  = 0.01    # ion diffusivity [m²/s]
E_app  = 1.0     # applied field V/d [V/m]
eps    = 1.0     # permittivity
```

Named top-level parameters allow changing the physics in one place. They are
referenced via `${mu_ion}`, `${D_ion}`, etc. throughout the file. The Peclet
number Pe = mu_ion * E_app * d / D_ion can be changed by adjusting D_ion (lower
D_ion = sharper front, harder to resolve; higher D_ion = smoother front, more
diffusion-dominated).

### `[Mesh]`

```
[gen]
  type = GeneratedMeshGenerator
  dim  = 2
  nx   = 60
  ny   = 4
  xmax = 1
  ymax = 0.067
```

A 60x4 structured quadrilateral mesh. The element size dx = 1/60 ≈ 0.017 m. The
charge front width due to diffusion at time t is proportional to sqrt(D_i * t) ≈
sqrt(0.01 * 1) = 0.1 m, which is resolved by approximately 6 elements — sufficient
for a qualitative result. For publication-quality results, nx = 200 would be appropriate.

The quasi-1D geometry (ny = 4, ymax = 0.067) keeps the computational cost low while
remaining formally 2D. All physics variables are uniform in y.

### `[Variables]`

Two primary variables are solved simultaneously:

- `rho_e` — free charge density, initial condition 0 everywhere.
- `phi` — electric potential, initial condition 0 everywhere. (The correct initial
  potential for the zero-charge state is the Laplace solution phi(x) = 1 - x, but
  starting from phi = 0 is acceptable because the solver converges to the correct
  solution within the first Newton iteration of the first timestep.)

### `[Kernels]`

Five kernels implement the two PDEs:

| Kernel | Type | Variable | Term implemented |
|--------|------|----------|-----------------|
| `rho_time` | ADTimeDerivative | rho_e | ∂ρ_e/∂t |
| `rho_diffusion` | ADMatDiffusion | rho_e | -D_i · ∇²ρ_e |
| `rho_advection` | ConservativeAdvection | rho_e | div(v · ρ_e) |
| `phi_laplacian` | ADHeatConduction | phi | -div(ε · grad(φ)) |
| `phi_source` | ADCoupledForce | phi | -ρ_e (RHS source) |

**`ConservativeAdvection`** (rho_advection):

This kernel implements the conservative divergence form:

```
R_i = ∫_Omega v · grad(ρ_e) · ψ_i dV + boundary terms
    = - ∫_Omega ρ_e · (v · grad(ψ_i)) dV + ∫_∂Omega (v · n) ρ_e ψ_i dS
```

The conservative form is critical for charge balance: the total charge in the domain
changes only due to boundary fluxes, never due to interior numerical error. With
`upwinding_type = full`, the interface value of ρ_e between neighbouring elements
is taken entirely from the upwind element (the one the flow is coming from), which
adds numerical diffusion equal to (1/2) * |v| * dx. This stabilises the sharp front
at the cost of some front smearing.

**`ADHeatConduction`** (phi_laplacian):

This kernel's standard use is for heat conduction: ∫ k · grad(T) · grad(ψ) dV.
For the Poisson equation, the mapping is exact: k → ε, T → φ. The material
property name must match the `thermal_conductivity` parameter.

**`ADCoupledForce`** (phi_source):

Adds -∫ ρ_e · ψ dV to the residual for φ. This places ρ_e as a body source on the
RHS of the Poisson equation, representing the charge density driving the potential.

### `[BCs]`

Four Dirichlet boundary conditions:

| Block | Variable | Boundary | Value | Physical meaning |
|-------|----------|----------|-------|-----------------|
| `rho_inject` | rho_e | left | 1.0 | Continuous charge injection |
| `rho_collect` | rho_e | right | 0.0 | Perfect charge collection |
| `phi_anode` | phi | left | 1.0 | Applied voltage (anode) |
| `phi_cathode` | phi | right | 0.0 | Ground reference (cathode) |

The top and bottom boundaries carry no explicit BCs. For `rho_e`, the natural Neumann
condition (zero normal flux) applies — this is correct because no charge leaves through
the insulating side walls. For `phi`, the natural Neumann condition (zero normal
field component) is also correct — the side walls are insulating.

### `[Materials]`

Three material objects:

**`ion_diffusivity`** (`ADGenericConstantMaterial`): provides scalar property `D_ion`
consumed by `ADMatDiffusion`. The AD prefix ensures derivative information is propagated
for exact Jacobian assembly.

**`permittivity_mat`** (`ADGenericConstantMaterial`): provides scalar property
`permittivity` consumed by `ADHeatConduction` (phi_laplacian). Naming the property
`permittivity` rather than `thermal_conductivity` is purely cosmetic — the kernel
accepts any name via the `thermal_conductivity` parameter.

**`drift_velocity_mat`** (`GenericConstantVectorMaterial`): provides a constant vector
material property `drift_velocity = (1.0, 0, 0)` consumed by `ConservativeAdvection`.
The vector components are `(v_x, v_y, v_z)`. Setting v_x = mu_ion = 1.0 encodes
rightward drift at the applied-field speed. Note that `GenericConstantVectorMaterial`
(without `AD`) is used here because `ConservativeAdvection` does not carry AD support
(it is not an AD kernel); mixing AD and non-AD objects in the same system is
permissible when they act on different variables, but requires care with the Jacobian.

### `[Postprocessors]`

| Name | Type | Variable | Meaning |
|------|------|----------|---------|
| `avg_rho` | ElementAverageValue | rho_e | Mean charge density (rises from 0 to ~0.5) |
| `max_rho` | ElementExtremeValue (max) | rho_e | Peak charge (should stay near 1.0) |
| `max_phi` | ElementExtremeValue (max) | phi | Maximum potential (should stay near 1.0) |
| `min_phi` | ElementExtremeValue (min) | phi | Minimum potential (should stay near 0.0) |

The `avg_rho` postprocessor is especially diagnostic: starting from 0, it rises as
the charge front traverses the gap. Its steady-state value depends on the profile shape.
For a linear rho_e profile (rho_e = 1 - x), the average is exactly 0.5.

### `[Executioner]`

```
type       = Transient
solve_type = NEWTON
dt         = 0.02
end_time   = 1.0
```

`Transient` with fixed dt = 0.02 s. The implicit BDF1 (backward Euler) time integrator
is unconditionally stable, so the Courant number C = v * dt / dx = 1 * 0.02 / 0.017 ≈ 1.2
does not limit stability (though it does affect temporal accuracy). Using 50 timesteps
to traverse the domain gives good resolution of the advancing front.

Newton with BoomerAMG algebraic multigrid preconditioner handles the coupled system
efficiently. The full SMP preconditioning (see below) captures the rho_e–phi coupling
in the Jacobian.

### `[Preconditioning]`

```
[SMP]
  type = SMP
  full = true
[]
```

`SMP` (Single Matrix Preconditioner) with `full = true` assembles all off-diagonal
Jacobian blocks, including the block coupling phi to rho_e (from the `phi_source`
kernel) and the block coupling rho_e to phi (from `rho_advection`, if phi-dependence
were included). Without `full = true`, the preconditioner ignores the cross-variable
coupling and requires more Newton iterations. With `full = true`, the coupling is
captured and convergence is typically 3-6 Newton iterations per timestep.

---

## Running the Simulation

This case uses only framework kernels and materials (no physics modules beyond the
base framework), so the framework test executable suffices:

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case24-drift-diffusion \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case24_drift_diffusion.i 2>&1 | tail -30'
```

The `combined-opt` executable includes all modules and is always safe to use. The
simulation runs 50 timesteps (dt = 0.02, end_time = 1.0), which typically completes
in under a minute.

Output files written to the same directory:

- `case24_drift_diffusion_out.e` — Exodus file with rho_e and phi at every timestep
- `case24_drift_diffusion_out.csv` — Postprocessor time history (avg_rho, max_rho,
  max_phi, min_phi vs. time)

Expected console output per timestep:

```
Time Step 1, time = 0.02
 0 Nonlinear |R| = ...
 1 Nonlinear |R| = ...  (drops several orders of magnitude)
 2 Nonlinear |R| = ~1e-10
  Solve Converged!
```

Newton should converge in 2-4 iterations per timestep for most of the run. Near
t = 1.0 when the front reaches the cathode, convergence may require 1-2 extra
iterations as the solution adjusts to the new steady state.

---

## Expected Results

### Charge Density Field

The charge front propagates from left to right at the drift speed v = 1.0 m/s:

| Time [s] | Approximate front position [m] | avg_rho | max_rho |
|----------|-------------------------------|---------|---------|
| 0.0      | 0.0 (just starting)           | ~0.0    | ~0.0    |
| 0.2      | ~0.2                          | ~0.18   | ~1.0    |
| 0.4      | ~0.4                          | ~0.37   | ~1.0    |
| 0.6      | ~0.6                          | ~0.48   | ~1.0    |
| 0.8      | ~0.8                          | ~0.53   | ~1.0    |
| 1.0      | reached cathode (x=1)         | ~0.50   | ~1.0    |

The front position is approximate due to diffusion spreading. At Pe = 100, the
diffuse front width is approximately sqrt(4 * D_i * t / v) ≈ sqrt(4 * 0.01 / 1)
= 0.2 m — about 20% of the gap width.

Once the front reaches the cathode (t ~ 1 s), the system transitions to the
quasi-steady state: rho_e decreasing monotonically from 1 at the anode to 0 at
the cathode, with the exact profile determined by the balance of drift and diffusion.

### Electric Potential Field

At t = 0, the potential is the Laplace solution: phi(x) = 1 - x (linear).

As charge accumulates, space charge bends the potential. Positive charge (rho_e > 0)
acts as a distributed source, pushing phi upward in the charged region. By t = 1 s,
the potential is convex-upward in the charged zone compared to the linear baseline.
The max_phi postprocessor should read ~1.0 throughout (anode BC is fixed); the
min_phi should read ~0.0 throughout (cathode BC is fixed).

The degree of bending is quantified by the Poisson equation. For a uniform charge
density rho_e = C across the full gap with phi = 0 at both ends, the solution is:

```
phi(x) = C/(2*eps) * x * (1 - x)     (parabolic correction)
```

With C = 1 and eps = 1, the peak correction at x = 0.5 is 1/8 = 0.125 V — a 12.5%
distortion above the Laplace solution value of 0.5 V at the midpoint.

### Steady-State Profile

The 1D steady-state drift-diffusion equation is:

```
v * d(rho_e)/dx  =  D_i * d²(rho_e)/dx²

Solution: rho_e(x) = (exp(Pe * x) - exp(Pe)) / (1 - exp(Pe))
                   where Pe = v * d / D_i = 100
```

At Pe = 100, exp(Pe) is extremely large, so the profile is essentially:

```
rho_e(x) ≈ exp(-Pe * (1 - x)) ≈ 0   for x < 1
rho_e(x) ≈ 1                         for x very close to 0
```

In practice this means the steady-state charge density is ~1 throughout most of the
gap except for a very thin boundary layer of width ~1/Pe = 0.01 m near the cathode.
At the resolution of this mesh (dx = 0.017 m), this thin layer is only marginally
resolved. This is the physical regime of space-charge-limited conduction: virtually
all the charge injected piles up throughout the gap.

---

## Key Takeaways

- **`ConservativeAdvection`** implements the conservative form div(v * rho_e) and
  is the correct MOOSE kernel for transport problems that must conserve total charge
  (or mass, or any scalar). The `upwinding_type = full` parameter is essential for
  stability at high Peclet numbers; without it, sharp fronts produce spurious
  oscillations.

- **Coupled Poisson + transport** is the prototype for many self-consistent field
  problems (electrostatics, semiconductor drift-diffusion, Nernst-Planck, Vlasov).
  MOOSE solves both equations simultaneously in a single Newton iteration, which is
  more robust than alternating between them (operator splitting).

- **Full SMP preconditioning** (`[Preconditioning] [SMP] full = true`) is essential
  when multiple coupled variables share off-diagonal Jacobian blocks. It ensures that
  the Newton update accounts for the cross-variable coupling, cutting iteration counts
  significantly compared to block-diagonal preconditioning.

- **`ADHeatConduction` as a general Laplacian operator**: the kernel is not limited
  to heat problems. Any equation of the form -div(k * grad(u)) = source can be
  implemented by mapping k → `thermal_conductivity` and u → the MOOSE variable.
  This reuse avoids writing new kernels for standard operators.

- **Prescribed vs. self-consistent coupling**: this case uses a prescribed drift
  velocity (E_applied = V/d = constant). Full self-consistency requires expressing v
  as a function of grad(phi) at every quadrature point. MOOSE supports this through
  coupled material properties or through the AD system, but introduces stronger
  nonlinearity that requires more Newton iterations and careful preconditioning.

- **Peclet number and upwinding**: at Pe >> 1, central differencing of the advection
  term produces oscillations that grow in amplitude and eventually dominate the
  solution. Full upwinding cures the oscillations at the cost of numerical diffusion
  proportional to (1/2) * v * dx. A practical rule: always use upwinding when
  Pe_element = v * dx / D > 2.

- **`GenericConstantVectorMaterial`** provides a constant 3-component vector property
  without AD support. It is the correct choice when the velocity is a prescribed
  constant and the consuming kernel is non-AD (as `ConservativeAdvection` is). For
  a spatially varying or time-varying drift velocity, a custom material or a
  `DerivativeParsedMaterial` would be needed.

---

## Experiments to Try

### Experiment 1: Vary the Peclet Number

Change `D_ion = 0.1` (Pe = 10) and re-run. The charge front will be noticeably wider
and more diffuse. The transition from injection-controlled (Pe >> 1) to diffusion-
controlled (Pe << 1) behaviour is one of the most important concepts in charge transport.

### Experiment 2: Reverse the Field (Push Charges Away)

Change the phi BCs to `phi_left = 0` and `phi_right = 1` (reverse polarity) and the
drift velocity to `prop_values = '-1.0 0 0'`. Ions injected at x=0 now face a field
that pushes them back toward the injection electrode. This models the space-charge-
blocking regime in Melcher §5.6.

### Experiment 3: Increase Injection Density

Change `rho_inject` value to 10.0. The space-charge distortion of the potential
becomes large (C_sc = 10). The phi field will show a strongly curved profile
deviating significantly from the linear Laplace solution.

### Experiment 4: Refine the Mesh Near the Cathode

Add a biased mesh with more refinement near x = 1 to resolve the thin boundary layer.
Use `bias_x = 0.9` in the GeneratedMeshGenerator to cluster elements toward the right:

```
[gen]
  type  = GeneratedMeshGenerator
  ...
  bias_x = 0.9    # elements shrink toward x=1
```

This resolves the Pe=100 steady-state boundary layer (width ~0.01 m) with several
elements instead of less than one.

### Experiment 5: Visualize in ParaView

Open `case24_drift_diffusion_out.e` in ParaView and:

1. Color by `rho_e` and step through timesteps to watch the charge front advance.
   Use a blue-white-red colormap with range [0, 1].

2. Color by `phi` to see the potential evolving from linear (t=0) to curved (t=1).
   The curvature is concave toward the injection electrode where space charge
   accumulates.

3. Use Filters > Plot Over Line along y = 0.033 (centerline) to extract the 1D
   rho_e(x) and phi(x) profiles at selected times and compare to the analytical
   estimates above.
