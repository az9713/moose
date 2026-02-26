# Case 27: MHD Hartmann Flow — Magnetic Braking of Channel Flow

## Overview

This case models **Hartmann flow** — the magnetohydrodynamic (MHD) analogue of
Poiseuille flow in which a uniform transverse magnetic field B0 brakes the motion
of an electrically conducting fluid flowing through a channel. As the conducting fluid
moves, it generates induced currents j = sigma*(v x B0). Those currents experience a
Lorentz force j x B0 that opposes the flow, redistributing momentum and flattening
the velocity profile from the familiar parabola of viscous Poiseuille flow into the
characteristic flat-topped shape known as the Hartmann profile.

The physics follows Melcher, *Continuum Electromechanics* (MIT Press, 1981),
Chapter 9, Sections 9.9-9.10, one of the canonical treatments of
low-magnetic-Reynolds-number MHD.

The key observation exploited in this input file is that the Lorentz drag
`-sigma*B0^2*v` is mathematically identical to Darcy friction in a porous medium
`-mu/(rho*K)*v`. Setting the permeability `K = mu/(rho*sigma*B0^2) = 1/Ha^2` makes
the two expressions identical. This allows the standard MOOSE `porous_medium_treatment`
path in the `NavierStokesFV` action to model the magnetic braking without any custom
kernels or `ADParsedFunctorMaterial`.

This case builds directly on the incompressible Navier-Stokes machinery established
in Case 15 (lid-driven cavity) and Case 16 (natural convection), replacing the closed
cavity setup with a proper inlet-outlet channel configuration.

**Connection to earlier cases:**

- **Case 15 (Lid-Driven Cavity)**: Introduced `[Modules/NavierStokesFV]` for
  incompressible Navier-Stokes on a closed domain with all no-slip walls. Case 27
  uses the same action with inlet and outlet boundary conditions instead.
- **Case 16 (Natural Convection)**: Demonstrated adding a body force (buoyancy) via
  the action's built-in parameters. Case 27 uses the built-in `friction_types = 'darcy'`
  parameter to add the Lorentz drag without any explicit `FVKernels` block.
- **Case 19 (Porous Flow)**: Introduced the porous medium treatment in the
  NavierStokesFV action. Case 27 borrows that mechanism for a non-porous physics
  purpose — MHD braking.

---

## The Physics

### Physical Setup

Consider a channel of length 5 and height 1 with electrically insulating walls at
y = 0 and y = 1. The fluid has electrical conductivity sigma and dynamic viscosity mu.
A uniform magnetic field B0 is applied in the z-direction (perpendicular to both the
flow direction and the wall-normal direction). Flow enters from the left with a uniform
velocity of 1.0 and exits freely at zero pressure on the right.

As the fluid moves with velocity v_x(y), it cuts through the magnetic field lines and
generates an induced electromotive force. For insulating walls, the resulting Lorentz
body force opposes the flow:

```
f_Lorentz = j x B = -sigma * B0^2 * v_x    [N/m^3, in the x-direction]
```

This force is linear in velocity: faster-moving fluid in the channel core experiences
stronger braking than the slower near-wall fluid. The equilibrium between the
driving inlet momentum, viscous diffusion, and Lorentz drag produces the Hartmann
profile.

### Darcy-Lorentz Equivalence

The Darcy friction term in a porous medium is:

```
f_Darcy = -mu / (rho * K) * v
```

Setting `K = mu / (rho * sigma * B0^2)` gives:

```
f_Darcy = -sigma * B0^2 * v = f_Lorentz
```

With non-dimensional parameters rho = 1, mu = 1, and Ha^2 = sigma*B0^2 = 25, the
Darcy coefficient 1/K = Ha^2 = 25. This equivalence means the MOOSE porous medium
module handles the Lorentz physics exactly, with no custom code required.

### Governing Equation

Steady momentum balance in x:

```
0  =  -dp/dx  +  mu * d^2v_x/dy^2  -  sigma*B0^2 * v_x
```

The Hartmann number:

```
Ha = B0 * L * sqrt(sigma/mu)
```

With Ha = 5, the analytical profile for fully developed flow (far from the inlet) with
no-slip conditions v_x(0) = v_x(1) = 0 and bulk average velocity v_avg is:

```
v(y) = v_max * [1 - cosh(Ha * (y - 0.5)) / cosh(Ha/2)]
```

where v_max is determined by the average velocity constraint. For an inlet velocity
of 1.0 and the domain geometry used here, the flow approaches this profile in the
developed region.

### The Hartmann Number

The Hartmann number Ha is the fundamental dimensionless parameter for MHD channel flow:

```
Ha  =  B0 * L * sqrt(sigma/mu)
```

It measures the ratio of electromagnetic braking to viscous diffusion. For large Ha,
the Lorentz force dominates everywhere except in thin layers near the walls of
thickness delta ~ L/Ha (the Hartmann layer). In the limit Ha = 0, the solution
recovers the parabolic Poiseuille profile.

| Ha    | Flow character                                              |
|-------|-------------------------------------------------------------|
| 0     | Parabolic Poiseuille profile (no magnetic field)            |
| 1     | Slight flattening, visible departure from parabola          |
| 5     | Clearly flat core, thin Hartmann layers (delta ~ 0.2*L)     |
| 10    | Very flat core, layers beginning to require mesh refinement |
| 50    | Plug flow in the core, extremely thin layers (delta ~ 0.02) |

This case uses Ha = 5 (Ha^2 = 25), giving a clearly flattened profile with Hartmann
boundary layer thickness delta ~ L/Ha = 0.2 — well resolved by the 25 cells in the
y-direction (dy = 0.04, giving ~5 cells inside the layer).

### ASCII Domain Sketch

```
y=1  +--inlet--+-------- no-slip top wall (Hartmann wall) --------+--outlet--+
     |         |                                                   |          |
     | v_x=1.0 |  -->  -->  -->  Hartmann channel flow  -->  -->  | p = 0    |
     |         |   flat core, v_x > 1 in the developed region     |          |
     | v_y=0   |                                                   |          |
y=0  +---------+-------- no-slip bot wall (Hartmann wall) --------+----------+
     x=0                                                           x=5

Magnetic field B0 in z-direction (out of page)
Hartmann layers: thin boundary layers at y=0 and y=1
Inlet: fixed velocity (v_x=1, v_y=0)       Outlet: fixed pressure (p=0)
Mesh: 50 x 25 uniform quadrilateral cells (finite volume, cell-centered)
```

---

## Input File Walkthrough

The input file is `case27_hartmann_flow.i`.

### Top-Level HIT Variable

```
Ha2 = 25.0   # Hartmann number squared Ha^2 = 25  (Ha = 5)
```

A single variable defined at the top of the file and referenced throughout as `${Ha2}`.
Changing this one value updates all Darcy coefficients consistently. The value Ha^2 = 25
corresponds to Ha = 5 and sets the Darcy coefficient 1/K = 25 in all coordinate
directions.

### Block: `[Mesh]`

A 50x25 uniform quadrilateral mesh on the rectangle [0, 5] x [0, 1]:

```
[gen]
  type = GeneratedMeshGenerator
  dim  = 2
  xmin = 0
  xmax = 5
  ymin = 0
  ymax = 1
  nx   = 50
  ny   = 25
```

The longer x-dimension (5 channel heights) ensures the flow has sufficient length
to develop from the inlet toward the Hartmann profile before reaching the outlet.
`GeneratedMeshGenerator` automatically creates named sidesets `left`, `right`,
`top`, and `bottom`, which are referenced in the NavierStokesFV block.

### Block: `[Modules/NavierStokesFV]`

This action block creates all FV variables, kernels, boundary conditions, and the
Rhie-Chow interpolator automatically.

**`porous_medium_treatment = true`**: Activates the Darcy friction term in the
momentum equations. When this flag is set, the variable names change from `vel_x`
and `vel_y` to `superficial_vel_x` and `superficial_vel_y` (the porous-medium
superficial, or Darcy, velocity). In this case the porosity is set to 1.0, so the
superficial velocity equals the physical velocity.

**`friction_types = 'darcy'` and `friction_coeffs = 'Darcy_coefficient'`**: These
parameters specify that the drag model is Darcy friction (linear in velocity) and
that the coefficient is read from the functor material property `Darcy_coefficient`
defined in `[FunctorMaterials]`. The Darcy coefficient vector holds `[Ha^2, Ha^2, Ha^2]`,
giving the isotropic drag `-Ha^2 * v` in each momentum component. The x-component
drag is the Lorentz drag; the y- and z-components are formally zero in a perfect
Hartmann problem but their presence is harmless in the 2D simulation.

**Inlet boundary (`left`)**: `momentum_inlet_types = 'fixed-velocity'` with
`momentum_inlet_functors = '1 0'` prescribes a uniform inlet profile
v_x = 1.0, v_y = 0.0. The flow enters with a flat (plug) profile and develops
the Hartmann shape downstream.

**Outlet boundary (`right`)**: `momentum_outlet_types = 'fixed-pressure'` with
`pressure_functors = '0'` sets the outlet gauge pressure to zero. This is the
standard pressure reference for a channel-flow problem. The pressure drop along
the channel is then determined by the momentum balance (viscous and Lorentz losses).

**Wall boundaries (`top` and `bottom`)**: `momentum_wall_types = 'noslip noslip'`
enforces zero velocity at the Hartmann walls. These are the walls where the
Hartmann boundary layers form.

**No `pin_pressure`**: Unlike the closed-cavity cases (15, 16, 26), this channel
has a defined pressure outlet. The pressure is anchored at the right boundary and
no pinning constraint is needed.

**`momentum_advection_interpolation = 'upwind'`**: Upwind differencing for the
convective term. The effective magnetic Reynolds number is small (Lorentz drag
suppresses inertia) so the convective term is relatively weak and upwind is both
stable and sufficiently accurate.

### Block: `[FunctorMaterials]`

Three functor materials provide the fluid properties:

**`fluid_props` (`ADGenericFunctorMaterial`)**: Declares `rho = 1.0` and `mu = 1.0`.
Unit density and viscosity make Ha the single controlling parameter.

**`darcy_coeff` (`ADGenericVectorFunctorMaterial`)**: Declares the vector property
`Darcy_coefficient` with value `[${Ha2}, ${Ha2}, ${Ha2}]` = [25, 25, 25]. The
porous medium treatment reads a vector coefficient to allow anisotropic drag;
setting all three components to Ha^2 gives isotropic braking.

**`porosity` (`ADGenericFunctorMaterial`)**: Declares `porosity = 1.0`. Required
by the porous medium treatment to convert between superficial (Darcy) and
interstitial (physical) velocities. Porosity = 1 means the two are identical —
this is not a real porous medium, only the drag mechanism is borrowed.

### Block: `[Postprocessors]`

Four postprocessors record the converged state:

- **`max_vel_x`**: Peak x superficial velocity via `ADElementExtremeFunctorValue`.
  For inlet velocity 1.0 and Ha = 5, the Hartmann profile in the developed region
  has a flat core with velocity above the inlet average, so max_vel_x > 1.0.

- **`avg_vel_x`**: Domain-average x velocity via `ElementAverageValue`. By mass
  conservation in the channel with uniform inlet, this is expected to be close to
  the inlet velocity of 1.0.

- **`max_vel_y`**: Maximum y velocity. This should be small (confined to the inlet
  development region), confirming the flow is predominantly streamwise.

- **`avg_pressure`**: Domain-average pressure. With the outlet pinned to zero and
  Lorentz drag along the full channel length, this will be a positive value
  reflecting the pressure drop.

### Block: `[Executioner]`

Standard steady-state Newton solve with LU factorisation:

```
type        = Steady
solve_type  = 'NEWTON'
petsc_options_iname = '-pc_type -pc_factor_shift_type'
petsc_options_value = 'lu       NONZERO'
nl_rel_tol  = 1e-8
nl_abs_tol  = 1e-10
```

The Darcy drag term adds a positive diagonal contribution Ha^2 = 25 to the velocity
Jacobian block, which improves the condition number compared to pure Navier-Stokes.
This is why convergence is fast — typically 2 Newton iterations.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case27-hartmann-flow \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case27_hartmann_flow.i 2>&1 | tail -30'
```

The `MSYS_NO_PATHCONV=1` environment variable prevents MSYS/Git Bash on Windows from
mangling the absolute path in the `-v` mount argument.

The simulation is steady-state and typically completes in under 10 seconds. Outputs:

- `case27_hartmann_flow_out.e` — Exodus file with `superficial_vel_x`,
  `superficial_vel_y`, and `pressure` fields
- `case27_hartmann_flow_out.csv` — CSV file with the four postprocessor values at
  convergence

---

## Expected Results

### Postprocessor Values

| Quantity       | Expected (50x25 mesh)  | Notes                                    |
|----------------|------------------------|------------------------------------------|
| `max_vel_x`    | ~1.378                 | Peak in the flat Hartmann core           |
| `avg_vel_x`    | ~1.0                   | Mass conservation from inlet = 1.0       |
| `max_vel_y`    | small (O(0.01))        | Transverse velocity in development region|
| `avg_pressure` | positive               | Pressure drop driven by Lorentz + viscous|

The simulation converges in approximately 2 Newton iterations, reflecting the
well-conditioned system created by the positive Darcy diagonal.

### Velocity Profile Shape

In the developed flow region (away from the inlet), a cross-section at fixed x
should show:

1. v_x = 0 at y = 0 and y = 1 (no-slip Hartmann walls)
2. A rapid rise from zero over approximately 0.2 L (the Hartmann layer, delta ~ L/Ha)
3. A nearly flat core with v_x above the average (approximately 1.4) across
   y in [0.2, 0.8]
4. A symmetric fall back to zero on the upper Hartmann layer

This flat-topped shape is visually distinct from the smooth parabola that pure
Poiseuille flow with the same mean velocity would produce.

Opening `case27_hartmann_flow_out.e` in ParaView and using the Plot Over Line
filter along a vertical line at x = 2.5 will show the developed Hartmann profile.

### ASCII Profile (Cross-Section in Developed Region)

```
y=1  +--wall (v=0)---+
     |                \
y=0.8|                 | Hartmann layer (rapid change)
     |     flat core   |
y=0.5|  v_x ~ 1.4     |
     |     flat core   |
y=0.2|                 | Hartmann layer (rapid change)
     |                /
y=0  +--wall (v=0)---+
     0    v_x -->  1.4

  Compare to Poiseuille: smooth parabola, same mean velocity
  Hartmann: strong Lorentz drag suppresses the parabolic peak
```

---

## Key Takeaways

| Concept | Where in Input |
|---------|----------------|
| Darcy friction as Lorentz drag analog | `friction_types = 'darcy'`, `friction_coeffs = 'Darcy_coefficient'` |
| `porous_medium_treatment = true` activates drag | `[Modules/NavierStokesFV]` |
| Superficial velocity variables (porous medium naming) | `superficial_vel_x`, `superficial_vel_y` |
| Isotropic Darcy coefficient vector | `ADGenericVectorFunctorMaterial` with `Ha^2` in all directions |
| Channel flow: inlet + outlet instead of closed cavity | `inlet_boundaries`, `outlet_boundaries` |
| No `pin_pressure` needed with pressure outlet | Outlet pressure fixes the reference |
| Ha^2 on Jacobian diagonal aids convergence | 2 Newton iterations to convergence |
| Porosity = 1 decouples porous-medium naming from actual porosity | `porosity = 1.0` in `[FunctorMaterials]` |

---

## Experiments to Try

### Experiment 1: Vary the Hartmann Number

Change `Ha2` at the top of the file:

```
Ha2 = 1.0    # Ha = 1: slight flattening, nearly Poiseuille
Ha2 = 25.0   # Ha = 5: clearly flat core (default)
Ha2 = 100.0  # Ha = 10: very flat core, thin layers (delta ~ 0.1)
```

For Ha = 10, the Hartmann layer is only 2-3 cells thick on the current mesh. Increase
`ny` to 50 for better resolution of the boundary layer at Ha = 10.

### Experiment 2: Compare to Poiseuille Flow

Set `Ha2 = 0.0` and remove the `friction_types` and `friction_coeffs` parameters
(or set the darcy coefficient to zero). With no magnetic braking and no-slip top/bottom
walls, the flow should converge to the parabolic Poiseuille profile in the developed
region:

```
v(y) = 6 * v_avg * y * (1 - y)
```

with peak velocity 1.5 at y = 0.5. This confirms the solver is correct before
magnetic braking is applied.

### Experiment 3: Shorten the Channel

Reduce `xmax` from 5 to 1 (keeping `nx = 10` or similar). The flow will not have
enough length to develop the Hartmann profile before reaching the outlet. The inlet
plug profile transitions to the Hartmann shape over a development length of
approximately L_dev ~ H * Re_eff; observe how far downstream the flat core extends.

### Experiment 4: Extract the Developed Profile

Add a `LineValueSampler` vector postprocessor to extract the velocity profile along
a vertical line in the developed region:

```
[VectorPostprocessors]
  [profile]
    type        = LineValueSampler
    variable    = superficial_vel_x
    start_point = '4.0 0 0'
    end_point   = '4.0 1 0'
    num_points  = 50
    sort_by     = y
  []
[]
```

The output CSV can be compared directly against the analytical Hartmann formula.
Note that the variable name is `superficial_vel_x` (not `vel_x`) because
`porous_medium_treatment = true` is active.

### Experiment 5: Increase Ha to the Plug-Flow Regime

Try Ha = 20 (`Ha2 = 400`, with `ny = 100` to resolve delta ~ 0.05). The Hartmann
layer shrinks to a very thin region near the walls and the channel core becomes
nearly perfectly uniform. This visually demonstrates why MHD is used in industrial
applications (aluminium casting, nuclear liquid-metal blankets) where a plug-flow
velocity profile suppresses turbulence and controls mixing.

### Experiment 6: Anisotropic Drag

Change the `Darcy_coefficient` to be anisotropic, for example:

```
prop_values = '${Ha2} 0 0'
```

This applies Lorentz drag only in the x-direction. Compare the resulting profile
to the isotropic case and to pure Poiseuille flow — the x-only drag still brakes
the streamwise velocity but the y-momentum equation is unaffected, altering the
pressure-velocity coupling in the development region.
