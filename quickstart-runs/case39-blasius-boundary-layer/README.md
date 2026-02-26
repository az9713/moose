# Case 39: Blasius Boundary Layer — Flat Plate Laminar Flow

## Overview

This case simulates **laminar boundary layer growth over a flat plate**, one of the
most important exact (similarity) solutions in fluid mechanics. Uniform flow at speed
U approaches a semi-infinite flat plate from the left. Viscosity causes the fluid
at the plate surface to stick (no-slip), and this retardation propagates outward into
the flow through viscous diffusion. The result is a growing **shear layer** whose
thickness increases as sqrt(x) along the plate.

The problem was solved analytically by Heinrich Blasius in 1908 by recognizing that the
velocity profile u(x, y) / U depends on position only through the dimensionless
**similarity variable** eta = y / delta(x), where delta(x) is the local boundary layer
thickness. This reduces the partial differential equations to a single nonlinear ordinary
differential equation that can be solved by numerical shooting. Every CFD code targeting
external aerodynamics is validated against the Blasius solution.

Key concepts introduced here:

- External flow with an inlet, outlet, no-slip wall, and slip (symmetry) boundary
- The boundary layer concept and the sqrt(x) growth law
- Blasius similarity solution and the skin friction coefficient C_f = 0.664 / sqrt(Re_x)
- Mesh biasing (clustering) toward the wall to resolve steep velocity gradients
- How incompressible FV Navier-Stokes handles open-domain external flow

**Reference**: Rieutord, M., *Fluid Dynamics* (Springer, 2015), Ch. 4, Sec. 4.3.
Original analysis: Blasius, H. (1908), *Grenzschichten in Fluessigkeiten mit kleiner Reibung*, Z. Math. Phys., 56, 1-37.

---

## The Physics

### Boundary Layer Growth on a Flat Plate

Consider uniform flow at speed U over a flat plate aligned with the x-axis. Far from
the plate, the flow is undisturbed (u = U, v = 0). At the plate surface (y = 0), the
no-slip condition requires u = v = 0. In between, viscosity communicates the effect of
the wall into the interior, creating a transition region — the boundary layer — where
u increases from 0 to approximately U.

The boundary layer thickness delta(x) is defined as the height at which u/U = 0.99.
The Blasius solution gives:

```
delta_99(x) = 5.0 * sqrt(mu * x / (rho * U))
            = 5.0 * x / sqrt(Re_x)

where  Re_x = rho * U * x / mu  (local Reynolds number based on x)
```

The factor 5.0 comes from numerically integrating the Blasius ODE. The boundary layer
grows as sqrt(x): thicker further downstream, because the wall has had more distance
(and therefore time) to decelerate the fluid.

### Governing Equations

Steady, incompressible, 2D Navier-Stokes:

```
Momentum (x):   rho * (u du/dx + v du/dy) = -dp/dx + mu * (d^2u/dx^2 + d^2u/dy^2)
Momentum (y):   rho * (u dv/dx + v dv/dy) = -dp/dy + mu * (d^2v/dx^2 + d^2v/dy^2)
Continuity:     du/dx + dv/dy = 0
```

For a flat plate with zero pressure gradient (dp/dx = 0, which holds for the Blasius
case), the boundary layer approximation reduces this to:

```
Boundary layer equation:  u du/dx + v du/dy = mu/rho * d^2u/dy^2
Continuity:               du/dx + dv/dy = 0
```

The wall-normal pressure gradient is negligible inside the boundary layer, so dp/dy = 0
and pressure is uniform across the layer thickness. This case solves the full Navier-Stokes
equations (not the boundary layer approximation), which is the correct approach for
a finite-volume CFD simulation.

### Blasius Similarity Solution

Blasius recognized that the velocity profile is **self-similar**: at every x location,
u(x, y) / U looks the same when plotted against the scaled variable

```
eta = y * sqrt(rho * U / (mu * x))  =  y / delta(x) * const
```

Introducing the stream function psi = sqrt(mu * U * x / rho) * f(eta), the Navier-Stokes
equations reduce to the **Blasius ODE**:

```
2 f''' + f * f'' = 0

Boundary conditions:
  f(0)   = 0    (no wall-normal velocity at y=0: v = 0 at plate... actually v|_wall is satisfied)
  f'(0)  = 0    (no-slip: u = 0 at y=0)
  f'(inf) = 1   (free stream: u/U = 1 far from plate)
```

The solution f'(eta) gives the universal velocity profile. Key values:

```
f''(0) = 0.4696  (determines wall shear stress)
f'(eta) = 0.99   is reached at eta_99 ~ 5.0
```

### Important Derived Quantities

**Skin friction coefficient** at position x:

```
C_f(x) = tau_w / (0.5 * rho * U^2)
        = 0.664 / sqrt(Re_x)
```

where tau_w = mu * (du/dy)|_{y=0} is the wall shear stress. The 0.664 comes from
the numerical value f''(0) * sqrt(2) from the Blasius solution.

**Displacement thickness** (how much the streamlines are displaced outward):

```
delta_1(x) = 1.72 * sqrt(mu * x / (rho * U))  =  1.72 * x / sqrt(Re_x)
```

**Momentum thickness** (related to the drag on the plate):

```
delta_2(x) = 0.664 * sqrt(mu * x / (rho * U))  =  0.664 * x / sqrt(Re_x)
```

**von Karman momentum integral** relates these: d(delta_2)/dx = C_f / 2, which is
satisfied by the Blasius solution.

### Parameters for This Simulation

```
rho    = 1.0       kg/m^3   (density)
U      = 1.0       m/s      (free-stream speed)
mu     = 0.005     Pa s     (dynamic viscosity)
L      = 2.0       m        (plate length)

Re_L   = rho * U * L / mu = 400   (well within laminar regime, transition ~ 5e5)

delta_99 at x = 1.0 m:  5.0 * sqrt(0.005 * 1.0 / 1.0)  = 0.354 m
delta_99 at x = 2.0 m:  5.0 * sqrt(0.005 * 2.0 / 1.0)  = 0.500 m
```

The flow is laminar throughout (Re_L = 400 << transition Re ~ 5e5).

### Boundary Conditions

```
x=0  (left, inlet):   u = U = 1,  v = 0       Uniform free stream
x=2  (right, outlet): p = 0                   Fixed pressure outlet
y=0  (bottom, plate): u = 0, v = 0            No-slip flat plate
y=1.5 (top):          slip (zero normal flux)  Approximates free stream at top
```

The top boundary uses a **slip condition** (zero wall-normal velocity, zero shear
stress) rather than no-slip, because the top represents the undisturbed free stream —
there is no physical wall there. Using a slip boundary avoids imposing an incorrect
no-slip condition that would add a spurious second boundary layer from above.

The outlet uses a fixed-pressure boundary (p = 0), which allows fluid to exit freely.
This fixes the pressure level and avoids the need for additional pressure pinning.

### Domain and Mesh

```
Domain: [0, 2] x [0, 1.5]  (length x height in meters)
Cells:  60 x 30 (streamwise x wall-normal)
Bias:   bias_y = 0.5 — cells cluster toward y = 0 (the plate)
```

The bias_y = 0.5 parameter makes cell heights decrease geometrically toward the plate,
providing much finer resolution where velocity gradients are steepest. Without wall
clustering, the coarse cells near the plate would fail to resolve the sharp velocity
gradient du/dy|_{y=0} that determines the wall shear stress.

At x = 1.0 with 30 biased cells over [0, 1.5], the smallest cells near the plate have
heights on the order of 0.02 m, which gives approximately 17 cells inside the boundary
layer at that location (delta_99 = 0.354 m). This is adequate for qualitative
verification.

### ASCII Domain Diagram

```
y=1.5  +------------------------------------------+
       |  slip (free stream — no shear)           |  top
       |  u ~ 1.0 throughout this region           |
       |                                           |
       |    ~~~ boundary layer grows ~~~           |
       |         (shear region)                    |
       |                                           |
y=0    +==========================================+
       |  no-slip flat plate (u=v=0)              |  bottom
x=0                                            x=2
(inlet: u=1, v=0)                   (outlet: p=0)

Mesh: 60 cells in x, 30 cells in y with bias_y=0.5 (fine near y=0)
```

---

## Input File Walkthrough

The input file is `case39_blasius_boundary_layer.i`.

### HIT Top-Level Variables

```
mu_val  = 0.005
rho_val = 1.0
U_inf   = 1.0
```

Three top-level scalar parameters that appear throughout the file via `${...}`
substitution. Changing `mu_val` alone rescales the Reynolds number without touching
any other block. For example, `mu_val = 0.002` gives Re_L = 1000.

### Block: `[Mesh]`

```
[Mesh]
  [gen]
    type   = GeneratedMeshGenerator
    dim    = 2
    xmin   = 0
    xmax   = 2
    ymin   = 0
    ymax   = 1.5
    nx     = 60
    ny     = 30
    bias_y = 0.5
  []
[]
```

The `bias_y = 0.5` parameter makes each successive row of cells (from top to bottom)
0.5 times the height of the row above it. This creates a geometric sequence of cell
heights with the smallest cells adjacent to the flat plate (y = 0). The ratio 0.5 is
aggressive; values between 0.7 and 0.9 are more common in engineering practice but
a stronger bias is acceptable here because the boundary layer is thick relative to
the domain height at Re_L = 400.

The domain extends to y = 1.5, well above the maximum expected delta_99 = 0.5 at
x = 2, ensuring the slip top boundary does not distort the boundary layer profile.

### Block: `[Modules/NavierStokesFV]`

The NSFV action handles all variable creation, kernel assembly, and boundary condition
setup. The key boundary specification for flat plate flow is:

**Inlet (`left`)**: `momentum_inlet_types = 'fixed-velocity'` with
`momentum_inlet_functors = '1 0'`. This prescribes (u, v) = (U, 0) on the left face
of every boundary cell. The profile is uniform (plug flow), which is the correct
far-upstream condition before any boundary layer has formed.

**Outlet (`right`)**: `momentum_outlet_types = 'fixed-pressure'` with
`pressure_functors = '0'`. The outlet pressure is fixed at zero, allowing the
velocity to adjust freely. This is the correct condition for a subsonic external flow
exit: no velocity is prescribed, only the pressure level. The convective terms in the
momentum equation handle the outflow automatically.

**Bottom wall (`bottom`)**: `momentum_wall_types = 'noslip'`. The flat plate imposes
(u, v) = (0, 0) through the no-slip condition. The NSFV action creates
`INSFVNoSlipWallBC` objects for both velocity components on this boundary.

**Top boundary (`top`)**: `momentum_wall_types = 'slip'`. The slip (free-slip or
symmetry) condition imposes zero normal velocity (v = 0 at y = 1.5) and zero tangential
stress (du/dy = 0 at y = 1.5). This is the correct model for the free-stream far field
where there is no physical wall.

**Advection scheme**: `momentum_advection_interpolation = 'upwind'`. Upwind (first-order)
differencing is more stable for advection-dominated flow. Near the plate leading edge
and inside the thin initial boundary layer, local cell Reynolds numbers can be O(10-100),
making upwind the safer choice. Second-order central differencing (`average`) would
also work here but is more susceptible to oscillations.

**No pressure pinning**: Because the outlet already fixes p = 0 via the Dirichlet
pressure boundary condition, the pressure is uniquely determined. Adding `pin_pressure`
would be redundant (and potentially over-constraining).

### Block: `[FunctorMaterials]`

```
[FunctorMaterials]
  [fluid_props]
    type        = ADGenericFunctorMaterial
    prop_names  = 'rho  mu'
    prop_values = '1.0  0.005'
  []
[]
```

`ADGenericFunctorMaterial` provides constant functor properties with automatic
differentiation support. The NSFV action requires AD-compatible materials to assemble
exact Jacobians for Newton iteration.

### Block: `[Postprocessors]`

Six postprocessors track the key flow features:

- `max_vel_x`: Maximum x-velocity. Should remain near U = 1.0 in the free stream.
  A value significantly above 1.0 would indicate flow acceleration from the displacement
  effect of the boundary layer (possible near the plate trailing edge).

- `max_vel_y`: Maximum wall-normal velocity. The boundary layer displaces fluid outward,
  inducing a small positive v in the free stream. Blasius theory gives
  v_e / U ~ 0.86 * sqrt(mu / (rho * U * x)) at the boundary layer edge.

- `avg_vel_x`: Domain-average x-velocity. Decreases slightly below 1.0 because the
  growing boundary layer (with u < U) occupies an increasing fraction of the domain.

- `u_at_x1_y035` and `u_at_x1_y010`: Point samples at x = 1. The Blasius solution
  predicts u/U = 0.99 at y = delta_99 = 0.354 and u/U ≈ 0.42 at y ≈ 0.10.

- `u_at_x2_y050` and `u_at_x2_y020`: Point samples at x = 2. Blasius: u/U = 0.99
  at y = 0.5 and u/U ≈ 0.63 at y = 0.20 (eta ≈ 2.0, f'(2.0) ≈ 0.63).

### Block: `[Executioner]`

```
[Executioner]
  type       = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 30
  automatic_scaling = true
[]
```

Newton with direct LU factorization. The 60x30 mesh produces approximately 5400 cells
with 3 degrees of freedom each (vel_x, vel_y, pressure), giving a ~16200-unknown
system. LU factorization is affordable at this size and provides reliable convergence
for the moderately nonlinear problem at Re_L = 400.

`automatic_scaling = true` balances the velocity and pressure unknowns in the Newton
system, which improves convergence when the advective terms are significant.

---

## Running the Simulation

```bash
cd quickstart-runs/case39-blasius-boundary-layer

# Run with the combined modules application (includes Navier-Stokes module)
/path/to/moose/modules/combined/combined-opt -i case39_blasius_boundary_layer.i
```

For a Docker-based environment:

```bash
docker run --rm -v $(pwd):/work -w /work \
  idaholab/moose:latest \
  /opt/moose/modules/combined/combined-opt -i case39_blasius_boundary_layer.i
```

The simulation produces:

- `case39_blasius_boundary_layer_out.e` — Exodus file with vel_x, vel_y, pressure
  fields for ParaView visualization
- `case39_blasius_boundary_layer_out.csv` — CSV with postprocessor values for
  direct comparison with Blasius theory

Console output shows Newton iteration residuals. Expect 10-20 Newton iterations
from the uniform initial condition, with residuals decreasing monotonically.

---

## Expected Results

### Velocity Field

The converged solution should show:

1. **Boundary layer growth**: vel_x contours near y = 0 show a thickening low-velocity
   region as x increases. The delta_99 contour (u = 0.99 U) should trace approximately
   the curve y = 5 * sqrt(0.005 * x), growing from 0 at the inlet to 0.5 m at x = 2.

2. **Wall-normal velocity**: vel_y is small but non-zero above the plate. The boundary
   layer displaces fluid outward, driving a weak upward velocity in the free stream.
   The maximum v/U is on the order of 1 / sqrt(Re_L) ~ 0.05, so about 5% of U_inf.

3. **Pressure field**: Nearly uniform (dp/dx ~ 0) for flow over a flat plate. There
   is a very slight favorable pressure gradient near the inlet (the inlet condition
   forces a sharp transition from uniform to developing flow) and the outlet fixes p = 0.

### Blasius Profile Comparison

The most important validation is that the velocity profiles at different x stations
collapse onto the Blasius similarity solution when plotted against eta = y / delta(x).

Blasius f'(eta) values (from numerical integration of 2f''' + f f'' = 0):

```
eta     f'(eta) = u/U   Description
------  ---------------  -----------------------------
0.0     0.000            Wall (no-slip)
0.5     0.166            Deep inside boundary layer
1.0     0.330            ~1/3 of boundary layer thickness
1.5     0.487
2.0     0.630            eta_50 ~ 1.95 (u/U = 0.50)
2.5     0.751
3.0     0.846
3.5     0.913
4.0     0.956
4.5     0.980
5.0     0.992            delta_99 definition
```

To check similarity collapse in ParaView:
1. Use the Plot Over Line filter at x = 0.5, 1.0, 1.5, 2.0
2. Extract u/U vs. y
3. Scale y by 1 / delta_99(x) = sqrt(rho * U / (mu * x)) / 5.0
4. All profiles should overlay on the Blasius curve above

### Expected Postprocessor Values

After convergence, the CSV file should show approximately:

```
Quantity          Expected value    Notes
max_vel_x         ~1.01 to 1.05     Slight acceleration from displacement effect
max_vel_y         ~0.03 to 0.05     Weak upward flow from boundary layer growth
avg_vel_x         ~0.90 to 0.95     Reduced by boundary layer mass deficit
u_at_x1_y035      ~0.99             At delta_99 of x=1 — should be near free stream
u_at_x1_y010      ~0.35 to 0.45     Inside boundary layer; Blasius gives ~0.42 at eta=1.4
u_at_x2_y050      ~0.99             At delta_99 of x=2 — should be near free stream
u_at_x2_y020      ~0.55 to 0.65     Inside boundary layer; Blasius gives ~0.63 at eta=2.0
```

The point values depend on how well the mesh resolves the profile at those locations.
Finer meshes (higher nx and ny) give results closer to the Blasius reference values.

### Skin Friction Coefficient

The Blasius skin friction coefficient is:

```
C_f(x) = tau_w / (0.5 * rho * U^2)  =  0.664 / sqrt(Re_x)
```

At x = 1 (Re_x = 200):  C_f = 0.664 / sqrt(200) = 0.0470
At x = 2 (Re_x = 400):  C_f = 0.664 / sqrt(400) = 0.0332

The wall shear stress tau_w = mu * (du/dy)|_{y=0} can be estimated from the
numerical solution by examining the velocity gradient in the first cell above the plate.

---

## Key Takeaways

- **Slip vs. no-slip boundary**: The top boundary uses `slip` rather than `noslip`.
  Slip (zero normal velocity, zero shear stress) is the correct model for a free-stream
  boundary or symmetry plane where no physical wall exists. Using no-slip on the top
  would create a spurious second boundary layer and corrupt the solution.

- **Wall-biased mesh**: The `bias_y` parameter clusters cells near the flat plate where
  velocity gradients are largest. Without wall clustering, a uniform mesh would need
  many more cells to resolve the boundary layer accurately. Mesh biasing is essential
  in boundary layer simulations.

- **Fixed-pressure outlet**: For external flow with a well-defined downstream exit,
  `fixed-pressure` outlet at p = 0 is the standard choice. It allows the velocity
  profile to develop naturally and exit without forcing any particular velocity
  distribution. This avoids the pressure pinning needed for closed cavities.

- **sqrt(x) growth law**: The Blasius boundary layer thickens as delta ~ sqrt(x),
  which follows from dimensional analysis: the only relevant length scale in the problem
  is sqrt(nu * x / U) where nu = mu/rho is the kinematic viscosity. This scaling
  is universal for all laminar flat plate boundary layers.

- **Self-similarity**: The Blasius solution is a **similarity solution** — profiles
  at different x locations are geometrically similar when the wall-normal coordinate is
  scaled by the local boundary layer thickness. Similarity solutions exist when the
  governing equations admit a group symmetry, here the scaling (x, y) -> (lambda^2 x,
  lambda y). Plotting u/U vs. eta collapses all profiles onto a single curve.

- **Validation significance**: The Blasius boundary layer is the foundation of external
  aerodynamics, skin friction drag prediction, heat transfer correlations (Pohlhausen
  solution for Pr ~ 1 fluids), and turbulent boundary layer models (which use the
  Blasius profile as the laminar sub-layer foundation). Every serious CFD code for
  external flows is validated against it.

- **Re_L = 400 regime**: The simulation uses Re_L = 400, which is comfortably laminar
  (transition typically occurs around Re_L = 5e5 for natural transition, lower with
  disturbances). The Blasius solution is exact for Re_x << 5e5.

---

## Experiments to Try

### Experiment 1: Change Reynolds Number

Reduce `mu_val` to change the Reynolds number and observe how the boundary layer thins:

```
mu_val = 0.005   =>  Re_L = 400,  delta_99(L) = 0.500 m  (current)
mu_val = 0.002   =>  Re_L = 1000, delta_99(L) = 0.316 m
mu_val = 0.001   =>  Re_L = 2000, delta_99(L) = 0.224 m
```

A thinner boundary layer at higher Re requires finer wall-normal resolution (increase ny
or reduce bias_y). Verify that the postprocessor values at delta_99 remain near 0.99.

### Experiment 2: Measure the Blasius Profile

Add a line of PointValue postprocessors at x = 1 to sample the u-velocity profile at
multiple y heights:

```
[Postprocessors]
  [u_y000]
    type  = PointValue
    variable = vel_x
    point = '1.0 0.00 0'
  []
  [u_y010]
    type  = PointValue
    variable = vel_x
    point = '1.0 0.10 0'
  []
  [u_y020]
    type  = PointValue
    variable = vel_x
    point = '1.0 0.20 0'
  []
  # ... continue at 0.30, 0.40, 0.50, 0.60 m
[]
```

Then compute eta = y * sqrt(rho * U / (mu * x)) for each point (at x=1 with
rho=1, U=1, mu=0.005: eta = y / 0.0707) and compare u/U to the Blasius f'(eta) table.

### Experiment 3: Visualize the Similarity Profile in ParaView

1. Open the Exodus output file in ParaView
2. Color by `vel_x` to see the boundary layer growing from left to right
3. Apply Plot Over Line at x = 0.5, 1.0, 1.5, 2.0 (lines from y=0 to y=1.5)
4. For each line, the velocity profile should look qualitatively similar but thicker
   at larger x
5. Rescale the y-axis by 1/delta_99(x) at each station to collapse the profiles

### Experiment 4: Apply a Streamwise Pressure Gradient

Modify the outlet pressure or add a body force in x to impose dp/dx != 0. This
produces the **Falkner-Skan** family of boundary layers. An adverse pressure gradient
(dp/dx > 0) causes the boundary layer to grow faster and can lead to **separation**
(flow reversal near the wall). A favorable pressure gradient (dp/dx < 0) thins the
boundary layer. To impose a favorable gradient, change outlet pressure to a negative
value, e.g.:

```
pressure_functors = '-0.1'
```

Observe how the boundary layer thickness changes compared to the zero-gradient Blasius case.

### Experiment 5: Refine the Mesh

Increase `nx = 120` and `ny = 60` (doubling in each direction). Run both the original
and refined meshes and compare the PointValue postprocessors. The refined mesh should
give values closer to the Blasius analytical results. The difference between the two
meshes quantifies the **discretization error** at the coarser resolution.

### Experiment 6: Compute Wall Shear Stress

The Blasius skin friction coefficient C_f = 0.664/sqrt(Re_x) can be compared to the
numerical wall shear stress. Add an `INSFVWallFunctionSideForce` or a boundary flux
postprocessor to estimate tau_w = mu * (du/dy)|_{y=0} at the plate. Compare:

```
tau_w_Blasius(x=1) = C_f(x=1) * 0.5 * rho * U^2 = 0.0470 * 0.5 * 1 * 1 = 0.0235 Pa
tau_w_Blasius(x=2) = C_f(x=2) * 0.5 * rho * U^2 = 0.0332 * 0.5 * 1 * 1 = 0.0166 Pa
```
