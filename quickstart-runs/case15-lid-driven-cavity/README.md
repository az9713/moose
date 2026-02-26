# Case 15: Lid-Driven Cavity — Incompressible Navier-Stokes

## Overview

This case solves the **lid-driven cavity problem**, the most widely used benchmark
in computational fluid dynamics (CFD). A square enclosure is filled with a viscous,
incompressible fluid. Three walls are stationary (no-slip), and the top wall (the
"lid") moves at a constant horizontal velocity of U = 1. The shearing motion of the
lid drives a recirculating vortex that fills the cavity.

Every CFD code in the world has been validated against this benchmark. Ghia, Ghia,
and Shin (1982) published stream-function / vorticity reference solutions at Reynolds
numbers from 100 to 10,000. At Re = 100, the flow is steady with a single dominant
vortex and two small corner eddies at the bottom corners.

This case introduces MOOSE's **finite-volume Navier-Stokes module** via the `NSFVAction`
action system. Rather than listing individual kernels, boundary conditions, and
variables by hand, a single `[Modules/NavierStokesFV]` block declares the complete
physics and the action infrastructure constructs everything automatically.

The key concepts introduced here are:

- The `[Modules/NavierStokesFV]` action block and why it simplifies incompressible
  Navier-Stokes input files significantly compared to listing every object by hand
- The incompressible constraint and pressure pinning for unique solutions
- Reynolds number as the single dimensionless parameter governing cavity flow
- Finite volume discretization of the momentum and continuity equations
- Newton's method with direct LU factorization for saddle-point systems

---

## The Physics

### Physical Problem in Plain English

Picture a square box, 1 meter on each side, filled with a viscous fluid (like light
oil). The bottom, left, and right walls are fixed and stationary. The top wall slides
to the right at speed U = 1 m/s, dragging the fluid near it along with it.

The dragged fluid cannot pile up at the right wall; it turns downward, flows along the
right wall toward the bottom, turns left along the bottom, and returns upward along the
left wall. This creates a large, roughly circular vortex that fills most of the cavity.
At the bottom corners, the sharp 90-degree turns create small secondary corner eddies
that rotate in the opposite direction to the primary vortex.

The relative importance of inertia to viscosity is measured by the Reynolds number.
At Re = 100, the flow is steady: after a transient startup, it settles into a
time-independent pattern. At Re around 8,000, the flow becomes unsteady.

### Governing Equations

Incompressible, isothermal, steady Navier-Stokes:

```
Momentum:   rho * (u . grad) u  =  -grad p  +  mu * Laplacian(u)
Continuity: div(u)  =  0
```

In component form (2D):

```
rho * (u * du/dx + v * du/dy)  =  -dp/dx  +  mu * (d^2u/dx^2 + d^2u/dy^2)
rho * (u * dv/dx + v * dv/dy)  =  -dp/dy  +  mu * (d^2v/dx^2 + d^2v/dy^2)
du/dx + dv/dy  =  0
```

Symbol definitions:

- `u = (u, v)` — velocity vector field [m/s]. u is the x-component, v is the y-component.
- `p` — pressure [Pa or dimensionless]. In incompressible flow, pressure is a
  Lagrange multiplier enforcing `div(u) = 0`, not a thermodynamic state variable.
- `rho = 1` — fluid density [kg/m^3].
- `mu = 0.01` — dynamic viscosity [Pa s].
- `rho * (u . grad) u` — the convective (inertia) term. Nonlinear in u.
- `mu * Laplacian(u)` — the diffusive (viscous) term. Linear in u.

The nonlinear inertia term is what makes Navier-Stokes hard compared to the Stokes
equations (which drop it). At Re = 100 it is ten times larger than the viscous term,
so the nonlinearity is significant.

### The Reynolds Number

The Reynolds number is the fundamental dimensionless parameter:

```
Re  =  rho * U * L / mu
```

where U is the characteristic velocity (lid speed = 1) and L is the characteristic
length (cavity width = 1):

```
Re  =  1 * 1 * 1 / 0.01  =  100
```

| Re range     | Flow character                                    |
|--------------|---------------------------------------------------|
| Re < 1       | Stokes (creeping) flow — inertia negligible        |
| Re ~ 100     | Steady laminar — single primary vortex             |
| Re ~ 1000    | Steady laminar — primary vortex + stronger corners |
| Re ~ 8000    | Onset of unsteadiness and oscillation              |
| Re > 10000   | Turbulent                                          |

At Re = 100, the primary vortex center is slightly above and to the right of the
geometric center of the cavity. As Re increases, the vortex center migrates downward
and the corner eddies grow.

### Boundary Conditions

```
y=1 (top):     u = 1, v = 0       Moving lid — no-slip, but moving
y=0 (bottom):  u = 0, v = 0       No-slip, stationary wall
x=0 (left):    u = 0, v = 0       No-slip, stationary wall
x=1 (right):   u = 0, v = 0       No-slip, stationary wall
```

The top boundary has a velocity discontinuity at the two upper corners: the lid moves
at U = 1 but the side walls are fixed at U = 0. This corner singularity is an inherent
feature of the lid-driven cavity problem. It is integrable (finite energy) and does
not prevent convergence of the numerical solution.

Pressure has no explicit boundary condition. The incompressible pressure is determined
only up to an additive constant (adding any constant C to p satisfies the equations).
The pressure level is fixed by pinning the domain-average pressure to zero.

### ASCII Domain Diagram

```
x=0                       x=1
 |  --> U=1 (lid moves right) -->  |  y=1
 +----------------------------------+
 |                                  |
 |     Primary vortex (CW)          |  u=0
u=0                                  u=0
 |   +------+                       |
 |   |      |  vortex center at    |
 |   | ~(0.6, 0.74) at Re=100     |
 |   +------+                       |
 |                                  |
 | [SW]               [SE]          |  y=0
 +----------------------------------+
         u=0, v=0 (bottom)

CW  = clockwise primary vortex
[SW] = small counter-clockwise secondary vortex at bottom-left corner
[SE] = small counter-clockwise secondary vortex at bottom-right corner

Mesh: 30 x 30 uniform quadrilateral cells (finite volume, cell-centered)
```

---

## Input File Walkthrough

The input file is `case15_lid_driven_cavity.i`.

### Block: `[Mesh]`

```
[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim  = 2
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
    nx   = 30
    ny   = 30
  []
[]
```

A 30x30 structured quadrilateral mesh on the unit square. This gives 900 cells and a
uniform cell width of dx = dy = 1/30 ≈ 0.033. At Re = 100, the boundary layers are
relatively thick and this resolution is adequate for a qualitatively correct solution.

The mesh uses `GeneratedMeshGenerator` (the generator-style syntax) rather than the
older `type = GeneratedMesh` syntax. The generator approach is required when using the
finite volume Navier-Stokes action, which expects to identify boundary sidesets by
name (`left`, `right`, `top`, `bottom`). The `GeneratedMeshGenerator` automatically
creates these named sidesets.

For higher Reynolds numbers (Re > 1000), finer meshes are needed near the walls where
velocity gradients are large. A 30x30 mesh is standard for Re = 100 validation.

### Block: `[Modules/NavierStokesFV]`

This is the central block of the input file. The `NSFVAction` action system reads
the parameters here and automatically constructs:

- FV variables: `vel_x`, `vel_y`, `pressure`
- FV kernels: advection, diffusion, pressure gradient for momentum; mass advection
  for continuity
- FV boundary conditions: no-slip walls, fixed-velocity inlet (the lid), pressure
  pinning constraint

Without the action, setting up incompressible Navier-Stokes requires approximately
15-25 individual block entries across `[Variables]`, `[FVKernels]`, `[FVBCs]`,
and additional objects. The action collapses all of this to one block.

Key parameters:

**`compressibility = 'incompressible'`**: Selects the incompressible Navier-Stokes
formulation. The continuity equation becomes `div(u) = 0`. Density does not appear
in the continuity equation, only in the momentum equation as a constant coefficient.

**`add_energy_equation = false`**: This is an isothermal simulation. No temperature
variable, no heat equation, no buoyancy coupling. The energy equation would be needed
for natural convection or heated cavity problems.

**`density = 'rho'`** and **`dynamic_viscosity = 'mu'`**: Names of functor material
properties. The actual values (rho = 1, mu = 0.01) are defined in `[FunctorMaterials]`.

**`initial_velocity = '0 0 0'`** and **`initial_pressure = 0.0`**: Starting point for
the Newton iteration. A quiescent (motionless) initial condition is standard. Newton's
method will converge from this starting point for Re = 100.

**`wall_boundaries` and `momentum_wall_types`**: The left, right, and bottom boundaries
are no-slip walls (velocity = 0 on these boundaries). The NSFVAction creates the
appropriate `INSFVNoSlipWallBC` objects automatically.

**`inlet_boundaries = 'top'`** and **`momentum_inlet_types = 'fixed-velocity'`**:
The moving lid is modeled as a fixed-velocity inlet. This is the standard MOOSE FV
approach: a wall with prescribed velocity is treated identically to an inlet with
that velocity. The `momentum_inlet_function = '1 0'` prescribes (u, v) = (1, 0).

**`momentum_advection_interpolation = 'average'`**: Selects central-difference (second-
order) interpolation for the convective flux. At Re = 100 on a 30x30 mesh, central
differencing is stable. For higher Re or coarser meshes, upwind schemes are safer.

**`pin_pressure = true`**: Activates the pressure pinning constraint. The average
pressure over the domain is constrained to `pinned_pressure_value = 0`. This removes
the pressure null space that exists in incompressible flow (pressure is only determined
up to an additive constant). Without pinning, the linear system is singular.

### Block: `[FunctorMaterials]`

```
[FunctorMaterials]
  [fluid_properties]
    type        = ADGenericFunctorMaterial
    prop_names  = 'rho mu'
    prop_values = '1   0.01'
  []
[]
```

`ADGenericFunctorMaterial` declares constant functor material properties with automatic
differentiation (AD) support. The FV Navier-Stokes system uses AD throughout to
assemble exact Jacobians for Newton's method, so AD-compatible material types are required.

The choice `rho = 1`, `mu = 0.01` gives exactly Re = 100 with U = 1 and L = 1.
To change the Reynolds number, adjust `mu`:

```
Re = 100:  mu = 0.01     (current setting)
Re = 400:  mu = 0.0025
Re = 1000: mu = 0.001
```

### Block: `[Postprocessors]`

Five postprocessors record key scalar quantities at the converged steady-state solution:

- `max_vel_x`, `min_vel_x` — the extreme x-velocity values. The maximum occurs
  near the lid; the minimum (negative) occurs in the return flow below the vortex center.
- `max_vel_y`, `min_vel_y` — the extreme y-velocity values. These occur near the
  right and left walls respectively, where the fluid turns the corners.
- `avg_pressure` — the domain-average pressure. Due to pressure pinning, this should
  be approximately zero. A non-zero value here would indicate the pinning is not
  working correctly.

### Block: `[Executioner]`

```
[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 50
[]
```

**`type = Steady`**: Directly solves the time-independent equations. No time-marching
is needed because Re = 100 cavity flow is steady. The steady executioner repeatedly
applies the nonlinear solver until convergence.

**`solve_type = 'NEWTON'`**: Full Newton method with exact Jacobian (provided by AD).
For the incompressible Navier-Stokes saddle-point system (velocity-pressure coupling),
Newton converges significantly faster than PJFNK because the coupling between momentum
and continuity is captured exactly in the Jacobian.

**`-pc_type lu -pc_factor_shift_type NONZERO`**: Direct LU factorization as the
preconditioner. For a 900-cell 2D mesh this is affordable (the system has ~2700
unknowns: 900 for vel_x, 900 for vel_y, 900 for pressure). The `NONZERO` shift
prevents factorization failure when zeros appear on the diagonal, which can happen
in pressure rows of the saddle-point system.

**`nl_rel_tol = 1e-8`** and **`nl_abs_tol = 1e-10`**: Tight tolerances appropriate
for a steady benchmark comparison. Newton typically converges in 5-15 iterations from
the quiescent initial condition.

---

## Running the Simulation

```bash
cd quickstart-runs/case15-lid-driven-cavity

# Run with the combined modules application (includes Navier-Stokes module)
/path/to/moose/modules/combined/combined-opt -i case15_lid_driven_cavity.i
```

The `combined-opt` executable includes all MOOSE modules, including the Navier-Stokes
finite volume module that provides `NSFVAction`. The framework-only `moose_test-opt`
executable does not include this module.

For a Docker-based environment:

```bash
docker run --rm -v $(pwd):/work -w /work \
  idaholab/moose:latest \
  /opt/moose/modules/combined/combined-opt -i case15_lid_driven_cavity.i
```

The simulation produces:
- `case15_lid_driven_cavity_out.e` — Exodus file with vel_x, vel_y, pressure fields
- `case15_lid_driven_cavity_out.csv` — CSV with postprocessor values at convergence

Console output will show Newton iteration progress. Expect 8-15 Newton iterations
with each iteration's linear residual dropping roughly one to two orders of magnitude.

---

## Expected Results

### Ghia et al. (1982) Benchmark Comparison

The canonical reference for lid-driven cavity flow is:

> Ghia, U., Ghia, K. N., & Shin, C. T. (1982). High-Re solutions for incompressible
> flow using the Navier-Stokes equations and a multigrid method.
> Journal of Computational Physics, 48(3), 387-411.

At Re = 100, Ghia et al. report the following u-velocity profile along the vertical
centerline (x = 0.5):

| y        | u (Ghia) | Description              |
|----------|----------|--------------------------|
| 1.000    | 1.0000   | Lid surface (boundary)   |
| 0.9766   | 0.8412   | Just below lid           |
| 0.9688   | 0.7887   |                          |
| 0.9609   | 0.7372   |                          |
| 0.9531   | 0.6872   |                          |
| 0.8516   | 0.2305   |                          |
| 0.7344   | 0.0033   | Approximate vortex center|
| 0.6172   | -0.1364  | Negative return flow     |
| 0.5000   | -0.2058  | Centerline minimum       |
| 0.4531   | -0.2109  | True minimum of u        |
| 0.2813   | -0.1566  |                          |
| 0.1719   | -0.1015  |                          |
| 0.1016   | -0.0643  |                          |
| 0.0625   | -0.0381  |                          |
| 0.0000   | 0.0000   | Bottom wall (boundary)   |

The 30x30 mesh used here is coarser than Ghia's 129x129 finite-difference grid, so
agreement will be qualitative. On the centerline, the MOOSE solution should show:
- u > 0 near the top (dragged by the lid)
- u passing through zero near y ≈ 0.74 (vortex center height)
- u < 0 in the lower half (return flow)
- u returning to zero at the bottom

### Key Flow Features

The converged solution should display:

1. **Primary vortex**: A large clockwise-rotating vortex filling most of the cavity.
   At Re = 100, the vortex center is near (x, y) ≈ (0.62, 0.74).

2. **High-velocity region near the lid**: The x-velocity peaks at U = 1 on the lid
   and decreases rapidly into the interior. On the 30x30 mesh, vel_x near y=1 should
   show values approaching 1.

3. **Pressure field**: Higher pressure at the right side of the cavity (where the lid
   drives fluid toward the right wall), lower pressure at the left side and vortex
   center. The domain-average is pinned to zero.

4. **Bottom corner eddies**: At Re = 100, there are small secondary vortices in the
   bottom-left and bottom-right corners. On a 30x30 mesh they may be only marginally
   resolved; a 64x64 or finer mesh makes them clearly visible.

### Expected Postprocessor Values

After convergence, the CSV file should show approximately:

| Quantity    | Expected value (30x30 mesh, Re=100) |
|-------------|--------------------------------------|
| max_vel_x   | ~1.0 (near the lid)                  |
| min_vel_x   | ~-0.21 (return flow, near center)    |
| max_vel_y   | ~0.17 (near the right wall)          |
| min_vel_y   | ~-0.25 (near the left wall)          |
| avg_pressure| ~0.0 (pinned by construction)        |

These values are approximate; the exact numbers depend on the mesh resolution and
interpolation scheme. Finer meshes and the Ghia grid converge toward the reference values.

---

## Key Takeaways

- **NSFVAction**: The `[Modules/NavierStokesFV]` block is a high-level action that
  constructs the complete incompressible Navier-Stokes system automatically — variables,
  kernels, and BCs — from a compact set of parameters. This is MOOSE's recommended
  approach for standard incompressible flow problems.

- **Finite Volume discretization**: Unlike earlier cases which used finite elements
  (nodal Lagrange variables), this case uses cell-centered FV. Velocities and pressure
  are stored as cell-average values, not nodal values. The `FVKernels` enforce the
  divergence theorem on each cell face.

- **Pressure pinning**: Incompressible flow determines pressure only to within an
  additive constant. The `pin_pressure = true` setting constrains the average pressure,
  making the linear system non-singular and the solution unique.

- **Moving wall via fixed-velocity inlet**: MOOSE models the moving lid by treating
  the top boundary as a fixed-velocity inlet with (u, v) = (1, 0). This is physically
  equivalent to a no-slip condition on a wall moving at that speed.

- **Newton + LU for saddle-point systems**: The incompressible Navier-Stokes equations
  have a saddle-point structure (velocity-pressure coupling). Full Newton with LU
  factorization handles this robustly. PJFNK can also work but requires a more
  sophisticated preconditioner.

- **ADGenericFunctorMaterial**: The FV Navier-Stokes action requires AD-compatible
  functor materials (`ADGenericFunctorMaterial`) rather than the standard
  `GenericConstantMaterial` used in FE cases. Functors are evaluated at cell centers
  (or faces) rather than at quadrature points.

- **Ghia et al. benchmark**: This benchmark is the standard reference for incompressible
  Navier-Stokes solvers. Comparing your simulation to Ghia's tabulated velocity profiles
  is the accepted way to verify a CFD code's incompressible flow capability.

---

## Experiments to Try

### Experiment 1: Increase the Reynolds Number

Change `mu = 0.01` to `mu = 0.001` in `[FunctorMaterials]`. This raises Re from 100
to 1000. The vortex center migrates downward and the corner eddies become larger and
more distinct. You may need to refine the mesh (increase nx and ny to 50 or 64) for
a well-resolved solution. Newton may require more iterations.

### Experiment 2: Refine the Mesh

Change `nx = 30, ny = 30` to `nx = 64, ny = 64`. Run again and compare the centerline
u-velocity profile to the Ghia table. The 64x64 result should agree with Ghia to within
a few percent. A 30x30 result gives the correct qualitative picture but not full
quantitative accuracy.

### Experiment 3: Switch Advection Interpolation

Change `momentum_advection_interpolation = 'average'` to
`momentum_advection_interpolation = 'upwind'`. Upwind differencing introduces more
numerical diffusion, which artificially reduces the effective Reynolds number. The
vortex structure will be qualitatively similar but quantitatively different (weaker
vortex, velocity profiles deviating further from Ghia). This illustrates the trade-off
between stability and accuracy in advection schemes.

### Experiment 4: Add a Point Value at the Vortex Center

Add a postprocessor to sample the pressure at the approximate vortex center:

```
[Postprocessors]
  [p_vortex_center]
    type    = PointValue
    variable = pressure
    point   = '0.62 0.74 0'
  []
[]
```

The pressure at the vortex center should be a local minimum (lower than the surrounding
region). Compare the pressure at the center to the high-pressure region near the
right wall.

### Experiment 5: Visualize in ParaView

Open `case15_lid_driven_cavity_out.e` in ParaView:

1. Color the domain by `vel_x` using a blue-white-red diverging colormap. You will
   see red (positive u) near the top/lid and blue (negative u, return flow) in the
   lower half.

2. Apply the Glyph filter (Filters > Common > Glyph) with Arrow glyphs scaled by
   the velocity magnitude. This shows the direction and speed of the flow everywhere
   in the cavity.

3. Apply the Stream Tracer filter (Filters > Common > Stream Tracer) with seed points
   along a line. You will see streamlines spiraling around the primary vortex center,
   with small tight spirals in the bottom corners.

4. Color by `pressure` to see the pressure field: high pressure near the upper-right
   corner where the lid drives fluid into the wall, low pressure at the vortex center.
