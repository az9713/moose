# Case 43: Ekman Spiral — Rotating Boundary Layer

## Overview

This case models the **Ekman spiral** — the classic solution for steady viscous
flow near a solid wall in a uniformly rotating frame. The Coriolis acceleration
couples the two horizontal velocity components, so momentum diffused downward from
the geostrophic interior is continuously rotated as it reaches the surface. The
result is that the velocity vector traces a spiral in the horizontal plane as the
wall is approached: pointing in the geostrophic direction far from the wall, rotating
clockwise (in the Northern Hemisphere convention) with decreasing height, and
vanishing at the no-slip surface.

The physics follows Rieutord, *Fluid Dynamics* (Springer, 2015), Chapter 8,
Section 8.4, one of the clearest modern derivations of the Ekman layer. The same
problem appears in Pedlosky, *Geophysical Fluid Dynamics* (Springer, 1987), Chapter 4,
as the canonical boundary-layer solution for geostrophic flow.

Two MOOSE features are central to this case:

1. **ADMatDiffusion kernel**: Implements the Laplacian term `nu * d²v/dz²` for each
   velocity component. Using `ADGenericConstantMaterial` for the viscosity property
   provides automatic differentiation support for the Jacobian.

2. **CoupledForce kernel**: Adds a term proportional to a *different* variable to
   the equation for the current variable. Here, the Coriolis force in the vx equation
   depends on vy and vice versa. Two `CoupledForce` kernels wire this cross-coupling,
   with carefully chosen signs derived from the residual convention.

**Connection to earlier cases:**

- **Case 09 (Coupled Reaction-Diffusion)**: Introduced `CoupledForce` for coupling
  two scalar equations. Case 43 uses the same mechanism but for a momentum equation
  where the coupling coefficient has physical meaning (Coriolis parameter).
- **Case 27 (Hartmann Flow)**: Another fluid mechanics case with a body force that
  modifies the standard Poiseuille profile. Case 43 is simpler (no pressure solver,
  no Navier-Stokes action) because the Ekman problem reduces to two coupled 1D ODEs.
- **Case 13 (Custom Kernel)**: Demonstrated how MOOSE kernels contribute to the
  weak-form residual. Case 43 is a study in composing the correct residual from
  standard MOOSE primitives rather than writing custom code.

---

## The Physics

### Physical Setup

Consider the lower boundary layer of a large-scale atmospheric or oceanic flow. Far
above the surface (the "geostrophic interior"), the flow is in balance between the
horizontal pressure gradient and the Coriolis force, giving a horizontal velocity
of magnitude U_g directed in the x-direction. Near the surface, viscous friction
slows the flow, breaking the geostrophic balance and allowing the Coriolis force to
deflect the velocity vector.

The problem is steady and one-dimensional: everything varies only with z (height
above the surface). The two horizontal velocity components vx and vy obey:

```
nu * d²vx/dz² + 2*Omega*vy = 0
nu * d²vy/dz² - 2*Omega*vx = -2*Omega*U_g
```

where:
- `nu` — kinematic viscosity [m²/s]
- `Omega` — rotation rate of the reference frame [rad/s] (Earth: 7.27e-5 rad/s;
  here set to 1.0 for a dimensionless problem)
- `U_g` — geostrophic velocity far from the wall [m/s], set to 1.0
- `2*Omega` — the Coriolis parameter f (here f = 2.0)

The Ekman layer thickness is:

```
delta_E = sqrt(nu / Omega)
```

With nu = 0.01 and Omega = 1.0, delta_E = 0.1. This is the length scale over which
the velocity vector rotates from zero at the wall to nearly (but not exactly) the
geostrophic value.

### Analytical Solution

The system of two coupled second-order ODEs has the exact solution:

```
vx(z) = U_g * (1 - exp(-z/delta_E) * cos(z/delta_E))
vy(z) = U_g * exp(-z/delta_E) * sin(z/delta_E)
```

Both solutions vanish at z = 0 (no-slip wall). As z → infinity, vx → U_g and
vy → 0, recovering the geostrophic flow. The exponential decay means the layer is
essentially fully established by z ~ 5*delta_E = 0.5; the domain extends to
z = 0.7 = 7*delta_E, well into the geostrophic region.

Key validation points:

| Quantity | Formula | Value |
|---|---|---|
| vx at z = delta_E | 1 - e^{-1}*cos(1) | 0.8011 |
| vy at z = delta_E | e^{-1}*sin(1) | 0.3096 |
| max(vy) | at z = pi/4 * delta_E = 0.0785 | 0.3224 |
| vx at z = 5*delta_E | 1 - e^{-5}*cos(5) | 0.9981 |

### The Ekman Spiral in the Hodograph Plane

Plotting (vx, vy) as z increases from 0 to infinity traces a spiral in velocity
space. At z = 0: (0, 0). As z increases, the tip of the velocity vector spirals
outward and clockwise, reaching peak vy = 0.322 at z = pi/4 * delta_E = 0.0785,
then turning back. The spiral passes above the geostrophic point (vx > 1, vy > 0)
before crossing the vx-axis at z = pi*delta_E = 0.314, where vy = 0 again and
vx overshoots U_g by about 4.3% (vx = 1.043). The spiral then decays inward,
converging asymptotically to (U_g, 0) = (1, 0) as z → infinity.

This spiral is the defining feature of the Ekman layer and the reason the solution
is named an "Ekman spiral" rather than simply an "Ekman profile." It does not appear
in any scalar diffusion problem — it is intrinsically a two-component phenomenon
arising from the Coriolis coupling.

### ASCII Domain Sketch

```
z = 0.7  +-------------------------------------------+  vx = U_g = 1, vy = 0
(7*dE)   |  geostrophic interior: vx ~ 1, vy ~ 0    |  (far-field Dirichlet BC)
         |                                           |
z = 0.5  |  Ekman spiral nearly complete             |
(5*dE)   |  vx ~ 0.998, vy ~ -0.006 (slight undershoot)|
         |                                           |
z = 0.3  |  vx OVERSHOOT: spiral past geostrophic    |
(3*dE)   |  vx ~ 1.049, vy ~ 0.007                  |
         |                                           |
z = 0.1  |  z = delta_E                              |
(1*dE)   |  vx ~ 0.801, vy ~ 0.310                  |
         |                                           |
z = 0    +-------------------------------------------+  vx = 0, vy = 0
         x-axis = MOOSE x-coordinate (represents z) (no-slip Dirichlet BC)

Domain:  [0, 0.7] x [0, 0.01]  (quasi-1D: 200 elements in x, 1 element in y)
Mesh:    GeneratedMesh, 200 x 1 uniform quadrilateral elements
Note:    MOOSE x-coordinate plays the role of the physical z (height) variable.
         The y-direction is a dummy dimension of width 0.01 (irrelevant physically).
         vx overshoots U_g = 1 between z ~ 0.16 and z ~ 0.5 — this is physical.
```

---

## Kernel Sign Derivation

Getting the signs right in a coupled system requires care. This section works through
the derivation explicitly so that the input file can be verified independently.

### MOOSE Residual Convention

MOOSE assembles the residual R = 0, where:
- `ADMatDiffusion` with diffusivity D contributes `+integral(D * grad(u) . grad(phi)) dV`
  to the residual (weak form). In strong form this is `-D * d²u/dz²`.
- `CoupledForce` with coupled variable v and coefficient coef contributes
  `-coef * integral(v * phi) dV` to the residual. In strong form this is `-coef * v`.
- `BodyForce` with value F contributes `-F * integral(phi) dV` to the residual.
  In strong form this is `-F`.

### Equation for vx

Strong-form equation: `nu * d²vx/dz² + 2*Omega*vy = 0`

Rewrite as residual (move all terms to one side): `-nu * d²vx/dz² - 2*Omega*vy = 0`

Map to kernel contributions:

```
-nu * d²vx/dz²          ->  ADMatDiffusion(diffusivity=nu)         contributes -nu*d²vx/dz²
-coef_vx * vy           ->  CoupledForce(v=vy, coef=coef_vx)       contributes -coef_vx*vy
```

Matching: `-coef_vx = -2*Omega`, so `coef_vx = 2*Omega = 2.0`. The `[vx_coriolis]`
kernel uses `coef = ${coriolis}` where `coriolis = 2.0`.

### Equation for vy

Strong-form equation: `nu * d²vy/dz² - 2*Omega*vx = -2*Omega*U_g`

Rewrite as residual: `-nu * d²vy/dz² + 2*Omega*vx - 2*Omega*U_g = 0`

Map to kernel contributions:

```
-nu * d²vy/dz²          ->  ADMatDiffusion(diffusivity=nu)              contributes -nu*d²vy/dz²
-coef_vy * vx           ->  CoupledForce(v=vx, coef=coef_vy)            contributes -coef_vy*vx
-body_val               ->  BodyForce(value=body_val)                   contributes -body_val
```

Matching: `-coef_vy = +2*Omega`, so `coef_vy = -2*Omega = -2.0`. The `[vy_coriolis]`
kernel uses `coef = '${fparse -2.0 * Omega}'`.

Matching: `-body_val = -2*Omega*U_g`, so `body_val = 2*Omega*U_g = 2.0`. The
`[vy_pressure]` kernel uses `value = ${body_force_val}` where `body_force_val = 2.0`.

---

## Input File Walkthrough

The input file is `case43_ekman_spiral.i`.

### Top-Level HIT Variables

```
nu = 0.01
Omega = 1.0
U_g = 1.0

coriolis = '${fparse 2.0 * Omega}'
body_force_val = '${fparse 2.0 * Omega * U_g}'
```

Five variables at the top of the file control all physics parameters. The `${fparse ...}`
expressions compute derived values at parse time. Changing `Omega` updates all coupling
coefficients consistently. Changing `nu` adjusts both the material viscosity and (through
delta_E = sqrt(nu/Omega)) the implicit length scale of the spiral.

### Block: `[Mesh]`

```
[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = 0
  xmax = 0.7
  ymin = 0
  ymax = 0.01
  nx = 200
  ny = 1
[]
```

The MOOSE x-coordinate plays the role of the physical vertical coordinate z. The
domain extends from z = 0 (wall) to z = 0.7 = 7*delta_E, placing the far-field
boundary deep inside the geostrophic region where the exponential terms are negligible
(exp(-7) = 0.0009). The y-direction is a dummy dimension of width 0.01 with a single
element — it exists only because MOOSE requires `dim = 2` for a structured mesh,
but all physics are independent of y.

With 200 elements over xmax = 0.7, the element size is dx = 0.0035. The Ekman layer
has thickness delta_E = 0.1, which is resolved by 0.1 / 0.0035 ≈ 29 elements — more
than sufficient for a polynomial FEM solution. The mesh is uniform because the
analytical solution decays smoothly without any sharp sub-layer features.

The `GeneratedMesh` automatically creates sidesets named `left` (x=0, the wall) and
`right` (x=0.7, the far field), used in the `[BCs]` block.

### Block: `[Variables]`

```
[Variables]
  [vx]
  []
  [vy]
  []
[]
```

Two scalar unknowns are declared. Each uses the default first-order Lagrange shape
functions with zero initial conditions. Both are solved simultaneously in a coupled
system: the Jacobian has four blocks (dvx/dvx, dvx/dvy, dvy/dvx, dvy/dvy), with the
off-diagonal blocks provided by the `CoupledForce` kernels.

Total DOFs: 2 variables × (200+1) nodes × (1+1) nodes = 2 × 201 × 2 = 804 DOFs.
The system is small enough that a direct LU solver (`-pc_type lu`) is efficient and
avoids any preconditioning tuning.

### Block: `[Kernels]`

Five kernels implement the two coupled ODEs.

**Kernels for the vx equation** (`nu*d²vx/dz² + 2*Omega*vy = 0`):

- `[vx_diff]`: `type = ADMatDiffusion, variable = vx, diffusivity = nu_mat` — viscous
  diffusion of vx. Uses the AD variant to provide exact Jacobian contributions from
  the diffusivity material property.

- `[vx_coriolis]`: `type = CoupledForce, variable = vx, v = vy, coef = ${coriolis}`
  — Coriolis coupling. The kernel adds `-coef * vy` to the residual for vx, which
  in strong form represents `+coef * vy = +2*Omega*vy` on the left-hand side, matching
  the Coriolis term in the governing equation.

**Kernels for the vy equation** (`nu*d²vy/dz² - 2*Omega*vx + 2*Omega*U_g = 0`):

- `[vy_diff]`: `type = ADMatDiffusion, variable = vy, diffusivity = nu_mat` — viscous
  diffusion of vy. Identical in structure to `vx_diff`.

- `[vy_coriolis]`: `type = CoupledForce, variable = vy, v = vx, coef = '${fparse -2.0 * Omega}'`
  — Coriolis coupling with negative coefficient. The kernel adds `-(-2*Omega)*vx =
  +2*Omega*vx` to the residual, representing `-2*Omega*vx` when moved to the strong-form
  left-hand side — consistent with the `-2*Omega*vx` term in the vy equation.

- `[vy_pressure]`: `type = BodyForce, variable = vy, value = ${body_force_val}` — the
  geostrophic pressure gradient forcing. This constant term (`2*Omega*U_g = 2.0`)
  represents the pressure gradient that drives the geostrophic interior flow. In the
  absence of this term (and with the far-field BC changed), the solution would
  trivially be vx = vy = 0.

### Block: `[BCs]`

```
[BCs]
  [vx_wall]   type = DirichletBC, variable = vx, boundary = left,  value = 0   []
  [vy_wall]   type = DirichletBC, variable = vy, boundary = left,  value = 0   []
  [vx_farfield] type = DirichletBC, variable = vx, boundary = right, value = ${U_g} []
  [vy_farfield] type = DirichletBC, variable = vy, boundary = right, value = 0  []
[]
```

Four Dirichlet conditions completely specify the problem:

- **Wall (x=0)**: No-slip — both velocity components vanish. This is the physical
  boundary condition at the solid rotating surface.

- **Far field (x=0.7)**: The geostrophic velocity — vx = U_g = 1, vy = 0. By placing
  this boundary at 7*delta_E, the exponential terms in the analytical solution contribute
  less than 0.1% to the field, so the Dirichlet far-field BC is essentially exact.

There are no Neumann conditions: the natural BC (zero normal flux) would be imposed
at the top and bottom of the 2D strip (y=0 and y=1), but since the physics are
uniform in y, these vanish automatically.

### Block: `[Materials]`

```
[Materials]
  [diffusivity]
    type = ADGenericConstantMaterial
    prop_names = 'nu_mat'
    prop_values = '${nu}'
  []
[]
```

A single material property `nu_mat = 0.01` provides the kinematic viscosity to both
`[vx_diff]` and `[vy_diff]`. Using `ADGenericConstantMaterial` (rather than the
non-AD `GenericConstantMaterial`) ensures that the `ADMatDiffusion` kernels receive
a material property with AD type information, allowing exact Jacobian assembly.

### Block: `[Postprocessors]`

```
[Postprocessors]
  [vx_at_delta]   type = PointValue, variable = vx, point = '0.1 0.005 0'  []
  [vy_at_delta]   type = PointValue, variable = vy, point = '0.1 0.005 0'  []
  [max_vy]        type = ElementExtremeValue, variable = vy, value_type = max  []
  [avg_vx]        type = ElementAverageValue, variable = vx  []
[]
```

Four postprocessors validate the solution against the analytical result:

- `[vx_at_delta]`: vx at z = delta_E = 0.1. Analytical value: 0.8011.
- `[vy_at_delta]`: vy at z = delta_E = 0.1. Analytical value: 0.3096.
- `[max_vy]`: Peak vy over the entire domain. Analytical value: 0.3224 at z = 0.0785 (= pi/4 * delta_E).
- `[avg_vx]`: Domain-average vx. Approaches U_g = 1 as most of the domain is in the
  geostrophic regime.

The point `'0.1 0.005 0'` lies at x = delta_E = 0.1, y = 0.005 (midpoint of the
0.01-wide strip), z = 0. MOOSE will find the element containing this point and
interpolate the solution.

### Block: `[Executioner]`

```
[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  nl_max_its = 20
[]
```

The problem is linear (all kernels are linear in vx and vy), so Newton's method
converges in a single iteration. LU direct factorisation (`-pc_type lu`) is ideal for
the small 804-DOF system and eliminates any concern about iterative solver convergence.
The tight tolerances (`1e-10`, `1e-12`) confirm that the solution is numerically exact
to nearly machine precision.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case43-ekman-spiral \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case43_ekman_spiral.i 2>&1 | tail -20'
```

The `MSYS_NO_PATHCONV=1` environment variable prevents MSYS/Git Bash on Windows from
mangling the absolute path in the `-v` mount argument.

The simulation is linear and steady-state: it completes in under 2 seconds. Outputs:

- `case43_ekman_spiral_out.e` — Exodus file with the `vx` and `vy` fields over the mesh
- `case43_ekman_spiral_out.csv` — CSV file with the four postprocessor values

---

## Expected Results

### Postprocessor Values

| Postprocessor | Numerical (200-element mesh) | Analytical | Notes |
|---|---|---|---|
| `vx_at_delta` | ~0.8011 | 0.8011 | vx at z = delta_E = 0.1 |
| `vy_at_delta` | ~0.3096 | 0.3096 | vy at z = delta_E = 0.1 |
| `max_vy` | ~0.3224 | 0.3224 | Peak vy in the spiral |
| `avg_vx` | ~0.929 | ~0.929 | Average over [0, 0.7]; includes overshoot region |

The numerical solution should match the analytical values to 4 or more significant
figures with the 200-element mesh. The problem is linear with smooth exponential
solutions, so the FEM error converges rapidly with mesh refinement.

### Velocity Profiles Along z

```
z/dE    vx(z)     vy(z)     comments
----    -----     -----     --------
0.0     0.000     0.000     no-slip wall
0.5     0.468     0.291     vy still rising toward maximum
pi/4    0.678     0.322     z = pi/4 * dE = 0.0785  (MAXIMUM vy)
1.0     0.801     0.310     z = delta_E
1.5     0.984     0.223     vx approaching geostrophic
2.0     1.056     0.123     vx OVERSHOOTS U_g (spiral past geostrophic direction)
3.0     1.049     0.007     spiral turning back toward geostrophic
pi      1.043     0.000     vy crosses zero, vx ~4% above geostrophic
5.0     0.998    -0.006     slight undershoot, essentially geostrophic
7.0     ~1.000    ~0.000    far-field (matches BC)
```

(dE = delta_E = 0.1; z/dE = x/0.1 in the MOOSE mesh)

### The Ekman Spiral (Hodograph)

Plotting vy vs vx as a parametric curve with z as the parameter traces the spiral:

```
vy
^
0.32 |   *  (max vy at z = pi/4*dE = 0.0785, vx ~ 0.644)
     | *   *
0.20 |*       *
     |           *
0.10 |               *
     |                    *
0.00 +*-----------+----*---->  vx
    0.0          1.0  1.05
    (wall)   (geostrophic) (overshoot)

The curve starts at (0,0), rises steeply to max vy ~ 0.322 (small z),
curves past the geostrophic point, overshoots to vx ~ 1.043 at z = pi*dE,
then spirals back to (1, 0) as z → infinity.
```

---

## Understanding the Plots

Running `python visualize_all.py` from the `quickstart-runs/` directory produces
plots for this case.

### Velocity Profiles (vx and vy vs z)

**What the plot shows.** Two curves plotted against z (MOOSE x-coordinate). One curve
(blue) shows vx rising from 0 at the wall, overshooting U_g = 1 around z ~ 0.2, then
settling back to 1.0 at the far-field. The second curve (orange) shows vy starting at 0,
rising steeply to a peak near z = 0.0785, then decaying back toward 0.

**How to judge correctness.** Both curves should start exactly at zero (no-slip BC).
vx should overshoot U_g = 1.0 before converging back (this overshoot is physical, not
a numerical artifact). vy should show a single clear peak at approximately z = 0.0785
(= pi/4 * delta_E) with maximum value approximately 0.322. The vx overshoot should
reach approximately 1.056 at z ~ 0.2. If vy has no peak or vx shows no overshoot,
the Coriolis coupling coefficient has the wrong sign or magnitude.

**What would indicate a problem.**
- vx not reaching 1.0 at the right boundary: the far-field BC is not applied.
- vy strongly negative (below -0.01): a Coriolis coupling sign is wrong. Note that a
  small negative vy (order -0.006 at z = 5*dE) is physically correct — the spiral
  slightly undershoots zero before converging. Large negative values are unphysical.
- Both profiles identical to the trivial solution (vx proportional to z, vy = 0): the
  coupling kernels are absent or have zero coefficient.
- No vx overshoot above 1.0: if vx rises monotonically without exceeding 1.0, the
  Coriolis coupling is too weak or the domain is too short to exhibit the overshoot.

### Hodograph (vy vs vx)

**What the plot shows.** A parametric curve in the (vx, vy) plane, with z increasing
from 0 to 0.7 along the curve. The curve starts at the origin, curves upward and to
the right, reaches a maximum vy, then turns back toward (1, 0).

**How to judge correctness.** The curve should form a clean inward spiral. The starting
point must be exactly (0, 0) and the endpoint (at z = 0.7 = 7*delta_E) should be very
close to (1, 0). The curve crosses the vx-axis at z = pi*delta_E = 0.314 (inside the
domain) with vx ~ 1.043 — this crossing is physical and expected. The spiral then
continues above vx = 1.0 briefly before converging back to (1, 0).

**What would indicate a problem.**
- A straight line (vy = 0 everywhere): the Coriolis coupling is missing.
- A loop that goes into negative vy: the sign of one Coriolis term is flipped.
- The curve not approaching (1, 0) at the far end: either the domain is too short
  (increase xmax) or the far-field BC is incorrect.

---

## Key Takeaways

| Concept | Where in Input |
|---|---|
| Two coupled 1D ODEs solved monolithically | `[Variables]` with two entries |
| ADMatDiffusion for viscous terms | `[vx_diff]`, `[vy_diff]` kernels |
| CoupledForce for cross-variable Coriolis coupling | `[vx_coriolis]`, `[vy_coriolis]` |
| Negative CoupledForce coefficient for a coupling that acts as a source | `coef = '${fparse -2.0 * Omega}'` |
| BodyForce for constant geostrophic pressure gradient | `[vy_pressure]` |
| ADGenericConstantMaterial required for ADMatDiffusion | `[diffusivity]` |
| Sign derivation: CoupledForce adds -coef*v to residual | See sign derivation section |
| Quasi-1D: 200 x 1 mesh, x represents the physical z-coordinate | `[Mesh]` |
| Direct LU solver for small linear system | `-pc_type lu` |
| PointValue postprocessor for analytical comparison | `[vx_at_delta]`, `[vy_at_delta]` |

---

## Experiments to Try

### Experiment 1: Vary the Rotation Rate (Ekman Number)

The Ekman number is `Ek = nu / (Omega * L^2)` where L is a reference length. In this
non-dimensional problem, varying `Omega` with fixed `nu` and `U_g` changes delta_E:

```
Omega = 0.25  =>  delta_E = sqrt(0.01/0.25) = 0.2  (thicker layer, slower rotation)
Omega = 1.0   =>  delta_E = sqrt(0.01/1.0)  = 0.1  (default)
Omega = 4.0   =>  delta_E = sqrt(0.01/4.0)  = 0.05 (thinner layer, faster rotation)
```

For Omega = 4.0, the Ekman layer is only 0.05 thick and requires `nx = 400` or more to
resolve. Observe how the spiral becomes tighter (more rotations fit in the same domain)
and the transition from wall to geostrophic is sharper.

### Experiment 2: Check Against the Analytical Solution

Add a `[Functions]` block with the analytical expressions and use
`ElementL2Error` postprocessors to compute the L2 norm of the error:

```
[Functions]
  [vx_analytical]
    type = ParsedFunction
    expression = 'Ug * (1 - exp(-x/dE) * cos(x/dE))'
    symbol_names = 'Ug dE'
    symbol_values = '1.0 0.1'
  []
  [vy_analytical]
    type = ParsedFunction
    expression = 'Ug * exp(-x/dE) * sin(x/dE)'
    symbol_names = 'Ug dE'
    symbol_values = '1.0 0.1'
  []
[]

[Postprocessors]
  [error_vx]
    type = ElementL2Error
    variable = vx
    function = vx_analytical
  []
  [error_vy]
    type = ElementL2Error
    variable = vy
    function = vy_analytical
  []
[]
```

Run with `nx = 50, 100, 200, 400` and record the errors. First-order Lagrange
elements on a uniform mesh give second-order convergence in L2: doubling the mesh
should reduce the error by a factor of 4.

### Experiment 3: Extend to the Transient Problem

Replace the `Steady` executioner with `Transient` and add `TimeDerivative` kernels
to both vx and vy. Start from zero initial conditions (both at rest) and watch the
Ekman spiral spin up from rest. The spin-up time scale is approximately `1/Omega = 1.0`.
Use `end_time = 5.0` and `dt = 0.05` to capture the full spin-up.

Expected behavior: initially vx grows faster than vy (the first effect of the
pressure gradient driving vx). As vx grows, the Coriolis force deflects momentum into
the vy direction. The system oscillates (inertial oscillations at frequency Omega)
while the Ekman layer develops, eventually damping to the steady spiral.

### Experiment 4: Two-Layer Ekman Problem (Surface and Bottom)

Modify the domain to have no-slip walls on both sides:

```
xmin = -0.5, xmax = 0.5
```

with Dirichlet BC `vx = vy = 0` at both x = -0.5 and x = 0.5, and a geostrophic
interior in the middle driven only by the BodyForce. This models the Ekman layers on
both the top and bottom boundaries of an ocean layer. The solution is the sum of two
mirror-image Ekman spirals that merge in the interior. If the domain is much thicker
than delta_E, the spirals are independent; if comparable, they interact.

### Experiment 5: Vary the Geostrophic Velocity Direction

The geostrophic flow need not be in the x-direction. Change the far-field BC to
`vx = U_g * cos(theta), vy = U_g * sin(theta)` and change the BodyForce to reflect
the new pressure gradient direction. The Ekman spiral will still develop but oriented
at angle theta. Since the equations are linear, the solution is simply the original
solution rotated by theta in the hodograph plane.

Expected outcome: the spiral in the (vx, vy) plane rotates rigidly by angle theta,
demonstrating the rotational symmetry of the Ekman problem.
