# Case 26: EHD Pumping — Coulomb Force Drives Fluid Flow

## Overview

This case demonstrates **electrohydrodynamic (EHD) pumping**: a mechanism by which
an electric body force (the Coulomb force) drives a viscous fluid into circulation
without any moving parts. The physical concept appears in Melcher's *Continuum
Electromechanics* (MIT Press, 1981), Chapter 9, Sections 9.11-9.12.

In an EHD pump, space charge (free ions or injected charge) in a dielectric liquid
experiences a Coulomb force when an electric field is applied:

```
f_Coulomb = rho_e * E
```

where rho_e is the volumetric free charge density and E is the electric field vector.
This body force enters the Navier-Stokes momentum equation in exactly the same
position as gravity or any other volume force. The fluid accelerates in response,
setting up a recirculating flow pattern determined by the force distribution.

In this case the force is **prescribed analytically** rather than computed
self-consistently from Poisson's equation and a charge-transport equation. This
simplification isolates the fluid-mechanics part of the problem — the coupling
between the Coulomb body force and the Navier-Stokes equations — and avoids the
additional complexity of solving the charge distribution simultaneously. The
prescribed-force approach is the standard starting point in EHD analysis and
directly parallels the way Melcher introduces the problem in Chapter 9.

**Connection to earlier cases:**

- **Case 15 (Lid-Driven Cavity)**: Introduced the `[Modules/NavierStokesFV]` action
  for incompressible Navier-Stokes. Case 26 uses the same action and adds a body
  force on top of it, demonstrating how to extend the action-based setup.
- **Case 16 (Natural Convection)**: Used the Boussinesq buoyancy body force, which
  is built into the action. Case 26 shows how to add an *external* body force that
  is not a built-in option of the action, using `INSFVBodyForce` explicitly.
- **Case 24 (Drift-Diffusion)**: Solved the charge transport equation in a prescribed
  field. Case 26 takes the complementary view: the charge distribution is prescribed
  and the flow response is computed.

---

## The Physics

### Governing Equations

The fluid satisfies the incompressible, isothermal Navier-Stokes equations augmented
with a prescribed body force:

```
Momentum:   rho * (v . grad) v  =  -grad p  +  mu * Laplacian(v)  +  f_Coulomb(x, y)
Continuity: div(v)  =  0
```

The Coulomb body force is prescribed as an analytical function:

```
f_x(x, y)  =  A * (1 - x) * sin(pi * y)
f_y(x, y)  =  0
```

with force amplitude A = 500 and fluid properties rho = 1, mu = 0.01.

### Physical Interpretation of the Force Shape

The factor `(1 - x)` reflects the typical spatial distribution of space charge in
an EHD pump with unipolar injection:

- At x = 0 (injection electrode): maximum force — space charge density is highest
  here, freshly injected from the electrode surface.
- At x = 1 (collecting electrode): zero force — charge has been swept away and
  collected, leaving the electrode region nearly charge-free.
- Between 0 and 1: linear decay represents the simplest approximation to the
  steady-state charge profile when charge recombination is negligible.

The factor `sin(pi * y)` sets the y-variation:

- At y = 0 (bottom wall): sin(0) = 0 — no force (consistent with no-slip boundary).
- At y = 0.5 (midheight): sin(pi/2) = 1 — maximum force drives strongest flow.
- At y = 1 (top wall): sin(pi) = 0 — no force (consistent with no-slip boundary).

Throughout the domain [0,1]^2, sin(pi * y) is non-negative, so the force is
uniformly rightward. This drives fluid from left to right across the cavity interior.
The return path runs along the top and bottom walls and down the right wall, producing
a single large recirculation loop — the EHD analog of the natural convection cell
seen in Case 16.

### EHD vs. Natural Convection

The analogy between EHD pumping and natural convection is direct:

| Natural Convection (Case 16) | EHD Pumping (Case 26) |
|------------------------------|-----------------------|
| Buoyancy body force: f = -rho*g*alpha*(T-T_ref) | Coulomb body force: f = rho_e * E |
| Force depends on temperature field | Force depends on charge density and E-field |
| Temperature field governed by energy equation | Charge field governed by charge transport |
| Both fields solved simultaneously | Force prescribed analytically (one-way) |
| Driven by heated/cooled walls | Driven by injection electrodes |

In both cases, a spatially varying body force drives recirculating flow in a
closed cavity. The key difference is that in natural convection the buoyancy
force is self-consistently coupled to the temperature field (two-way coupling),
while here the Coulomb force is prescribed (one-way coupling).

### Dimensionless Parameters

The ratio of body force to viscous force defines an EHD Reynolds number:

```
Re_EHD  ~  rho * U * L / mu
```

where U ~ A * L^2 / mu is the Stokes-flow velocity scale for a body force of
magnitude A. With A = 500, L = 1, rho = 1, mu = 0.01:

```
U_Stokes  ~  A * L^2 / mu  =  500 * 1 / 0.01  =  50000
Re_EHD    ~  rho * U * L / mu  =  1 * 50000 * 1 / 0.01  =  5,000,000
```

This extremely large estimate assumes purely viscous Stokes flow; in reality,
inertia (the convective term) limits the velocity to far smaller values. The
actual steady-state velocity magnitude will be on the order of 1-10 in
dimensionless units due to inertia. The effective Reynolds number of the
actual flow is moderate (~50-200), making Newton's method convergent.

### Boundary Conditions

```
All four walls: no-slip (vel_x = vel_y = 0)
Pressure: average value pinned to zero (removes pressure null space)
```

There are no inlet or outlet boundaries — the cavity is fully enclosed. The
EHD body force provides the only source of momentum. The pressure pinning
removes the free constant that incompressible pressure always has in a
closed cavity.

### Domain Diagram

```
x=0                                x=1
 |                                  |  y=1
 +----------------------------------+
 |  no-slip top wall (vel=0)        |
 |                                  |
 |  --> --> f_x = A(1-x)sin(pi*y)   |
 |  --> --> (rightward body force)  |
 |  --> --> -->  (return via walls) |
 |                                  |
 |  no-slip bottom wall (vel=0)     |
 +----------------------------------+  y=0

Force is strongest at x=0 (left), zero at x=1 (right).
Force is zero at y=0 and y=1, maximum at y=0.5.
Drives one large clockwise recirculation loop.

Mesh: 25 x 25 uniform quadrilateral cells (finite volume, cell-centered)
```

---

## Input File Walkthrough

The input file is `case26_ehd_pumping.i`.

### Top-Level Variable

```
force_amp = 500.0   # Coulomb force amplitude
```

HIT allows scalar variables to be defined at the file top level and referenced
anywhere using `${force_amp}`. Changing this single value scales the entire body
force, making it easy to explore different force magnitudes without editing the
mathematical expression.

### Block: `[Mesh]`

A 25x25 structured quadrilateral mesh on the unit square. This gives 625 cells
and a uniform cell width of dx = dy = 1/25 = 0.04. The resolution is adequate to
capture the single recirculation vortex for the present force amplitude.
`GeneratedMeshGenerator` creates the named sidesets `left`, `right`, `top`,
`bottom` required by the NavierStokesFV action.

### Block: `[Modules/NavierStokesFV]`

This action block creates all FV variables (`vel_x`, `vel_y`, `pressure`), all
FV kernels (momentum advection, momentum diffusion, pressure gradient, mass
conservation), Rhie-Chow interpolation, and all wall boundary conditions
automatically.

Key parameters for this case:

**`compressibility = 'incompressible'`**: Selects the incompressible formulation.
The continuity equation is div(v) = 0 with no density variation.

**`add_energy_equation = false`**: This is an isothermal case. No temperature
variable is needed; the EHD force drives the flow without any thermal coupling.

**`wall_boundaries` and `momentum_wall_types`**: All four walls are no-slip.
There is no moving lid or inlet; all fluid motion is driven entirely by the
body force.

**`momentum_advection_interpolation = 'upwind'`**: Upwind differencing for the
convective term. More dissipative than central differencing but more stable when
the body force creates strong velocity gradients. At moderate resolution the
extra dissipation is acceptable.

**`pin_pressure = true`**: Required in a fully enclosed cavity. Without pinning,
the pressure is determined only up to a constant and the linear system is singular.

### Block: `[FVKernels]` — The EHD Body Force

```
[FVKernels]
  [ehd_force_x]
    type                  = INSFVBodyForce
    variable              = vel_x
    functor               = coulomb_force_x
    momentum_component    = 'x'
    rhie_chow_user_object = 'ins_rhie_chow_interpolator'
  []
[]
```

`INSFVBodyForce` adds a prescribed volumetric body force to the finite-volume
momentum equation. The key parameters:

**`variable = vel_x`**: This kernel contributes to the x-momentum equation.
A separate kernel would be needed for `variable = vel_y` if there were a y-component
of the force (there is none in this case).

**`functor = coulomb_force_x`**: The name of a functor property that provides
the force value at each cell face or cell center where it is evaluated. The functor
is defined in `[FunctorMaterials]` as an `ADParsedFunctorMaterial`.

**`momentum_component = 'x'`**: Tells `INSFVBodyForce` which momentum component
it is acting on. This is needed for the Rhie-Chow interpolation correction.

**`rhie_chow_user_object = 'ins_rhie_chow_interpolator'`**: The Rhie-Chow
interpolation user object created automatically by the `NavierStokesFV` action.
`INSFVBodyForce` must register with this user object so that the body force is
included in the Rhie-Chow velocity correction (which prevents pressure-velocity
decoupling on collocated meshes). The fixed name `ins_rhie_chow_interpolator` is
always used by the `NavierStokesFV` action.

This pattern — using the `NavierStokesFV` action for the base NS physics and then
adding `INSFVBodyForce` kernels in a separate `[FVKernels]` block — is the standard
MOOSE approach for any body force that is not a built-in option of the action (such
as a custom electric force, magnetic force, or any other spatially varying source).

### Block: `[FunctorMaterials]` — Fluid Properties and Force Field

Two functor materials are defined:

**`fluid_props` (`ADGenericFunctorMaterial`)**: Declares the constant material
properties `rho = 1.0` and `mu = 0.01`. These are AD-compatible functor materials
required by the `NavierStokesFV` action.

**`ehd_force` (`ADParsedFunctorMaterial`)**: Defines the analytical Coulomb force
functor. The key parameter is:

```
expression = '${force_amp} * (1 - x) * sin(3.14159265358979 * y)'
```

`ADParsedFunctorMaterial` evaluates an algebraic expression at every point where
the functor is queried (cell centers, face centers, quadrature points). The spatial
coordinates `x`, `y`, `z` and the time `t` are always available as built-in symbols.
The `${force_amp}` reference substitutes the top-level HIT variable at parse time.

The result is an AD-compatible functor named `coulomb_force_x` that `INSFVBodyForce`
can consume. No additional `functor_names` or `functor_symbols` are needed because
the expression only uses the built-in spatial coordinates and the substituted
constant — no other functor properties are referenced within the expression.

### Block: `[Postprocessors]`

Four postprocessors record key scalar quantities at the converged steady state:

- **`max_vel_x`**: Maximum x-velocity — confirms rightward flow is being driven.
- **`min_vel_x`**: Most negative x-velocity — measures the return flow magnitude.
- **`max_vel_y`**: Maximum y-velocity — measures the turning flow near the walls.
- **`avg_pressure`**: Domain-average pressure — should remain near zero due to
  pinning; confirms the constraint is active.

`ADElementExtremeFunctorValue` is used for velocity because `vel_x` and `vel_y`
are FV functor quantities (not standard nodal FE variables). For the FE pressure
variable, `ElementAverageValue` is used directly.

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

**`type = Steady`**: Finds the steady-state solution directly. EHD-driven flow at
the chosen force amplitude is steady (no time oscillations), making this efficient.

**`solve_type = NEWTON`**: Full Newton method with exact AD Jacobians. For the
incompressible NS saddle-point system with an added body force, Newton converges
reliably when starting from a quiescent initial condition.

**LU factorization**: Direct solve appropriate for the 25x25 system (~1875 cells,
~5625 unknowns). The NONZERO diagonal shift prevents factorization failures in the
pressure rows of the saddle-point block.

---

## Running the Simulation

Run using the MOOSE `combined-opt` executable, which includes the Navier-Stokes
finite-volume module:

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case26-ehd-pumping \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case26_ehd_pumping.i 2>&1 | tail -30'
```

The `combined-opt` binary includes all MOOSE modules. The `MSYS_NO_PATHCONV=1`
environment variable prevents MSYS/Git Bash on Windows from mangling the
absolute path in the `-v` mount argument.

The simulation produces:
- `case26_ehd_pumping_out.e` — Exodus file with `vel_x`, `vel_y`, `pressure` fields
- `case26_ehd_pumping_out.csv` — CSV file with postprocessor values at convergence

Console output will show Newton iteration progress. Expect 10-30 Newton iterations
before the residuals satisfy `nl_rel_tol = 1e-8`.

---

## Expected Results

### Flow Pattern

The converged solution should show a single large clockwise recirculation loop:

```
x=0                                x=1
 |                                  |  y=1
 +------------ top wall ------------+
 |                                  |
 | -->-->-->-->-->-->  (rightward)  |
 |                           |      |
 |                    return |      |
 |                       flow|      |
 | <--<--<--<--<--<--        |      |
 |                           v      |
 +------------ bot wall ------------+  y=0
```

Because the force `f_x = A*(1-x)*sin(pi*y)` drives fluid rightward throughout the
domain (it is always positive in [0,1]^2), fluid moves right through the interior.
The no-slip right wall deflects it upward and downward; return flow runs back along
the top and bottom walls. The result is one large recirculation cell (as opposed to
natural convection which has flow turning near the hot/cold walls).

### Velocity Field

- **`vel_x`**: Positive (rightward) in the interior. Strongest near x = 0 where
  the Coulomb force is maximum, weakening toward x = 1.
- **`vel_y`**: Positive (upward) near the right half and negative (downward) near
  the left half as fluid turns at the walls.
- The velocity magnitude will be much smaller than the Stokes estimate because
  inertia redistributes and limits the driven flow. Typical values of max |vel_x|
  in the range 1-10 (dimensionless) are expected.

### Pressure Field

The pressure field develops from the incompressibility constraint and the body force.
Higher pressure will appear near the right wall (where fluid driven rightward
stagnates) and lower pressure near the upper-left and lower-left corners.

### Expected Postprocessor Values

After convergence, the CSV file should show approximately:

| Quantity     | Expected value (25x25, A=500) |
|--------------|-------------------------------|
| max_vel_x    | ~1-10 (rightward flow peak)   |
| min_vel_x    | ~0 to slightly negative       |
| max_vel_y    | ~1-5 (turning flow at walls)  |
| avg_pressure | ~0.0 (pinned by construction) |

The exact values depend on the force amplitude and mesh resolution. Doubling
`force_amp` roughly doubles the velocities (in the Stokes-flow limit).

---

## Key Takeaways

- **`INSFVBodyForce` for external body forces**: When the `NavierStokesFV` action
  does not natively support a particular body force (e.g., Coulomb, Lorentz, or
  other electromagnetic forces), add it as an explicit `INSFVBodyForce` kernel in
  `[FVKernels]`. The action creates the base NS system; the manual kernel adds the
  extra physics.

- **`rhie_chow_user_object` is mandatory**: Any `INSFVBodyForce` kernel added
  outside the action must reference the Rhie-Chow interpolation user object
  (`ins_rhie_chow_interpolator`) created by the action. Omitting this reference
  causes the body force to be excluded from the pressure-velocity coupling
  correction, leading to incorrect results or solver divergence.

- **`ADParsedFunctorMaterial` for analytical force fields**: Analytical body forces
  expressed as functions of position can be defined concisely using
  `ADParsedFunctorMaterial`. The spatial coordinates `x`, `y`, `z` and time `t` are
  always available as built-in symbols. No custom C++ object is needed for simple
  prescribed functions.

- **EHD as the electromagnetic analog of natural convection**: In natural convection,
  a buoyancy body force (driven by temperature gradients) drives recirculating flow
  in a closed cavity. In EHD pumping, a Coulomb body force (driven by charge density
  and electric field) does the same thing. The Navier-Stokes equations and the MOOSE
  input file structure are identical; only the source of the body force differs.

- **Prescribed vs. self-consistent forcing**: This case uses a prescribed force
  (one-way coupling: force drives flow, but flow does not modify the charge or field).
  The fully coupled EHD problem requires solving Poisson's equation for the electric
  potential and a charge transport equation simultaneously with the Navier-Stokes
  equations — a natural extension of Case 26 combining elements from Cases 15, 16,
  and 24.

---

## Experiments to Try

### Experiment 1: Vary the Force Amplitude

Change `force_amp` at the top of the input file:

```
force_amp = 100.0   # Weaker force, slower flow
force_amp = 1000.0  # Stronger force, faster flow
```

At small amplitudes, the flow is in the Stokes regime (linear in force amplitude).
At large amplitudes, inertia becomes important and the relationship between force
amplitude and maximum velocity becomes sublinear. Track the change in `max_vel_x`
to observe this transition.

### Experiment 2: Increase Viscosity (Reduce Effective Reynolds Number)

Change `mu = 0.01` to `mu = 0.1` in `[FunctorMaterials]`. This decreases the
effective Reynolds number by ten and pushes the flow toward the Stokes (creeping
flow) limit. In this limit, the velocity should scale linearly with A/mu, and the
streamline pattern becomes smoother and more symmetric.

### Experiment 3: Add a Y-Component Body Force

Add a second `INSFVBodyForce` kernel for the y-momentum:

```
[FVKernels]
  [ehd_force_y]
    type                  = INSFVBodyForce
    variable              = vel_y
    functor               = coulomb_force_y
    momentum_component    = 'y'
    rhie_chow_user_object = 'ins_rhie_chow_interpolator'
  []
[]
```

and define a y-force in `[FunctorMaterials]`:

```
[ehd_force_y]
  type          = ADParsedFunctorMaterial
  property_name = coulomb_force_y
  expression    = '${force_amp} * sin(3.14159265 * x) * (1 - y)'
[]
```

This creates a y-force analogous to the x-force (strongest near the bottom, decaying
upward, sinusoidal in x). With both x and y forces present, the flow pattern becomes
more complex with two interacting vortices.

### Experiment 4: Refine the Mesh

Change `nx = 25, ny = 25` to `nx = 50, ny = 50`. The finer mesh resolves sharper
velocity gradients near the walls more accurately. Compare the `max_vel_x` values
between the two resolutions; when they agree to within ~5%, the solution is
mesh-converged.

### Experiment 5: Visualize in ParaView

Open `case26_ehd_pumping_out.e` in ParaView:

1. Color the domain by `vel_x`. You will see the rightward-driven interior flow
   (red/warm colors) and the return flow near the walls (blue/cool colors).

2. Apply a Glyph filter (Arrows) to visualize the velocity direction. The arrows
   should clearly show the single large recirculation loop.

3. Apply a Stream Tracer filter with seed points distributed over the domain.
   Streamlines will spiral around the single vortex center.

4. Color by `pressure`. Higher pressure near the right wall (where the EHD force
   pushes fluid into the wall) and lower pressure near the left wall (where fluid
   is drawn away from the wall by the force).
