# Case 14: Thermoelasticity — Heated Plate with Thermal Stress

## Overview

This case couples two distinct physics: heat conduction and linear elasticity.
A 2D steel plate is held at different temperatures on its left and right edges.
The resulting temperature gradient causes the material to expand unevenly, which
generates internal stress and visible displacement — even though no mechanical
force is applied.

This is a one-way coupling: the temperature field drives deformation, but the
deformation does not feed back to change the temperature. This is the standard
"thermal stress analysis" pattern used throughout mechanical engineering.

New concepts introduced in this case:

- **`Physics/SolidMechanics/QuasiStatic` action**: the modern MOOSE block that
  sets up all solid mechanics kernels, variables, and strain materials in one place.
- **Eigenstrains**: a general mechanism for "stress-free" strains — deformations
  that would occur freely if the material were unconstrained, but that generate
  stress when the material is held by boundaries or by mechanical compatibility.
- **`ADHeatConduction` + `ADComputeThermalExpansionEigenstrain`**: the AD-aware
  objects that connect the thermal and mechanical physics.
- **`combined-opt` executable**: required because this case uses both the
  `heat_transfer` and `solid_mechanics` physics modules.

---

## The Physics

### The Physical Problem in Plain English

Imagine a one-metre square steel plate lying flat. The left edge is clamped to a
hot pipe (500 K); the right edge is in contact with a cold reservoir (300 K). The
top and bottom edges are insulated (no heat flows through them). After a long time,
the temperature reaches a steady linear gradient — hotter on the left, cooler on
the right.

Steel expands when it gets hot. The hot left side wants to be longer than the cold
right side. But the plate is a single connected solid: every element must remain
in contact with its neighbours. This mechanical compatibility constraint is what
generates stress. Where the material is forced to be shorter than it wants (hot
regions restrained from expanding), compressive stress builds up. Where it is
forced to be longer than it would be if free, tensile stress appears.

The bottom edge is pinned (zero displacement), which prevents rigid-body motion
but still allows the plate to deform freely above it.

### Governing Equations

**Heat conduction (steady state):**

```
-div( k * grad(T) ) = 0    in Omega = [0,1] x [0,1]

T = 500 K    on the left boundary  (x = 0)
T = 300 K    on the right boundary (x = 1)
```

With insulated top and bottom (natural Neumann condition: zero flux), the steady
solution is a linear temperature ramp:

```
T(x) = 500 - 200*x    (K)
```

**Linear elasticity with thermal eigenstrain:**

```
-div( sigma ) = 0    in Omega

sigma = C : (epsilon - epsilon_th)

epsilon_th = alpha * (T - T_ref) * I
```

where:
- `sigma` — Cauchy stress tensor (Pa)
- `C` — fourth-order isotropic elasticity tensor (from E = 200 GPa, nu = 0.3)
- `epsilon` — total small strain tensor = sym(grad(u))
- `epsilon_th` — thermal eigenstrain (stress-free expansion)
- `alpha = 12e-6 /K` — coefficient of thermal expansion
- `T_ref = 300 K` — reference (stress-free) temperature
- `I` — identity tensor (isotropic expansion)
- `u = (disp_x, disp_y)` — displacement vector

The eigenstrain represents the deformation the material would undergo if it were
unconstrained. The actual mechanical strain is `epsilon - epsilon_th`; stress is
proportional to that mechanical strain via Hooke's law.

### Domain Diagram

```
T=500 K (hot)             T=300 K (cold)
    |                         |
    |    Steel Plate           |
    |    E = 200 GPa           |
    |    nu = 0.3              |
    |    alpha = 12e-6 /K      |
    |    k = 50 W/(m K)        |
    |                         |
    +-------- bottom --------+
         disp_x = disp_y = 0   (pinned)

Temperature at steady state:
  T(x) = 500 - 200*x  (linear ramp from left to right)

Thermal expansion:
  Hot left  (T=500, delta_T = +200 K) wants to expand most
  Cold right (T=300, delta_T =   0 K) has no thermal expansion
  --> left side pushes right; plate bends upward (net disp_y > 0)
```

---

## Input File Walkthrough

The input file is `case14_thermoelasticity.i`.

### Header Parameters

```
E     = 200e9
nu    = 0.3
alpha = 12e-6
k_th  = 50
cp    = 500
T_ref = 300
```

Named parameters at the top allow changing material properties in one place.
They are referenced throughout as `${E}`, `${nu}`, etc.

### `[GlobalParams]`

```
displacements = 'disp_x disp_y'
```

Declares the displacement variable names once so that every block that needs them
(the action, BCs, postprocessors) picks them up automatically.

### `[Variables]`

Only `T` is declared here. The displacement variables `disp_x` and `disp_y` are
generated automatically by the `Physics/SolidMechanics/QuasiStatic` action
(`add_variables = true`), so they do not need to appear in `[Variables]`.

The `FunctionIC` for `T` sets an initial guess of `300 + 200*(1-x)`, which
matches the expected linear ramp and helps the Newton solver start close to the
answer.

### `[Kernels]`

```
[heat_conduction]
  type     = ADHeatConduction
  variable = T
[]
```

`ADHeatConduction` implements the weak form of `-div(k*grad(T))`. The `AD` prefix
means automatic differentiation is used to build the Jacobian exactly. It reads the
`thermal_conductivity` material property provided by `ADHeatConductionMaterial`.

There is no mechanics kernel here because the `Physics/SolidMechanics/QuasiStatic`
action generates those automatically.

### `[Physics/SolidMechanics/QuasiStatic]`

```
[solid]
  strain             = SMALL
  add_variables      = true
  eigenstrain_names  = 'thermal_eigenstrain'
  generate_output    = 'vonmises_stress'
  use_automatic_differentiation = true
[]
```

This single block replaces what would otherwise be dozens of lines of manual
kernel, variable, and strain material declarations:

- `strain = SMALL` — uses the linearised small-strain tensor `(grad(u) + grad(u)^T)/2`.
  Appropriate for metals at moderate temperatures where strains are well below 1%.
- `add_variables = true` — the action creates `disp_x` and `disp_y` variables.
- `eigenstrain_names` — the list of eigenstrains to subtract from total strain before
  computing stress. Must match the `eigenstrain_name` in the material block below.
- `generate_output = 'vonmises_stress'` — the action creates an `AuxVariable` and
  `AuxKernel` to compute and store the von Mises equivalent stress at every node.
- `use_automatic_differentiation = true` — all residuals and Jacobians are computed
  via AD, enabling the `NEWTON` solver to use an exact Jacobian.

### `[BCs]`

Four boundary conditions are set:

| Name | Variable | Boundary | Value | Purpose |
|------|----------|----------|-------|---------|
| `T_hot` | T | left | 500 K | Hot wall |
| `T_cold` | T | right | 300 K | Cold wall |
| `pin_bottom_x` | disp_x | bottom | 0 | No sliding |
| `pin_bottom_y` | disp_y | bottom | 0 | No lifting |

The top boundary has no displacement BC, so it is free to move. The left and right
boundaries also have no displacement BCs — the plate can expand laterally at those
edges, which partially relieves the thermal stress.

### `[Materials]`

Four material objects work together:

1. **`ADHeatConductionMaterial`**: provides `thermal_conductivity = 50` and
   `specific_heat = 500`. The conductivity is consumed by `ADHeatConduction`.

2. **`ADComputeIsotropicElasticityTensor`**: builds the 4th-order stiffness tensor
   `C` from `E = 200 GPa` and `nu = 0.3`.

3. **`ADComputeThermalExpansionEigenstrain`**: computes `epsilon_th = alpha*(T-T_ref)*I`
   at every quadrature point. The `temperature = T` parameter couples it to the
   solution variable. The result is named `thermal_eigenstrain`, which must appear
   in the action's `eigenstrain_names` list.

4. **`ADComputeLinearElasticStress`**: computes `sigma = C : (epsilon - epsilon_th)`.
   Depends on the strain computed by the action's auto-generated `ADComputeSmallStrain`
   material and on the eigenstrain from (3).

### `[Postprocessors]`

Five postprocessors track key scalar results:

| Name | Type | Variable | Meaning |
|------|------|----------|---------|
| `max_temperature` | ElementExtremeValue (max) | T | Peak temperature |
| `avg_temperature` | ElementAverageValue | T | Spatial average T |
| `max_disp_x` | ElementExtremeValue (max) | disp_x | Max horizontal expansion |
| `max_disp_y` | ElementExtremeValue (max) | disp_y | Max vertical displacement |
| `max_vonmises` | ElementExtremeValue (max) | vonmises_stress | Peak stress |

### `[Executioner]`

```
type       = Steady
solve_type = 'NEWTON'
petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
petsc_options_value = 'lu       mumps'
```

`Steady` solves the fully coupled `(T, disp_x, disp_y)` system once. `NEWTON` with
an LU preconditioner is appropriate because:
- The problem is linear (linear elasticity + linear heat conduction), so Newton
  converges in one iteration.
- The system is small enough (about 3 x 21 x 21 = 1323 DOFs) that a direct solver
  is fast and robust.

---

## Running the Simulation

This case requires the `combined-opt` executable because it uses the `heat_transfer`
and `solid_mechanics` modules. From a Docker container:

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/you/moose-next/quickstart-runs/case14-thermoelasticity:/work" \
  -w /work \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case14_thermoelasticity.i'
```

Expected console output:

```
Mesh Information:
  Spatial dimension: 2
  Mesh: 400 elements, 441 nodes

Nonlinear System:
  Num DOFs: 1323    (441 nodes * 3 variables: T, disp_x, disp_y)

 0 Nonlinear |R| = ...
 1 Nonlinear |R| = ~1e-12

Solve Converged!
```

Newton should converge in 1-2 iterations because the physics is linear.

Output files written to the same directory:
- `case14_thermoelasticity_out.e` — Exodus file with all field variables
- `case14_thermoelasticity_out.csv` — Postprocessor scalar values

---

## Expected Results

### Temperature Field

The steady temperature follows the exact linear ramp (insulated top/bottom enforce
zero y-gradient):

```
T(x) = 500 - 200*x    K
```

- Left edge: 500 K
- Centre (x=0.5): 400 K
- Right edge: 300 K

Postprocessor `max_temperature` should read 500 K; `avg_temperature` should read
approximately 400 K.

### Displacement Field

The dominant displacement is horizontal (`disp_x`), driven directly by thermal
expansion. The hot left side expands more than the cold right side:

- The plate expands predominantly to the right (positive `disp_x` increases with x
  from the pinned bottom).
- The thermal gradient also generates a bending moment that lifts the plate upward
  (`disp_y > 0` away from the pinned bottom).

Approximate magnitude for `max_disp_x`:

```
delta ~ alpha * delta_T * L = 12e-6 /K * 200 K * 1 m ~ 2.4 mm
```

The actual value will be somewhat smaller because boundary constraints partially
restrain the expansion, but the order of magnitude is correct.

### Stress Field

The von Mises stress is highest near the pinned bottom edge and near the hot left
corner, where geometric constraints most strongly resist the thermal expansion.

Order-of-magnitude estimate for peak stress:

```
sigma ~ E * alpha * delta_T ~ 200e9 * 12e-6 * 200 ~ 480 MPa
```

This is near the yield stress of steel (~250-500 MPa depending on grade), which is
physically consistent — large thermal gradients in restrained steel structures do
generate near-yield stresses. In a real design analysis you would compare this to
the allowable stress for the grade and add safety factors.

---

## Key Takeaways

- **One-way multiphysics coupling** in MOOSE: solve temperature first (or
  simultaneously), then use it as input to mechanics via an eigenstrain.
- **`Physics/SolidMechanics/QuasiStatic` action**: one block replaces many manual
  kernel/variable/material declarations; always use it for solid mechanics problems.
- **Eigenstrains** are the correct MOOSE mechanism for any "stress-free" strain:
  thermal expansion, swelling, creep pre-strain, phase transformation volume change.
- **`combined-opt`** is required whenever physics modules beyond the base framework
  are needed. Use it in place of `moose_test-opt` for all module-based problems.
- The linear temperature ramp and approximate `delta ~ alpha * delta_T * L` formula
  provide an instant sanity check on the numerical result before opening ParaView.
