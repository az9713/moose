# Case 05: Spatially Varying Conductivity

## Overview

Cases 01 through 04 all assumed **uniform material properties** — the diffusion
coefficient was the same everywhere in the domain. Real materials are rarely uniform.
The thermal conductivity of a composite material varies by position. The permeability
of a geological formation changes from point to point. The diffusivity in a biological
tissue depends on local structure.

This case introduces **spatially varying coefficients**: the conductivity is no longer a
constant but a function of position, `k(x) = 1 + x`. This single change makes the
problem physically richer and requires a fundamentally different approach in MOOSE:
using the **Material system** and **GenericFunctionMaterial** to provide position-
dependent properties to kernels.

The problem solved is:

```
-div( k(x) * grad(u) ) = 0,    u=0 at x=0,  u=1 at x=1
```

with `k(x) = 1 + x`. This is a steady-state heat conduction problem where the material
gets more conductive from left to right. Heat flows from the hot right boundary (u=1)
to the cold left boundary (u=0), but the varying conductivity distorts the temperature
profile away from a simple straight line.

This case teaches how MOOSE decouples **material property definitions** from **PDE
operators**, and how both can be independently varied to model different physical scenarios.

---

## The Physics

### The Physical Problem

Imagine a slab of material — like a layered composite board or a piece of graded alloy —
that sits between two temperature-controlled walls. The left wall is cold (u = 0, think
of this as 0 degrees), the right wall is hot (u = 1, think of this as 100 degrees). The
top and bottom are insulated (no heat flux through them).

In a uniform material, heat flows in straight horizontal lines and the temperature varies
linearly from left to right: `u = x`. But this material is not uniform. The conductivity
increases from left to right as `k(x) = 1 + x`. Near the left wall, k = 1. Near the
right wall, k = 2.

A more conductive region carries heat more easily. Intuitively, the high-conductivity
right half "prefers" to be at a uniform temperature — it conducts away any gradient
rapidly. The low-conductivity left half must carry the same total heat flux but with a
smaller k, so it needs a steeper gradient. This means the temperature profile curves:
it rises steeply on the left and flattens out on the right.

### The Governing Equation

```
-div( k(x) * grad(u) ) = 0    on [0,1] x [0,1]
```

Expanding in 2D:

```
  d           du         d           du
- -- [ k(x) * -- ]  -   -- [ k(x) * -- ]  =  0
  dx          dx         dy          dy
```

Since `k(x) = 1 + x` depends only on x (not y), and there is no y-direction forcing,
the solution does not vary in y. The effective 1D equation governs the behavior:

```
  d           du
- -- [ (1+x) * -- ]  =  0
  dx           dx
```

This is the **variable-coefficient diffusion equation** (or variable-coefficient heat
equation in steady state).

Symbol definitions:
- `u(x, y)` — the scalar field (temperature, concentration, potential)
- `k(x)` — the **conductivity** (or diffusivity): how easily the quantity moves through
  the material at location x. Units: [quantity * length / (flux * time)]
- `grad(u)` — spatial gradient of u: vector pointing in the direction of greatest increase
- `k(x) * grad(u)` — the **flux**: how much of `u` passes through a unit area per unit time
- `div(...)` — divergence of the flux: net flux leaving a small volume
- `-div(k * grad u) = 0` — divergence of flux is zero, meaning what flows in equals
  what flows out: **conservation** of the transported quantity

Setting the divergence of flux to zero means the system is in steady state with no sources.

### Boundary Conditions

```
u = 0    on x=0 (left boundary, cold wall)
u = 1    on x=1 (right boundary, hot wall)
```

No boundary condition is applied on the top and bottom. This is a **Neumann zero** (or
**natural**) boundary condition — it arises automatically from the weak form and means
zero flux through those boundaries. Heat only enters and exits through the left and right walls.

### Analytical Solution

For the 1D version (ignoring y), the steady-state equation with `k(x) = 1+x` has an
analytical solution. Setting the heat flux `q = -(1+x) * du/dx` to a constant C (the
same flux passes through every cross-section at steady state):

```
-(1+x) * du/dx = C
du/dx = -C / (1+x)
u(x) = -C * ln(1+x) + D
```

Applying boundary conditions `u(0) = 0` and `u(1) = 1`:

```
u(0) = 0 = D,               so D = 0
u(1) = 1 = -C * ln(2),      so C = -1/ln(2)

u(x) = ln(1+x) / ln(2)
```

This is a **logarithmic profile**: it rises steeply near x=0 (low conductivity) and
flattens near x=1 (high conductivity). At x=0.5, the exact temperature is:

```
u(0.5) = ln(1.5)/ln(2) ≈ 0.405/0.693 ≈ 0.585
```

Compare this with the linear profile (uniform k): `u(0.5) = 0.5`. The bias toward
values above 0.5 confirms the physical intuition: high conductivity on the right
compresses the gradient there, pushing the midpoint temperature upward.

### What Spatially Varying Conductivity Means Physically

Think of it this way: the heat flux (amount of heat flowing per unit area per unit time)
must be continuous everywhere in a steady state. At every x-plane, the same amount of
heat passes through. This is conservation of energy.

But the flux is `q = -k * du/dx`. If k is large (right side), du/dx must be small to
carry the same q. If k is small (left side), du/dx must be large. The temperature
gradient is inversely proportional to conductivity.

Materials with spatially varying conductivity appear in:
- **Functionally graded materials (FGM)**: deliberately engineered composition gradients
- **Geological formations**: permeability varies with rock type
- **Biological tissue**: diffusivity varies with cell density
- **Temperature-dependent properties**: if k depends on u, the problem becomes nonlinear

### ASCII Domain Diagram

```
  y=1  (insulated: zero flux)
       +----------------------------------------+
       |  steep gradient  |  shallow gradient   |
       |  (k=1, low k)    |  (k~2, high k)      |
  u=0  |                  |                     |  u=1
 (cold)|  u rises fast    |  u rises slowly     | (hot)
       |  from 0          |  toward 1           |
       |                  |                     |
       +----------------------------------------+
  y=0  (insulated: zero flux)
       x=0              x=0.5                 x=1
       k=1                                   k=2

  Temperature profile along x (1D cross-section):
  1.0 |                                    *
  0.8 |                              *
  0.6 |                        *
  0.5 |. . . . . . . .  . (linear profile would be here)
  0.4 |              *   <- actual profile is above linear
  0.2 |       *
  0.0 +-----------------------------------
      x=0                              x=1

  Exact: u(x) = ln(1+x) / ln(2)
  Profile is concave (curves upward), above the linear solution
```

---

## Input File Walkthrough

The input file is `case05_varying_k.i`. Let us go through every block and parameter.

### Block: `[Mesh]`

```
[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 30
  ny   = 30
[]
```

- `type = GeneratedMesh` — built-in structured rectangular mesh, domain [0,1] x [0,1].
- `dim = 2` — two-dimensional mesh, even though the physics is effectively 1D.
  The 2D mesh is used to demonstrate that material properties work in 2D just as
  well as 1D.
- `nx = 30` — 30 elements in x. Element width h_x = 1/30 ≈ 0.033.
- `ny = 30` — 30 elements in y. The y-direction is included for completeness; the
  solution will be constant in y (no y-variation in either BCs or conductivity).

Why 30x30 instead of 20x20 from case 04? The finer mesh provides a smoother
visualization and better resolution of the curved logarithmic profile.

Total: 900 elements, 961 nodes.

### Block: `[Variables]`

```
[Variables]
  [u]
  []
[]
```

- Single scalar unknown `u`, first-order Lagrange (default).
- 961 degrees of freedom.
- Initial condition: zero everywhere (default).

### Block: `[Functions]`

```
[Functions]
  [k_fn]
    type       = ParsedFunction
    expression = '1 + x'
  []
[]
```

- `[k_fn]` — defines the conductivity as a function of position.
- `type = ParsedFunction` — evaluates the string expression at any (x, y, z, t).
- `expression = '1 + x'` — the conductivity field. At x=0: k=1. At x=1: k=2.
  The conductivity increases linearly from left to right.

This is a **Function object** — not a material property yet. The Function just knows
how to evaluate `1 + x` at any point. The Material block below connects this Function
to a named material property that the Kernel can use.

The distinction matters: Functions and Materials serve different roles in MOOSE's
architecture. Functions are pure mathematical evaluators (any (x,y,z,t) to a value).
Materials are physics-layer properties that are evaluated on elements at quadrature
points and cached for use by Kernels. GenericFunctionMaterial is the bridge between
the two.

### Block: `[Materials]`

```
[Materials]
  [conductivity]
    type        = GenericFunctionMaterial
    prop_names  = 'k'
    prop_values = 'k_fn'
  []
[]
```

This is the key new concept of Case 05.

- `[conductivity]` — the user-chosen name for this material object (any name works).
- `type = GenericFunctionMaterial` — a built-in MOOSE material class that creates
  material properties by evaluating Function objects at each quadrature point.
- `prop_names = 'k'` — declares a material property named `'k'`. This is the name
  that kernels, BCs, and other objects will use to request the conductivity.
- `prop_values = 'k_fn'` — specifies that the value of property `'k'` comes from
  the Function object named `k_fn`.

**How it works internally:** During the finite element assembly loop, MOOSE visits each
element. For each element, it computes the values at Gauss quadrature points. When
`GenericFunctionMaterial` is asked to compute the property `k`, it calls:

```
k_value_at_qp = k_fn.value(qp_x, qp_y, qp_z, t)
```

where `(qp_x, qp_y)` is the physical coordinate of the quadrature point. This value
is stored in the material property system and made available to kernels.

If you had a constant conductivity instead, you would use:

```
[Materials]
  [conductivity]
    type        = GenericConstantMaterial
    prop_names  = 'k'
    prop_values = '2.5'
  []
[]
```

`GenericFunctionMaterial` is the generalization of this for spatially (and temporally)
varying properties.

**Why use the Material system at all?** Why not just evaluate `1 + x` directly in the
kernel? Two reasons:

1. **Separation of concerns**: The physics (diffusion) and the material (conductivity)
   are independent. You can swap materials without touching kernels.
2. **Multiple clients**: If multiple kernels need the same property (e.g., both a heat
   equation kernel and a stress kernel need Young's modulus), the Material evaluates it
   once and all kernels read the cached value.

### Block: `[Kernels]`

```
[Kernels]
  [diffusion]
    type        = MatDiffusion
    variable    = u
    diffusivity = k
  []
[]
```

- `type = MatDiffusion` — implements the variable-coefficient diffusion term:
  ```
  integral( k(x) * grad(phi_i) . grad(u) dV )
  ```
  This is the weak form of `-div(k * grad u)`.
- `variable = u` — the unknown field.
- `diffusivity = k` — the **name of the material property** to use as the coefficient.
  This must match exactly the `prop_names` string declared in the `[Materials]` block.

Contrast with Case 04's kernel:

```
[diffusion]
  type     = Diffusion      # constant coefficient (k=1 implicitly)
  variable = u
[]
```

vs Case 05's kernel:

```
[diffusion]
  type        = MatDiffusion   # variable coefficient from material system
  variable    = u
  diffusivity = k              # reads material property 'k' at each quad point
[]
```

`MatDiffusion` internally calls `getMaterialProperty<Real>("k")` to obtain a reference
to the material property array. At each quadrature point during assembly, it reads the
spatially varying value of `k` from that array.

The weak form integral for MatDiffusion is:

```
R_i = integral_over_element ( k(x) * grad(phi_i) . grad(u_h) ) dV
```

where `k(x)` is evaluated at each quadrature point separately. This is critical: the
coefficient must vary *within* an element, not just between elements. Because fparser
expressions can be evaluated at any continuous point, MOOSE captures sub-element
variation correctly.

### Block: `[BCs]`

```
[BCs]
  [left]
    type     = DirichletBC
    variable = u
    boundary = left
    value    = 0
  []
  [right]
    type     = DirichletBC
    variable = u
    boundary = right
    value    = 1
  []
[]
```

- Two constant Dirichlet BCs: zero on the left wall, one on the right wall.
- No BCs on the top and bottom: these are natural (zero-flux, insulated) boundaries.
- `type = DirichletBC` — the simplest BC: pins the solution to a fixed value at all
  boundary nodes on the specified boundary.
- `boundary = left` / `boundary = right` — named boundaries assigned automatically by
  `GeneratedMesh`. The left boundary is all nodes at x=0; the right is all at x=1.

### Block: `[Postprocessors]`

```
[Postprocessors]
  [avg_u]
    type     = ElementAverageValue
    variable = u
  []
[]
```

- `type = ElementAverageValue` — computes the volume-weighted average of `u` over the
  entire domain.
- `variable = u` — the field to average.

For the uniform-k linear case, the average would be exactly 0.5. For `k = 1+x`, the
profile is biased upward (temperature is higher than linear in the interior), so the
average should be above 0.5. From the exact solution:

```
avg_u = integral from 0 to 1 of ln(1+x)/ln(2) dx
      = [  (1+x)*ln(1+x) - (1+x) ]_0^1 / ln(2)
      = [ (2*ln(2) - 2) - (0 - 1) ] / ln(2)
      = [ 2*ln(2) - 1 ] / ln(2)
      = 2 - 1/ln(2)
      ≈ 2 - 1.4427
      ≈ 0.5573
```

This sanity check — average above 0.5 as expected — is a quick indicator that the
spatially varying conductivity is having the correct physical effect.

### Block: `[Executioner]`

```
[Executioner]
  type       = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]
```

- Same steady-state PJFNK solver as Case 04.
- For the linear problem `MatDiffusion` with a *position-dependent but u-independent*
  coefficient, the problem is still linear in u. Newton converges in one iteration.
- BoomerAMG is highly effective for this type of problem.

### Block: `[Outputs]`

```
[Outputs]
  exodus = true
  csv    = true
[]
```

- `exodus = true` — writes the solution `u` to `case05_varying_k_out.e`.
- `csv = true` — writes the `avg_u` postprocessor to `case05_varying_k_out.csv`.

---

## What Happens When You Run This

### Invocation

```bash
cd quickstart-runs/case05
mpirun -n 1 moose-app-opt -i case05_varying_k.i
```

### Material Property Evaluation

The key difference from earlier cases is in the assembly phase. When MOOSE visits an
element and calls `GenericFunctionMaterial::computeQpProperties()`:

1. For each quadrature point `qp` with physical coordinate `(x_qp, y_qp)`:
2. Evaluates `k_fn.value(x_qp, y_qp, 0, 0) = 1 + x_qp`
3. Stores this in `_prop_k[qp]` (the material property array for this element)

Then when `MatDiffusion::computeQpResidual()` is called for the same quadrature point:

1. Reads `_diffusivity[qp]` which is the same `1 + x_qp` value
2. Uses it to weight the gradient term

This per-quadrature-point evaluation is what makes the property truly spatially varying.

### Convergence

The problem is linear (k does not depend on u), so Newton converges in one step:

```
 0 Nonlinear |R| = 1.000000e+00
      0 Linear |R| = 8.23e-15
 1 Nonlinear |R| = 8.23e-15

Converged!

Postprocessor Values:
  avg_u = 5.5730e-01
```

The `avg_u` value of approximately 0.557 confirms the logarithmic profile is being
computed correctly (matches the analytical average computed above).

---

## Output Files

### `case05_varying_k_out.e`

The Exodus file contains the solution field `u` at all 961 nodes. In ParaView:
1. Open the file, click Apply
2. Color by `u`
3. Select a vertical cut plane at y=0.5 (or use a "Plot Over Line" filter along y=0.5)
4. The 1D profile should match `ln(1+x)/ln(2)` closely

Characteristics to observe:
- The color gradient is NOT uniform left-to-right (unlike case 02)
- The midpoint x=0.5 has temperature near 0.585, not 0.5
- The solution is constant in the y-direction (all horizontal lines have identical color)
- The gradient (color change per unit length) is steeper near x=0 than near x=1

### `case05_varying_k_out.csv`

```
time,avg_u
0,5.5730e-01
```

The average of approximately 0.557 is above 0.5, consistent with the logarithmic exact
solution. For a constant k, the solution would be linear and the average would be exactly 0.5.

---

## Understanding the Plots

Running `python visualize_all.py` from the `quickstart-runs/` directory produces the
following two plots saved into this directory.

### `case05_contour_2d.png`

**What the plot shows.** A 2D filled-contour of the solution `u(x,y)` on the unit
square with conductivity `k(x) = 1 + x`. The viridis colormap is used with contour
bands filled by u value.

**Physical quantities.** The color encodes the scalar field u, which represents
temperature in a material whose conductivity increases from left (k=1) to right (k=2).

**How to judge correctness.** The contour bands should still be roughly vertical
(since BCs are applied at the left and right walls and there is no y-variation in
the problem), but they should be more closely spaced on the left (low-k region needs
a steeper gradient to carry the same flux) and more widely spaced on the right
(high-k region carries the flux easily with a shallow gradient). This gives the
characteristic "crowded on the left" appearance compared to Case 02's uniform spacing.

**What would indicate a problem.**
- Uniformly spaced vertical bands (as in Case 02): the varying conductivity k(x) was
  not applied — the material is using constant k=1 instead.
- Curved or tilted bands: the top/bottom boundary conditions are incorrect.
- Bands more closely spaced on the right (opposite of expected): the conductivity
  function has the wrong sign or direction.

### `case05_line_exact.png`

**What the plot shows.** A line plot comparing two curves along the horizontal
midline y ≈ 0.5: blue circles connected by lines for the MOOSE solution extracted at
nodes near y=0.5, and a dashed red line for the exact analytical solution
`u = ln(1 + x) / ln(2)`.

**Physical quantities.** The horizontal axis is x (0 to 1). The vertical axis is u
(0 to 1). The exact solution is a logarithmic function that is concave-up: it rises
steeply near x=0 (low k, steep gradient) and flattens near x=1 (high k, shallow
gradient), in contrast to the straight-line solution of Case 02.

**How to judge correctness.** The MOOSE blue dots should closely follow the red dashed
curve. The curve should be visibly concave-up (bowing below the straight line y=x).
The endpoints should be exactly 0 at x=0 and 1 at x=1. For a 20-element mesh the
MOOSE dots may deviate slightly from the exact curve, but the shape should match.

**What would indicate a problem.**
- Blue dots lying on the straight line y=x instead of the logarithmic curve: constant
  conductivity was used instead of k(x) = 1+x.
- Blue dots with a concave-down shape: conductivity is decreasing instead of increasing.
- Large gaps between blue dots and the red curve beyond what mesh discretization
  would explain: the exact solution formula or material property is incorrect.

---

## Interpreting the Results

### The Solution Profile

The temperature field `u(x)` follows a logarithmic curve:

```
u(x) = ln(1+x) / ln(2)
```

Key values:
- u(0.0) = 0.000  (left BC, exactly enforced)
- u(0.5) = 0.585  (above the linear midpoint of 0.500)
- u(1.0) = 1.000  (right BC, exactly enforced)

The profile is concave upward (like a logarithm) meaning it rises steeply at first
(where k is small) and flattens out (where k is large).

### Verifying Correctness

1. **Average value**: should be ~0.557 — confirmed by postprocessor.
2. **Heat flux**: the product `k(x) * du/dx` should be constant at steady state.
   - At x=0: k=1, du/dx = 1/ln(2) ≈ 1.443, flux = 1.443
   - At x=1: k=2, du/dx = 0.5/ln(2) ≈ 0.721, flux = 2 * 0.721 = 1.443  (same)
   - This confirms energy conservation (same flux through every cross-section).
3. **Compare to case 02** (constant k): the profile should visibly curve relative to
   the straight-line temperature distribution of the uniform case.

### Physical Insight

The key physical lesson: **heat flux is conserved, but temperature gradient is not**.
In steady state, the same amount of heat energy passes through every vertical cross-
section. But how steep the temperature gradient must be to carry that heat depends on
the local conductivity. Low-k regions need a steeper gradient; high-k regions can carry
the heat with a shallower gradient.

This principle underlies the design of thermal barrier coatings, the analysis of
geothermal gradients in layered rock, and the modeling of heterogeneous materials.

---

## Key Concepts Learned

- **Variable-coefficient PDEs**: the diffusion coefficient is a function of position
  rather than a constant, leading to curved (non-linear) solution profiles
- **ParsedFunction**: evaluate any mathematical expression at spatial points — used here
  to define the conductivity field k(x) = 1+x
- **GenericFunctionMaterial**: the MOOSE bridge from Function objects to material
  properties; evaluates a Function at each quadrature point and stores the result as
  a named material property
- **Material property system**: separates material definitions from PDE operators,
  enabling clean model-material separation
- **MatDiffusion kernel**: the variable-coefficient version of the Diffusion kernel,
  reads its coefficient from the material property system by name
- **Natural (Neumann zero) boundary conditions**: top and bottom boundaries have no
  explicit BC, which automatically means zero flux through those boundaries
- **ElementAverageValue postprocessor**: computes a domain-averaged quantity for a
  quick sanity check against an analytical estimate
- **Logarithmic temperature profile**: the analytical solution for linear k(x) = 1+x,
  and how to verify the simulation matches it
- **Heat flux conservation**: at steady state, the same flux passes through every
  cross-section regardless of local conductivity

---

## Experiments to Try

### Experiment 1: Reverse the Conductivity Gradient

Change `k_fn` expression to `'2 - x'` so that k=2 on the left and k=1 on the right.
Now the high-conductivity side is on the left (cold side). Predict: the profile will
curve the other way — below the linear solution, with avg_u < 0.5.

The analytical solution becomes `u(x) = (1 - ln(2-x)/ln(2))`, which is concave downward.
Expected avg_u ≈ 0.443.

### Experiment 2: Stronger Variation

Try `k_fn` expression `'0.1 + x'` (k ranges from 0.1 to 1.1, a 10x variation instead
of 2x). The logarithmic profile becomes much more pronounced. The midpoint temperature
rises to approximately:

```
u(0.5) = ln(1.1/0.1) / ln(1.1/0.1 + 1) ... (rederive using the integral formula)
```

The solution becomes much more curved — nearly all the temperature drop is concentrated
near the low-k left boundary. This simulates a thermal boundary layer.

### Experiment 3: Two-Direction Variation

Try `k_fn` expression `'(1 + x) * (1 + y)'`. Now k varies in both x and y, ranging
from 1 at (0,0) to 4 at (1,1). The solution is no longer constant in y. This creates
a genuinely 2D problem with no simple analytical solution. Use ParaView to explore
the 2D temperature field.

### Experiment 4: Add a Source Term

Add a `BodyForce` kernel with a constant function `f = 1.0`:

```
[Kernels]
  [source]
    type     = BodyForce
    variable = u
    value    = 1.0
  []
[]
```

This represents uniform heat generation inside the material (like electrical resistance
heating). The solution is no longer zero-divergence — it must now accommodate both the
source and the boundary conditions. The effect of varying k on the source distribution
is more subtle and requires numerical exploration.

### Experiment 5: Check Mesh Convergence

Run with nx=ny=15, 30, 60. Extract the midpoint temperature (use an additional
postprocessor of type `PointValue` at point '0.5 0.5 0') and compare to the exact value
of `ln(1.5)/ln(2) ≈ 0.5850`. Confirm second-order convergence.
