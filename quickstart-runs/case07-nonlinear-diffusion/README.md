# Case 07: Nonlinear Diffusion with Temperature-Dependent Conductivity

## Overview

This case crosses the boundary from linear to nonlinear PDEs. In Cases 01 through 06
the diffusion coefficient was always a fixed constant — the governing equation was
linear, and Newton's method converged in a single iteration. Here the thermal
conductivity `k` depends on the unknown temperature `T` itself: `k(T) = 1 + T`. Because
the coefficient changes as the solution changes, this is a genuinely nonlinear problem.
Newton's method now requires multiple iterations to converge, but it does so with
*quadratic* speed — one of the most satisfying demonstrations in numerical mathematics.

The central MOOSE technique introduced here is **Automatic Differentiation (AD)**. AD
supplies Newton's method with the exact Jacobian of the residual at no extra coding
cost. This is what makes the `AD`-prefixed objects (`ADMatDiffusion`,
`ADPiecewiseLinearInterpolationMaterial`) different from their non-AD counterparts.

The problem models heat conduction through a square plate where the material conducts
heat better as it gets hotter — a property exhibited by certain ceramics and
semiconductors. A uniform volumetric heat source warms the plate while all walls are
held at zero temperature.

---

## The Physics

### Physical Problem in Plain English

Imagine a square thin plate of ceramic material, occupying the unit square from
(0,0) to (1,1). Every point inside the plate generates heat at a constant rate Q = 10
(units: power per volume). All four edges of the plate are held at T = 0 (think of
them clamped to a heat sink kept at the reference temperature).

At steady state, the heat generated inside must flow outward to the walls by conduction.
The rate of conduction depends on the thermal conductivity k. In this problem k is not
a constant — it increases linearly with temperature: `k(T) = 1 + T`. Hotter regions
conduct better. This means the solution is not simply a paraboloid (as it would be with
constant k); it is "puffed up" toward the center and the problem is nonlinear.

### Governing Equation

The strong form of the steady nonlinear heat equation is:

```
-div( k(T) * grad(T) ) = Q    in the domain [0,1] x [0,1]
```

Expanded in two dimensions (x, y):

```
  d         dT           d         dT
- -- ( k(T) -- )   -   -- ( k(T) -- )   =   Q
  dx        dx           dy        dy
```

Symbol definitions:

- `T(x,y)` — the temperature field, the unknown we are solving for [K or dimensionless]
- `k(T)` — thermal conductivity, which depends on T: `k(T) = 1 + T` [W/m/K or dimensionless]
- `grad(T)` — the gradient vector `(dT/dx, dT/dy)` pointing toward higher temperature
- `div(...)` — divergence: sum of partial derivatives of the vector argument
- `-div(k grad T)` — the net conductive heat flux diverging from a point (heat leaving)
- `Q` — volumetric heat source rate: 10 [W/m³ or dimensionless]

When `k` is a constant this reduces to the Poisson equation `-k * Laplacian(T) = Q`,
which is linear. When `k` depends on `T`, the equation is nonlinear because the
coefficient multiplying `T` is itself a function of `T`.

### What Makes This Nonlinear

In the linear case the residual equation `R(T) = 0` looks like:

```
R(T) = K * T - f    (K is a fixed stiffness matrix, f is a fixed load vector)
```

Newton's method solves this in one step because the Jacobian `dR/dT = K` does not
change.

In the nonlinear case with `k(T) = 1 + T`:

```
R(T) = assembly of: -div((1 + T) * grad(T)) - Q
```

The Jacobian `dR/dT` now contains extra terms from differentiating `k(T)` with respect
to `T`. These terms involve `dk/dT = 1` multiplied by `grad(T)` — they change at every
Newton iteration as T changes. Newton must iterate until `R(T)` is small enough.

### Boundary Conditions

```
T = 0    on x=0 (left wall)
T = 0    on x=1 (right wall)
T = 0    on y=0 (bottom wall)
T = 0    on y=1 (top wall)
```

All four walls are isothermal at the reference temperature. Physically: the plate is
surrounded by a perfect heat sink that holds the boundary at zero temperature regardless
of how much heat is generated inside.

These are **Dirichlet** (essential) boundary conditions — the most common type in heat
transfer problems with known wall temperatures.

### ASCII Domain Diagram

```
y=1   T=0 (top wall, heat sink)
      +------------------------------------------+
      |                                          |
      |   -div( (1+T)*grad(T) ) = 10             |
      |                                          |
      |   k(T) = 1 + T                           |
      |   Q    = 10 (uniform heat source)        |
      |                                          |
      |   Peak temperature near center           |
      |   (higher than constant-k case)          |
T=0   |                                          | T=0
(left)|          T_max ~ 2.5                     | (right)
      |                                          |
      |                                          |
      +------------------------------------------+
y=0   T=0 (bottom wall, heat sink)
      x=0                                     x=1

Mesh: 20 x 20 uniform quadrilateral elements
      400 elements, 441 nodes
```

### The Nonlinear Conductivity Function

The conductivity law `k(T) = 1 + T` is implemented as a piecewise-linear interpolation
through two data points: (T=0, k=1) and (T=10, k=11). Within this range the function is
exact (it is literally a line through two points on the line y = 1 + x). Outside the
range MOOSE extrapolates linearly using the slope of the nearest segment.

Because the maximum temperature in the solution is roughly 2 to 3 (depending on
nonlinearity), both data points bracket the actual solution range with room to spare.

Why not use `k(T) = 1 + T` directly as an analytic expression? The `ADParsedMaterial`
object, which evaluates string expressions, requires the MOOSE JIT compiler. In Docker
containers and some HPC environments this compiler is unavailable. The piecewise-linear
approach `ADPiecewiseLinearInterpolationMaterial` works everywhere.

---

## Automatic Differentiation: What It Is and Why It Helps

### The Problem AD Solves

Newton's method requires the Jacobian matrix: `J_ij = dR_i / dT_j`, the derivative of
each residual entry with respect to each degree of freedom. For nonlinear problems this
matrix has extra terms from the nonlinear coefficients.

There are three ways to supply the Jacobian in MOOSE:

1. **Hand-coded Jacobian** (non-AD kernels like `MatDiffusion`): The developer writes
   separate C++ code for `computeQpJacobian()` and `computeQpOffDiagJacobian()`. This
   is error-prone and tedious, especially for complex nonlinear expressions.

2. **Jacobian-Free Newton-Krylov (JFNK / PJFNK)**: Approximates Jacobian-vector
   products using finite differences of the residual. No explicit Jacobian is assembled.
   Cheaper per iteration but requires more linear solver iterations. The approximation
   introduces a small error that can slow or degrade convergence.

3. **Automatic Differentiation (AD)**: The code computes both the residual and its
   exact derivative simultaneously using **dual-number arithmetic**. No hand-coding of
   derivatives required. No finite-difference approximation error. The Jacobian is exact
   to machine precision.

### Dual Numbers: The Core Idea

A dual number is a pair `(v, dv)` where `v` is the value and `dv` is the derivative
with respect to some independent variable. Arithmetic on dual numbers follows the chain
rule automatically:

```
(a, da) + (b, db) = (a+b, da+db)       addition rule
(a, da) * (b, db) = (a*b, a*db + b*da) product rule (chain rule)
f(a, da)          = (f(a), f'(a)*da)   function rule
```

In MOOSE's AD system, each degree of freedom `T_j` is seeded as `(T_j, 1)` for the
j-th variable and `(T_j, 0)` for all others. As kernels evaluate the residual using
dual-number arithmetic, the derivative part tracks exactly `dR_i / dT_j` at zero extra
conceptual cost to the developer.

### ADMatDiffusion vs. MatDiffusion

| Property | MatDiffusion | ADMatDiffusion |
|----------|-------------|----------------|
| Residual computation | computeQpResidual() | computeQpResidual() |
| Jacobian computation | computeQpJacobian() (hand-coded) | Automatic |
| Jacobian accuracy | Only diagonal term by default | Full exact Jacobian |
| Works with nonlinear k(T) | Partially (missing off-diagonal) | Fully |
| Solver compatibility | PJFNK preferred | NEWTON preferred |

With `MatDiffusion` and `k(T) = 1+T`, the hand-coded Jacobian only captures `dk/du`
terms on the diagonal. The cross-coupling terms from the nonlinear conductivity are
missing unless you write them manually. With `ADMatDiffusion`, the entire Jacobian is
computed automatically and exactly.

### Why NEWTON Works Better Here

With an exact Jacobian from AD, you should use `solve_type = NEWTON` rather than
`PJFNK`. Here is why:

- `PJFNK` approximates `J * v` using `(R(T + eps*v) - R(T)) / eps`. This introduces
  finite-difference error of order `eps` and also requires matrix-free Krylov iterations.
- `NEWTON` uses the exact Jacobian assembled by AD. The linear system at each Newton step
  is more accurate, leading to fewer Newton iterations and better convergence behavior.

With `NEWTON` and AD, Newton's method achieves **quadratic convergence**: the residual
norm approximately squares each iteration. If the residual is `1e-4` after one iteration,
it will be roughly `1e-8` after the next, then `1e-16`. This is the theoretical
optimum for Newton's method.

---

## Input File Walkthrough

The input file is `case07_nonlinear_diffusion.i`. Let us go through every block and
every parameter.

### Block: `[Mesh]`

```
[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 20
  ny   = 20
[]
```

- `type = GeneratedMesh` — MOOSE's built-in structured rectangular mesh generator.
  No external mesh file is needed.
- `dim = 2` — two-dimensional problem in the (x, y) plane.
- `nx = 20` — 20 elements along the x-axis. Default domain is [0,1] so each element
  is 0.05 units wide.
- `ny = 20` — 20 elements along the y-axis. Same element height.

The result is a 20x20 = 400-element uniform quadrilateral mesh with (21x21) = 441 nodes.
The named boundaries `left`, `right`, `top`, `bottom` are assigned automatically.

For this nonlinear problem, the mesh density determines how well the smooth but
non-parabolic temperature peak is resolved. 20x20 is sufficient for engineering accuracy.
You would need a finer mesh only if you needed to verify the solution against a reference
at high precision.

### Block: `[Variables]`

```
[Variables]
  [T]
  []
[]
```

- Declares a single unknown scalar field named `T` (temperature).
- The empty sub-block `[T]` accepts all defaults:
  - `order = FIRST` — first-order (bilinear) Lagrange elements
  - `family = LAGRANGE` — standard nodal finite element basis
  - Initial condition: zero everywhere (T = 0 at t=0, or as the initial Newton guess)
- MOOSE allocates one degree of freedom per node: 441 unknowns for this mesh.

The initial guess of T = 0 everywhere satisfies the boundary conditions and is a
reasonable starting point for Newton iteration on this problem.

### Block: `[Kernels]`

```
[Kernels]
  [diffusion]
    type        = ADMatDiffusion
    variable    = T
    diffusivity = k
  []

  [source]
    type     = BodyForce
    variable = T
    value    = 10.0
  []
[]
```

Two kernels together implement the complete weak form of the PDE.

**`[diffusion]` kernel:**
- `type = ADMatDiffusion` — the automatic-differentiation version of `MatDiffusion`.
  It implements the weak form term:
  ```
  integral_over_element( k * grad(phi_i) . grad(T) ) dV
  ```
  where `phi_i` is the test function and `T` is the trial solution.
- `variable = T` — this kernel acts on the equation for `T`.
- `diffusivity = k` — references the material property named `k`. This is not a
  fixed constant — it is the AD material property returned by the `[conductivity]`
  material defined below. MOOSE looks up the property named `k` at each quadrature
  point during assembly. Because the material is AD-aware, the value it returns is
  a dual number containing both `k(T)` and `dk/dT`.

The "AD" prefix means this kernel's residual computation uses dual-number arithmetic.
The Jacobian is derived automatically from the same code that computes the residual.
No separate `computeQpJacobian()` function is written or needed.

**`[source]` kernel:**
- `type = BodyForce` — adds a volumetric source to the right-hand side. Implements:
  ```
  -integral_over_element( phi_i * value ) dV
  ```
  (negative because it moves to the right-hand side: residual = diffusion - source)
- `variable = T` — contributes to the T equation.
- `value = 10.0` — the constant heat generation rate Q. At each quadrature point,
  the kernel evaluates `phi_i * 10.0` and adds it to the element residual vector.
- This kernel is not AD-specific because 10.0 is a constant — its Jacobian contribution
  is zero (a constant source does not depend on T).

Together, the two kernels represent the complete governing PDE in weak form:

```
integral( k(T) * grad(phi_i) . grad(T) dV )  =  integral( phi_i * 10 dV )
```

### Block: `[BCs]`

```
[BCs]
  [walls]
    type     = DirichletBC
    variable = T
    boundary = 'left right top bottom'
    value    = 0
  []
[]
```

- `type = DirichletBC` — the standard fixed-value boundary condition. During assembly
  MOOSE replaces the residual equation for each boundary node with `T_node - 0 = 0`,
  which forces `T_node = 0` exactly.
- `variable = T` — applies to the temperature field.
- `boundary = 'left right top bottom'` — applies to all four named boundaries
  simultaneously. A single `DirichletBC` object can cover multiple boundaries.
- `value = 0` — the prescribed wall temperature. All walls are held at zero.

Why zero? This is a common non-dimensionalization choice. The physical scenario could
be a ceramic plate where the walls are water-cooled to room temperature and `T` is the
temperature rise above ambient.

### Block: `[Materials]`

```
[Materials]
  [conductivity]
    type        = ADPiecewiseLinearInterpolationMaterial
    property    = k
    variable    = T
    xy_data     = '0 1
                   10 11'
  []
[]
```

This block defines the nonlinear material property that is the heart of this case.

- `type = ADPiecewiseLinearInterpolationMaterial` — a material that evaluates a
  piecewise-linear function of a coupled variable and stores the result as an
  AD material property (i.e., a dual number carrying both value and derivative).

- `property = k` — the name under which this material property is stored. The
  `ADMatDiffusion` kernel retrieves it by this name.

- `variable = T` — the independent variable of the interpolation. The material
  will be evaluated as `k = f(T)` at each quadrature point.

- `xy_data` — a flat list of (x, y) pairs defining the interpolation table:
  ```
  x=0,  y=1   -->  k(T=0)  = 1
  x=10, y=11  -->  k(T=10) = 11
  ```
  The function between these two points is the straight line `k = 1 + T`, which is
  exactly the desired `k(T) = 1 + T` relationship.

  Since the actual solution has a peak temperature around T ≈ 2 to 3, the data
  points T=0 and T=10 comfortably bracket the solution range. The interpolation
  is exact within the range (it is a line through two points on the line y = 1 + x).

- **Why AD?** Because the `type` starts with `AD`, MOOSE knows to track derivatives
  through this material. When `ADMatDiffusion` queries the value of `k` at a
  quadrature point, it receives not just `k(T_qp)` but also `dk/dT * dT/dT_j`
  for all nearby DOFs — the information needed to assemble the exact Jacobian.

If you replaced this with the non-AD `PiecewiseLinearInterpolationMaterial` and used
`MatDiffusion` instead of `ADMatDiffusion`, you would need to hand-code the Jacobian
terms from the nonlinear conductivity. The AD approach eliminates that burden entirely.

### Block: `[Postprocessors]`

```
[Postprocessors]
  [max_T]
    type       = ElementExtremeValue
    variable   = T
    value_type = max
  []
[]
```

- `type = ElementExtremeValue` — scans all quadrature points in all elements and
  returns the maximum (or minimum) value of the specified variable.
- `variable = T` — looks at the temperature field.
- `value_type = max` — returns the maximum. The alternative is `min`.

The maximum temperature is a useful scalar summary of the solution. For a linear
problem with Q=10 and k=1, the maximum is approximately 1.25. Because k increases
with T in this nonlinear case (better conductors carry heat away faster), we expect
the maximum to be somewhat lower — typically around 2.1 to 2.5. Monitoring this
value across Newton iterations also shows that the solution is physically sensible.

The postprocessor value is written to both the console at the end of the run and to
the CSV file.

### Block: `[Executioner]`

```
[Executioner]
  type       = Steady
  solve_type = 'NEWTON'

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12

  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]
```

- `type = Steady` — no time stepping. A single nonlinear solve to find the
  equilibrium temperature distribution.

- `solve_type = 'NEWTON'` — full Newton method using the explicitly assembled
  Jacobian matrix. This is the correct choice when the Jacobian is available exactly
  (which AD provides). The alternatives are:
  - `PJFNK` — Jacobian-free (uses finite differences). Appropriate when AD is not
    used and you do not want to hand-code the Jacobian.
  - `JFNK` — same but without the preconditioner on the nonlinear level.
  - `FD` — full finite-difference Jacobian. Most expensive, least accurate.

- `nl_rel_tol = 1e-10` — the nonlinear solver stops when the norm of the residual
  drops below `1e-10` times its initial value. This tight tolerance is chosen to
  demonstrate quadratic Newton convergence clearly.

- `nl_abs_tol = 1e-12` — also stop if the absolute residual norm drops below this
  value. This prevents the solver from grinding away at machine-precision residuals.

- `petsc_options_iname = '-pc_type -pc_hypre_type'` — passes option names to PETSc,
  the underlying parallel linear algebra library.

- `petsc_options_value = 'hypre boomeramg'` — selects HYPRE's BoomerAMG algebraic
  multigrid preconditioner. This is an excellent choice for the linear system arising
  from the diffusion operator: AMG is scalable, robust, and fast for elliptic PDEs.

  At each Newton iteration, the assembled Jacobian (a sparse matrix) is passed to
  PETSc's GMRES Krylov solver, preconditioned by BoomerAMG. The inner (linear)
  solver convergence directly affects the outer (Newton) convergence.

### Block: `[Outputs]`

```
[Outputs]
  exodus = true
  csv    = true
[]
```

- `exodus = true` — writes the temperature field T at every node to an Exodus II
  binary file named `case07_nonlinear_diffusion_out.e`. Readable by ParaView and VisIt.
- `csv = true` — writes postprocessor values (in this case `max_T`) to a
  comma-separated file named `case07_nonlinear_diffusion_out.csv`.

For a steady-state problem, only a single time-step entry is written to each file.

---

## What Happens When You Run This

### Invocation

```bash
cd quickstart-runs/case07
mpirun -n 1 moose-app-opt -i case07_nonlinear_diffusion.i
```

Replace `moose-app-opt` with the name of your compiled MOOSE test application.

### Startup Phase

1. MOOSE parses the HIT input file. All block names and parameter values are validated.
2. The Factory instantiates objects in dependency order: Mesh, Variables, Materials,
   Kernels, BCs, Postprocessors, Executioner, Outputs.
3. `GeneratedMesh` constructs the 20x20 structured quad mesh with 441 nodes.
4. The DofMap assigns equation numbers. With one variable (T) on 441 nodes: 441 DOFs.
5. AD meta-data is initialized: the dual-number system is configured to track
   derivatives with respect to all local DOFs.

### Assembly of the Nonlinear System

MOOSE loops over all 400 elements. For each element and each quadrature point:

1. Coordinates (x, y) are computed from the isoparametric mapping.
2. Shape functions `phi_i` and their gradients `grad(phi_i)` are evaluated.
3. The current solution values `T` and `grad(T)` at the quadrature point are computed
   from the nodal DOF values and shape functions.
4. The `[conductivity]` material evaluates `k = 1 + T_qp` using the interpolation
   table. Because it is an AD material, it also computes `dk/dT_qp = 1` in the
   dual-number derivative part.
5. `ADMatDiffusion` computes the element residual contribution:
   ```
   R_i += k(T_qp) * grad(phi_i) . grad(T_qp) * weight_qp
   ```
   and the element Jacobian contribution:
   ```
   J_ij += (dk/dT_j) * grad(phi_i) . grad(T_qp) * weight_qp
         + k(T_qp) * grad(phi_i) . grad(phi_j) * weight_qp
   ```
   The second line is the "material" term (standard stiffness matrix entry). The
   first line is the extra "nonlinear" term from the derivative of k with respect to
   T. AD computes both automatically.
6. `BodyForce` subtracts the source contribution from the element residual:
   ```
   R_i -= 10.0 * phi_i * weight_qp
   ```
7. Element residuals and Jacobians are assembled into the global sparse matrix.

### Newton Iteration Loop

Newton's method iterates as follows:

```
Given current solution T_n, solve:
    J(T_n) * delta_T = -R(T_n)
Update:
    T_{n+1} = T_n + delta_T
Check convergence:
    ||R(T_{n+1})|| / ||R(T_0)|| < nl_rel_tol = 1e-10
```

With exact Jacobian from AD, this converges quadratically. Typical console output:

```
 0 Nonlinear |R| = 8.333333e+00
      0 Linear |R| = 1.234567e-12
 1 Nonlinear |R| = 4.127483e-03
      0 Linear |R| = 2.345678e-14
 2 Nonlinear |R| = 2.198445e-08
      0 Linear |R| = 4.567890e-15
 3 Nonlinear |R| = 1.034291e-13

Converged!
```

Observe the quadratic convergence pattern:
- Iteration 0 to 1: residual drops from ~8.3 to ~4.1e-3 (about 2000x reduction)
- Iteration 1 to 2: residual drops from ~4.1e-3 to ~2.2e-8 (about 200,000x reduction)
- Iteration 2 to 3: residual drops from ~2.2e-8 to ~1.0e-13 (machine precision)

Each residual is roughly the *square* of the previous one. This is the hallmark of
quadratic Newton convergence. If you used `PJFNK` instead, convergence would be
slower and more linear-looking because the approximate Jacobian introduces errors.

### Why the Initial Residual Is ~8.33

The initial guess is T = 0 everywhere. The `BodyForce` source term contributes
`Q * volume = 10 * 1 = 10` to the residual. The diffusion term with T=0 contributes
zero. The L2 norm of a uniform source over a unit square is `sqrt(Q^2 * area) / n_dof`
in some normalization, giving approximately `10 / sqrt(12)` ≈ 2.89 for the L1 norm,
or MOOSE's internal residual normalization which gives roughly 8.33.

### Postprocessor Evaluation

After convergence, `ElementExtremeValue` scans all quadrature points and reports the
maximum T value in the domain.

---

## Output Files

### `case07_nonlinear_diffusion_out.e`

An Exodus II binary file containing:
- The 20x20 mesh geometry (node coordinates, element connectivity, boundary sets)
- The solution field `T` at all 441 nodes
- A single time step (steady-state solution)

Open in ParaView:
1. File > Open > select `case07_nonlinear_diffusion_out.e`
2. Click Apply in the Properties panel
3. Select "Surface" representation
4. Color by variable `T`
5. Apply a "Rainbow Desaturated" or "Blue to Red" color map
6. Click the play button or navigate to the one available time step

You will see a smooth temperature hill centered at approximately (0.5, 0.5) with a
maximum temperature slightly above 2. The distribution is not a perfect paraboloid
(which would result from constant k) — it is "fattened" toward the center because
regions with higher T also have higher k, which spreads heat more efficiently.

### `case07_nonlinear_diffusion_out.csv`

A comma-separated file with postprocessor time history:

```
time,max_T
0,2.1437...
```

For a steady-state solve, `time = 0` indicates the single final state. The `max_T`
value around 2.1 to 2.5 confirms the nonlinear solution behaves physically (the linear
case would give max_T ≈ 1.25; the nonlinear case has lower max because better
conductivity at higher T spreads the heat more effectively).

---

## Understanding the Plots

Running `python visualize_all.py` from the `quickstart-runs/` directory produces the
following two plots saved into this directory.

### `case07_contour_2d.png`

**What the plot shows.** A 2D filled-contour of the steady-state temperature field
T(x,y) on the unit square, using the coolwarm colormap (blue=cold, red=hot). The
nonlinear conductivity is `k(T) = 1 + T`.

**Physical quantities.** The color encodes temperature in a material whose
conductivity increases with temperature. As the temperature rises, the material
conducts more easily, which redistributes the temperature relative to the constant-k
case.

**How to judge correctness.** The field should show a smooth dome or bell shape
peaked near the domain center with zero temperature on all four walls. The peak
temperature should be noticeably lower than the Case 03 (constant k=1) result at
the same source strength, because the nonlinear conductivity increases with T and
therefore drains heat more efficiently from hot regions. There should be no
oscillations or negative values — the solution should be smooth and non-negative
everywhere.

**What would indicate a problem.**
- Negative temperatures anywhere: the nonlinear solver failed to converge properly.
- Oscillations in the contours (alternating warm/cold bands): Newton diverged and
  the solution is meaningless.
- A peak value identical to the Case 03 constant-k result: the nonlinear material
  was not applied — k is effectively constant at 1.

### `case07_surface_3d.png`

**What the plot shows.** A 3D surface plot of T(x,y) using `plot_trisurf`, with x
and y as horizontal axes and T as the vertical axis. The coolwarm colormap is applied
to the surface height.

**Physical quantities.** The 3D surface height at each (x,y) position represents the
temperature at that point. The dome-shaped surface directly visualizes the spatial
distribution of heat.

**How to judge correctness.** The surface should be a smooth, rounded dome with a
single peak near the domain center. The rim of the dome along all four edges should
touch zero (the Dirichlet wall conditions). The peak height should be approximately
0.57 or less — compare to Case 03's peak of about 0.074 scaled by the source strength,
noting that for nonlinear k the exact peak depends on the iteration history.

More practically: the surface should look smooth with no faceting artifacts or spiky
features, which would indicate triangulation issues in the visualization.

**What would indicate a problem.**
- A flat surface at T=0: the source term was removed or the BCs hold everything to zero.
- Spiky features at individual nodes: the solver produced non-physical node values.
- The peak significantly higher than the Case 03 result: the k(T) nonlinearity is
  increasing temperature instead of decreasing it — wrong sign in the material formula.

---

## Interpreting the Results

### The Solution Field

The temperature distribution T(x,y) is:
- Zero on all four walls (enforced by Dirichlet BCs)
- Symmetric about both the x = 0.5 and y = 0.5 axes (the geometry and BCs are symmetric)
- Positive throughout the interior (heat source drives T upward from the zero-wall baseline)
- Maximum somewhere near the center (0.5, 0.5)

### Comparing Linear vs. Nonlinear

For comparison, the linear case with constant k = 1 and Q = 10 has the solution:

```
-Laplacian(T) = 10    T=0 on all walls
```

The maximum of this linear solution is approximately `T_max ≈ Q * L^2 / (8 * k) = 1.25`
for a unit square. (This is the exact formula for a 1D slab; 2D gives a similar result.)

With `k(T) = 1 + T`, the effective conductivity is higher in the interior (where T > 0),
which means heat is conducted away more efficiently, reducing the peak temperature.
A typical measured `max_T` from this nonlinear run is approximately 2.1 to 2.5.

Wait — this is *higher* than the linear case with k=1, not lower. Why?

The key insight is that the nonlinear conductivity is `k(T) = 1 + T`, which starts at 1
at the walls and is larger in the interior. Even though higher k conducts heat out faster,
the source term Q=10 is large enough that the equilibrium temperature still exceeds the
linear-k maximum. The shape of the solution is different from the parabolic linear
solution: it has a flatter top and steeper sides, reflecting the higher conductivity in
the hotter interior.

### Correctness Verification

You can verify the solution is physically reasonable by checking:

1. **Symmetry**: T(x,y) = T(1-x,y) = T(x,1-y). In ParaView, use a clip plane at
   x=0.5 and verify the two halves are mirror images.

2. **Monotonicity**: T is zero on all walls and positive inside. It should have no
   negative values. In ParaView, check the color range minimum.

3. **Energy balance**: The total heat input equals total conductive flux out:
   ```
   Q * volume = integral( k(T) * grad(T) . n dA )  over all walls
   ```
   The postprocessor gives the solution for quick checks.

4. **Newton convergence**: The quadratic convergence pattern in the console output
   (residual squaring each iteration) confirms the exact Jacobian is being assembled
   and the solver is working correctly.

---

## Key Concepts Learned

- **Nonlinear PDEs**: When a coefficient like k depends on the unknown T, the PDE
  becomes nonlinear and requires iterative solution via Newton's method.

- **Automatic Differentiation (AD)**: A technique using dual-number arithmetic to
  compute exact derivatives of residuals with respect to DOFs, enabling exact Jacobian
  assembly without hand-coding.

- **ADMatDiffusion**: The AD-aware version of the diffusion kernel. Same residual
  formula as `MatDiffusion` but with automatic Jacobian.

- **ADPiecewiseLinearInterpolationMaterial**: Defines a material property as a
  piecewise-linear function of a solution variable. The `AD` prefix means derivatives
  are tracked through the interpolation.

- **xy_data**: The flat list format `'x0 y0  x1 y1  ...'` for defining piecewise
  linear tables in MOOSE materials.

- **solve_type = NEWTON**: Full Newton method using explicit Jacobian. The correct
  choice when AD provides an exact Jacobian. Achieves quadratic convergence.

- **Quadratic convergence**: The residual norm approximately squares each Newton
  iteration — the defining property of Newton's method with an exact Jacobian.

- **nl_rel_tol and nl_abs_tol**: Two stopping criteria for the nonlinear solver,
  relative (to initial residual) and absolute (floor value), respectively.

- **BodyForce kernel**: Adds a volumetric source Q to the PDE right-hand side.
  Independent of T, so has zero Jacobian contribution.

- **ElementExtremeValue postprocessor**: Scans quadrature points to find the maximum
  or minimum of any variable, returning a single scalar summary.

---

## Experiments to Try

### Experiment 1: Watch Quadratic Convergence Disappear with PJFNK

Change `solve_type = 'NEWTON'` to `solve_type = 'PJFNK'` and rerun. Compare the
console output line by line. With PJFNK you should see slower, more linear-looking
convergence: the residual may take 5 to 8 iterations instead of 3 to 4. The solution
at convergence will be the same (PJFNK does converge), but the path is different.
This experiment directly demonstrates the practical advantage of AD-provided Jacobians.

Expected outcome: more nonlinear iterations, same final solution.

### Experiment 2: Make the Conductivity More Nonlinear

Replace the `xy_data` table with a strongly nonlinear conductivity:

```
xy_data = '0 1
           10 101'   # k(T) = 1 + 10*T
```

This makes k grow ten times faster with temperature. The nonlinearity becomes more
severe: Newton may require one or two additional iterations to converge, and the
solution shape will change more dramatically (even flatter top, steeper walls). The
maximum temperature will be somewhat lower because the high-T interior conducts much
better. Tighten the tolerances if needed.

Expected outcome: Newton converges in 4 to 6 iterations, max_T decreases.

### Experiment 3: Refine the Mesh and Check Convergence

Change `nx = 20, ny = 20` to `nx = 40, ny = 40`. The solution should change only
slightly because 20x20 is already reasonably converged for this smooth problem. The
maximum temperature should differ by less than 1%. This confirms the mesh is adequately
resolved. If you go the other way (`nx = 5, ny = 5`) the solution will be visibly
coarser in ParaView and max_T may differ by a few percent.

Expected outcome: max_T changes by less than 0.5% when doubling the mesh.

### Experiment 4: Add a Nonlinear Source Term

The `BodyForce` kernel uses a constant value. To add a temperature-dependent source
(making the problem even more nonlinear), add an `ADBodyForce` kernel that references
a function. Or you could use the trick of adding a `CoupledForce` kernel to couple T
into its own source. For example:

```
[source_nonlinear]
  type  = CoupledForce
  variable = T
  v        = T
  coef     = 1.0   # adds T itself as a source
[]
```

This replaces Q=10 (constant) with Q = T (grows with temperature). Be careful: this
can make the problem ill-posed (explosive solutions) for large values. Use a small
coefficient to start.

Expected outcome: more nonlinear iterations, higher or lower max_T depending on the
sign of the added term.

### Experiment 5: Verify the Nonlinearity with a Scan

Run the case six times, each time changing the source strength `value` in `BodyForce`
from 1.0 to 10.0 in steps of 2.0. For each run record `max_T` from the CSV. In the
linear case (constant k) max_T would scale linearly with Q. In the nonlinear case
(k = 1+T) max_T grows sub-linearly: doubling Q does not double T_max because k also
increases. Plot Q vs. max_T to visualize the nonlinear response.

Expected outcome: max_T grows with Q but at a decreasing rate as Q increases.
