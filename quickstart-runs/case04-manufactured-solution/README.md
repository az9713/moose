# Case 04: Poisson Equation with Method of Manufactured Solutions

## Overview

This case introduces one of the most powerful and widely-used techniques in scientific
computing for verifying simulation codes: the **Method of Manufactured Solutions (MMS)**.

Instead of solving a problem and hoping the answer is right, MMS works in reverse. You
*start* with a known answer, derive the source term that would produce that answer, then
solve the problem and check how closely the numerical result matches the known answer.
The gap between the two — measured with a mathematical norm called the L2 error — tells
you whether your solver is working correctly and with the expected order of accuracy.

This case solves the **Poisson equation** (steady-state diffusion with a source term) on
a unit square domain using a sinusoidal exact solution. Because the exact answer is known,
the test is self-contained and fully verifiable without any experimental data.

Why does this matter? In real engineering simulations the true answer is never known. MMS
lets developers *prove* their code is correct on simple, controllable problems before
trusting it on complicated ones.

---

## The Physics

### The Physical Problem

Imagine a thin, flat plate occupying the unit square (x from 0 to 1, y from 0 to 1).
The plate has some physical quantity `u` distributed across it — think of it as temperature,
concentration of a chemical, or electrostatic potential. Heat (or mass, or charge) diffuses
through the plate according to Fick's law or Fourier's law.

At every interior point, diffusion carries the quantity away from high-concentration regions.
If an internal source `f` continuously adds the quantity, a steady state is reached where
diffusion away from each point exactly balances production at that point.

The governing equation states: the rate at which `u` spreads out via diffusion equals the
source that replenishes it.

### The Governing Equation

```
-div(grad u) = f    on the domain [0,1] x [0,1]
```

In two dimensions (x, y) this expands to:

```
  d²u     d²u
- --- - ------- = f(x, y)
  dx²     dy²
```

Symbol definitions:
- `u` — the scalar field being solved for (temperature, concentration, potential, etc.)
- `div` — the divergence operator: takes a vector field and returns a scalar
- `grad u` — the gradient of u: a vector pointing in the direction of steepest increase of u
- `-div(grad u)` — the Laplacian of u with a minus sign; measures how much u at a point
  differs from the average of its neighbors. Positive means u is locally low.
- `f(x, y)` — the source term: how much `u` is injected per unit area per unit time

When `f = 0` this is Laplace's equation. With `f` nonzero it is Poisson's equation.

### Boundary Conditions

On all four sides of the square, `u` is fixed to a prescribed value. This is called a
**Dirichlet boundary condition** — the most intuitive type, equivalent to saying "the
temperature at this wall is exactly X."

For this case, the boundary values are determined by evaluating the exact solution on the
boundary. Because `sin(pi*x)` is zero at `x=0` and `x=1`, and `sin(pi*y)` is zero at
`y=0` and `y=1`, the boundary value is zero on all four sides.

```
u = 0    on x=0 (left),  x=1 (right),  y=0 (bottom),  y=1 (top)
```

### The Manufactured Solution

The key idea: we *choose* the exact solution first, then derive `f` from it.

**Chosen exact solution:**

```
u_exact(x, y) = sin(pi*x) * sin(pi*y)
```

This function:
- Is zero on all four boundaries (as computed above — satisfies our BCs naturally)
- Is smooth and positive in the interior
- Has a single peak at (0.5, 0.5) with value 1.0

**Deriving the source term:**

Compute `-div(grad u_exact)` analytically:

```
d/dx [sin(pi*x)*sin(pi*y)] =  pi*cos(pi*x)*sin(pi*y)
d²/dx² [...]               = -pi²*sin(pi*x)*sin(pi*y)

d/dy [sin(pi*x)*sin(pi*y)] =  pi*sin(pi*x)*cos(pi*y)
d²/dy² [...]               = -pi²*sin(pi*x)*sin(pi*y)

-div(grad u_exact) = -(-pi² - pi²)*sin(pi*x)*sin(pi*y)
                   =  2*pi²*sin(pi*x)*sin(pi*y)
```

**Therefore the manufactured source term is:**

```
f(x, y) = 2*pi²*sin(pi*x)*sin(pi*y)
```

If we solve `-div(grad u) = f` with this `f` and zero Dirichlet BCs, the exact answer
is `sin(pi*x)*sin(pi*y)`. Any finite element solution that is not exactly this (and it
won't be, due to discretization error) is in error by a quantifiable amount.

### ASCII Domain Diagram

```
  y=1  u=0 (top boundary)
       +------------------------------------+
       |            *                       |
       |         *     *                    |
       |       *         *                  |  u=0
  u=0  |      *   PEAK    *                 |  (right)
(left) |     *   u~1 at    *                |
       |      *  (0.5,0.5) *                |
       |       *         *                  |
       |         *     *                    |
       |            *                       |
       +------------------------------------+
  y=0  u=0 (bottom boundary)
       x=0                              x=1

  Exact solution: u = sin(pi*x)*sin(pi*y)
  Single smooth peak at center, zero on all walls
```

### The L2 Error Norm

After solving, the L2 error measures how far the numerical solution `u_h` is from the
exact solution `u_exact` over the entire domain:

```
||u_h - u_exact||_L2 = sqrt( integral of (u_h - u_exact)^2 dV )
```

For a well-implemented second-order finite element method with mesh spacing `h`, theory
predicts the L2 error scales as `h²`. This means halving the mesh spacing should reduce
the error by a factor of 4. This is called **second-order convergence**.

With a 20x20 mesh (as used here), `h = 0.05`, so the L2 error should be on the order of
`h² * C` where `C` is a problem-dependent constant.

---

## Input File Walkthrough

The input file is `case04_manufactured.i`. Let us examine every block and parameter.

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
  Requires no external mesh file.
- `dim = 2` — two-dimensional problem. The domain is a rectangle in (x, y).
- `nx = 20` — number of elements in the x-direction. With the default domain [0,1],
  each element is 0.05 units wide.
- `ny = 20` — number of elements in the y-direction. Same spacing vertically.

The mesh has 20 x 20 = 400 quadrilateral elements and (21 x 21) = 441 nodes. Default
domain is [0,1] x [0,1] unless `xmin`, `xmax`, `ymin`, `ymax` are specified.

What happens if you change `nx` and `ny`? A finer mesh (larger values) gives a smaller
L2 error; a coarser mesh gives a larger error. Doubling both from 20 to 40 should
reduce the L2 error by approximately a factor of 4 for a second-order method.

### Block: `[Variables]`

```
[Variables]
  [u]
  []
[]
```

- Declares the single unknown field `u` to be solved for.
- The empty sub-block `[u]` uses all defaults: first-order Lagrange shape functions,
  initial condition of zero.
- MOOSE automatically allocates degrees of freedom (DOFs) at every mesh node for `u`.
- There are 441 DOFs for this 21x21 node mesh.

First-order Lagrange elements are the simplest: `u` is approximated as a piecewise
bilinear function over the mesh, continuous across element boundaries.

### Block: `[Functions]`

```
[Functions]
  [exact_fn]
    type       = ParsedFunction
    expression = 'sin(pi*x)*sin(pi*y)'
  []

  [forcing_fn]
    type       = ParsedFunction
    expression = '2*pi^2*sin(pi*x)*sin(pi*y)'
  []
[]
```

This block defines mathematical functions that can be evaluated at any point (x, y, z)
and time t. Functions are a central concept in MOOSE — they are not just constants but
callable objects that other MOOSE objects reference by name.

**`[exact_fn]`:**
- `type = ParsedFunction` — evaluates a mathematical expression string. MOOSE uses
  the `fparser` library to parse and JIT-compile the expression at startup.
- `expression = 'sin(pi*x)*sin(pi*y)'` — the manufactured exact solution.
- Built-in variables: `x`, `y`, `z`, `t`. Built-in constants: `pi`, `e`.
- This function is used in two places: as the Dirichlet BC (boundary condition) and
  as the reference solution in the L2 error postprocessor.

**`[forcing_fn]`:**
- `expression = '2*pi^2*sin(pi*x)*sin(pi*y)'` — the analytically derived source term.
- `pi^2` means pi squared. MOOSE's expression parser supports `^` for exponentiation.
- This function is evaluated at every quadrature point inside every element during
  assembly of the right-hand side vector.
- At the domain center (x=0.5, y=0.5): f = 2*pi²*1*1 ≈ 19.74. At the corners: f = 0.

### Block: `[Kernels]`

```
[Kernels]
  [diffusion]
    type     = Diffusion
    variable = u
  []

  [source]
    type     = BodyForce
    variable = u
    function = forcing_fn
  []
[]
```

Kernels are the core of MOOSE — they implement the terms of the weak form of the PDE.
Each kernel contributes to the global residual and Jacobian.

**Weak form background:** The finite element method converts the strong-form PDE into
a weak form by multiplying by a test function `phi_i` and integrating over the domain.
Integration by parts turns second derivatives into first derivatives (making the
problem solvable with piecewise linear functions).

**`[diffusion]` kernel:**
- `type = Diffusion` — implements the term `integral( grad(phi_i) . grad(u) dV )`.
- This is the diffusion operator `-div(grad u)` expressed in weak form.
- No coefficient: this assumes unit diffusivity (k = 1).
- Every element contributes a 4x4 (for quad elements) element stiffness matrix.

**`[source]` kernel:**
- `type = BodyForce` — implements the term `integral( phi_i * f dV )`.
- `function = forcing_fn` — references the `forcing_fn` function defined above.
- At each quadrature point, MOOSE calls `forcing_fn.value(x, y, z, t)` to get `f`.
- This adds the right-hand-side load vector contribution from the source term.
- Without this kernel, you would be solving `-div(grad u) = 0`, which with zero
  boundary conditions has only the trivial solution u = 0.

Together, these two kernels represent the complete weak form:
```
integral( grad(phi_i) . grad(u) dV ) = integral( phi_i * f dV )
```

### Block: `[BCs]`

```
[BCs]
  [all_walls]
    type     = FunctionDirichletBC
    variable = u
    boundary = 'left right top bottom'
    function = exact_fn
  []
[]
```

- `type = FunctionDirichletBC` — applies a Dirichlet (fixed-value) boundary condition
  by evaluating a Function at each boundary node.
- `boundary = 'left right top bottom'` — applies to all four named boundaries. These
  names are automatically assigned by `GeneratedMesh`.
- `function = exact_fn` — references our exact solution function.

Since `sin(pi*x)*sin(pi*y)` evaluates to zero on all four boundaries (as shown in the
physics section), this is equivalent to `DirichletBC` with `value = 0`. However, using
`FunctionDirichletBC` with the exact solution makes the intent explicit: the boundary
condition is *derived from* the manufactured solution.

This design choice is important for MMS in general. If you change the manufactured
solution, you only need to update the `exact_fn` expression — the boundary condition
automatically updates too.

### Block: `[Postprocessors]`

```
[Postprocessors]
  [L2_error]
    type     = ElementL2Error
    variable = u
    function = exact_fn
  []
[]
```

- `type = ElementL2Error` — computes the L2 norm of the difference between the
  numerical solution and a reference function, integrated over all elements.
- `variable = u` — the numerical solution to compare.
- `function = exact_fn` — the reference (exact) solution.

The formula computed is:
```
L2_error = sqrt( sum over all elements of: integral( (u_h - u_exact)^2 dV ) )
```

This single number summarizes the total error in the entire solution. For a correct
second-order FEM implementation, it should be approximately `h² * C` where `h` is the
element size. For our 20x20 mesh with h=0.05, the expected L2 error is on the order
of `0.05² * C = 0.0025 * C`.

The result is written to both the console at the end of the run and to the CSV file.

### Block: `[Executioner]`

```
[Executioner]
  type       = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]
```

- `type = Steady` — steady-state solve. No time stepping is performed.
- `solve_type = 'PJFNK'` — **Preconditioned Jacobian-Free Newton-Krylov**. This is
  MOOSE's default nonlinear solver strategy.
  - Newton's method linearizes the residual to find the correction step.
  - JFNK approximates matrix-vector products with the Jacobian using finite differences
    of the residual, avoiding explicit Jacobian assembly.
  - Krylov refers to the iterative linear solver (GMRES) used at each Newton step.
  - The preconditioner (PC) improves convergence of the Krylov solver.
- `petsc_options_iname = '-pc_type -pc_hypre_type'` — passes options to PETSc.
- `petsc_options_value = 'hypre boomeramg'` — selects HYPRE's BoomerAMG algebraic
  multigrid as the preconditioner. AMG is highly effective for elliptic PDEs like this.

For a linear problem like Poisson's equation, PJFNK converges in one Newton iteration
because the Jacobian is the same as the stiffness matrix (no nonlinearity). The nonlinear
solver is overkill here but is used consistently across cases.

### Block: `[Outputs]`

```
[Outputs]
  exodus = true
  csv    = true
[]
```

- `exodus = true` — writes the solution field `u` to an Exodus II file. This binary
  format is readable by ParaView, VisIt, and Cubit.
- `csv = true` — writes all postprocessor values (in this case `L2_error`) to a
  comma-separated values file.

Output filenames are auto-derived from the input filename: `case04_manufactured_out.e`
and `case04_manufactured_out.csv`.

---

## What Happens When You Run This

### Invocation

```bash
cd quickstart-runs/case04
mpirun -n 1 moose-app-opt -i case04_manufactured.i
```

(Replace `moose-app-opt` with the actual executable name for your MOOSE application.)

### Startup Phase

1. MOOSE parses `case04_manufactured.i` using the HIT parser.
2. The Factory instantiates all objects: Mesh, Variables, Functions, Kernels, BCs,
   Postprocessors, Executioner, and Outputs.
3. `GeneratedMesh` creates the 20x20 structured quad mesh.
4. The two `ParsedFunction` expressions are compiled by `fparser`.
5. The DofMap assigns equation numbers to all 441 nodes.

### Assembly Phase

MOOSE loops over all 400 elements. For each element:
1. Maps element geometry (Jacobian of isoparametric transformation).
2. Evaluates shape functions and their gradients at Gauss quadrature points.
3. Calls each Kernel to add its contributions to the element residual and Jacobian.
4. For the `Diffusion` kernel: fills 4x4 element stiffness matrix.
5. For the `BodyForce` kernel: evaluates `forcing_fn` at quadrature points, fills
   element load vector.
6. Assembles element contributions into the global sparse matrix and right-hand side.

### Boundary Condition Application

Before solving, MOOSE enforces the Dirichlet BCs. For `FunctionDirichletBC`, it
evaluates `exact_fn` at each boundary node to get the prescribed value (in this case,
all zeros), then modifies rows of the system matrix to enforce `u_node = 0`.

### Linear Solve

Since Poisson's equation is linear, Newton's method converges in **one iteration**:
- The linear system `K * u = f` is assembled (K is the stiffness matrix, f the load).
- PETSc's GMRES + BoomerAMG preconditioner solves the system.
- Convergence typically occurs in a few Krylov iterations for this small problem.

### Postprocessor Evaluation

After the solve, `ElementL2Error` integrates `(u_h - exact_fn)^2` over the mesh using
Gauss quadrature, takes the square root, and stores the result.

### Typical Console Output

```
Framework Information:
  MOOSE Version:          git commit XXXXXXX

Mesh Information:
  Spatial dimension:      2
  Mesh:                   400 elements, 441 nodes

Nonlinear System:
  Num DOFs:               441

 0 Nonlinear |R| = 9.869604e+00
      0 Linear |R| = 2.049390e-13
 1 Nonlinear |R| = 3.127912e-13

Converged!

Postprocessor Values:
  L2_error = 2.4168e-04
```

- The nonlinear residual starts at ~9.87 (the norm of the right-hand side).
- One Newton step drives it to machine precision (~3e-13).
- The L2 error of ~2.4e-4 is the quantitative verification metric.

---

## Output Files

### `case04_manufactured_out.e`

The Exodus II file contains:
- The mesh geometry (node coordinates, element connectivity)
- The solution field `u` at every node
- Can be opened in ParaView: File > Open, select the `.e` file, click Apply

In ParaView, apply a color map to `u` to see the smooth sinusoidal hill centered at (0.5, 0.5).

### `case04_manufactured_out.csv`

A comma-separated file with postprocessor values:

```
time,L2_error
0,2.4168e-04
```

- `time = 0` — steady-state problems are output at "time" 0.
- `L2_error` — the computed L2 error norm between the FEM solution and the exact solution.

This file is the primary verification artifact. To perform a convergence study, run
with different mesh sizes and plot L2_error vs h on a log-log scale.

---

## Understanding the Plots

Running `python visualize_all.py` from the `quickstart-runs/` directory produces the
following three plots saved into this directory.

### `case04_numerical.png`

**What the plot shows.** A 2D filled-contour of the MOOSE numerical solution `u_h(x,y)`
on the unit square. The color axis and contour levels are shared with the exact-solution
plot below so that both plots can be compared directly.

**Physical quantities.** The color encodes the value of the solution field at each
point. For this manufactured solution `u = sin(pi*x) * sin(pi*y)`, the peak value
(approximately 1.0) appears at the domain center (0.5, 0.5) and the field decays to
exactly 0 on all four boundaries.

**How to judge correctness.** The contours should form a smooth, concentric dome
pattern — roughly circular near the center, becoming squarish near the corners. The
peak should be at the exact center. The colormap range should be 0 to 1.

**What would indicate a problem.**
- Flat contours or a uniform color: the source term f was not applied or was computed
  incorrectly.
- Peak shifted away from center: asymmetric BCs or incorrect function expression.
- Negative values: the manufactured solution formula was computed incorrectly.

### `case04_exact.png`

**What the plot shows.** A 2D filled-contour of the exact manufactured solution
`u_exact = sin(pi*x) * sin(pi*y)`, evaluated at the same node positions as the
numerical solution and plotted with the same color scale.

**Physical quantities.** Same as the numerical plot. The exact solution is computed
analytically in Python using `np.sin(np.pi * x) * np.sin(np.pi * y)`.

**How to judge correctness.** This plot should look visually identical to
`case04_numerical.png` at the scale of the plot. Both share the same color axis,
so any visible color difference between them directly reveals numerical error.

**What would indicate a problem.** If this plot looks noticeably different from the
numerical plot, there is a large numerical error — which would also show up in
the error plot below.

### `case04_error.png`

**What the plot shows.** A 2D filled-contour of the pointwise absolute error
`|u_h - u_exact|`, plotted with the magma colormap (dark=small, bright=large).

**Physical quantities.** The color encodes how far the MOOSE solution deviates from
the exact manufactured solution at each mesh node. This is the primary correctness
metric for MMS.

**How to judge correctness.** For a 20x20 mesh with first-order Lagrange elements,
the error should be small — on the order of `O(h^2)` where `h = 1/20 = 0.05`. This
gives expected errors of roughly `h^2 = 0.0025`, so the colormap maximum should be
in the range 1e-3 to 1e-2. The error is typically largest near the domain boundaries
(where the solution gradient is large) and smallest at the center (where it is smooth).

**What would indicate a problem.**
- Errors close to 1 (same order as the solution itself): the solver failed or the
  source term is wrong.
- Errors below 1e-14: this is machine precision — if you see this, you have
  accidentally tested a problem that is exactly representable by first-order elements
  (which should not happen for a sinusoidal manufactured solution).
- Error concentrated in one corner: a BC is missing or incorrect on one edge.

---

## Interpreting the Results

### The Solution Field

The numerical solution `u_h` approximates `sin(pi*x)*sin(pi*y)` over the 20x20 mesh.
It shows:
- A smooth dome shape with a single maximum near (0.5, 0.5)
- Value of approximately 1.0 at the peak
- Zero on all four boundaries
- Perfect fourfold symmetry about both axes (the physics is symmetric)

### Verifying Correctness

For second-order Lagrange elements (which is what default MOOSE uses), the theoretical
L2 error estimate is:

```
||u_h - u_exact||_L2 <= C * h^2
```

With h = 1/20 = 0.05:

```
expected order: L2_error ~ C * 0.0025
```

The actual value near 2.4e-4 is consistent with this estimate. To confirm second-order
convergence, run with nx=ny=10, 20, 40 and verify that each halving of h reduces the
error by a factor of ~4.

| nx=ny | h     | L2 error   | Ratio |
|-------|-------|------------|-------|
| 10    | 0.10  | ~9.7e-4    |  —    |
| 20    | 0.05  | ~2.4e-4    | ~4.0  |
| 40    | 0.025 | ~6.1e-5    | ~4.0  |

This table confirms that the solver is correctly implementing a second-order FEM scheme.

### Physical Insight

The fact that MMS *works* at all — that the solver finds the correct answer — validates:
1. The weak form implementation in the `Diffusion` kernel
2. The `BodyForce` kernel's function evaluation
3. The Dirichlet BC enforcement
4. The linear algebra assembly and solve
5. The L2 error postprocessor integration

If any of these were buggy, the L2 error would not converge to zero at the expected rate.

---

## Key Concepts Learned

- **Method of Manufactured Solutions**: reverse-engineer source terms from a known
  exact solution to create a fully verifiable test problem
- **ParsedFunction**: evaluate arbitrary mathematical expressions as MOOSE Function
  objects, callable at any (x, y, z, t)
- **FunctionDirichletBC**: apply Dirichlet boundary conditions from a Function object
  rather than a fixed constant
- **BodyForce kernel**: add a distributed source term (right-hand side) to the PDE
- **ElementL2Error postprocessor**: compute the L2 norm of the error between the
  numerical solution and an exact reference function
- **Convergence rate verification**: use mesh refinement studies to confirm the solver
  achieves the theoretical order of accuracy
- **PJFNK solver**: understand why Newton converges in one iteration for linear problems
- **Second-order FEM convergence**: L2 error scales as h² for standard linear elements

---

## Experiments to Try

### Experiment 1: Mesh Refinement Convergence Study

Change `nx` and `ny` in the `[Mesh]` block to 10, 20, 40, and 80. Record the `L2_error`
from each run. Plot log(L2_error) vs log(h) where h = 1/nx. The slope should be
approximately 2.0. If it is not, there is a bug somewhere in the implementation.

Expected outcome: each doubling of resolution reduces the error by a factor of ~4.

### Experiment 2: Change the Manufactured Solution

Replace `exact_fn` with a polynomial:

```
expression = 'x*(1-x)*y*(1-y)'
```

The new forcing function is:
```
f = -div(grad u) = 2*y*(1-y) + 2*x*(1-x)
```

Update `forcing_fn` accordingly. The exact solution is zero on all boundaries (since
it contains `x*(1-x)` and `y*(1-y)` as factors). This polynomial is exactly representable
by second-order elements, so the L2 error should be near machine precision (< 1e-14)
for any mesh. This is a useful sanity check that the integration is exact for polynomials.

### Experiment 3: Non-Zero Boundary Conditions

Replace `exact_fn` with:

```
expression = 'x + sin(pi*x)*sin(pi*y)'
```

This solution is not zero on all boundaries. Update `forcing_fn` to match and observe
that `FunctionDirichletBC` correctly applies the non-zero boundary values. The L2 error
should still converge at second order.

### Experiment 4: Coarsen the Mesh Until the Solution Degrades

Set `nx = ny = 4` (only 16 elements). The L2 error will be large and the solution will
look "blocky" in ParaView. This illustrates why mesh refinement matters in practice.
The solution still converges to the right answer as the mesh is refined — this is the
mathematical guarantee that FEM provides.

### Experiment 5: Use Second-Order Elements

Add `order = SECOND` to the `[Variables]` block:

```
[Variables]
  [u]
    order = SECOND
  []
[]
```

Second-order (quadratic) Lagrange elements add midside nodes to each element. The
theoretical L2 convergence rate for smooth problems becomes h³ (third order). With
the sinusoidal exact solution, each halving of the mesh should now reduce the error
by a factor of ~8. The mesh size is now (2*nx+1)*(2*ny+1) nodes instead of (nx+1)*(ny+1).
