# Case 01: Steady-State 1D Diffusion

## Overview

This is the first and simplest MOOSE problem: a one-dimensional steady-state diffusion
equation on a line segment. There is no source term, no material property, no time
dependence, and no coupling to other physics. The only ingredients are the domain, one
unknown field, one differential operator, and two boundary conditions.

Despite its simplicity, Case 01 contains every structural element that every MOOSE input
file must have. Every block introduced here reappears in every subsequent case. Master
this case and the scaffolding of MOOSE will feel familiar from this point forward.

What you will learn:
- The minimum required blocks in a MOOSE input file
- How to build a structured mesh programmatically (no external mesh file needed)
- How to declare an unknown variable that MOOSE will solve for
- How to apply the diffusion (Laplace) operator using the `Diffusion` kernel
- How to fix boundary values with `DirichletBC`
- How to run a steady-state solve using the `Steady` executioner
- How to verify the result against an exact analytical solution

The exact solution to this problem is `u(x) = x` — a perfect straight line from 0 to 1.
This means after running the case you can immediately verify that MOOSE is working
correctly without any post-processing tools.

---

## The Physics

### The Physical Problem in Plain English

Imagine a thin metal rod, one meter long, lying along the x-axis from x=0 to x=1.
The left end is held at 0 degrees (perhaps immersed in ice water) and the right end is
held at 1 degree (perhaps by a small heater). No heat is generated inside the rod itself.

After a long time, the rod reaches a steady state: the temperature at every point stops
changing. What does the temperature distribution look like at steady state?

Because no heat is generated or consumed inside the rod, and the rod conducts heat
uniformly, the temperature must vary linearly from 0 at the left end to 1 at the right
end. There is no reason for the temperature profile to bulge up or dip down — it simply
ramps. The answer is `u(x) = x`.

This same mathematics describes many other physical situations:
- Electric potential along a resistive wire with fixed voltages at both ends
- Concentration of a dissolved substance between two reservoirs at fixed concentrations
- Pressure along a permeable channel between two reservoirs

The governing equation is the same in all cases: the **Laplace equation** in one dimension.

### The Governing Equation

The PDE being solved is:

```
  d²u
- ---- = 0    on the domain 0 < x < 1
  dx²
```

Every symbol:
- `u` — the scalar field we are solving for (temperature, concentration, potential, etc.)
- `x` — the spatial coordinate along the rod, ranging from 0 to 1
- `d²u/dx²` — the second derivative of u with respect to x; measures the curvature of u
- The minus sign — convention inherited from the general diffusion operator `-div(grad u)`;
  ensures the operator is positive definite (physically: diffusion smooths out peaks)
- `= 0` — there is no source or sink anywhere inside the domain

The equation says: the second derivative of u is everywhere zero. A function whose second
derivative is zero everywhere is a straight line (or, in higher dimensions, a linear
function). So we know analytically that u must be of the form `u = a*x + b` for some
constants `a` and `b`. The boundary conditions determine `a` and `b`.

In the general multi-dimensional form (which Case 02 uses), the operator is written:

```
-div(grad u) = 0
```

where `grad u` is the gradient vector and `div` is the divergence. In 1D this reduces
exactly to `-d²u/dx² = 0`.

### Boundary Conditions

Two **Dirichlet boundary conditions** pin the value of u at both ends of the rod:

```
u(0) = 0     (left boundary: x = 0)
u(1) = 1     (right boundary: x = 1)
```

A Dirichlet boundary condition specifies the value of the unknown directly. It is the
most intuitive type: "at this boundary, u equals this number." No derivative is involved.

Physical interpretation:
- `u(0) = 0` — the left end of the rod is held at exactly 0 degrees
- `u(1) = 1` — the right end is held at exactly 1 degree

Together these two conditions uniquely determine the solution. There is exactly one
straight-line function that passes through (0, 0) and (1, 1): the line `u(x) = x`.

### Exact Solution

The Laplace equation with these boundary conditions has a known analytical solution:

```
u(x) = x
```

This is a linear ramp from 0 to 1. The MOOSE numerical solution should reproduce this
exactly (to within machine precision) because linear elements can represent linear
functions exactly on any mesh.

### ASCII Domain Diagram

```
u=0                                                             u=1
 |                                                               |
 |   Rod (1D domain, x = 0 to x = 1)                           |
 +---o---o---o---o---o---o---o---o---o---o---o---o---o---o---o---+
     |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
  element element element element element element element element
  (20 total)

Solution u(x) = x:
 1.0                                                         *
 0.8                                                    *
 0.6                                               *
 0.4                                          *
 0.2                                     *
 0.0   *
      x=0                                                   x=1

 o = mesh node (21 nodes for 20 elements)
 * = value of u at corresponding x position
```

The solution is a perfectly straight ramp. Every mesh node's computed `u` value equals
its x-coordinate.

---

## Input File Walkthrough

The input file is `case01_diffusion_1d.i`. Every line is explained below.

### Header Comments

```
# ============================================================
# Case 1: Steady-State 1-D Diffusion
# Solves -d^2u/dx^2 = 0, u(0)=0, u(1)=1
# Exact solution: u(x) = x
# ============================================================
```

Comments in MOOSE input files begin with `#` and run to the end of the line. These
header comments are purely informational — MOOSE ignores them. They state the PDE being
solved, the boundary conditions, and the expected exact solution. Writing the exact
solution in comments is good practice: it lets anyone reading the file immediately check
whether the output is correct.

---

### Block: `[Mesh]`

```
[Mesh]
  type = GeneratedMesh
  dim  = 1        # 1-D problem: a line segment
  nx   = 20       # 20 uniform elements along x
  xmin = 0        # left endpoint
  xmax = 1        # right endpoint
[]
```

The `[Mesh]` block defines the geometry of the problem domain and its discretization.
Every MOOSE problem must have a mesh — there is no default.

**`type = GeneratedMesh`**

MOOSE's built-in structured mesh generator. It creates a rectangular (or line segment,
or box) mesh entirely from parameters — no external mesh file is needed. This is the
simplest way to create a mesh and is ideal for problems on simple domains.

Alternative mesh types include `FileMesh` (reads from an Exodus or other mesh file),
`StitchedMesh`, and various mesh generator compositions. For anything other than simple
rectangles you would use the `[MeshGenerators]` system instead, but `GeneratedMesh`
is perfect for the introductory cases.

**`dim = 1`**

Declares this a one-dimensional problem. The mesh lives on a line segment.
Changing this to `dim = 2` would require adding `ny` and domain y-bounds, which is
exactly what Case 02 does. Changing to `dim = 3` would add a third dimension.

The dimensionality setting affects which basis functions are used, how integrals are
computed, and how boundary names are assigned.

**`nx = 20`**

Number of elements in the x-direction. With `xmin=0` and `xmax=1`, each element spans
`1/20 = 0.05` units. The mesh has:
- 20 elements (each is a line segment in 1D)
- 21 nodes (one at each element endpoint; the interior nodes are shared between elements)

Why 20? It is a reasonable starting resolution for a 1D problem — fine enough that the
linear solution is captured exactly, coarse enough that the run is instant. For nonlinear
or higher-dimensional problems a much finer mesh may be needed.

What happens if you change `nx`? For this particular problem (which has a linear exact
solution), any value of `nx >= 1` gives the exact answer, because linear finite elements
represent linear functions exactly. For nonlinear problems, larger `nx` means smaller
discretization error.

**`xmin = 0` and `xmax = 1`**

The left and right endpoints of the domain. The domain is the interval [0, 1].

On a 1D `GeneratedMesh`, the named boundaries are:
- `left` — the node at `x = xmin`
- `right` — the node at `x = xmax`

These names are used in the `[BCs]` block below.

---

### Block: `[Variables]`

```
[Variables]
  [u]
  []
[]
```

The `[Variables]` block declares the unknown fields that MOOSE will solve for. There
is exactly one variable here: `u`.

**`[u]`**

This sub-block declares a variable named `u`. The name is arbitrary — you could call it
`temperature`, `concentration`, or `phi` — but `u` is conventional for a generic scalar.

The sub-block is empty (between `[u]` and `[]` there are no parameters). This means
all defaults apply:
- **Family: LAGRANGE** — the finite element basis function family. Lagrange basis
  functions are the standard choice: they are nodal (each basis function is 1 at one
  node and 0 at all others), making the degrees of freedom (DOFs) directly interpretable
  as values of `u` at mesh nodes.
- **Order: FIRST** — first-order (linear) Lagrange elements. In 1D these are hat
  functions: triangles centered at each node that are 1 at the node and 0 at its
  neighbors. The approximation of `u` is piecewise linear and continuous.
- **Initial condition: 0** — MOOSE starts the Newton iteration from `u = 0` everywhere.
  For a linear problem this does not matter (Newton converges in one step regardless
  of the initial guess). For nonlinear problems the initial condition can matter greatly.

MOOSE automatically creates one DOF at every mesh node for variable `u`. With 21 nodes,
there are 21 DOFs in this problem. The DOF vector `u_h = [u_0, u_1, ..., u_20]` is what
the solver computes.

---

### Block: `[Kernels]`

```
[Kernels]
  [diffusion]
    type     = Diffusion
    variable = u
  []
[]
```

The `[Kernels]` block defines the volume integrals that make up the weak form of the
PDE. This is the heart of the input file — it tells MOOSE what PDE you are solving.

**Background: The Weak Form**

The finite element method does not solve the PDE directly in its strong form
(`-d²u/dx² = 0`). Instead it multiplies by a test function `phi_i`, integrates over
the domain, and applies integration by parts to reduce the derivative order. The result
is the **weak form**:

```
integral( d(phi_i)/dx * du/dx  dx ) = 0     for all test functions phi_i
```

The left side is what the `Diffusion` kernel computes. The right side is zero because
there is no source term.

Integration by parts is not just a mathematical trick — it has physical meaning. It
transfers a derivative from `u` to the test function, which means both `u` and the test
functions only need to be piecewise differentiable (C0 continuity). This allows the use
of simple polynomial basis functions.

**`type = Diffusion`**

The built-in `Diffusion` kernel implements the weak-form term:

```
integral( grad(phi_i) . grad(u)  dV )
```

In 1D this is:

```
integral( d(phi_i)/dx * du/dx  dx )
```

This term is evaluated element by element using Gauss quadrature. For each element:
1. MOOSE evaluates the shape function gradients at each quadrature point
2. Multiplies them together
3. Weights by the quadrature weight and element Jacobian
4. Accumulates into the element stiffness matrix

The `Diffusion` kernel assumes unit diffusivity (no coefficient). If the problem had
a diffusivity `k` (like thermal conductivity in heat conduction), you would use
`MatDiffusion` instead, as Case 03 demonstrates.

**`variable = u`**

Associates this kernel with the variable `u`. MOOSE uses this association to:
- Know which column of the Jacobian this kernel contributes to
- Know which DOFs are involved in this kernel's residual
- Automatically loop over the right elements and DOFs during assembly

---

### Block: `[BCs]`

```
[BCs]
  [pin_left]
    type     = DirichletBC
    variable = u
    boundary = left
    value    = 0
  []
  [pin_right]
    type     = DirichletBC
    variable = u
    boundary = right
    value    = 1
  []
[]
```

The `[BCs]` block defines how the PDE solution is constrained at the domain boundary.
Without boundary conditions the problem is underdetermined (infinitely many solutions).

**Two `DirichletBC` objects are created — one for each end of the rod.**

**`type = DirichletBC`**

A `DirichletBC` fixes the value of a variable at a boundary. It is the simplest and most
common boundary condition type. Mathematically it modifies the system of equations:
instead of the assembled residual equation at a boundary node, MOOSE replaces it with
the equation `u_node = value`.

`DirichletBC` is called a "strong" BC because it directly constrains the DOF value,
as opposed to "weak" BCs (like Neumann conditions) that are enforced via the residual.

**`[pin_left]`**

- `boundary = left` — applies to the boundary named `left`. In a 1D `GeneratedMesh`,
  `left` refers to the single node at `x = xmin = 0`.
- `value = 0` — sets `u = 0` at this node. This is our boundary condition `u(0) = 0`.

**`[pin_right]`**

- `boundary = right` — applies to the boundary named `right`, the node at `x = xmax = 1`.
- `value = 1` — sets `u = 1` at this node. This is our boundary condition `u(1) = 1`.

**Why are these called "Dirichlet" conditions?**

Johann Peter Gustav Lejeune Dirichlet (1805-1859) was a German mathematician who
studied potential theory. His name is attached to the type of boundary condition that
fixes the value of a field (as opposed to Neumann conditions that fix the flux, or Robin
conditions that mix both).

**What happens without BCs?**

If you removed the `[BCs]` block entirely, the Laplace equation `-d²u/dx²=0` would
have infinitely many solutions (any linear function `u=ax+b` satisfies it). The system
matrix would be singular and the solver would fail. Boundary conditions are not optional
for elliptic PDEs.

---

### Block: `[Executioner]`

```
[Executioner]
  type = Steady

  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]
```

The `[Executioner]` block controls how MOOSE drives the solve — the outer loop strategy.

**`type = Steady`**

A `Steady` executioner solves `F(u) = 0` exactly once, then stops. There is no time
loop, no timestep, no history. This is appropriate because we want the equilibrium
solution, not the transient approach to equilibrium.

For steady problems, MOOSE:
1. Assembles the global residual vector `F` and Jacobian matrix `J`
2. Calls the nonlinear solver to find `u` such that `F(u) = 0`
3. Writes output
4. Exits

**`solve_type = 'PJFNK'`**

PJFNK stands for **Preconditioned Jacobian-Free Newton-Krylov**. Let us unpack this:

- **Newton's method**: an iterative scheme to solve `F(u) = 0`. Starting from an initial
  guess `u_0`, each iteration solves the linear system `J * delta_u = -F(u)` where `J`
  is the Jacobian matrix `dF/du`, then updates `u <- u + delta_u`. This converges
  quadratically near the solution.

- **Jacobian-Free**: instead of explicitly assembling the full Jacobian matrix `J`
  (which can be expensive and complex to implement), JFNK approximates matrix-vector
  products `J*v` using a finite difference:
  ```
  J*v ≈ (F(u + eps*v) - F(u)) / eps
  ```
  This means MOOSE only needs to evaluate the residual function `F`, not its derivative.
  The Jacobian is never formed explicitly.

- **Krylov**: the linear system at each Newton step is solved with a Krylov subspace
  method (typically GMRES — Generalized Minimum Residual). Krylov methods work purely
  with matrix-vector products, which is why JFNK can use the finite-difference
  approximation above.

- **Preconditioned**: the Krylov solver converges much faster if the system is
  preconditioned — multiplied by an approximate inverse of `J`. This is where
  BoomerAMG comes in.

For the linear Laplace equation, Newton converges in **one iteration** because the
residual is linear in `u`. The PJFNK framework is overkill here but is used consistently
across all cases so the pattern is familiar when nonlinear problems appear.

**`petsc_options_iname` and `petsc_options_value`**

These two parameters pass options directly to PETSc, the underlying solver library.
They work as paired lists: each name in `iname` gets the corresponding value from
`petsc_options_value`.

- `-pc_type hypre` — select HYPRE as the preconditioner package
- `-pc_hypre_type boomeramg` — within HYPRE, use the BoomerAMG algorithm

**BoomerAMG** (Bootstrap AMG) is an algebraic multigrid preconditioner. Multigrid
methods are highly effective for elliptic PDEs (like the Laplace equation) because they
eliminate error at all spatial scales simultaneously. "Algebraic" means it works
directly from the matrix entries, without needing to know the mesh geometry explicitly.

For this tiny 21-DOF problem, the preconditioner choice has no visible effect on speed.
It matters enormously for large 3D problems with millions of DOFs.

---

### Block: `[Outputs]`

```
[Outputs]
  exodus = true
  csv    = true
[]
```

The `[Outputs]` block controls what files MOOSE writes after the solve.

**`exodus = true`**

Writes the solution to an **Exodus II** file named `case01_diffusion_1d_out.e`.
Exodus II is a binary format based on HDF5/NetCDF developed at Sandia National
Laboratories specifically for finite element data. It stores:
- The mesh (node coordinates, element connectivity)
- Solution fields at nodes (the variable `u` at all 21 nodes)
- Metadata (time, variable names, etc.)

Exodus files can be opened in ParaView and VisIt for visualization.

The output filename is automatically derived from the input filename by stripping the
`.i` extension and appending `_out.e`.

**`csv = true`**

Writes a comma-separated values file named `case01_diffusion_1d_out.csv`. For this
case with no postprocessors, the CSV file contains minimal content (just a time column
for the steady-state output at t=0). It becomes more useful in Case 03 where
postprocessors provide scalar summary quantities.

---

## What Happens When You Run This

### Invocation

From inside the `case01` directory, assuming your MOOSE test application is at
`../../test/moose_test-opt`:

```bash
cd quickstart-runs/case01
../../test/moose_test-opt -i case01_diffusion_1d.i
```

Or from any directory using an absolute path:

```bash
/path/to/moose/test/moose_test-opt -i /path/to/case01/case01_diffusion_1d.i
```

### Step-by-Step Internally

**1. Startup and Parsing**

MOOSE reads `case01_diffusion_1d.i` using the HIT parser. Every block and parameter is
validated against the registered object list. If you misspell `DirichletBC` or omit a
required parameter, MOOSE reports an error here before doing any computation.

**2. Object Construction**

The Factory pattern instantiates all objects declared in the input file:
- One `GeneratedMesh` object
- One `MooseVariable` object for `u`
- One `Diffusion` kernel object
- Two `DirichletBC` objects
- One `Steady` executioner object
- One `Exodus` output object, one `CSV` output object

**3. Mesh Generation**

`GeneratedMesh` creates the 1D mesh: 21 nodes at x = 0, 0.05, 0.10, ..., 0.95, 1.00,
connected by 20 line elements. The boundaries `left` and `right` are tagged.

**4. DOF Numbering**

The DofMap assigns global equation numbers to every node. With one variable and 21 nodes,
there are 21 equations (DOFs numbered 0 through 20).

**5. Initial Condition**

All DOFs are initialized to 0.0 (the default initial condition for variable `u`).

**6. Nonlinear Solve**

The `Steady` executioner calls the nonlinear solver (Newton-Krylov):

  **Newton iteration 0:**
  - MOOSE loops over all 20 elements, calling the `Diffusion` kernel's `computeQpResidual()`
    at each quadrature point (2 Gauss points per 1D element)
  - Assembles the 21x21 global stiffness matrix `K` and right-hand side `f`
  - `DirichletBC` modifies rows 0 and 20 of the system to enforce `u=0` and `u=1`
  - Computes the initial residual norm `|R|`

  **Linear solve:**
  - GMRES + BoomerAMG solves `K * delta_u = -R` for the correction
  - For this linear problem, GMRES converges in 1-2 iterations

  **Newton update:**
  - `u <- u + delta_u`
  - The residual `|R|` drops to machine precision (around 1e-15)
  - Newton declares convergence (one iteration for a linear problem)

**7. Output**

MOOSE writes `case01_diffusion_1d_out.e` and `case01_diffusion_1d_out.csv`, then exits.

### Reading the Console Output

Typical console output looks like this (line numbers added for reference):

```
Framework Information:
  MOOSE Version:          git commit abc1234
  PETSc Version:          3.19.0
  LibMesh Version:        ...

Mesh Information:
  Spatial dimension:      1
  Mesh:                   20 elements, 21 nodes

Nonlinear System:
  Num DOFs:               21

Executing initial conditions...

 0 Nonlinear |R| = 1.000000e+00
      0 Linear |R| = 1.000000e+00
      1 Linear |R| = 2.204860e-16
 1 Nonlinear |R| = 1.234567e-15

Solve Converged!
```

Line-by-line explanation:
- `Mesh: 20 elements, 21 nodes` — confirms the mesh was built correctly
- `Num DOFs: 21` — 21 unknowns to solve for
- `0 Nonlinear |R| = 1.000000e+00` — the residual norm at Newton iteration 0 (before any
  correction). The value 1.0 reflects the initial residual from the boundary conditions.
- `0 Linear |R| = 1.000000e+00` — the Krylov (GMRES) linear residual at iteration 0
- `1 Linear |R| = 2.2e-16` — after 1-2 Krylov iterations, the linear solve is at
  machine precision
- `1 Nonlinear |R| = 1.2e-15` — after applying the Newton correction, the nonlinear
  residual is at machine precision. Newton has converged.
- `Solve Converged!` — the executioner declares success

**What does convergence mean?** The solver has found `u` such that the residual
`F(u) = K*u - f` has norm below a tolerance (default: `nl_rel_tol = 1e-8` relative,
`nl_abs_tol = 1e-50` absolute). For a linear problem this happens in exactly one Newton
step.

---

## Output Files

### `case01_diffusion_1d_out.e`

This Exodus II binary file is created in the same directory as the input file. It contains:

| Data | Description |
|------|-------------|
| Node coordinates | x positions of all 21 nodes |
| Element connectivity | Which nodes belong to each of the 20 elements |
| Variable `u` | Solution value at every node |
| Time | 0.0 (steady-state is treated as "time 0") |

**How to visualize with ParaView:**

1. Open ParaView
2. File > Open, navigate to `case01_diffusion_1d_out.e`, click OK
3. In the Properties panel, click Apply
4. The solution appears. For a 1D result, use Filters > Data Analysis > Plot Over Line
5. Set Point1 to (0, 0, 0) and Point2 to (1, 0, 0), click Apply
6. The plot shows u vs x — it should be a perfectly straight line from (0,0) to (1,1)

**How to read nodal values programmatically:**

```python
# Using the meshio Python library
import meshio
mesh = meshio.read("case01_diffusion_1d_out.e")
# or use the ExodusII Python API, or read with numpy
```

### `case01_diffusion_1d_out.csv`

A plain-text file with one row per output time step. For this steady-state case with no
postprocessors declared, the CSV has minimal content:

```
time
0
```

The file grows more useful in Case 03 where postprocessors are added.

---

## Understanding the Plots

Running `python visualize_all.py` from the `quickstart-runs/` directory produces the
following plot saved into this directory.

### `case01_diffusion_1d.png`

**What the plot shows.** A line plot with x on the horizontal axis (0 to 1) and u on
the vertical axis (0 to 1). Two curves are drawn: blue circles connected by lines
representing the MOOSE numerical solution at each mesh node, and a dashed red line
representing the exact solution `u = x`.

**Physical quantities.** The horizontal axis is the spatial coordinate along the 1D
rod. The vertical axis is the solution field `u`, which could represent temperature,
concentration, or electrostatic potential depending on context.

**How to judge correctness.** The MOOSE dots should lie exactly on the red dashed
line — not approximately, but exactly. Because first-order Lagrange elements represent
linear functions without any approximation error, the maximum pointwise difference
between the two curves should be below 1e-14 (floating-point machine precision).

**What would indicate a problem.**
- Any visible gap between the blue dots and the red line signals a solver failure or
  a wrong boundary condition setup.
- A curved or non-monotone blue line means the solver did not converge or the mesh is
  degenerate.
- A flat horizontal line at u=0 means the right Dirichlet BC was not applied.
- A flat line at u=1 means the left BC was not applied.

---

## Interpreting the Results

### What the Solution Looks Like

The numerical solution field `u` should be exactly the linear ramp `u(x) = x`. At every
mesh node the computed value equals the node's x-coordinate:

```
Node 0  (x=0.00):  u = 0.00
Node 1  (x=0.05):  u = 0.05
Node 2  (x=0.10):  u = 0.10
...
Node 10 (x=0.50):  u = 0.50
...
Node 20 (x=1.00):  u = 1.00
```

This result is exact (not merely an approximation) because:
1. The true solution `u = x` is a degree-1 polynomial
2. First-order Lagrange elements represent degree-1 polynomials exactly
3. The Gauss quadrature used by MOOSE integrates degree-3 polynomials exactly
4. There is no rounding error in representing a linear function on a uniform mesh

The maximum pointwise error between the MOOSE solution and the exact solution should be
below 1e-14 (floating-point machine precision), not just some small finite number.

### Verifying Correctness

Three ways to confirm the answer is right:

1. **Visual check**: Open in ParaView, plot over line — it should be a perfectly straight
   line from (0,0) to (1,1). Any deviation would be visible.

2. **Nodal value check**: Read the Exodus file and compare node values to their x
   coordinates. Every difference should be below 1e-14.

3. **Convergence is reported**: The solver converges in 1 Newton iteration. If it takes
   more iterations or fails to converge, something is wrong with the setup.

### Physical Insight

The linear temperature profile in a 1D rod with no internal heat generation and fixed
endpoint temperatures tells us:
- Heat flows at a **constant rate** through the rod (the heat flux is `-k * du/dx = -1`,
  uniform everywhere)
- There are **no hot spots or cold spots** in the interior
- The rod acts like a pure resistor: voltage (temperature) drops linearly from one end
  to the other

This is the thermal equivalent of Ohm's law. The "resistance" is the length of the rod
divided by its conductivity. More resistance (longer rod, lower conductivity) means a
steeper temperature gradient for the same fixed temperature difference.

---

## Key Concepts Learned

- **HIT format**: MOOSE input files use block-based Hierarchical Input Text format with
  `[BlockName]` / `[]` delimiters and `parameter = value` assignments
- **GeneratedMesh**: creates structured rectangular meshes from parameters, no external
  file required; assigns named boundaries (`left`, `right` in 1D)
- **Variables block**: declares the unknown fields to be solved; defaults to first-order
  Lagrange elements with zero initial condition
- **Diffusion kernel**: implements the weak form of `-div(grad u)` via the integral
  `integral( grad(phi_i) . grad(u) dV )`; no coefficient means unit diffusivity
- **DirichletBC**: fixes the value of a variable at a named boundary by directly
  replacing the DOF equation with `u = value`
- **Steady executioner**: solves `F(u) = 0` once with no time loop
- **PJFNK**: Preconditioned Jacobian-Free Newton-Krylov — MOOSE's default nonlinear
  solver; converges in one iteration for linear problems
- **BoomerAMG**: algebraic multigrid preconditioner from HYPRE; ideal for elliptic PDEs
- **Exodus output**: binary file format for FEM data, readable by ParaView and VisIt
- **Exact solvability**: linear FEM reproduces linear solutions exactly, providing a
  built-in correctness check for the simplest case

---

## Experiments to Try

### Experiment 1: Change Mesh Resolution

Modify `nx` in the `[Mesh]` block:

```
nx = 5      # very coarse
nx = 100    # fine
nx = 1000   # very fine
```

For this linear problem, all three give the same exact answer. Observe that:
- The run time barely changes (even 1000 DOFs solves instantly in 1D)
- The solution values at x = 0.5 are all exactly 0.5
- Newton always converges in 1 iteration regardless of mesh size

This demonstrates that for problems with exact polynomial solutions, mesh refinement
does not improve accuracy because accuracy is already perfect.

### Experiment 2: Change the Boundary Values

Try different Dirichlet values:

```
# Make the rod go from -1 to 1:
[pin_left]   value = -1
[pin_right]  value =  1
# Exact solution: u(x) = 2x - 1

# Make both ends the same:
[pin_left]   value = 5
[pin_right]  value = 5
# Exact solution: u(x) = 5 everywhere (trivial uniform field)

# Reverse the gradient:
[pin_left]   value = 1
[pin_right]  value = 0
# Exact solution: u(x) = 1 - x
```

Each modification changes the solution predictably. This builds intuition for how
Dirichlet BCs control the solution.

### Experiment 3: Add a CSV Postprocessor to Quantify the Solution

Add a `[Postprocessors]` block to compute the average value of u over the domain:

```
[Postprocessors]
  [avg_u]
    type     = ElementAverageValue
    variable = u
  []
[]
```

For the linear solution `u = x` on [0,1], the spatial average is exactly 0.5. The
postprocessor value will appear in the console output and the CSV file. Try changing
the boundary values and predicting the average before running.

### Experiment 4: Remove a Boundary Condition

Comment out the `[pin_right]` sub-block:

```
[BCs]
  [pin_left]
    type     = DirichletBC
    variable = u
    boundary = left
    value    = 0
  []
  # [pin_right] removed
[]
```

Without a Dirichlet condition on the right, the `right` boundary has a **natural
Neumann condition**: `du/dn = 0`, which means zero flux through the right end. Combined
with `u(0) = 0` and no source, the only solution satisfying both `-d²u/dx²=0` and zero
flux on the right is `u = 0` everywhere (a trivial flat solution). Observe this in
the output.

### Experiment 5: Extend the Domain

Change `xmin` and `xmax`:

```
xmin = -5
xmax =  5
```

With the same boundary conditions (`u=0` on left, `u=1` on right), the exact solution
becomes `u(x) = (x - (-5)) / (5 - (-5)) = (x + 5) / 10`. The solution still ramps
linearly from 0 to 1, but now over a domain of length 10 instead of 1. The slope
`du/dx = 0.1` is shallower. This demonstrates that the Laplace equation is scale-invariant:
the solution shape depends only on boundary values, not on domain size.
