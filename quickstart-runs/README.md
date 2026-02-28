# MOOSE Quickstart Tutorial Cases: Complete Reference

This directory contains 93 self-contained simulation cases for learning MOOSE from zero.
Each case has its own subdirectory with an input file (`.i`) and pre-run output files.
You do not need to build or install anything to study the input files, understand the physics,
and read the results. If you want to run the simulations yourself, see Section 5.

This document is designed so that someone who has never used MOOSE, never written a finite
element simulation, and is not familiar with scientific computing file formats can read it
from top to bottom and understand everything in these 93 cases.

Read every section. Do not skip ahead. The later cases build directly on concepts introduced
in the earlier ones.

---

## Table of Contents

1. [What is MOOSE?](#1-what-is-moose)
2. [Understanding the HIT Input File Format](#2-understanding-the-hit-input-file-format)
3. [Understanding Output Files](#3-understanding-output-files)
4. [How MOOSE Solves Problems](#4-how-moose-solves-problems-conceptual)
5. [Running Simulations](#5-running-simulations)
6. [The 93 Cases at a Glance](#6-the-93-cases-at-a-glance)
7. [Creating Your Own Simulations](#7-creating-your-own-simulations)
8. [Glossary](#8-glossary)

---

## 1. What is MOOSE?

MOOSE stands for **Multiphysics Object-Oriented Simulation Environment**. It is an open-source
software framework developed by Idaho National Laboratory for writing physics simulation programs.

### What does "simulation" mean here?

A simulation answers the question: "How does a physical quantity vary across space and time?"

Examples of physical questions a simulation can answer:
- How hot does the center of a nuclear fuel rod get, given a certain power level?
- How does a chemical concentration spread through a porous rock over 1000 years?
- How does stress distribute through a metal beam when a load is applied?
- Where does fluid pile up when flow enters a pipe with a bend?

In each case, the answer is a **field** — a number (temperature, concentration, stress,
pressure) at every point in a physical domain (the rod, the rock, the beam, the pipe).

### What is a PDE, and why do we need to solve one?

Physical laws governing how quantities spread through space and change over time are written
as **partial differential equations** (PDEs). You have probably seen the heat equation written
in a physics or engineering class:

```
dT/dt = alpha * (d²T/dx² + d²T/dy² + d²T/dz²)
```

This says: the rate of change of temperature in time equals the thermal diffusivity times
the spatial "curvature" of the temperature field (roughly: how much a point differs from its
neighbors). Where temperature is locally higher than its surroundings, it decreases. Where it
is lower, it increases. The system naturally smooths itself out over time.

For a real three-dimensional object with complex geometry and non-constant material
properties, you cannot solve this equation by hand. You need a computer.

### What does MOOSE do?

MOOSE breaks the process into manageable steps:

1. **You describe the physics**: You write an input file stating which PDE to solve, what
   the material properties are, what the boundary conditions are (temperatures, fluxes, etc.
   at the edges of the domain), and what the initial conditions are.

2. **MOOSE discretizes the domain**: The continuous region of space (your fuel rod, rock layer,
   beam) is divided into a large number of small, simple shapes called **elements** — typically
   triangles or quadrilaterals in 2D, tetrahedra or hexahedra in 3D. This grid of shapes is
   called the **mesh**.

3. **MOOSE sets up the equations**: Inside each element, the unknown field (temperature,
   concentration, etc.) is approximated using simple polynomial functions called **shape
   functions**. Substituting these approximations into the PDE and requiring that the equation
   is satisfied in an averaged sense (the **weak form**) turns the PDE into a large system of
   algebraic equations — one equation per unknown at each mesh node.

4. **MOOSE solves the equations**: The algebraic system (often millions of equations) is solved
   using numerical linear algebra libraries, primarily **PETSc** (Portable, Extensible Toolkit
   for Scientific Computation) and **libMesh** (the finite element library underneath MOOSE).

5. **MOOSE writes results**: The solution (the field values at every node in the mesh, at
   every time step) is written to output files that can be visualized with tools like
   **ParaView** or **VisIt**.

### Why is MOOSE useful?

Writing a finite element code from scratch is a research-level software engineering project.
MOOSE provides a framework so that a scientist or engineer who understands their physics
can focus on writing the PDE terms ("kernels") and material models, without reimplementing
mesh handling, linear algebra, parallel computing, and I/O from scratch.

MOOSE also handles:
- **Parallel computing**: automatically distributes work across multiple CPU cores
- **Adaptivity**: automatically refines the mesh where the solution changes rapidly
- **Multiphysics coupling**: allows multiple physics problems to communicate through
  shared variables and data transfers
- **Nonlinear problems**: handles problems where material properties depend on the solution
  itself (like a material that gets stiffer as it heats up)

### What are the 93 cases in this directory?

These cases form a progressive tutorial starting from the simplest possible problem
(1D steady-state diffusion with an exact solution of u = x) and building up to
multi-application coupled solves with adaptive time stepping. By the end you will have
seen every major MOOSE feature used in everyday simulation work.

---

## 2. Understanding the HIT Input File Format

Every MOOSE simulation is controlled by a plain-text **input file** with a `.i` extension.
These files use a format called **HIT** (Hierarchical Input Text). You need to be able to
read and write these files fluently to use MOOSE.

### Basic Syntax Rules

**Blocks**: HIT organizes settings into named blocks. A block opens with `[BlockName]` and
closes with `[]`. Blocks can be nested.

```text
[Mesh]              # opens the Mesh block
  type = GeneratedMesh
  dim  = 2
[]                  # closes the Mesh block
```

**Sub-blocks**: Blocks can contain named sub-blocks. The sub-block name is user-defined
(you choose it). The `type` parameter inside tells MOOSE which C++ class to instantiate.

```text
[Kernels]
  [my_diffusion_term]   # user-chosen name for this kernel
    type     = Diffusion
    variable = u
  []
[]
```

**Key-value pairs**: Inside a block, settings are written as `parameter = value`. No quotes
needed for single-word values. Use single quotes for values with spaces:

```text
boundary = 'left right top bottom'
```

**Comments**: Anything after `#` on a line is a comment. MOOSE ignores it.

```text
nx = 20    # 20 elements in the x direction — this is a comment
```

**Top-level variables**: You can define scalar variables at the top of the file (before any
block) and reference them using `${varname}` syntax. This is useful for avoiding
repetition:

```text
k_value = 1.5   # define once at the top

[Materials]
  [mat]
    type        = GenericConstantMaterial
    prop_names  = 'k'
    prop_values = '${k_value}'   # reference the top-level variable
  []
[]
```

### A Complete Annotated Mini-Example

```text
# -------------------------------------------------------
# Mini-example: 1D steady diffusion
# This comment block explains the whole file.
# -------------------------------------------------------

# Top-level variable: define the domain length once
L = 1.0

[Mesh]
  # GeneratedMesh creates a structured grid without needing an external file.
  type = GeneratedMesh
  dim  = 1          # one-dimensional problem (a line segment)
  nx   = 10         # number of elements along x
  xmin = 0          # left endpoint of the line
  xmax = ${L}       # right endpoint, uses the top-level variable L = 1.0
[]

[Variables]
  # Declare the unknown field to solve for.
  # 'u' is just a name — you can call it anything.
  [u]
    # No extra settings needed for simple cases.
    # Defaults to first-order Lagrange shape functions.
  []
[]

[Kernels]
  # Each kernel implements one term of the PDE.
  [diffusion]
    type     = Diffusion   # implements -div(grad(u))
    variable = u           # which unknown does this kernel act on?
  []
[]

[BCs]
  # Boundary conditions fix the solution at the domain edges.
  [left_end]
    type     = DirichletBC   # fix the value (not the flux)
    variable = u
    boundary = left          # named boundary on the generated mesh
    value    = 0             # u = 0 at x = 0
  []
  [right_end]
    type     = DirichletBC
    variable = u
    boundary = right         # x = xmax
    value    = 1             # u = 1 at x = 1
  []
[]

[Executioner]
  # How to march toward the solution.
  type = Steady   # no time stepping: just find F(u) = 0

  # The nonlinear solver strategy.
  solve_type = 'PJFNK'

  # Preconditioner options passed to PETSc.
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]

[Outputs]
  exodus = true   # write case_out.e for ParaView
  csv    = true   # write case_out.csv for postprocessors
[]
```

### All Standard Block Types Explained

The following sections explain every block type you will encounter in the 93 cases.
Each explanation defines what the block does, what parameters mean, and gives a
realistic example.

---

#### `[Mesh]` — The Computational Grid

The mesh block defines the geometry of the physical domain and how it is divided into
elements. Think of the mesh as a grid of tiny cells that tiles your physical region.
The solution is computed at the corners (nodes) and interior quadrature points of these cells.

MOOSE supports two main approaches: building a mesh automatically from parameters
(`GeneratedMesh`, `GeneratedMeshGenerator`) or reading a mesh from an external file
created by a meshing tool like Cubit or Gmsh.

**For `GeneratedMesh`** (simple rectangular/box domains):

```text
[Mesh]
  type = GeneratedMesh
  dim  = 2        # dimensionality: 1, 2, or 3
  nx   = 20       # elements in x direction
  ny   = 20       # elements in y direction
  xmin = 0.0      # x coordinate of the left face
  xmax = 1.0      # x coordinate of the right face
  ymin = 0.0      # y coordinate of the bottom face
  ymax = 1.0      # y coordinate of the top face
[]
```

The `GeneratedMesh` automatically names its boundaries:
- In 1D: `left` (x = xmin) and `right` (x = xmax)
- In 2D: `left`, `right`, `bottom` (y = ymin), `top` (y = ymax)
- In 3D: `left`, `right`, `bottom`, `top`, `back` (z = zmin), `front` (z = zmax)

**For mesh generators** (more complex geometry, multiple subdomains):

```text
[Mesh]
  [gen]
    type = GeneratedMeshGenerator   # same as above, but generator syntax
    dim  = 2
    nx   = 40
    ny   = 20
  []
  [split_right]
    type        = SubdomainBoundingBoxGenerator
    input       = gen        # which mesh generator to start from
    bottom_left = '0.5 0 0'  # lower-left corner of the bounding box
    top_right   = '1.0 1 0'  # upper-right corner
    block_id    = 2          # elements inside the box get this subdomain ID
  []
[]
```

The mesh generator system processes generators in dependency order. Later generators
can reference earlier ones by name (the `input` parameter). This lets you build complex
meshes by chaining simple operations.

**Key parameters summary:**

| Parameter    | Meaning                                          |
|--------------|--------------------------------------------------|
| `type`       | Which mesh generator or mesh type to use         |
| `dim`        | Spatial dimension: 1, 2, or 3                    |
| `nx/ny/nz`   | Number of elements in each direction             |
| `xmin/xmax`  | Domain extents in x (similarly ymin/ymax, zmin/zmax) |
| `block_id`   | Integer ID for a mesh subdomain (used in [Materials]) |

---

#### `[Variables]` — The Unknowns Being Solved

The `[Variables]` block declares the primary unknown fields — the quantities the solver
is computing at each mesh node. Examples: temperature `T`, concentration `c`, displacement
components `ux`, `uy`, `uz`.

Each sub-block declares one variable. The sub-block name is the variable name:

```text
[Variables]
  [T]
    # MOOSE defaults: first-order Lagrange, initial value 0
  []
  [c]
    order  = SECOND    # use quadratic (higher-accuracy) shape functions
    family = LAGRANGE  # the mathematical family of the shape functions
    initial_condition = 0.5   # start with c = 0.5 everywhere
  []
[]
```

**Understanding order and family:**
- `order = FIRST` (default): piecewise linear approximation. Values are stored at mesh
  nodes (element corners). Simple, fast, sufficient for most problems.
- `order = SECOND`: piecewise quadratic. Adds midside nodes. More accurate for the same
  mesh but more memory and work.
- `family = LAGRANGE` (default): nodal shape functions. The solution is continuous
  (no jumps across element boundaries). Appropriate for temperature, concentration, etc.

For the 21 tutorial cases, you will always see `[Variables]` with simple empty sub-blocks
(using all defaults) because first-order Lagrange is standard.

---

#### `[Kernels]` — The PDE Terms

Kernels are the heart of MOOSE. Each kernel object implements one term of the governing
partial differential equation. If your PDE has four terms, you put four kernels in this block.

To understand what a kernel does, you need to understand the **weak form** briefly:

**The weak form in plain English**: Instead of requiring the PDE to hold exactly at every
single point in space (which is hard to enforce with simple polynomials), the finite element
method requires the PDE to hold in an *averaged* sense. You multiply the PDE by a test
function (one for each mesh node), integrate over the whole domain, and require the
integral (called the **residual**) to be zero. Integration by parts converts second-order
terms into first-order terms, which polynomials can represent accurately.

The result is that each term in the PDE maps to an integral that a kernel computes.

**Common kernels:**

`Diffusion` — implements the operator -div(grad u) with no coefficient:
```text
[diffusion]
  type     = Diffusion
  variable = u
[]
```
This contributes the integral: `integral( grad(phi_i) * grad(u) dV )`
where `phi_i` is the test function for node `i`.

`MatDiffusion` — same as Diffusion but multiplied by a material property (conductivity):
```text
[heat_conduction]
  type        = MatDiffusion
  variable    = T
  diffusivity = k    # name of a material property declared in [Materials]
[]
```
This contributes: `integral( k * grad(phi_i) * grad(T) dV )`

`TimeDerivative` — implements the time rate of change dT/dt:
```text
[time_deriv]
  type     = TimeDerivative
  variable = T
[]
```
This contributes: `integral( phi_i * dT/dt dV )`
Required for any transient (time-evolving) problem.

`BodyForce` — implements a volumetric source term Q:
```text
[heat_source]
  type     = BodyForce
  variable = T
  value    = 1.0        # constant source term (can also use 'function =')
[]
```
This contributes: `integral( phi_i * Q dV )` (adds Q to the right-hand side)

`CoupledForce` — adds a term proportional to another variable (coupling between equations):
```text
[coupling_v_into_u]
  type     = CoupledForce
  variable = u
  v        = v      # the variable that drives the source
  coef     = 0.1    # coefficient; default 1.0, can be negative
[]
```
This contributes: `integral( phi_i * coef * v dV )` to the equation for `u`.

`ConservativeAdvection` — implements advection (transport by a velocity field):
```text
[advection]
  type              = ConservativeAdvection
  variable          = c
  velocity_material = velocity_vec   # vector material property
[]
```

`ADMatDiffusion` — the Automatic Differentiation version of `MatDiffusion`. Used when
the material property itself depends on the solution (nonlinear problem). The `AD` prefix
means MOOSE automatically computes the exact derivative (Jacobian) using dual-number
arithmetic rather than finite differences.

---

#### `[BCs]` — Boundary Conditions

The `[BCs]` block specifies what happens at the boundaries of the domain. Without boundary
conditions, the PDE has infinitely many solutions (think of a ramp u = x + C for any C).
Boundary conditions make the solution unique.

**Dirichlet boundary conditions** fix the value of the unknown at a boundary:

```text
[BCs]
  [pin_left]
    type     = DirichletBC
    variable = T
    boundary = left      # boundary name (from mesh)
    value    = 100.0     # T = 100 at the left boundary
  []
  [pin_right]
    type     = DirichletBC
    variable = T
    boundary = right
    value    = 0.0       # T = 0 at the right boundary
  []
[]
```

**Neumann boundary conditions** fix the flux (rate of flow) at a boundary. A zero-flux
Neumann condition (meaning nothing flows through the boundary) is the **natural** or
**default** boundary condition in MOOSE. You get it automatically for any boundary you
do not list in `[BCs]`. You never need to write it explicitly.

A non-zero flux Neumann BC is written with `NeumannBC`:
```text
[heat_flux_in]
  type     = NeumannBC
  variable = T
  boundary = left
  value    = 500.0   # 500 W/m^2 flowing into the domain
[]
```

**FunctionDirichletBC** applies a Dirichlet condition from a mathematical expression
(from the `[Functions]` block), useful when the boundary value varies in space or time:
```text
[time_varying_left]
  type     = FunctionDirichletBC
  variable = T
  boundary = left
  function = my_function_name   # evaluated at each node's coordinates and time
[]
```

**When NO boundary condition is listed for a boundary**: MOOSE applies zero-flux
(homogeneous Neumann). This means: no heat flows across that boundary, no concentration
diffuses through it. This is physically appropriate for an insulated wall, a symmetry plane,
or an outflow boundary.

---

#### `[ICs]` — Initial Conditions

Initial conditions are required for transient (time-evolving) problems. They set the
value of the unknown fields at time `t = 0`.

```text
[ICs]
  [T_start]
    type     = ConstantIC
    variable = T
    value    = 20.0    # uniform initial temperature of 20 degrees
  []

  [concentration_blob]
    type     = FunctionIC
    variable = c
    function = blob_fn    # initial condition from a [Functions] entry
  []
[]
```

For `FunctionIC`, the function is evaluated at each mesh node's coordinates at t=0.
A common use is to start with a Gaussian distribution:
```text
[Functions]
  [blob_fn]
    type       = ParsedFunction
    expression = 'exp(-((x-0.5)^2 + (y-0.5)^2) / 0.01)'
  []
[]
```

If you do not specify an IC for a variable, MOOSE assumes `ConstantIC` with `value = 0`.
For the heat equation starting from a uniform cold state, this default is often what you want.

---

#### `[Materials]` — Material Properties

The `[Materials]` block defines the physical properties of the material(s) in the domain.
These are numbers like thermal conductivity, density, specific heat, diffusivity, and
viscosity. Kernels and BCs look up these values at each quadrature point during assembly.

Material properties can be:
- Constant numbers
- Spatially varying functions
- Temperature-dependent (or dependent on any solution variable)

**GenericConstantMaterial** — declare one or more scalar properties with fixed values:
```text
[thermal_props]
  type        = GenericConstantMaterial
  prop_names  = 'k    rho  cp'
  prop_values = '50.0 7800 500'   # k=50 W/(mK), rho=7800 kg/m^3, cp=500 J/(kgK)
[]
```
The `prop_names` and `prop_values` lists must have the same number of entries.

**GenericFunctionMaterial** — material property from a `[Functions]` entry:
```text
[conductivity]
  type        = GenericFunctionMaterial
  prop_names  = 'k'
  prop_values = 'k_function'   # name of a Function in [Functions]
[]
```

**Block restriction**: When you have different materials in different mesh regions (subdomains),
add a `block` parameter:
```text
[mat_left]
  type        = GenericConstantMaterial
  block       = 0       # only applies to elements in subdomain 0
  prop_names  = 'k'
  prop_values = '1.0'
[]
[mat_right]
  type        = GenericConstantMaterial
  block       = 2       # only applies to elements in subdomain 2
  prop_names  = 'k'
  prop_values = '5.0'
[]
```

**ADPiecewiseLinearInterpolationMaterial** — material property that depends on a
solution variable (nonlinear), using automatic differentiation for the Jacobian:
```text
[conductivity]
  type     = ADPiecewiseLinearInterpolationMaterial
  property = k
  variable = T          # k depends on temperature T
  xy_data  = '0  1.0
              200 1.5
              400 2.0'  # piecewise linear: k(0)=1, k(200)=1.5, k(400)=2
[]
```

---

#### `[Functions]` — Mathematical Expressions

Functions define mathematical expressions that can be evaluated at any point in space
and time. They are referenced by name from `[BCs]`, `[ICs]`, `[Kernels]`, `[Materials]`,
and `[Postprocessors]`.

**ParsedFunction** — the most common type, evaluates a mathematical string:
```text
[my_func]
  type       = ParsedFunction
  expression = 'sin(pi*x) * cos(2*pi*t) + 1.0'
[]
```

Built-in variables in the expression: `x`, `y`, `z`, `t`
Built-in constants: `pi` (3.14159...), `e` (2.71828...)
Supported operators: `+`, `-`, `*`, `/`, `^` (power), unary `-`
Supported functions: `sin`, `cos`, `tan`, `exp`, `log`, `sqrt`, `abs`, `if(cond,a,b)`,
and many more (the expression is parsed by the `fparser` library).

**PiecewiseLinear** — define a function by a table of (x, y) data points:
```text
[ramp_up]
  type = PiecewiseLinear
  x = '0.0  0.5  1.0'   # time or x-coordinate values
  y = '0.0  500  500'   # corresponding function values
[]
```

For the 21 tutorial cases, `ParsedFunction` is used to define exact solutions for
verification (Method of Manufactured Solutions), initial concentration distributions,
and spatially varying material properties.

---

#### `[AuxVariables]` — Additional Fields (Not Solved For)

Auxiliary variables hold field data that is not part of the primary solve. They can be:
- Computed from the primary solution (strain from displacement, flux from temperature)
- Received from another MOOSE application in a MultiApp simulation
- Used to visualize intermediate quantities

```text
[AuxVariables]
  [heat_flux_x]
    order  = CONSTANT
    family = MONOMIAL    # piecewise constant, one value per element
  []
  [T_from_parent]
    order  = FIRST
    family = LAGRANGE    # nodal, same as a primary variable
  []
[]
```

Auxiliary variables do not have their own equations (no kernels act on them). They are
updated by `[AuxKernels]` objects or populated by `[Transfers]` in MultiApp setups.

---

#### `[Postprocessors]` — Scalar Diagnostic Quantities

Postprocessors compute a single number from the solution field at each output step
(each time step for transient problems, once for steady). They are written to the CSV
output file and printed to the console.

Common postprocessors:

**ElementAverageValue** — spatial average of a variable over the domain (or a block):
```text
[avg_temp]
  type     = ElementAverageValue
  variable = T
  # block = 0   # optional: restrict to one subdomain
[]
```

**ElementExtremeValue** — maximum or minimum of a variable:
```text
[max_temp]
  type       = ElementExtremeValue
  variable   = T
  value_type = max   # or 'min'
[]
```

**ElementIntegralVariablePostprocessor** — integral of a variable over the domain:
```text
[total_energy]
  type     = ElementIntegralVariablePostprocessor
  variable = T
[]
```
Computes `integral( T dV )`. For a unit square domain, this equals the average value.

**ElementL2Error** — L2 norm of the difference between solution and an exact function:
```text
[L2_error]
  type     = ElementL2Error
  variable = u
  function = exact_fn   # reference to a [Functions] entry
[]
```
Computes `sqrt( integral( (u - u_exact)^2 dV ) )`. Used for verification.

**ElementL2Norm** — L2 norm of the solution itself:
```text
[T_norm]
  type     = ElementL2Norm
  variable = T
[]
```
Computes `sqrt( integral( T^2 dV ) )`.

**TimestepSize** — reports the current time step size:
```text
[dt_pp]
  type = TimestepSize
[]
```
Useful for monitoring adaptive time stepping.

---

#### `[Executioner]` — The Solver Strategy

The executioner controls how MOOSE marches toward the solution: steady-state vs.
time-stepping, nonlinear solver settings, and convergence tolerances.

**Steady executioner** — find the solution where F(u) = 0 (no time dependence):
```text
[Executioner]
  type = Steady

  solve_type = 'PJFNK'

  # PETSc preconditioner options:
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]
```

**Transient executioner** — march forward in time from initial conditions:
```text
[Executioner]
  type = Transient

  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'

  start_time = 0.0
  end_time   = 10.0

  nl_rel_tol = 1e-8   # nonlinear relative convergence tolerance
  nl_abs_tol = 1e-10  # nonlinear absolute convergence tolerance
  nl_max_its = 20     # maximum Newton iterations per step

  [TimeStepper]
    type = ConstantDT
    dt   = 0.1
  []
[]
```

**`solve_type` options:**
- `PJFNK` (Preconditioned Jacobian-Free Newton-Krylov): the default. Approximates the
  Jacobian using finite differences of the residual. Works for most problems.
- `NEWTON`: uses the exact (or AD-computed) Jacobian. Faster convergence per iteration
  but requires explicitly coded or AD-generated Jacobians. Used in Case 7.
- `JFNK`: like PJFNK but without a preconditioner. Slower, rarely used.
- `LINEAR`: for problems that are truly linear (no nonlinear terms). Skips Newton iteration.

**Convergence tolerances:**
- `nl_rel_tol`: converged when the nonlinear residual drops to this fraction of the initial
  residual. Default 1e-8.
- `nl_abs_tol`: converged when the nonlinear residual is below this absolute value.
  Default 1e-50 (effectively disabled; relative tolerance dominates).
- `l_tol`: tolerance for the inner linear solver at each Newton step. Default 1e-5.

---

#### `[TimeStepper]` — Controlling the Time Step Size

The `[TimeStepper]` sub-block (inside `[Executioner]`) controls how the time step size
changes over the course of a simulation.

**ConstantDT** — fixed step size:
```text
[TimeStepper]
  type = ConstantDT
  dt   = 0.01
[]
```

**IterationAdaptiveDT** — adjusts the step size based on how hard the previous step was
(measured by the number of Newton iterations needed):
```text
[TimeStepper]
  type               = IterationAdaptiveDT
  dt                 = 0.001         # initial step size
  optimal_iterations = 5             # target Newton iteration count
  growth_factor      = 2.0           # multiply dt by this if iterations < optimal
  cutback_factor     = 0.5           # multiply dt by this if iterations > optimal
  iteration_window   = 2             # acceptable range around optimal
[]
```

Logic: if the previous step took fewer than `optimal_iterations` Newton iterations, the
problem is "easy" so MOOSE multiplies dt by `growth_factor`. If it took many iterations,
the problem is "hard" so dt is multiplied by `cutback_factor`. This automatically uses
large steps during slow transients and small steps during fast transients. Case 11 is
entirely dedicated to demonstrating this feature.

---

#### `[Outputs]` — What Files to Write

```text
[Outputs]
  exodus = true    # write the full solution field to an Exodus II file (.e)
  csv    = true    # write postprocessor values to a CSV file
[]
```

**Exodus output** creates a file named `<input_file_base>_out.e`. This binary file contains
the mesh geometry and the solution fields at every time step. Open it with ParaView.

**CSV output** creates a file named `<input_file_base>_out.csv`. This is a plain-text table
with one row per time step and one column per postprocessor, plus a `time` column.

Optional output control:
```text
[Outputs]
  [my_exodus]
    type     = Exodus
    interval = 5    # write only every 5th time step (reduces file size)
  []
  [my_csv]
    type = CSV
  []
[]
```

---

#### `[Adaptivity]` — Automatic Mesh Refinement

For problems where the solution varies greatly from smooth regions to sharp features,
automatic mesh refinement (AMR) places small elements where they are needed and
large elements where the solution is smooth. This uses computational resources
efficiently instead of applying a uniformly fine mesh everywhere.

```text
[Adaptivity]
  marker         = err_marker   # which Marker decides what to refine
  initial_steps  = 4            # refinement cycles before the first solve
  initial_marker = err_marker
  steps          = 1            # additional refinement steps after the solve
  max_h_level    = 5            # maximum number of refinements relative to base mesh

  [Indicators]
    # Indicators estimate the local error on each element.
    [jump_indicator]
      type     = GradientJumpIndicator
      variable = u
    []
  []

  [Markers]
    # Markers use indicator values to decide which elements to refine/coarsen.
    [err_marker]
      type      = ErrorFractionMarker
      indicator = jump_indicator
      refine    = 0.5    # refine the worst 50% of elements by error
      coarsen   = 0.05   # coarsen the best 5% of elements
    []
  []
[]
```

**Indicators** estimate how large the error is on each element. The
`GradientJumpIndicator` measures how sharply the gradient of the solution jumps
at element faces — large jumps indicate that the element is too coarse to resolve
the solution accurately.

**Markers** use indicator values to classify each element as "refine", "coarsen",
or "do nothing". The `ErrorFractionMarker` sorts elements by their estimated error
and refines the worst-performing fraction.

---

#### `[MultiApps]` — Running Multiple Coupled Simulations

The MultiApp system allows one MOOSE simulation (the parent) to launch and control
one or more sub-simulations (sub-apps). This is used for:
- Operator splitting (solve one physics first, then another)
- Scale-bridging (microscale model feeding a macroscale model)
- Domain decomposition across non-overlapping meshes

```text
[MultiApps]
  [my_sub_app]
    type        = FullSolveMultiApp   # run sub-app to convergence each time
    input_files = sub_app.i          # path to the sub-app's input file
    execute_on  = timestep_end       # when to run: 'initial', 'timestep_end', etc.
    positions   = '0 0 0'           # where in the parent's coordinate system the sub-app sits
  []
[]
```

`FullSolveMultiApp` runs the sub-application to complete convergence at the specified
execution point. The sub-application is a completely independent MOOSE solve with its
own mesh, variables, kernels, and BCs. It just receives and sends data via `[Transfers]`.

---

#### `[Transfers]` — Moving Data Between Applications

Transfers move data between the parent and sub-applications. Without transfers, the
two solves are completely decoupled.

```text
[Transfers]
  [send_temperature]
    type            = MultiAppCopyTransfer
    to_multi_app    = my_sub_app       # direction: parent -> sub
    source_variable = T                # variable in the parent
    variable        = T_received       # AuxVariable in the sub-app
  []
[]
```

`MultiAppCopyTransfer` copies nodal values from one application to another.
It requires the meshes to have the same topology (same nodes at the same positions).
For more complex cases where meshes differ, `MultiAppInterpolationTransfer` interpolates
between meshes.

The data transferred can go either direction:
- `to_multi_app = name` — parent sends to sub-app
- `from_multi_app = name` — sub-app sends to parent (bidirectional coupling)

---

## 3. Understanding Output Files

### Exodus II Files (.e)

**What is Exodus II?** Exodus II is a binary file format created by Sandia National
Laboratories for storing finite element analysis results. It is the standard output format
for MOOSE and most other FEA codes. The binary format is compact and efficient for large
meshes with many time steps.

**What's stored in an Exodus file?**
- The complete mesh: node coordinates (x, y, z for every mesh node), element connectivity
  (which nodes form each element), boundary names and their associated nodes/faces
- Subdomain (block) assignments for elements
- Solution fields: for every variable (`T`, `u`, `c`, etc.), the value at every mesh node
  at every stored time step
- Auxiliary variable fields
- Time step values

**How to view Exodus files: ParaView**

ParaView is a free, open-source visualization application. It runs on Windows, macOS, and
Linux. It is the standard tool for viewing MOOSE results.

Download: https://www.paraview.org/download/

Step-by-step for viewing a MOOSE result:

1. Start ParaView
2. Go to File > Open
3. Navigate to the `.e` file (e.g., `case01_diffusion_1d_out.e`)
4. Click OK
5. In the Pipeline Browser on the left, the file will appear. Click the "Apply" button
   in the Properties panel below the pipeline browser.
6. The mesh will appear in the 3D viewport, probably showing as a solid gray color.
7. To color it by the solution field: in the toolbar near the top of the 3D viewport,
   find a dropdown that says "Solid Color" or "vtkBlockColors". Click it and change it
   to the variable name (`u`, `T`, or whatever your variable is called).
8. The mesh will now show a color-coded representation of the field. A color bar appears
   on the side showing the mapping from color to value.
9. Click the "Rescale to Data Range" button (looks like a rescale symbol) to make the
   color bar span the actual min and max values of your solution.

**Navigating time steps (for transient results):**

Below the toolbar is a time control bar with a time slider and play/pause buttons.
Click "Play" to animate the time series and watch the solution evolve. You can also
type a specific time value in the box to jump to that step.

**Using "Plot Over Line" to extract 1D data:**

For 1D problems or 1D cuts through 2D/3D results:
1. Select the data in the Pipeline Browser
2. Go to Filters > Data Analysis > Plot Over Line
3. Set the start point and end point of the line
4. Click Apply
5. A 2D plot appears showing the field value along the line

**Using "Plot Data Over Time" (formerly "Plot Selection Over Time"):**

To see how a value at a specific point changes over time:
1. Select the data
2. Filters > Data Analysis > Plot Data Over Time
3. Use the Selection Inspector to pick a point
4. Click Apply for a time series plot

**Alternative: VisIt**

VisIt is another free visualization tool from Lawrence Livermore National Laboratory.
Download: https://visit-dav.github.io/visit-website/
It supports Exodus II and most other FEA formats.

**Programmatic access with Python:**

You can read Exodus files in Python using the `netCDF4` or `exodus` library (Exodus II
is built on top of netCDF4):

```python
import netCDF4 as nc

ds = nc.Dataset('case03_heat_transient_out.e')

# Print all available variables
print(list(ds.variables.keys()))

# Read x-coordinates of all nodes
x_coords = ds.variables['coordx'][:]

# Read the time values stored in the file
times = ds.variables['time_whole'][:]

# Read the solution field 'T' (named 'vals_nod_var1' if T is the first variable)
# The naming convention varies; inspect the variable names first.
T_field = ds.variables['vals_nod_var1'][:]  # shape: (n_timesteps, n_nodes)
```

The `mooseutils` Python library (in `python/mooseutils/`) provides higher-level access
to MOOSE output files if you have the MOOSE Python environment available.

---

### CSV Files (.csv)

CSV (Comma-Separated Values) is a plain text format. Every modern spreadsheet and data
analysis tool can read it.

**Format**: The first line is a header with column names. Subsequent lines contain
numerical values, separated by commas. The `time` column always comes first.

Example CSV from Case 3 (transient heat):
```
time,avg_temperature,max_temperature
0.01,0.000247,0.000393
0.02,0.000464,0.000742
0.05,0.001134,0.001813
0.10,0.002198,0.003515
0.20,0.004052,0.006481
0.50,0.073300,0.117130
```

**Opening in a spreadsheet**: Open the file directly in Microsoft Excel, LibreOffice Calc,
or Google Sheets. Select the columns you want and create a chart.

**Reading with Python:**

```python
import csv

# Simple approach using the built-in csv module
with open('case03_heat_transient_out.csv', newline='') as f:
    reader = csv.DictReader(f)
    rows = list(reader)

time        = [float(r['time'])            for r in rows]
avg_temp    = [float(r['avg_temperature']) for r in rows]
max_temp    = [float(r['max_temperature']) for r in rows]

print(f"Final average temperature: {avg_temp[-1]:.5f}")
print(f"Final maximum temperature: {max_temp[-1]:.5f}")
```

Or using `pandas` (a more powerful data library):
```python
import pandas as pd

df = pd.read_csv('case03_heat_transient_out.csv')
print(df.describe())   # summary statistics
print(df.tail())       # last few rows

# Plot with matplotlib
import matplotlib.pyplot as plt
plt.plot(df['time'], df['avg_temperature'], label='Average T')
plt.plot(df['time'], df['max_temperature'], label='Maximum T')
plt.xlabel('Time (s)')
plt.ylabel('Temperature')
plt.legend()
plt.savefig('temperature_history.png')
```

---

## 4. How MOOSE Solves Problems (Conceptual)

This section explains the computational process in plain language. You do not need to
follow every mathematical detail to use MOOSE, but understanding the big picture helps
you interpret console output and diagnose problems.

### The Mesh and Degrees of Freedom

The continuous physical domain (a line, square, cube, or complex geometry) is divided
into a mesh of small elements. In 2D, these are typically quadrilaterals (four-sided)
or triangles. In 3D, they are hexahedra (bricks) or tetrahedra.

At each mesh **node** (a corner of an element), MOOSE allocates one or more unknowns
called **degrees of freedom** (DOFs). For a temperature problem with 100 nodes and one
variable `T`, there are 100 DOFs. The entire solution is described by specifying all 100
temperature values simultaneously.

Within each element, the temperature is approximated using **shape functions** — smooth
polynomial functions that are 1 at one node and 0 at all other nodes. The temperature
at any point inside an element is a weighted sum of the nodal temperatures. For first-order
(bilinear) elements, this makes the temperature field piecewise linear — a flat plane
on each element, connected smoothly across elements.

### The Residual

The **residual** is the measure of how poorly the current solution satisfies the governing
PDE. If the solution were perfect, the residual would be exactly zero at every point.

For a steady-state problem `F(u) = 0`:
- The residual is the vector `R = F(u_current)` evaluated at the current guess for `u`.
- If `R = 0`, the solution has converged.
- If `R` is large, the current guess is far from the true solution.

MOOSE assembles the residual by looping over all mesh elements, calling each kernel to
compute its contribution to the element residual, and assembling these contributions into
a global vector of length equal to the number of DOFs.

### Newton's Method

Newton's method finds the solution to `F(u) = 0` by iteratively improving a guess:

1. Start with an initial guess `u_0` (usually all zeros, or the previous time step's solution)
2. Compute the residual `R = F(u_current)`
3. Check if `||R||` is small enough (convergence check). If yes, stop.
4. Compute the **Jacobian** matrix `J = dF/du` (how the residual changes as u changes)
5. Solve the linear system: `J * delta_u = -R`
6. Update: `u_new = u_current + delta_u`
7. Go to step 2 with `u_new`

For a linear problem (like pure diffusion with constant conductivity), Newton converges
in exactly **one iteration** because the Jacobian is the same as the stiffness matrix
and the correction `delta_u` is exact.

For a nonlinear problem (like diffusion with temperature-dependent conductivity), multiple
iterations are needed. The residual should decrease rapidly — ideally **quadratically** (each
iteration squares the error), which is the hallmark of Newton's method working correctly.

### The Linear Solve

Each Newton step requires solving the linear system `J * delta_u = -R`. For a mesh with
100,000 DOFs, this is a 100,000 x 100,000 system of equations.

Direct solution methods (Gaussian elimination) are impractical for systems this large.
MOOSE uses iterative solvers from **PETSc**:

- **GMRES** (Generalized Minimum Residual method): an iterative Krylov subspace method
  that finds increasingly good approximations to the solution. Requires a preconditioner.

- **Preconditioner**: a transformation applied to the system to make GMRES converge faster.
  The `boomeramg` preconditioner (algebraic multigrid) works by solving a sequence of
  coarser and coarser versions of the problem, using the coarse solutions to accelerate
  convergence on the fine mesh. It is very effective for elliptic problems (diffusion,
  elasticity, Poisson).

The linear solve is considered converged when the linear residual `||J * x - (-R)||`
falls below a tolerance (typically `l_tol = 1e-5`).

### The Transient Loop

For time-dependent problems, MOOSE adds an outer time-marching loop around the Newton
solve:

```
for each time step:
    guess: u_new = u_old (extrapolate from previous step)
    Newton solve to find u_new satisfying F(u_new, u_old, dt) = 0
    if converged:
        u_old = u_new
        advance time by dt
        (possibly adjust dt based on iteration count)
    if not converged:
        reduce dt and retry
```

MOOSE uses implicit time integration by default (**Backward Euler**, a first-order method).
This is numerically stable even for large time steps, unlike explicit methods. The tradeoff
is that each step requires solving a nonlinear system, but stability is not restricted by
any time step size limit.

### Interpreting Console Output

When you run a MOOSE simulation, the console prints a progress log. Here is a
line-by-line interpretation:

```
Framework Information:
  MOOSE Version:          git commit abc1234
```
The MOOSE framework version being used.

```
Mesh Information:
  Spatial dimension:      2
  Mesh:                   400 elements, 441 nodes
```
Basic mesh statistics. 441 nodes for a 20x20 quad mesh = 21x21 nodes.

```
Nonlinear System:
  Num DOFs:               441
```
Total number of unknowns to solve for (one temperature at each node).

```
 0 Nonlinear |R| = 1.000000e+00
```
Before the first Newton iteration, the residual norm is 1.0 (normalized to 1 by convention
for the initial step).

```
      0 Linear |R| = 1.000000e+00
      1 Linear |R| = 4.234567e-06
      ...
      9 Linear |R| = 2.345678e-13
```
The inner linear solver (GMRES) iterations, showing convergence of the linear solve at
this Newton step. `Linear |R|` is the norm of the residual of the linear system
`J * delta_u = -R`.

```
 1 Nonlinear |R| = 2.345678e-13
```
After one Newton step, the nonlinear residual dropped to machine precision. For a linear
problem, this is expected: one step is exact.

```
Converged in 1 its.
```
The nonlinear solve converged. "its." means "iterations."

For a nonlinear problem, you would see multiple nonlinear iterations, each with its own
set of linear iterations. The nonlinear residual should decrease rapidly (quadratically
if the Jacobian is exact, linearly or super-linearly for PJFNK with approximate Jacobian).

**Signs of trouble:**
- Nonlinear residual not decreasing: the Jacobian is wrong or the problem is ill-conditioned
- Linear solve not converging: bad preconditioner choice for this problem type
- "Nonlinear solve did not converge": need smaller time step, better initial guess,
  or looser tolerances

---

## 5. Running Simulations

### Prerequisites

To run the simulations, you need a compiled MOOSE test executable. These cases all use
the `moose_test-opt` executable from the MOOSE test application.

**Building from source (on a Linux/macOS system with MOOSE installed):**

```bash
cd /path/to/moose/test
make -j8
# Creates: test/moose_test-opt
```

**Using Docker (no local build required):**

Docker lets you run the MOOSE executable inside a pre-built container without compiling
anything yourself. Docker Desktop can be downloaded from https://www.docker.com/

```bash
# Pull the official MOOSE Docker image
docker pull idaholab/moose:latest

# Run case01 using Docker
docker run --rm \
  -v /full/path/to/quickstart-runs/case01:/work \
  -w /work \
  idaholab/moose:latest \
  /opt/moose/bin/moose_test-opt -i case01_diffusion_1d.i
```

Replace `/full/path/to/quickstart-runs/case01` with the actual absolute path to the
case directory on your machine.

### Running a Case

The basic command pattern is:
```bash
/path/to/moose_test-opt -i input_file.i
```

Run from the directory containing the input file so that output files are created there:

```bash
cd quickstart-runs/case01
/path/to/moose_test-opt -i case01_diffusion_1d.i
```

For Case 12 (MultiApp), run only the parent:
```bash
cd quickstart-runs/case12
/path/to/moose_test-opt -i case12_parent.i
```
Both `case12_parent.i` and `case12_sub.i` must be in the same directory.

### Useful Command-Line Options

| Option                    | Meaning                                               |
|---------------------------|-------------------------------------------------------|
| `-i input.i`              | Specify input file (required)                         |
| `--mesh-only`             | Generate and output the mesh, then exit (no solve)    |
| `-o output_prefix`        | Override the output file base name                    |
| `--n-threads=4`           | Use 4 CPU threads (shared memory parallelism)         |
| `--distributed-mesh`      | Distribute the mesh across MPI processes              |
| `--show-input`            | Print the parsed input file and exit                  |
| `--dump`                  | Dump all available object types and their parameters  |
| `-h` / `--help`           | Show all command-line options                         |
| `--no-color`              | Disable colored console output                        |

For MPI (distributed memory parallel runs):
```bash
mpirun -n 4 /path/to/moose_test-opt -i input.i
```
This distributes the mesh and computation across 4 MPI processes.

### Output File Naming

By default, output files are named `<input_file_base>_out.<extension>`.
For `case01_diffusion_1d.i`:
- Exodus: `case01_diffusion_1d_out.e`
- CSV: `case01_diffusion_1d_out.csv`

To override: add `file_base = my_custom_name` to the `[Outputs]` block.

---

## 6. The 93 Cases at a Glance

The cases are ordered from simplest to most complex. Cases 01-13 use only the MOOSE
framework (Diffusion, BodyForce, MatDiffusion, etc.). Cases 14-21 use physics **modules**
(heat_transfer, solid_mechanics, navier_stokes, phase_field, porous_flow) and demonstrate
genuine multi-physics coupling. Cases 22-29 cover continuum electromechanics — charge
transport, magnetic diffusion, electrohydrodynamic flows, and MHD — inspired by Melcher's
*Continuum Electromechanics* (MIT, 1981). Cases 30-36 cover waveguides, resonators, wave
scattering, coupled modes, noise, and solitons — inspired by **Professor Herman A. Haus**'s
*Electromagnetic Noise and Quantum Optical Measurements* (Springer, 2000), drawing from the
classical chapters (Chs 1-5, 10) only. Cases 37-44 cover advanced fluid dynamics —
instabilities, boundary layers, turbulence, compressible flow, rotating fluids, and MHD
waves — inspired by **Michel Rieutord**'s *Fluid Dynamics: An Introduction* (Springer, 2015),
spanning Chapters 4-10. Cases 45-48 cover stochastic tools and optimization — Monte Carlo UQ,
polynomial chaos surrogates, adjoint-based optimization, and Latin Hypercube parameter studies.

**Cases 49-73 (Batches A-E)** were generated autonomously by [Claude Code](https://claude.ai/claude-code)
using the `moose-simulation` skill. The entire 25-case expansion — covering nonlinear solid
mechanics, nuclear reactor physics, geomechanics, reactive transport, and advanced multiphysics —
was launched with a single user prompt:

> *"is it possible for you to complete each of Batch A, B, C, D, E using moose-simulation skill.
> In other words, apply moose-simulation skill to finish Batch A (push to GitHub), then apply
> moose-simulation skill to finish Batch B (push to GitHub), ...., and do the same to Batch E ?"*

Claude Code then executed the full 9-step workflow for each batch autonomously: researching
MOOSE module APIs, authoring input files, running simulations in Docker, validating convergence,
generating matplotlib plots, writing README documentation, and pushing to GitHub. All 25 cases
converge in Docker with `combined-opt` in under 2 minutes each.

Cases 49-53 cover nonlinear solid mechanics — J2 plasticity,
finite-strain kinematics, power-law creep, phase-field fracture, and pressure-vessel
verification. Cases 54-58 cover nuclear reactor physics — neutron diffusion eigenvalue
problems, fuel-pin heat transfer, xenon poisoning transients, and control-rod worth
calculations. Cases 59-63 cover geomechanics and porous flow — Terzaghi consolidation,
wellbore drawdown, unsaturated Richards' equation, Biot poroelasticity, and gravity dam
structural analysis. Cases 64-68 cover chemical reactions and transport — reaction-diffusion
with analytical verification, contaminant advection-diffusion-reaction driven by Darcy flow,
mineral precipitation via Arrhenius kinetics, aqueous equilibrium speciation with pH tracking,
and calcite dissolution combining equilibrium and kinetic reactions. Cases 69-73 cover advanced
multiphysics — MultiApp Picard coupling, mortar contact mechanics, XFEM heat conduction with
cracks, THM pipe flow using the Component DSL, and level-set bubble advection with SUPG
stabilization.

**Cases 74-83 (Batch F — Advanced Electromagnetism)** cover topics from MIT 6.635 *Advanced
Electromagnetism* (Spring 2003). Cases 74-75 demonstrate left-handed materials and Drude-model
skin depth in 1D. Case 76 is a full 3D rectangular waveguide TE10 mode using NEDELEC_ONE
(edge) elements on HEX20 mesh. Case 77 solves 2D EM scattering from a dielectric cylinder.
Case 78 is a time-domain Gaussian pulse reflection using NewmarkBeta integration. Case 79
demonstrates Snell's law and the scattered-field formulation for oblique incidence. Case 80
models a 5-period Bragg mirror achieving |R|^2 = 1 at the design frequency. Case 81 computes
photonic crystal eigenfrequencies at the Gamma point of a square lattice of alumina rods.
Case 82 extends the 2D waveguide eigenvalue problem to a 3D PEC rectangular cavity. Case 83
demonstrates Veselago flat-lens focusing of a point source through a negative-index slab.

**Cases 84-93 (Batch G — Receivers, Antennas & Sensing)** are based on MIT 6.661
*Electromagnetics and Applications* (Prof. D. H. Staelin). Case 84 solves a lossy
TEM cavity eigenvalue problem to extract the quality factor Q. Cases 85-86 model
Hertzian and half-wave dipole radiation patterns in 2D axisymmetric (RZ) coordinates.
Case 87 demonstrates phased array beam steering. Case 88 computes single-slit
Fraunhofer diffraction. Case 89 finds guided TE modes in a dielectric slab waveguide.
Case 90 models parabolic reflector focusing. Case 91 computes the radar cross section
of a flat conducting plate. Case 92 creates interference fringes from two coherent
sources (aperture synthesis). Case 93 uses MOOSE's optimization module to recover an
unknown source amplitude from sparse measurements using the adjoint method.

Study them in order.

| Case | Subdirectory | Title | Physics | Key Concepts Introduced | Difficulty |
|------|--------------|-------|---------|-------------------------|------------|
| 01 | `case01-1d-steady-diffusion/` | Steady-State 1D Diffusion | Laplace equation in 1D: -d²u/dx² = 0, exact solution u=x | Mesh, Variables, Kernels (Diffusion), BCs (DirichletBC), Executioner (Steady), Outputs | Beginner |
| 02 | `case02-2d-steady-diffusion/` | Steady-State 2D Diffusion | Laplace equation on unit square, same exact solution u=x (no y variation) | GeneratedMesh in 2D, boundary names (left/right/top/bottom), natural (zero-flux) Neumann BCs | Beginner |
| 03 | `case03-transient-heat/` | Transient Heat Equation | Heat equation with source: rho*cp*dT/dt = div(k*grad T) + Q | TimeDerivative, MatDiffusion, BodyForce, Materials, Postprocessors, Transient executioner, ConstantDT | Beginner |
| 04 | `case04-manufactured-solution/` | Poisson with Manufactured Solution | Poisson equation with known exact solution for code verification | Functions (ParsedFunction), FunctionDirichletBC, ElementL2Error, Method of Manufactured Solutions | Beginner |
| 05 | `case05-neumann-bc/` | Spatially Varying Conductivity | -div((1+x)*grad u) = 0, conductivity varies with x | GenericFunctionMaterial, spatially varying material properties, non-trivial exact solution | Beginner |
| 06 | `case06-two-material-domain/` | Two-Material Domain | Composite domain: left half k=1, right half k=5 | SubdomainBoundingBoxGenerator, block-restricted materials, flux continuity at interfaces | Intermediate |
| 07 | `case07-nonlinear-diffusion/` | Nonlinear Diffusion | k(T) = 1 + T, solution-dependent conductivity | ADMatDiffusion, ADPiecewiseLinearInterpolationMaterial, NEWTON solver, quadratic convergence, AD | Intermediate |
| 08 | `case08-advection-diffusion/` | Advection-Diffusion | dc/dt + v*grad(c) - D*div(grad c) = 0, moving Gaussian blob | ConservativeAdvection, ICs (FunctionIC), GenericConstantVectorMaterial, top-level HIT variables | Intermediate |
| 09 | `case09-coupled-reaction-diffusion/` | Coupled Two-Variable System | Coupled PDEs: du/dt = Du*Lap(u)+v, dv/dt = Dv*Lap(v)-u | Multiple variables, CoupledForce, IterationAdaptiveDT, off-diagonal coupling | Intermediate |
| 10 | `case10-adaptive-mesh-refinement/` | Adaptive Mesh Refinement | Laplace with boundary step discontinuity, locally refined mesh | Adaptivity block, GradientJumpIndicator, ErrorFractionMarker, non-uniform meshes | Intermediate |
| 11 | `case11-adaptive-timestepping/` | Adaptive Time Stepping | Transient heat with step-change source, dt adapts to difficulty | IterationAdaptiveDT (standalone case study), TimestepSize postprocessor | Intermediate |
| 12 | `case12-multiapp-coupling/` | MultiApp Coupling | Parent thermal solve + sub-app that uses parent temperature | MultiApps, FullSolveMultiApp, Transfers, MultiAppCopyTransfer, AuxVariables | Advanced |
| 13 | `case13-custom-kernel/` | Full Postprocessor Analysis | Transient heat with five postprocessors, Python plotting | Multiple postprocessors (ElementAverageValue, ElementExtremeValue, ElementIntegral, ElementL2Norm, TimestepSize), CSV analysis, Python CSV reader, matplotlib | Advanced |
| 14 | `case14-thermoelasticity/` | Thermoelasticity — Heated Plate | Heat conduction + thermal stress: T→eigenstrain→displacement | ADHeatConduction, Physics/SolidMechanics/QuasiStatic action, ADComputeThermalExpansionEigenstrain, vonmises_stress, one-way coupling | Advanced |
| 15 | `case15-lid-driven-cavity/` | Lid-Driven Cavity (Re=100) | Incompressible Navier-Stokes FV, classic CFD benchmark | Modules/NavierStokesFV action, finite volume, ADGenericFunctorMaterial, pressure pinning, moving wall | Advanced |
| 16 | `case16-natural-convection/` | Natural Convection (Ra=10⁴) | Buoyancy-driven flow with Boussinesq approximation, two-way fluid-thermal coupling | NavierStokesFV + energy equation, boussinesq_approximation, Rayleigh/Prandtl number | Advanced |
| 17 | `case17-joule-heating/` | Joule Heating | Electric potential (Laplace) + heat equation with Joule source Q=σ\|∇V\|² | ADJouleHeatingSource, ElectromagneticHeatingMaterial, electro-thermal coupling | Advanced |
| 18 | `case18-cahn-hilliard/` | Cahn-Hilliard Spinodal Decomposition | Phase separation via split Cahn-Hilliard equation, periodic BCs | SplitCHParsed, SplitCHWRes, CoupledTimeDerivative, DerivativeParsedMaterial, phase_field module | Advanced |
| 19 | `case19-porous-flow/` | Darcy Flow + Heat in Porous Media | Single-phase saturated thermo-hydro flow, thermal plume advection | PorousFlowBasicTHM action, SimpleFluidProperties, PorousFlowPermeabilityConst, variable scaling | Advanced |
| 20 | `case20-elastic-wave/` | Elastic Wave Propagation | Dynamic solid mechanics — stress wave in a bar with Newmark-beta | Physics/SolidMechanics/Dynamic action, NewmarkBeta, InertialForce, Pressure BC, wave reflection | Advanced |
| 21 | `case21-bimetallic-strip/` | Bimetallic Strip Bending | Two metals with different thermal expansion heated uniformly → bending | Multi-material (block-restricted), ComputeThermalExpansionEigenstrain, SubdomainBoundingBoxGenerator | Advanced |
| 22 | `case22-charge-relaxation/` | Charge Relaxation in Ohmic Medium | Free charge decays as exp(-t/tau), tau=eps/sigma; Poisson for potential | ADReaction (linear decay), ADCoupledForce (Poisson source), charge relaxation time | Intermediate |
| 23 | `case23-magnetic-diffusion/` | Magnetic Diffusion into Conductor | Step-applied B-field penetrates diffusively: erfc profile, skin depth | ADMatDiffusion reused for magnetic diffusion (D_m = 1/(mu0*sigma)), analytical verification | Intermediate |
| 24 | `case24-drift-diffusion/` | Charge Drift-Diffusion | Unipolar ion injection: drift under E-field + diffusion, self-consistent Poisson coupling | ConservativeAdvection with velocity_as_variable_gradient, upwinding, SMP preconditioning | Advanced |
| 25 | `case25-induction-heating/` | Induction Heating | Oscillating B → eddy currents → Joule heating at skin depth | VariableGradientComponent AuxKernel, ParsedAux for Q=D_m*(dB/dx)^2, CoupledForce source | Advanced |
| 26 | `case26-ehd-pumping/` | EHD Pumping — Coulomb Force | Prescribed Coulomb body force drives recirculating cavity flow | INSFVBodyForce, ParsedFunctorMaterial, body-force-driven FV Navier-Stokes | Advanced |
| 27 | `case27-hartmann-flow/` | MHD Hartmann Flow | Lorentz drag flattens channel flow profile; Darcy friction as MHD analog | porous_medium_treatment, Darcy friction = Lorentz drag, Hartmann profile verification | Advanced |
| 28 | `case28-twoway-joule-heating/` | Two-Way Joule Heating | Extends Case 17: sigma(T) decreases with T (metallic negative feedback) | ADPiecewiseLinearInterpolationMaterial, tabulated T-dependent conductivity, two-way coupling | Advanced |
| 29 | `case29-electroconvection/` | Electroconvection — EHD-Enhanced Convection | Natural convection + EHD body force via effective thermal expansion | Boussinesq + EHD combined as alpha_eff, parameter study (Fe controls enhancement/suppression) | Advanced |
| 30 | `case30-waveguide-cutoff/` | Rectangular Waveguide Cutoff Frequencies | Helmholtz eigenvalue ∇²ψ + k_c²ψ = 0 on 2:1 rectangle | Eigenvalue executioner, CoefReaction(eigen tag), EigenDirichletBC, PotentialToFieldAux | Intermediate |
| 31 | `case31-driven-cavity/` | Driven Resonant Cavity | Helmholtz ∇²E + k²E = −J, near resonance → large field | CoefReaction(negative coeff), BodyForce as source, resonance vs off-resonance | Intermediate |
| 32 | `case32-dielectric-slab/` | EM Wave Reflection from Dielectric Slab | 1D Helmholtz with εᵣ(x), real/imag field splitting | EMRobinBC port condition, ADMatReaction, ReflectionCoefficient, ADGenericFunctionMaterial | Advanced |
| 33 | `case33-coupled-resonators/` | Coupled Resonator Beating | Two coupled modes: ∂u/∂t = D∇²u − γu + κv (mirrored) | ADReaction (decay), CoupledForce (coupling), beating and energy exchange | Intermediate |
| 34 | `case34-thermal-noise/` | Thermal Noise Relaxation | Diffusion from random IC: ∂T/∂t = D∇²T, spectral mode decay | RandomIC, mode-dependent decay rates, fluctuation-dissipation | Beginner |
| 35 | `case35-dispersive-pulse/` | Dispersive Pulse Broadening | Advection-diffusion as GVD fiber: ∂A/∂t + v_g∂A/∂x = D∂²A/∂x² | ADConservativeAdvection, ADMatDiffusion, pulse broadening physics | Intermediate |
| 36 | `case36-soliton-pulse/` | Soliton Pulse Propagation | Nonlinear pulse: ∂A/∂t + v_g∂A/∂x = D∂²A/∂x² − αA³ | ADMatReaction + DerivativeParsedMaterial, soliton balance, JIT workaround | Advanced |
| 37 | `case37-rayleigh-benard/` | Rayleigh-Benard Convection Onset | Boussinesq NS heated from below, Ra=2000 > Ra_c=1708 | NavierStokesFV + Boussinesq, transient onset, convective rolls | Advanced |
| 38 | `case38-kelvin-helmholtz/` | Kelvin-Helmholtz Instability | Shear layer rollup, tanh velocity profile, passive scalar | NavierStokesFV + energy as passive scalar, inlet/outlet BCs, upwind | Advanced |
| 39 | `case39-blasius-boundary-layer/` | Blasius Boundary Layer | Flat plate laminar flow, Re_L=400, Blasius similarity | NavierStokesFV steady, biased mesh (bias_y), PointValue sampling | Advanced |
| 40 | `case40-turbulent-channel/` | Turbulent Channel k-epsilon | RANS k-epsilon with wall functions, Re=13700, log-law validation | SIMPLE segregated, MooseLinearVariableFVReal, k-epsilon, wall functions | Expert |
| 41 | `case41-rayleigh-taylor/` | Rayleigh-Taylor Instability | Heavy-over-light Boussinesq, mushroom fingers, RT growth | NavierStokesFV + Boussinesq, tanh interface IC, sinusoidal perturbation | Advanced |
| 42 | `case42-sod-shock-tube/` | Sod Shock Tube — 1D Riemann | Compressible Euler, shock/contact/rarefaction at t=0.2 | CNS HLLC flux splitting, ExplicitSSPRungeKutta, IdealGasFluidProperties | Expert |
| 43 | `case43-ekman-spiral/` | Ekman Spiral — Rotating BL | Coriolis-coupled vx/vy in rotating frame, exact solution | CoupledForce, BodyForce, ADMatDiffusion, analytical Ekman spiral | Intermediate |
| 44 | `case44-alfven-wave/` | Alfven Wave — MHD Elsasser | Elsasser decomposition: d+ rightward, d- zero, diffusive decay | ADConservativeAdvection, ADMatDiffusion, opposite velocities | Advanced |
| 45 | `case45-monte-carlo-uq/` | Monte Carlo UQ | 2D heat conduction with uncertain k, 30 MC samples | `SamplerFullSolveMultiApp`, `MonteCarlo` sampler | Advanced |
| 46 | `case46-polynomial-chaos/` | Polynomial Chaos | 1D diffusion-reaction surrogate, 36 training + 100 eval | `PolynomialChaosTrainer`, `EvaluateSurrogate` | Advanced |
| 47 | `case47-heat-source-inversion/` | Heat Source Inversion | Adjoint-based recovery of unknown heat source | `Optimize`, `GeneralOptimization`, adjoint method | Expert |
| 48 | `case48-parameter-study/` | Parameter Study | 3-parameter LHS with ParameterStudy action | `[ParameterStudy]` action, `LatinHypercube` | Advanced |
| 49 | `case49-j2-plasticity/` | J2 Plasticity — Uniaxial Tension | Elastic-plastic uniaxial tension with isotropic hardening; bilinear stress-strain | `ComputeMultipleInelasticStress`, `IsotropicPlasticityStressUpdate`, J2 return mapping | Advanced |
| 50 | `case50-finite-strain-compression/` | Finite Strain — Large Deformation Compression | Multiplicative F = Fe Fp kinematics; 50% compressive strain, correct Cauchy stress | `ComputeFiniteStrain`, `ComputeFiniteStrainElasticStress`, `RankTwoScalarAux` | Advanced |
| 51 | `case51-power-law-creep/` | Power-Law Creep — Column Under Sustained Load | Steady-state power-law creep rate = A*sigma^n*exp(-Q/RT); implicit time integration | `PowerLawCreepStressUpdate`, `ComputeMultipleInelasticStress`, `IterationAdaptiveDT` | Advanced |
| 52 | `case52-phase-field-fracture/` | Phase-Field Fracture — Notched Specimen | Coupled damage-mechanics: (1-d)^2 stiffness degradation, crack nucleation and propagation | `ADPFFracture`, `ComputeLinearElasticPFFractureStress`, `PhaseFieldFractureMechanicsOffDiag` | Expert |
| 53 | `case53-pressure-vessel/` | Pressure Vessel — Thick-Walled Cylinder | Lame analytic solution: sigma_r and sigma_theta verified to machine precision | `Physics/SolidMechanics/QuasiStatic` (axisymmetric), `ADComputeSmallStrain`, `ElementL2Error` | Intermediate |
| 54 | `case54-neutron-diffusion-bare-slab/` | 1-Group Neutron Diffusion — Bare Slab Criticality | One-group diffusion eigenvalue: -D*div(grad phi) + Sigma_a*phi = (nu*Sigma_f/k_eff)*phi; cosine flux shape; k_eff = 1.0000 | `Eigenvalue` executioner, `CoefReaction` (eigen tag), `EigenDirichletBC`, `Diffusion` | Intermediate |
| 55 | `case55-two-group-diffusion/` | 2-Group Neutron Diffusion — Fast/Thermal with Fission | Two coupled flux variables (fast + thermal); down-scattering Sigma_12 couples groups; k_eff = 1.342 | `Eigenvalue` executioner, two-variable `CoupledForce` (down-scatter and fission source), `EigenDirichletBC` | Advanced |
| 56 | `case56-fuel-pin-heat-transfer/` | Fuel Pin Heat Transfer — Radial Temperature Profile | Steady axisymmetric (RZ) heat conduction in fuel + cladding with volumetric source, contact resistance, and convective outer BC; T_center ~ 1200 °C | `GeneratedMesh` (coord_type = RZ), `HeatConduction`, `HeatSource`, `ConvectiveHeatFluxBC`, two subdomains | Intermediate |
| 57 | `case57-xenon-poisoning/` | Xenon-135 Poisoning Transient — Reactor Flux Decay | Coupled I-135/Xe-135/flux ODEs over 24 h; xenon peak at ~7 h suppresses flux by ~50% before recovery | `TimeDerivative` (×3), `CoupledForce` (iodine→xenon decay and xenon absorption), `IterationAdaptiveDT`, `ParsedMaterial` | Advanced |
| 58 | `case58-control-rod-worth/` | Control Rod Worth — Eigenvalue Shift from Absorber | Spatially varying Sigma_a in rod subdomain; Delta_k = k_unrodded - k_rodded quantifies rod worth; flux depression near rod tip | `SubdomainBoundingBoxGenerator` (rod region), block-restricted `GenericConstantMaterial`, `Eigenvalue` executioner, `ElementIntegral` | Advanced |
| 59 | `case59-terzaghi-consolidation/` | Terzaghi 1D Consolidation | Saturated soil column loaded at top; excess pore pressure dissipates upward; pore-pressure profile verified against Fourier-series analytic solution | `PorousFlowBasicTHM` action, `SimpleFluidProperties`, `PorousFlowPermeabilityConst`, `SideAverageValue`, `PointValue` (settlement) | Intermediate |
| 60 | `case60-wellbore-drawdown/` | Wellbore Drawdown — Theis Solution in RZ | Radial flow toward pumping well in RZ geometry; drawdown profile verified against Theis well function W(u) at t = 1, 10, 100 hours | `GeneratedMesh` (coord_type = RZ), `PorousFlowBasicTHM` (hydro only), `NeumannBC` (pump flux), `IterationAdaptiveDT`, `PointValue` | Intermediate |
| 61 | `case61-unsaturated-flow/` | Unsaturated Flow — Richards' Equation with van Genuchten | 1D infiltration into initially dry soil; nonlinear wetting front via van Genuchten retention and Mualem relative permeability | `PorousFlowUnsaturated` action, `PorousFlowCapillaryPressureVG`, `PorousFlowRelativePermeabilityVG`, `IterationAdaptiveDT` | Advanced |
| 62 | `case62-biot-poroelasticity/` | Biot Poroelasticity — Explicit Coupling | Full Biot THM using explicit PorousFlow + SolidMechanics kernels; fluid-solid coupling terms visible in input; 2D plane-strain block verified against 1D Terzaghi | `Physics/SolidMechanics/QuasiStatic`, `PorousFlowFullySaturated`, `PorousFlowEffectiveFluidPressure`, `PorousFlowVolumetricStrain`, `ElementL2Error` | Advanced |
| 63 | `case63-gravity-dam/` | Gravity Dam — Multi-Material Structural Analysis | Plain concrete dam + rock foundation under self-weight and hydrostatic load; von Mises stress, principal stress, and heel tensile stress design check | `SubdomainBoundingBoxGenerator`, block-restricted `ADComputeIsotropicElasticityTensor`, `Gravity`, `Pressure` BC (hydrostatic), `ADRankTwoScalarAux`, `NodalExtremeValue` | Advanced |
| 64 | `case64-reaction-diffusion/` | Reaction-Diffusion — First-Order Decay + Diffusion | dc/dt = D*Lap(c) - lambda*c, Gaussian initial condition, analytical solution c(x,t)=c0*exp(-lambda*t)*G(x,t); L2 error verification | `TimeDerivative`, `MatDiffusion`, `Reaction` (first-order decay), `FunctionIC` (Gaussian blob), `ElementL2Error` | Intermediate |
| 65 | `case65-contaminant-transport/` | Contaminant Transport — ADR with Darcy Convection | Advection-diffusion-reaction in porous medium: phi*dc/dt + div(v*c) - div(D*grad c) + lambda*c = 0; Darcy velocity from upstream pressure gradient | `chemical_reactions` module `PrimaryConvection`, `PrimaryDiffusion`, `PrimaryTimeDerivative`, `PrimaryDecay`; `DarcyVelocity` AuxKernel | Advanced |
| 66 | `case66-mineral-precipitation/` | Mineral Precipitation — A+B→Mineral Arrhenius Kinetics | Bimolecular A + B → mineral at rate r = A_freq*exp(-Ea/RT)*[A][B]; coupled PDEs for [A], [B], and solid volume fraction; temperature-dependent rate | `SolidKineticReactions` action, `TimeDerivative`, `Diffusion`, `CoupledForce`, `ParsedMaterial` (Arrhenius rate) | Advanced |
| 67 | `case67-aqueous-equilibrium/` | Aqueous Equilibrium — CO2-H2O Speciation and pH | Batch equilibrium: CO2 + H2O ⇌ H2CO3 ⇌ HCO3⁻ + H⁺ ⇌ CO3²⁻ + H⁺; pH computed from [H⁺] via PHAux; charge balance enforced | `geochemistry` module `GeochemicalModelInterrogator`, `PHAux`, `GeochemistryTimeDependentReactor`; equilibrium constants from LLNL database | Advanced |
| 68 | `case68-calcite-dissolution/` | Calcite Dissolution — Combined Equilibrium + Kinetic Reactions | CaCO3 dissolves: equilibrium speciation + TST kinetic rate r = k*(1 - Q/K); Ca²⁺ and CO3²⁻ concentrations evolve; saturation index tracked as postprocessor | `GeochemistryTimeDependentReactor`, `GeochemistrySpatialReactor`, `KineticRate` object, `SaturationIndex` postprocessor | Expert |
| 69 | `case69-multiapp-coupling/` | MultiApp Coupled Diffusion — Picard Iteration | Bidirectional coupling: main (heat + source) ↔ sub (source + feedback); Picard fixed-point iteration converges in 2 iterations per step | `TransientMultiApp`, `MultiAppNearestNodeTransfer`, `CoupledForce`, fixed-point iteration, `fixed_point_max_its` | Advanced |
| 70 | `case70-contact-mechanics/` | Mortar Frictionless Contact | Two elastic blocks pushed together; mortar-based frictionless contact prevents interpenetration; von Mises stress at interface | `MeshCollectionGenerator`, `SubdomainIDGenerator`, `RenameBlockGenerator`, `[Contact]` action (mortar formulation), block-restricted postprocessors | Advanced |
| 71 | `case71-xfem-crack/` | XFEM Heat Conduction — Stationary Crack | Transient heat conduction with vertical edge crack as perfect insulator; XFEM enrichment without remeshing | `[XFEM]` block, `LineSegmentCutUserObject`, `qrule = volfrac`, enriched elements, no Constraints block needed | Advanced |
| 72 | `case72-thm-pipe-flow/` | THM 1D Pipe Flow — Compressible Gas | 1D compressible ideal gas flow using THM Component DSL; temperature front propagation and pressure equilibration | `[Components]` DSL, `FlowChannel1Phase`, `InletMassFlowRateTemperature1Phase`, `Outlet1Phase`, `scaling_factor_1phase`, `IdealGasFluidProperties` | Advanced |
| 73 | `case73-level-set-advection/` | Level Set Bubble Advection — SUPG | Smooth bubble in solid-body rotation; SUPG-stabilized advection; mass conservation monitoring | `LevelSetAdvection`, `LevelSetAdvectionSUPG`, `LevelSetTimeDerivativeSUPG`, `LAGRANGE_VEC`, `VectorFunctionIC`, BDF2 | Advanced |
| 74 | `case74-left-handed-material/` | Left-Handed Material — Reversed Phase Velocity | 1D plane wave through vacuum/LHM/vacuum slab; phase reversal inside n = -1.5 material | `ADMatDiffusion` (1/mu_r weighted), `ADMatReaction`, `EMRobinBC` (port + absorbing), piecewise `ADGenericFunctionMaterial` | Advanced |
| 75 | `case75-drude-slab/` | Lossy Drude Metal Slab — Skin Depth | 1D wave through Drude metal; coupled E_real/E_imag; exponential attenuation and skin depth | `ADMatCoupledForce` (Im(eps_r) cross-coupling), `ReflectionCoefficient`, complex Helmholtz splitting | Advanced |
| 76 | `case76-3d-waveguide/` | 3D Rectangular Waveguide — TE10 Mode | Full 3D TE10 mode propagation in PEC waveguide; NEDELEC_ONE edge elements on HEX20 | `CurlCurlField`, `VectorFunctionReaction`, `VectorEMRobinBC`, `VectorCurlPenaltyDirichletBC`, NEDELEC_ONE/HEX20 | Expert |
| 77 | `case77-cylinder-scattering/` | EM Scattering from a Dielectric Cylinder | 2D TE plane wave scatters from eps_r = 4 cylinder; diffraction pattern, shadow, forward lobe | Circular inhomogeneity via `if()` in ParsedFunction, `EMRobinBC` port injection, scattered-field | Advanced |
| 78 | `case78-pulse-reflection/` | Time-Domain Pulse Reflection from a Slab | Gaussian EM pulse hits dielectric slab; splits into reflected + transmitted pulses | `VectorSecondTimeDerivative`, `VectorTransientAbsorbingBC`, `VectorFunctionIC`, NewmarkBeta, NEDELEC_ONE/QUAD9 | Expert |
| 79 | `case79-snell-law-tir/` | Oblique Incidence — Snell's Law | 2D oblique plane wave at dielectric interface; scattered-field formulation; Fresnel verification | Scattered-field decomposition, oblique `BodyForce` source, `ParsedAux` with `use_xyzt = true` | Advanced |
| 80 | `case80-bragg-mirror/` | Multilayer Dielectric Stack — Bragg Mirror | 5-period quarter-wave stack; |R|^2 = 1 at design frequency; stopband demonstration | Nested-if `ParsedFunction` for 10 layers, `ReflectionCoefficient`, transfer-matrix verification | Advanced |
| 81 | `case81-photonic-crystal/` | 2D Photonic Crystal Band Gap — Gamma Point Eigenvalue | TM eigenfrequencies of square lattice alumina rods; eps_r-weighted mass matrix | `MatReaction` (non-AD) with `extra_vector_tags = 'eigen'`, `GenericFunctionMaterial`, Neumann (PMC) BCs | Expert |
| 82 | `case82-3d-cavity-resonance/` | 3D Rectangular Cavity Resonance — TM Eigenvalues | PEC cavity eigenvalue problem on HEX8 mesh; analytical k^2(m,n,p) verification | 3D `GeneratedMesh` (dim=3), `EigenDirichletBC` + `DirichletBC`, `PotentialToFieldAux` (Ex, Ey, Ez), KRYLOVSCHUR | Expert |
| 83 | `case83-veselago-lens/` | Veselago Flat Lens — Point Source Focusing | 2D point source through n = -1 LHM slab; three-peak focusing pattern at source, slab centre, image | `ADMatDiffusion` (negative diffusivity in slab), `BodyForce` (Gaussian source), loss regularisation, MUMPS LU | Expert |
| 84 | `case84-lossy-tem-cavity/` | Lossy TEM Cavity — Q Factor | 1D coupled eigenvalue for complex ε_r: coupled E_real/E_imag; resonant modes and quality factor Q = Re(ω)/2·Im(ω) | Eigenvalue solver (KRYLOVSCHUR), `ADMatReaction`, `ADMatCoupledForce`, shift-invert targeting | Expert |
| 85 | `case85-hertzian-dipole/` | Hertzian Dipole — Radiation Pattern | 2D RZ Helmholtz with Gaussian point source at origin; Robin absorbing BC; far-field sin²θ directivity | `coord_type = RZ`, `DiracKernel`-like BodyForce, `EMRobinBC` (absorbing), axisymmetric radiation | Advanced |
| 86 | `case86-half-wave-dipole/` | Half-Wave Dipole — Gain Pattern | 2D RZ Helmholtz with distributed current J(z) = I₀·cos(πz/L); gain = 1.64 (2.15 dBi) | `coord_type = RZ`, `BodyForce` with cosine profile, far-field gain extraction | Advanced |
| 87 | `case87-phased-array/` | Phased Array — Beam Steering | 2D Cartesian Helmholtz with two phased Gaussian sources; beam steering by phase offset δ | Multiple `BodyForce` sources with complex phase, `EMRobinBC`, interference pattern | Advanced |
| 88 | `case88-single-slit-diffraction/` | Single-Slit Diffraction — Resolution Limit | 2D Helmholtz: plane wave through slit in conducting screen; Fraunhofer sinc² pattern | `FunctionDirichletBC` (screen/slit), `EMRobinBC`, first null at θ = λ/a | Advanced |
| 89 | `case89-dielectric-waveguide/` | Dielectric Waveguide — Guided TE Modes | 1D eigenvalue for TE modes in symmetric slab: core n₁=1.5, cladding n₂=1.0; V-number mode cutoff | Eigenvalue solver, piecewise `GenericFunctionMaterial`, total internal reflection | Advanced |
| 90 | `case90-parabolic-reflector/` | Parabolic Reflector — Focal Spot | 2D Helmholtz: plane wave reflects off parabolic PEC surface; focal intensity gain ~6× | Lossy material inside reflector (PEC approximation), `EMRobinBC`, curved boundary | Advanced |
| 91 | `case91-radar-cross-section/` | Radar Cross Section — Flat Plate | 2D scattered-field: plane wave incident on conducting plate; σ ≈ 4πL²/λ² at normal incidence | Scattered-field formulation, `EMRobinBC` (port + absorbing), RCS postprocessor | Advanced |
| 92 | `case92-interferometer/` | Aperture Synthesis Interferometer | 2D Helmholtz with two coherent point sources; fringe pattern shows angular resolution ≈ λ/B | Multiple `BodyForce` sources, correlation extraction, interferometric fringes | Advanced |
| 93 | `case93-inverse-source-recovery/` | Inverse Source Recovery — PDE-Constrained Estimation | 3-file optimization: recover source amplitude p₁ in screened Poisson from 9 sensors; TAO L-BFGS converges in 2 iterations | `Optimize` executioner, `GeneralOptimization`, `OptimizationData`, `ParsedOptimizationFunction`, `ElementOptimizationSourceFunctionInnerProduct`, adjoint method | Expert |

### What each case produces

| Case | Output Files |
|------|-------------|
| 01 | `case01_diffusion_1d_out.e`, `case01_diffusion_1d_out.csv` |
| 02 | `case02_diffusion_2d_out.e` |
| 03 | `case03_heat_transient_out.e`, `case03_heat_transient_out.csv` |
| 04 | `case04_manufactured_out.e`, `case04_manufactured_out.csv` |
| 05 | `case05_varying_k_out.e`, `case05_varying_k_out.csv` |
| 06 | `case06_two_materials_out.e`, `case06_two_materials_out.csv` |
| 07 | `case07_nonlinear_diffusion_out.e`, `case07_nonlinear_diffusion_out.csv` |
| 08 | `case08_advection_diffusion_out.e`, `case08_advection_diffusion_out.csv` |
| 09 | `case09_coupled_system_out.e`, `case09_coupled_system_out.csv` |
| 10 | `case10_adaptive_refinement_out.e` |
| 11 | `case11_adaptive_dt_out.e`, `case11_adaptive_dt_out.csv` |
| 12 | `case12_parent_out.e`, `case12_parent_out_thermal_sub0.e` |
| 13 | `case13_postprocessors_out.e`, `case13_postprocessors_out.csv` |
| 14 | `case14_thermoelasticity_out.e`, `case14_thermoelasticity_out.csv` |
| 15 | `case15_lid_driven_cavity_out.e`, `case15_lid_driven_cavity_out.csv` |
| 16 | `case16_natural_convection_out.e`, `case16_natural_convection_out.csv` |
| 17 | `case17_joule_heating_out.e`, `case17_joule_heating_out.csv` |
| 18 | `case18_cahn_hilliard_out.e`, `case18_cahn_hilliard_out.csv` |
| 19 | `case19_porous_flow_out.e`, `case19_porous_flow_out.csv` |
| 20 | `case20_elastic_wave_exodus.e`, `case20_elastic_wave_out.csv` |
| 21 | `case21_bimetallic_strip_out.e`, `case21_bimetallic_strip_out.csv` |
| 22 | `case22_charge_relaxation_out.e`, `case22_charge_relaxation_out.csv` |
| 23 | `case23_magnetic_diffusion_out.e`, `case23_magnetic_diffusion_out.csv` |
| 24 | `case24_drift_diffusion_out.e`, `case24_drift_diffusion_out.csv` |
| 25 | `case25_induction_heating_out.e`, `case25_induction_heating_out.csv` |
| 26 | `case26_ehd_pumping_out.e`, `case26_ehd_pumping_out.csv` |
| 27 | `case27_hartmann_flow_out.e`, `case27_hartmann_flow_out.csv` |
| 28 | `case28_twoway_joule_heating_out.e`, `case28_twoway_joule_heating_out.csv` |
| 29 | `case29_electroconvection_out.e`, `case29_electroconvection_out.csv` |
| 30 | `case30_waveguide_cutoff_out.e`, `case30_waveguide_cutoff_out.csv` |
| 31 | `case31_driven_cavity_out.e`, `case31_driven_cavity_out.csv` |
| 32 | `case32_dielectric_slab_out.e`, `case32_dielectric_slab_out.csv` |
| 33 | `case33_coupled_resonators_out.e`, `case33_coupled_resonators_out.csv` |
| 34 | `case34_thermal_noise_out.e`, `case34_thermal_noise_out.csv` |
| 35 | `case35_dispersive_pulse_out.e`, `case35_dispersive_pulse_out.csv` |
| 36 | `case36_soliton_pulse_out.e`, `case36_soliton_pulse_out.csv` |
| 37 | `case37_rayleigh_benard_out.e`, `case37_rayleigh_benard_out.csv` |
| 38 | `case38_kelvin_helmholtz_out.e`, `case38_kelvin_helmholtz_out.csv` |
| 39 | `case39_blasius_boundary_layer_out.e`, `case39_blasius_boundary_layer_out.csv` |
| 40 | `case40_turbulent_channel_out.e`, `case40_turbulent_channel_out.csv` |
| 41 | `case41_rayleigh_taylor_out.e`, `case41_rayleigh_taylor_out.csv` |
| 42 | `case42_sod_shock_tube_out.e`, `case42_sod_shock_tube_out.csv` |
| 43 | `case43_ekman_spiral_out.e`, `case43_ekman_spiral_out.csv` |
| 44 | `case44_alfven_wave_out.e`, `case44_alfven_wave_out.csv` |
| 45 | `case45_monte_carlo_uq_out_storage_0001.csv.0` | Sampled avg_T and max_T for 30 MC samples |
| 46 | `case46_polynomial_chaos_out_eval_0002.csv` | 100 PCE surrogate predictions |
| 47 | `case47_main_out.csv`, `case47_main_out_OptimizationReporter_0001.csv` | Objective convergence and recovered parameter |
| 48 | `case48_parameter_study_csv_study_results_0001.csv` | 50 LHS parameter samples and QoI values |
| 49 | `case49_j2_plasticity_out.e`, `case49_j2_plasticity_out.csv` |
| 50 | `case50_finite_strain_compression_out.e`, `case50_finite_strain_compression_out.csv` |
| 51 | `case51_power_law_creep_out.e`, `case51_power_law_creep_out.csv` |
| 52 | `case52_phase_field_fracture_out.e`, `case52_phase_field_fracture_out.csv` |
| 53 | `case53_pressure_vessel_out.e`, `case53_pressure_vessel_out.csv` |
| 54 | `case54_neutron_diffusion_bare_slab_out.e`, `case54_neutron_diffusion_bare_slab_out.csv` |
| 55 | `case55_two_group_diffusion_out.e`, `case55_two_group_diffusion_out.csv` |
| 56 | `case56_fuel_pin_heat_transfer_out.e`, `case56_fuel_pin_heat_transfer_out.csv` |
| 57 | `case57_xenon_poisoning_out.e`, `case57_xenon_poisoning_out.csv` |
| 58 | `case58_control_rod_worth_out.e`, `case58_control_rod_worth_out.csv` |
| 59 | `case59_terzaghi_consolidation_out.e`, `case59_terzaghi_consolidation_out.csv` |
| 60 | `case60_wellbore_drawdown_out.e`, `case60_wellbore_drawdown_out.csv` |
| 61 | `case61_unsaturated_flow_out.e`, `case61_unsaturated_flow_out.csv` |
| 62 | `case62_biot_poroelasticity_out.e`, `case62_biot_poroelasticity_out.csv` |
| 63 | `case63_gravity_dam_out.e`, `case63_gravity_dam_out.csv` |
| 64 | `case64_reaction_diffusion_out.e`, `case64_reaction_diffusion_out.csv` |
| 65 | `case65_contaminant_transport_out.e`, `case65_contaminant_transport_out.csv` |
| 66 | `case66_mineral_precipitation_out.e`, `case66_mineral_precipitation_out.csv` |
| 67 | `case67_aqueous_equilibrium_out.e`, `case67_aqueous_equilibrium_out.csv` |
| 68 | `case68_calcite_dissolution_out.e`, `case68_calcite_dissolution_out.csv` |
| 69 | `case69_multiapp_coupling_out.e`, `case69_multiapp_coupling_out.csv`, `case69_multiapp_coupling_out_sub0.e` |
| 70 | `case70_contact_mechanics_out.e`, `case70_contact_mechanics_out.csv` |
| 71 | `case71_xfem_crack_out.e`, `case71_xfem_crack_out.csv` |
| 72 | `case72_thm_pipe_flow_out.e`, `case72_thm_pipe_flow_out.csv` |
| 73 | `case73_level_set_advection_out.e`, `case73_level_set_advection_out.csv` |
| 74 | `case74_left_handed_material_out.e`, `case74_left_handed_material_out.csv` |
| 75 | `case75_drude_slab_out.e`, `case75_drude_slab_out.csv` |
| 76 | `case76_3d_waveguide_out.e`, `case76_3d_waveguide_out.csv` |
| 77 | `case77_cylinder_scattering_out.e`, `case77_cylinder_scattering_out.csv` |
| 78 | `case78_pulse_reflection_out.e` |
| 79 | `case79_snell_law_tir_out.e`, `case79_snell_law_tir_out.csv` |
| 80 | `case80_bragg_mirror_out.e`, `case80_bragg_mirror_out.csv` |
| 81 | `case81_photonic_crystal_out.e`, `case81_photonic_crystal_out.csv`, `case81_photonic_crystal_out_eigenvalues_0002.csv` |
| 82 | `case82_3d_cavity_resonance_out.e`, `case82_3d_cavity_resonance_out_eigenvalues_0002.csv` |
| 83 | `case83_veselago_lens_out.e`, `case83_veselago_lens_out.csv` |
| 84 | `case84_lossy_tem_cavity_out.e`, `case84_lossy_tem_cavity_out.csv`, `case84_lossy_tem_cavity_out_eigenvalues_0002.csv` |
| 85 | `case85_hertzian_dipole_out.e`, `case85_hertzian_dipole_out.csv` |
| 86 | `case86_half_wave_dipole_out.e`, `case86_half_wave_dipole_out.csv` |
| 87 | `case87_phased_array_out.e`, `case87_phased_array_out.csv` |
| 88 | `case88_single_slit_diffraction_out.e`, `case88_single_slit_diffraction_out.csv` |
| 89 | `case89_dielectric_waveguide_out.e`, `case89_dielectric_waveguide_out.csv`, `case89_dielectric_waveguide_out_eigenvalues_0002.csv` |
| 90 | `case90_parabolic_reflector_out.e`, `case90_parabolic_reflector_out.csv` |
| 91 | `case91_radar_cross_section_out.e`, `case91_radar_cross_section_out.csv` |
| 92 | `case92_interferometer_out.e`, `case92_interferometer_out.csv` |
| 93 | `case93_inverse_source_recovery_out.csv`, `case93_inverse_source_recovery_out_OptimizationReporter_0001.csv`, `case93_inverse_source_recovery_out_forward0.e` |

All pre-run output files for cases 01-13 are included in this directory so you can
examine them without running anything. Cases 14-93 require `combined-opt` (all modules)
to run — see each case's README for Docker instructions.

---

## 7. Creating Your Own Simulations

Once you understand the 93 cases, you will want to adapt them or build new simulations
from scratch. This section walks through the process systematically.

### Step 1: Define Your Physics

Ask: what PDE governs my problem?

Common choices:
- **Steady-state diffusion**: `-div(k * grad u) = f`. Use `MatDiffusion` or `Diffusion`.
- **Transient diffusion/heat equation**: `rho*cp * du/dt = div(k * grad u) + Q`. Add `TimeDerivative`.
- **Convection-diffusion**: `du/dt + div(v*u) - div(D*grad u) = 0`. Add `ConservativeAdvection`.
- **Poisson equation**: `-div(grad u) = f`. Use `Diffusion` + `BodyForce`.
- **Coupled system**: write one equation per variable, use `CoupledForce` for coupling terms.

If you do not immediately know the PDE, consult a textbook on the physical process
you are modeling. MOOSE ships with modules covering heat conduction
(`modules/heat_transfer`), solid mechanics (`modules/solid_mechanics`), and fluid flow
(`modules/navier_stokes`) that provide pre-built kernels for these equations.

### Step 2: Choose the Domain and Mesh

Simple domains (rectangles, boxes): use `GeneratedMesh` or `GeneratedMeshGenerator`.
Complex geometries: create a mesh with Cubit, Gmsh, or Salome and read it with `FileMesh`.

Decide on dimensionality (1D, 2D, or 3D) and a mesh resolution. Start coarse and refine
to check that results do not change significantly (mesh independence study).

For problems with localized features (thin layers, sharp corners), use `[Adaptivity]` to
let MOOSE refine the mesh automatically.

### Step 3: Identify Boundary Conditions

For every boundary of your domain (every external face in 3D, every edge in 2D):
- **Known value** (temperature, concentration): use `DirichletBC` or `FunctionDirichletBC`
- **Known flux** (heat flux, mass flux): use `NeumannBC`
- **Insulated / no-flow / symmetry**: leave unspecified (zero-flux is the default)
- **Convection to ambient**: use `ConvectiveHeatFluxBC` (from the heat_transfer module)

Write down all boundary conditions before writing the input file. A missing or wrong
boundary condition is the most common source of physically wrong results.

### Step 4: Select Appropriate MOOSE Objects

Match each term of your PDE to a kernel:
- `du/dt` term: `TimeDerivative`
- `-div(grad u)` or `-div(k * grad u)`: `Diffusion` or `MatDiffusion` (or AD versions)
- Source/sink `f(x,t)`: `BodyForce` with a constant value or a `function = my_fn`
- Coupling to another variable `v`: `CoupledForce`
- Advection `div(v * u)`: `ConservativeAdvection`

For nonlinear problems (material properties that depend on the solution), use the `AD`
(Automatic Differentiation) versions: `ADMatDiffusion`, `ADParsedMaterial`, etc.
These give you the exact Jacobian for free and enable `solve_type = NEWTON` for
quadratic convergence.

### Step 5: Set Solver Parameters

For most problems, the default PJFNK + BoomerAMG preconditioner works well:
```text
[Executioner]
  type       = Steady       # or Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]
```

For transient problems, start with `ConstantDT` with a small time step, then switch to
`IterationAdaptiveDT` once your problem runs successfully to completion.

If the solve does not converge:
- Try tighter convergence tolerances (`nl_rel_tol = 1e-10`)
- Try smaller initial time step
- Try a different preconditioner (`-pc_type lu` for a small direct solve, useful for debugging)
- Check that all material properties are positive and physical
- Check that boundary conditions make physical sense

### Step 6: Add Postprocessors

Postprocessors let you monitor the simulation without loading the full Exodus file.
Good defaults for most problems:

```text
[Postprocessors]
  [max_value]
    type       = ElementExtremeValue
    variable   = T
    value_type = max
  []
  [avg_value]
    type     = ElementAverageValue
    variable = T
  []
[]
```

For verification purposes, add `ElementL2Error` if you have an exact solution.
For transient problems, add `TimestepSize` to monitor the time step.

### Step 7: Template Input File

Copy this template and fill in the blanks:

```text
# ============================================================
# [Your problem name here]
# [Brief description of the PDE]
# ============================================================

# Top-level constants (optional, but good practice)
k_value = 1.0
Q_value = 0.0

[Mesh]
  type = GeneratedMesh
  dim  = 2      # change to 1 or 3 as needed
  nx   = 20
  ny   = 20
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
[]

[Variables]
  [u]           # rename to T, c, phi, etc.
  []
[]

# Uncomment for transient problems:
# [ICs]
#   [u_init]
#     type     = ConstantIC
#     variable = u
#     value    = 0.0
#   []
# []

[Kernels]
  # Add your PDE terms here.
  [diffusion]
    type     = Diffusion
    variable = u
  []
  # [time_deriv]     # uncomment for transient
  #   type     = TimeDerivative
  #   variable = u
  # []
  # [source]
  #   type     = BodyForce
  #   variable = u
  #   value    = ${Q_value}
  # []
[]

[BCs]
  # Add boundary conditions for each boundary.
  [left]
    type     = DirichletBC
    variable = u
    boundary = left
    value    = 0.0
  []
  [right]
    type     = DirichletBC
    variable = u
    boundary = right
    value    = 1.0
  []
[]

# [Materials]    # uncomment if using MatDiffusion or material properties
#   [props]
#     type        = GenericConstantMaterial
#     prop_names  = 'k'
#     prop_values = '${k_value}'
#   []
# []

[Postprocessors]
  [avg_u]
    type     = ElementAverageValue
    variable = u
  []
  [max_u]
    type       = ElementExtremeValue
    variable   = u
    value_type = max
  []
[]

[Executioner]
  type = Steady        # change to Transient for time-stepping

  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  # Uncomment for transient:
  # [TimeStepper]
  #   type = ConstantDT
  #   dt   = 0.01
  # []
  # start_time = 0.0
  # end_time   = 1.0
[]

[Outputs]
  exodus = true
  csv    = true
[]
```

---

## 8. Glossary

This glossary defines every technical term that appears in the 21 tutorial cases.
If you encounter an unfamiliar word while reading an input file or console output,
look it up here.

---

**Advection (also: convection)**
Transport of a quantity (temperature, concentration) by a flowing fluid or moving medium.
The advective flux of a scalar `c` carried by velocity `v` is `v * c`. Advection moves
the quantity without changing its value — like a leaf carried by a river current. In a PDE
it appears as `div(v * c)` (the divergence of the advective flux). Advection is distinct
from diffusion: advection transports by bulk flow, diffusion transports from high to low
concentration regardless of any mean flow.

---

**AMG (Algebraic Multigrid)**
A class of iterative linear solver preconditioners that work at multiple scales. AMG
builds a hierarchy of progressively coarser approximations to the matrix, uses the
coarse-grid solutions to accelerate convergence on the fine grid, and reverses the process.
It is very effective for elliptic PDEs (diffusion, elasticity, Poisson). In MOOSE,
`boomeramg` from the HYPRE library is the most commonly used AMG preconditioner.

---

**Assembly**
The process of building the global stiffness matrix `K` and right-hand side vector `f`
from element-level contributions. MOOSE loops over all elements, calls each kernel to
compute its contribution to the element stiffness matrix and load vector, and adds these
to the global sparse matrix. Boundary condition contributions are added similarly.
The assembled system `K*u = f` is then passed to the linear solver.

---

**Automatic Differentiation (AD)**
A computational technique for computing exact derivatives of a function automatically,
without hand-coding the Jacobian. MOOSE uses dual-number arithmetic: every floating-point
variable carries not just its value but also its derivative with respect to the degrees
of freedom. When you use `AD` variants of kernels and materials (e.g., `ADMatDiffusion`,
`ADParsedMaterial`), MOOSE computes the exact Jacobian automatically. This enables
`solve_type = NEWTON` and achieves quadratic convergence even for nonlinear problems.

---

**Backward Euler**
A first-order implicit time integration method. For the ODE `du/dt = f(u)`, Backward Euler
approximates: `(u_new - u_old) / dt = f(u_new)`. Note that `f` is evaluated at the
**new** time level — this makes the method **implicit** (u_new appears on both sides)
and requires solving a nonlinear system at each time step. The advantage is
unconditional stability: no matter how large the time step, the method does not blow up
(though accuracy degrades for large steps). MOOSE uses Backward Euler by default.

---

**Block (mesh subdomain)**
A named region of the mesh containing a subset of elements. Used to apply different
materials to different regions of the domain. In `GeneratedMesh`, all elements are in
block 0 by default. Additional blocks are created with `SubdomainBoundingBoxGenerator`
or similar mesh generators. The `block` parameter in `[Materials]` restricts a material
to specific blocks.

---

**Body force**
A volumetric source term in the PDE, distributed throughout the interior of the domain.
In heat transfer, a body force corresponds to a volumetric heat source `Q` (W/m³) such
as nuclear fission heating or Joule heating. In MOOSE, the `BodyForce` kernel adds
the term `integral( phi_i * Q ) dV` to the right-hand side.

---

**Boundary condition (BC)**
A constraint applied at the boundary (edge or face) of the computational domain.
Required to make the PDE well-posed (to have a unique solution). The two main types are:
- Dirichlet: specifies the value of the unknown directly
- Neumann: specifies the flux (derivative) of the unknown normal to the boundary
See also: Dirichlet BC, Neumann BC.

---

**Convergence (nonlinear)**
The nonlinear solve is said to have converged when the residual `||R||` has been reduced
to below a specified tolerance. MOOSE checks both relative convergence (`||R|| / ||R_0|| <
nl_rel_tol`) and absolute convergence (`||R|| < nl_abs_tol`). If Newton's method fails
to converge within `nl_max_its` iterations, the solve fails.

---

**Convergence (order of accuracy)**
For a numerical method, convergence order refers to how fast the numerical error decreases
as the mesh is refined. For standard bilinear finite elements (the MOOSE default), the L2
error in the solution decreases as `O(h^2)`: halving the mesh spacing reduces the error
by a factor of 4. This is called second-order spatial convergence. Higher-order elements
or post-processing techniques can achieve higher orders.

---

**Degrees of Freedom (DOFs)**
The unknowns in the discrete system. For a scalar field like temperature with first-order
Lagrange elements, there is one DOF at each mesh node. For a mesh with `N` nodes and one
variable, there are `N` DOFs. For multiple variables or higher-order elements, the DOF
count is higher. The total DOF count determines the size of the linear system that must
be solved.

---

**Diffusion**
The process by which a quantity spreads from regions of high concentration to regions of
low concentration due to random molecular motion (Fick's law for mass, Fourier's law for
heat). Mathematically described by the term `-div(k * grad u)` where `k` is the
diffusivity or conductivity. In steady state, diffusion drives the solution toward the
average of its neighboring values.

---

**Dirichlet boundary condition**
A type of boundary condition that specifies the value of the unknown directly at the
boundary. Named after the mathematician Peter Gustav Lejeune Dirichlet. Example:
fixing the temperature at a wall to 100°C is a Dirichlet BC. In MOOSE: `DirichletBC`
with a constant value, `FunctionDirichletBC` with a spatially or temporally varying value.

---

**DOF**
Abbreviation for Degree of Freedom. See: Degrees of Freedom.

---

**Element**
One small cell in the computational mesh. In 2D, elements are typically quadrilaterals
(quads) or triangles. In 3D, they are hexahedra (bricks), tetrahedra, prisms, or
pyramids. MOOSE uses linear (first-order) elements by default, meaning the solution
is approximated by a piecewise bilinear (2D) or trilinear (3D) polynomial within
each element.

---

**Exodus II**
A binary file format for storing finite element mesh and solution data. Developed by
Sandia National Laboratories. Used by MOOSE as the primary output format. Files have
the extension `.e`. Viewable with ParaView, VisIt, and Cubit. Built on top of the netCDF
data format.

---

**FEM (Finite Element Method)**
A numerical technique for solving PDEs. The domain is divided into elements, the solution
is approximated by polynomials within each element, and the PDE is required to hold in a
weighted-average (weak) sense. The result is a large sparse linear or nonlinear algebraic
system. The name "finite element" refers to the finite-size elements that partition the domain.

---

**Finite Element**
A small region (element) in the mesh together with its associated shape functions. The
"finite element" in the name of the method refers to the fact that the domain is divided
into a finite number of small regions, contrasting with the infinitely many points in the
continuous domain.

---

**GMRES (Generalized Minimum Residual method)**
An iterative linear algebra solver for non-symmetric sparse systems. At each iteration,
it finds the best solution in an expanding Krylov subspace (span of `{b, Ab, A²b, ...}`
where `A` is the system matrix and `b` is the right-hand side). Requires a preconditioner
for fast convergence on ill-conditioned systems. The default Krylov solver in MOOSE/PETSc.

---

**Implicit method**
A time integration method where the right-hand side is evaluated at the new time level
(`t + dt`), requiring the solution of a nonlinear system at each step. Implicit methods
are unconditionally stable (no time step size limit for stability) but require more work
per step. MOOSE's default Backward Euler method is implicit. Contrasted with explicit
methods (like Forward Euler) which evaluate at the old time level and are conditionally
stable.

---

**Initial condition (IC)**
The value of the unknown field(s) at time `t = 0`. Required for transient problems.
Set via the `[ICs]` block. Without an explicit IC, MOOSE defaults to zero.

---

**Jacobian**
The matrix of partial derivatives of the residual vector with respect to the unknowns:
`J_ij = dR_i / du_j`. In Newton's method, the Jacobian describes how the residual
changes when the solution changes, and the linear system `J * delta_u = -R` determines
the Newton correction. For nonlinear problems, the Jacobian changes at every Newton
iteration (it depends on the current solution). MOOSE can compute the Jacobian by finite
differences (PJFNK), from hand-coded formulas, or automatically via AD.

---

**Kernel**
In MOOSE, a kernel is a C++ object that implements one term of the governing PDE. Each
kernel computes a contribution to the residual vector (and optionally the Jacobian matrix)
by integrating over all mesh elements. The name comes from the mathematical concept of
an integral kernel. The `[Kernels]` block in the input file instantiates and configures
these objects.

---

**Krylov method**
A class of iterative linear solvers (including GMRES and CG) that work by building a
subspace spanned by matrix-vector products `{b, Ab, A^2b, ...}` (the Krylov subspace)
and finding the best approximate solution in that subspace. Krylov methods are the standard
approach for large sparse systems because they only require matrix-vector products (not
matrix inversions or factorizations). They require a preconditioner to converge in a
reasonable number of iterations.

---

**L2 norm / L2 error**
The L2 norm of a function `f` over a domain is `sqrt( integral( f^2 ) dV )`. It measures
the "size" of a function in a global, averaged sense (unlike the maximum norm which just
looks at the largest value). The L2 error between a numerical solution `u_h` and an exact
solution `u_exact` is `sqrt( integral( (u_h - u_exact)^2 ) dV )`. For second-order
finite elements, the L2 error scales as `O(h^2)` where `h` is the element size.

---

**libMesh**
The finite element library that MOOSE is built on. libMesh provides: mesh data structures
(storage and manipulation of nodes, elements, and connectivity), shape function
computation, quadrature rules, degree-of-freedom management, and parallel mesh partitioning
using ParMETIS. MOOSE adds the physics layer on top of libMesh.

---

**Material property**
A physical quantity associated with the material in the domain (conductivity, density,
Young's modulus, etc.) that is needed by kernels and BCs. In MOOSE, material properties
are declared in the `[Materials]` block and looked up by name at quadrature points during
assembly. Material properties can be scalars, vectors, or tensors, and can depend on
position, time, or the current solution values.

---

**Mesh**
The subdivision of the computational domain into a collection of small, simple shapes
(elements). The mesh determines the computational grid on which the solution is approximated.
A finer mesh (more, smaller elements) gives more accurate results but takes more
computational time and memory.

---

**Method of Manufactured Solutions (MMS)**
A verification technique for numerical codes. The idea: choose a desired exact solution
first, derive the source term `f` that makes it satisfy the PDE, solve the PDE with that
source term, and compare the numerical result to the known exact solution. If the error
does not decrease at the theoretically expected rate as the mesh is refined, the code has
a bug. MMS does not require the manufactured solution to be physically motivated —
it just needs to be smooth and satisfy any boundary conditions you apply.

---

**Neumann boundary condition**
A type of boundary condition that specifies the flux (the derivative of the unknown in
the direction normal to the boundary) at the boundary. Named after Carl Neumann. Example:
specifying that 500 W/m² of heat flows into the domain through a surface. In MOOSE:
`NeumannBC`. The natural (default) Neumann condition is zero flux — no heat flows through
unspecified boundaries. This is why you do not need to explicitly write a zero-flux BC in MOOSE.

---

**Newton's method (Newton-Raphson)**
An iterative root-finding algorithm for solving nonlinear equations `F(u) = 0`. At each
step: (1) compute the Jacobian `J = dF/du` at the current guess, (2) solve `J * delta = -F`,
(3) update `u += delta`. For well-behaved problems, convergence is quadratic: the number
of correct decimal places doubles each iteration. For a linear problem (where `F` is linear
in `u`), Newton converges in exactly one iteration.

---

**Node**
A point in the mesh. Nodes are the corners of elements (and midside points for
second-order elements). Primary variable DOFs are stored at nodes for Lagrange elements.
A 20x20 mesh in 2D has 21x21 = 441 nodes (one more than the number of elements in each
direction because the mesh is divided into edges, and both endpoints of each edge are nodes).

---

**Nonlinear residual**
See: Residual. The "nonlinear" qualifier emphasizes that this is the residual of the full
nonlinear system `F(u) = 0`, as opposed to the linear residual of the inner linear solve.

---

**Parallel computing**
Running a simulation on multiple CPU cores or multiple computers simultaneously, with
each processor handling a portion of the work. MOOSE supports two styles:
- **MPI** (distributed memory): the mesh is partitioned across processes, each with its
  own memory; communication occurs through message-passing. Run with `mpirun -n N`.
- **Threads** (shared memory): multiple threads share the same memory; controlled with
  `--n-threads=N`. Limited to one node.

---

**ParaView**
A free, open-source scientific visualization application for analyzing and visualizing
large datasets. The standard tool for viewing MOOSE Exodus output files. Available at
https://www.paraview.org/. Supports 2D/3D rendering, streamlines, isosurfaces, slices,
time animation, and plotting.

---

**PDE (Partial Differential Equation)**
An equation relating a function to its partial derivatives with respect to two or more
independent variables (usually space and time). The governing equations of most physical
phenomena (heat conduction, fluid flow, structural mechanics, electromagnetism) are PDEs.
Example: the 2D heat equation `dT/dt = k*(d²T/dx² + d²T/dy²)` is a PDE because it
involves partial derivatives with respect to `t`, `x`, and `y`.

---

**PETSc (Portable Extensible Toolkit for Scientific Computation)**
A large open-source library of data structures and algorithms for scientific computing,
particularly for solving large systems of PDEs in parallel. MOOSE uses PETSc for its
linear algebra (sparse matrices, vectors) and linear/nonlinear solvers (GMRES, Newton).
PETSc provides the Krylov solvers, preconditioners (including the interface to HYPRE/AMG),
and nonlinear solver infrastructure.

---

**PJFNK (Preconditioned Jacobian-Free Newton-Krylov)**
The default nonlinear solver strategy in MOOSE. It combines:
- **Newton**: outer nonlinear iteration (Newton's method) to drive `F(u) = 0`
- **Jacobian-free**: the Jacobian-vector product is approximated by finite differences of
  the residual (`J*v ≈ (F(u + ε*v) - F(u)) / ε`), avoiding explicit Jacobian assembly
- **Krylov**: GMRES iterative solver for the linear system at each Newton step
- **Preconditioned**: a preconditioner (typically AMG) is used to accelerate GMRES

PJFNK works for almost any problem and requires no hand-coded Jacobian. The tradeoff is
that the Jacobian approximation is slightly inexact, so Newton convergence is superlinear
rather than quadratic.

---

**Postprocessor**
A MOOSE object that computes a single scalar number from the solution at each output step.
Examples: maximum temperature, domain-averaged concentration, L2 error, current time step
size. Results are written to the CSV output file and printed on the console. Postprocessors
are defined in the `[Postprocessors]` block.

---

**Preconditioner**
A transformation applied to a linear system `A*x = b` to make it easier for an iterative
solver (like GMRES) to solve. The preconditioned system `M_inv * A * x = M_inv * b` has
a better-conditioned matrix (eigenvalues more tightly clustered) so GMRES converges in
fewer iterations. Choosing a good preconditioner is problem-specific. For elliptic PDEs
(diffusion, elasticity), algebraic multigrid (AMG) is typically the best choice.

---

**Quadrature point**
A location inside an element where integrals are evaluated by Gaussian quadrature.
Instead of integrating exactly (which is complicated for arbitrary element geometries),
MOOSE evaluates the integrand at a small number of carefully chosen points and takes
a weighted sum. For bilinear quad elements, a 2×2 arrangement of four quadrature points
per element is standard and integrates polynomials of degree up to 3 exactly. Material
properties are evaluated at quadrature points during assembly.

---

**Residual**
The vector `R = F(u)` that measures how poorly the current solution `u` satisfies the
governing equations. In the finite element context, `R_i` is the integral of the PDE
error weighted by test function `phi_i`. When `R = 0`, the discrete equations are
satisfied exactly. Newton's method drives `R` toward zero by iteratively improving `u`.
The residual norm `||R||` is printed at each Newton iteration and is the primary
convergence metric in MOOSE console output.

---

**Shape function**
A polynomial function associated with a mesh node used to approximate the solution within
elements. For first-order Lagrange elements, the shape function `phi_i` at node `i` equals
1 at node `i`, 0 at all other nodes, and varies linearly (bilinearly in 2D) in between.
The solution field is written as `u(x) = sum_i u_i * phi_i(x)` where `u_i` are the
nodal values (the DOFs). Shape functions are the basis for the finite element approximation.

---

**Steady state / Steady-state problem**
A problem where the solution does not change with time — the system has reached
equilibrium. Mathematically: `du/dt = 0`, so the time-derivative term vanishes and the
PDE reduces to a system in space only. MOOSE solves steady-state problems with
`type = Steady` in `[Executioner]`.

---

**Stiffness matrix**
The global sparse matrix `K` assembled from element-level stiffness matrices. For a linear
elliptic PDE like `-div(k * grad u) = f`, the stiffness matrix encodes how the residual
changes with respect to the solution: `K_ij = integral( k * grad(phi_i) . grad(phi_j) dV )`.
Solving `K * u = f` gives the finite element solution. For nonlinear problems, the
"stiffness matrix" at each Newton step is the Jacobian evaluated at the current iterate.

---

**Sub-application (sub-app)**
In MOOSE's MultiApp system, a sub-application is a separate MOOSE simulation launched
and controlled by a parent application. The sub-app has its own input file, mesh,
variables, and kernels. Data is exchanged between parent and sub-app via `[Transfers]`.

---

**Subdomain**
A named region of elements in the mesh (also called a "block"). Elements in different
subdomains can have different material properties. Subdomains are identified by integer
IDs. In `GeneratedMesh`, all elements are in subdomain (block) 0. Additional subdomains
are created with mesh generators like `SubdomainBoundingBoxGenerator`.

---

**Test function**
In the finite element method, the test function (or weight function) `phi_i` is a shape
function used to convert the PDE from a pointwise condition into an integral (weak form)
condition. The PDE is required to hold in the sense that its integral against every test
function is zero. The number of test functions equals the number of DOFs, giving a system
of equations of the same size as the number of unknowns.

---

**Time integration**
The method used to advance the solution in time in a transient simulation. MOOSE defaults
to Backward Euler (first-order implicit). Other options include Crank-Nicolson (second-order
implicit), BDF2 (second-order implicit multi-step), and others. Implicit methods are
preferred because they allow large time steps without numerical instability.

---

**Time step (dt)**
The size of one increment in time during a transient simulation. Smaller time steps give
more accurate time resolution but require more total computation. Larger time steps are
faster but may miss rapid transients or cause the nonlinear solver to fail. Adaptive time
stepping (`IterationAdaptiveDT`) automatically adjusts `dt` to balance accuracy and cost.

---

**Transient problem**
A time-dependent problem where the solution evolves from initial conditions as time advances.
Governed by a PDE with a time-derivative term (`du/dt`). MOOSE solves transient problems
with `type = Transient` in `[Executioner]`. The solution at each time step is found by
Newton iteration. Results at each stored time step appear as separate frames in the Exodus file.

---

**Upwinding**
A numerical technique for advection-dominated problems that adds artificial diffusion in the
direction of flow to prevent unphysical oscillations (Gibbs phenomenon). Without upwinding,
standard Galerkin finite elements can produce wiggly solutions near sharp concentration
fronts when advection dominates diffusion. In MOOSE's `ConservativeAdvection` kernel, set
`upwinding_type = full` to enable full upwinding. Case 8 demonstrates the behavior.

---

**Variable (MOOSE variable)**
A scalar or vector field that MOOSE solves for (primary variable) or computes as a
by-product (auxiliary variable). Primary variables are declared in `[Variables]` and have
their own PDE equations implemented by kernels. Auxiliary variables are declared in
`[AuxVariables]` and are computed by `[AuxKernels]` or populated by transfers.

---

**Weak form**
The foundation of the finite element method. Instead of requiring the PDE to hold pointwise
(the "strong form"), the weak form requires the PDE residual to be zero when integrated
against every test function. This lower regularity requirement allows approximate piecewise
polynomial solutions. Integration by parts in the weak form converts second-order spatial
derivatives into first-order, enabling the use of the piecewise linear (first-order) shape
functions that MOOSE uses by default.

For example, the strong form `-div(grad u) = f` becomes the weak form:
`integral( grad(phi_i) . grad(u) dV ) = integral( phi_i * f dV )` for every test function `phi_i`.

---

**VisIt**
A free, open-source scientific visualization tool from Lawrence Livermore National Laboratory.
An alternative to ParaView for viewing Exodus II files. Available at
https://visit-dav.github.io/visit-website/. Particularly well-regarded for large-scale
parallel simulation data.

---

*End of README — return to the individual case READMEs in each subdirectory for
detailed walkthroughs of each simulation.*
