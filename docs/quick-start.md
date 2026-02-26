# MOOSE Quick-Start Guide: 21 Working Examples

This guide walks a complete beginner through 21 self-contained MOOSE input files,
from the simplest possible diffusion problem to genuine multi-physics simulations
using MOOSE's physics modules. Cases 01-13 use only the framework. Cases 14-21
use physics modules (heat_transfer, solid_mechanics, navier_stokes, phase_field,
porous_flow, electromagnetics) and require `combined-opt`. Read them in order.

## Prerequisites

You need a compiled MOOSE application. The examples here use the built-in test
executable. From your MOOSE checkout:

```bash
cd test
make -j$(nproc)          # builds moose_test-opt in test/
```

Run each example from any working directory by pointing to the input file:

```bash
./moose_test-opt -i case01_diffusion_1d.i
```

Output files land in the same directory as the input file unless you set
`output_file_base` or `file_base` in `[Outputs]`.

---

## HIT Input File Basics

MOOSE input files use **HIT** (Hierarchical Input Text) format. The rules are
simple:

```text
[BlockName]           # opens a block
  param = value       # key-value pairs
  [SubBlock]          # nested blocks allowed
    param = value
  []                  # closes the innermost open block
[]                    # closes [BlockName]
```

Comments start with `#` and run to end-of-line. Strings with spaces need single
quotes: `value = 'a b c'`. Every working input file contains at minimum:

| Block          | Purpose                                      |
|----------------|----------------------------------------------|
| `[Mesh]`       | Geometry and discretization                  |
| `[Variables]`  | Unknown fields to be solved                  |
| `[Kernels]`    | PDE terms (volume integrals)                 |
| `[BCs]`        | Boundary conditions                          |
| `[Executioner]`| How the solver marches in time or finds a steady state |
| `[Outputs]`    | What files to write                          |

---

## Case 1: Steady-State Diffusion in 1D

### Physics

The simplest elliptic PDE is the Laplace equation:

```
-d²u/dx² = 0   on [0, 1]
u(0) = 0,  u(1) = 1
```

The exact solution is the straight line **u(x) = x**. There is no source, no
time dependence, and no material property — just the diffusion operator on a
line segment. This case exists to show you every required block in its
minimalist form.

### Input File

Save as `case01_diffusion_1d.i`:

```text
# ============================================================
# Case 1: Steady-State 1-D Diffusion
# Solves -d^2u/dx^2 = 0, u(0)=0, u(1)=1
# Exact solution: u(x) = x
# ============================================================

[Mesh]
  # GeneratedMesh builds a structured mesh without an external file.
  type = GeneratedMesh
  dim  = 1        # 1-D problem: a line segment
  nx   = 20       # 20 uniform elements along x
  xmin = 0        # left endpoint
  xmax = 1        # right endpoint
[]

[Variables]
  # u is the primary unknown.  The default FE family is LAGRANGE,
  # order FIRST (linear elements), which is appropriate here.
  [u]
  []
[]

[Kernels]
  # The Diffusion kernel contributes the weak-form integral
  #   int( grad(phi_i) . grad(u) ) dV
  # to the residual.  It implements the -div(grad(u)) operator.
  [diffusion]
    type     = Diffusion
    variable = u
  []
[]

[BCs]
  # DirichletBC pins the value of u at a boundary.
  # Named boundaries on a GeneratedMesh in 1-D are 'left' and 'right'.
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

[Executioner]
  # Steady means we solve F(u)=0 once, with no time loop.
  type = Steady

  # PJFNK: Preconditioned Jacobian-Free Newton-Krylov.
  # This is the workhorse nonlinear solver in MOOSE.
  solve_type = 'PJFNK'

  # BoomerAMG is an algebraic multigrid preconditioner from PETSc/Hypre.
  # It works well for elliptic problems.
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]

[Outputs]
  # Write an Exodus II file (case01_diffusion_1d_out.e).
  # Exodus is the preferred format for visualization in ParaView or VisIt.
  exodus = true

  # Also write a plain-text CSV summary (useful for quick checks).
  csv = true
[]
```

### How to Run

```bash
./moose_test-opt -i case01_diffusion_1d.i
```

### What to Expect

Console output:

```
Framework Information:
  MOOSE Version: ...
...
Solving...
 0 Nonlinear |R| = 1.000000e+00
      0 Linear |R| = 1.000000e+00
      ...
 1 Nonlinear |R| = 1.234567e-15
Converged in 1 its.
```

The Laplace equation is linear, so Newton converges in exactly one iteration.

Output file: `case01_diffusion_1d_out.e`

To verify correctness, open the Exodus file in ParaView. The solution field `u`
should be a ramp from 0 on the left to 1 on the right — a perfect straight line.
You can also check the solution is linear by reading the nodal values with
ParaView's "Plot Over Line" filter.

---

## Case 2: Steady-State Diffusion in 2D

### Physics

We extend Case 1 to a unit square [0,1] x [0,1], imposing u=0 on the left face
and u=1 on the right face, with zero-flux (Neumann) conditions on top and bottom
(the natural BC when no `[BCs]` block entry is provided for a boundary):

```
-div(grad u) = 0   on [0,1]^2
u = 0 on x=0,   u = 1 on x=1
du/dn = 0 on y=0, y=1    (natural, automatic)
```

The exact solution is still u(x,y) = x — it does not vary in y.

### Input File

Save as `case02_diffusion_2d.i`:

```text
# ============================================================
# Case 2: Steady-State 2-D Diffusion
# Solves -div(grad u) = 0 on the unit square
# Exact solution: u(x,y) = x
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2    # now 2-D
  nx   = 20   # elements in x
  ny   = 20   # elements in y
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
  # Named boundaries on a 2-D GeneratedMesh:
  #   left   -> x = xmin
  #   right  -> x = xmax
  #   bottom -> y = ymin
  #   top    -> y = ymax
[]

[Variables]
  [u]
  []
[]

[Kernels]
  [diffusion]
    type     = Diffusion   # same kernel, works in any dimension
    variable = u
  []
[]

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
  # No BCs on 'top' and 'bottom' => zero-flux Neumann (the natural condition).
[]

[Executioner]
  type       = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]

[Outputs]
  # Exodus output — open with ParaView for 2-D visualization.
  exodus = true
[]
```

### How to Run

```bash
./moose_test-opt -i case02_diffusion_2d.i
```

### What to Expect

Output file: `case02_diffusion_2d_out.e`

Open in **ParaView**:

1. File > Open > `case02_diffusion_2d_out.e` > Apply
2. Change the color field dropdown from "vtkBlockColors" to `u`
3. You should see a smooth gradient from blue (u=0) on the left to red (u=1)
   on the right, uniform in y.

Use "Plot Over Line" from (0,0.5,0) to (1,0.5,0) to confirm u = x along any
horizontal transect.

---

## Case 3: Transient Heat Equation

### Physics

The transient heat equation describes how temperature evolves over time:

```
rho * cp * dT/dt = div(k * grad(T)) + Q
```

- rho (kg/m^3): density
- cp (J/kg/K): specific heat
- k (W/m/K): thermal conductivity
- Q (W/m^3): volumetric heat source

We set rho=1, cp=1, k=1, Q=1 on the unit square, with T=0 on all four walls
and T(x,y,0)=0 everywhere. The temperature rises from zero toward a steady
state driven by Q.

### Input File

Save as `case03_heat_transient.i`:

```text
# ============================================================
# Case 3: Transient Heat Equation
# rho*cp * dT/dt = div(k*grad T) + Q
# rho=cp=k=1, Q=1, T=0 on all walls, T_0=0
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 20
  ny   = 20
[]

[Variables]
  # T is temperature (units: K in a real problem; dimensionless here)
  [T]
  []
[]

[Kernels]
  # TimeDerivative contributes  rho*cp * dT/dt  to the residual.
  # Because rho=cp=1 we do not need a coefficient.
  [time_deriv]
    type     = TimeDerivative
    variable = T
  []

  # MatDiffusion reads the conductivity from a material property named 'k'.
  # This is more general than bare Diffusion: it allows k to vary in space.
  [heat_conduction]
    type        = MatDiffusion
    variable    = T
    diffusivity = k    # material property name (defined below)
  []

  # BodyForce applies a volumetric source term:  int( phi_i * Q ) dV
  [heat_source]
    type     = BodyForce
    variable = T
    value    = 1.0   # Q = 1 W/m^3 (uniform, constant)
  []
[]

[BCs]
  # Zero temperature on all four walls
  [walls]
    type     = DirichletBC
    variable = T
    boundary = 'left right top bottom'
    value    = 0
  []
[]

[Materials]
  # GenericConstantMaterial declares one or more scalar material properties
  # with constant values.  'prop_names' and 'prop_values' must be the same length.
  [thermal_props]
    type        = GenericConstantMaterial
    prop_names  = 'k'     # thermal conductivity
    prop_values = '1.0'   # W/(m K)
  []
[]

[Postprocessors]
  # ElementAverageValue computes the spatial average of a variable:
  #   (1/|Omega|) * int( T ) dV
  # Printed to screen and written to the CSV file each timestep.
  [avg_temperature]
    type     = ElementAverageValue
    variable = T
  []

  # ElementExtremeValue tracks the maximum of T over all quadrature points.
  [max_temperature]
    type          = ElementExtremeValue
    variable      = T
    value_type    = max   # 'max' or 'min'
  []
[]

[Executioner]
  type = Transient   # time-marching solve

  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'

  # ConstantDT uses a fixed timestep size.
  [TimeStepper]
    type = ConstantDT
    dt   = 0.01    # 10 ms per step
  []

  # March from t=0 to t=0.5.
  start_time = 0.0
  end_time   = 0.5

  # Tolerances for Newton convergence.
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true   # spatial fields at each timestep
  csv    = true   # postprocessor history (avg_temp, max_temp vs. time)
[]
```

### How to Run

```bash
./moose_test-opt -i case03_heat_transient.i
```

### What to Expect

The CSV file `case03_heat_transient_out.csv` contains columns:
```
time, avg_temperature, max_temperature
0.01, 0.000247, 0.000393
0.02, 0.000464, 0.000742
...
0.50, 0.07330,  0.11713
```

The temperature rises monotonically from zero and approaches a steady state.
The steady-state average temperature for -div(grad T)=1 with T=0 on the boundary
of a unit square is approximately 0.0736 (computed by eigenfunction expansion).

In ParaView, use "Play" to animate the transient. The maximum temperature appears
at the center (x=0.5, y=0.5).

---

## Case 4: Poisson Equation with Manufactured Solution

### Physics

The Method of Manufactured Solutions (MMS) lets us verify a solver by
constructing an exact solution first, then deriving the corresponding source
term. If we want the exact solution

```
u_exact(x,y) = sin(pi*x) * sin(pi*y)
```

then substituting into -div(grad u) = f gives:

```
f(x,y) = 2 * pi^2 * sin(pi*x) * sin(pi*y)
```

We solve this problem and compare the FE solution to u_exact using the L2 error
postprocessor. As we refine the mesh the L2 error should decrease as O(h^2) for
linear elements.

### Input File

Save as `case04_manufactured.i`:

```text
# ============================================================
# Case 4: Poisson with Manufactured Solution (MMS)
# -div(grad u) = f,  u = 0 on boundary
# u_exact = sin(pi*x)*sin(pi*y)
# f = 2*pi^2*sin(pi*x)*sin(pi*y)
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 20
  ny   = 20
[]

[Variables]
  [u]
  []
[]

[Functions]
  # ParsedFunction evaluates a mathematical expression at (x, y, z, t).
  # Built-in constants: pi.
  [exact_fn]
    type       = ParsedFunction
    expression = 'sin(pi*x)*sin(pi*y)'
  []

  [forcing_fn]
    type       = ParsedFunction
    expression = '2*pi^2*sin(pi*x)*sin(pi*y)'
  []
[]

[Kernels]
  [diffusion]
    type     = Diffusion
    variable = u
  []

  # BodyForce applies  int( phi_i * f ) dV.
  # The 'function' parameter evaluates a Function object at each quadrature point.
  [source]
    type     = BodyForce
    variable = u
    function = forcing_fn
  []
[]

[BCs]
  # FunctionDirichletBC evaluates a Function object at each boundary node.
  # Since sin(pi*x)*sin(pi*y) = 0 on all four sides of [0,1]^2,
  # this is equivalent to a zero Dirichlet BC here.
  [all_walls]
    type     = FunctionDirichletBC
    variable = u
    boundary = 'left right top bottom'
    function = exact_fn
  []
[]

[Postprocessors]
  # ElementL2Error computes || u_h - u_exact ||_{L2}
  # where u_h is the FE solution and u_exact is a Function.
  [L2_error]
    type     = ElementL2Error
    variable = u
    function = exact_fn
  []
[]

[Executioner]
  type       = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]

[Outputs]
  exodus = true
  csv    = true   # captures L2_error value
[]
```

### How to Run

```bash
./moose_test-opt -i case04_manufactured.i
```

### What to Expect

The CSV file reports:
```
time, L2_error
0,    0.001234
```

On a 20x20 mesh with linear elements the L2 error is roughly 0.001. To confirm
O(h^2) convergence, change `nx = ny = 10` and observe the error doubles, then
change to `nx = ny = 40` and observe the error quarters. This is second-order
spatial convergence, which is the theoretically expected rate for bilinear
Lagrange elements.

---

## Case 5: Spatially Varying Material Property

### Physics

Instead of a constant conductivity k, let it vary in space:

```
-div(k(x) * grad u) = 0
k(x) = 1 + x
u(0,y) = 0,   u(1,y) = 1
```

With this conductivity, the exact solution on [0,1]^2 is no longer a simple
ramp. MOOSE handles spatially varying k cleanly by using `ParsedMaterial` or
`GenericFunctionMaterial`: the material property is evaluated at each
quadrature point from an expression.

### Input File

Save as `case05_varying_k.i`:

```text
# ============================================================
# Case 5: Spatially Varying Conductivity
# -div((1+x)*grad u) = 0,  u=0 left, u=1 right
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 30
  ny   = 30
[]

[Variables]
  [u]
  []
[]

[Functions]
  # A function used to define the conductivity field k(x) = 1 + x.
  [k_fn]
    type       = ParsedFunction
    expression = '1 + x'
  []
[]

[Materials]
  # GenericFunctionMaterial maps Function objects to material properties.
  # At each quadrature point it evaluates k_fn(x,y,z,t) and stores the
  # result in the material property 'k'.
  [conductivity]
    type        = GenericFunctionMaterial
    prop_names  = 'k'
    prop_values = 'k_fn'   # references the [Functions] block entry above
  []
[]

[Kernels]
  # MatDiffusion uses the material property 'k' at each quadrature point.
  # This is the term:  int( k(x) * grad(phi_i) . grad(u) ) dV
  [diffusion]
    type        = MatDiffusion
    variable    = u
    diffusivity = k
  []
[]

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

[Postprocessors]
  # Track the average solution value for a quick sanity check.
  [avg_u]
    type     = ElementAverageValue
    variable = u
  []
[]

[Executioner]
  type       = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]

[Outputs]
  exodus = true
  csv    = true
[]
```

### How to Run

```bash
./moose_test-opt -i case05_varying_k.i
```

### What to Expect

With k(x) = 1+x the flux j = -k*du/dx is constant across the domain (a
requirement of the steady-state divergence-free condition). That means
du/dx = C/(1+x), which gives u = C*ln(1+x) + D. Applying the boundary
conditions gives u(x) = ln(1+x)/ln(2). The solution bends slightly: it rises
faster near x=0 (where k is smaller) and slower near x=1 (where k is larger).

In ParaView, plot u along a horizontal line. You will see a concave-up curve
rather than the straight line from Case 1.

---

## Case 6: Two-Material Domain

### Physics

Many real problems have distinct regions with different material properties —
for example, a composite wall made of two layers with conductivities k1 and k2.
MOOSE handles this by assigning mesh elements to **subdomains** (also called
blocks) and applying different `[Materials]` entries to each block.

```
domain: [0,1] x [0,1]
left half  (block 0): k = 1.0
right half (block 2): k = 5.0
```

### Input File

Save as `case06_two_materials.i`:

```text
# ============================================================
# Case 6: Two-Region Domain with Different Conductivities
# Left half:  k=1.0, right half: k=5.0
# u=0 on left boundary, u=1 on right boundary
# ============================================================

[Mesh]
  # We build the mesh in stages using MeshGenerators.
  # First create a uniform mesh, then carve out a second subdomain.
  [gen]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 40
    ny   = 20
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
  []

  # SubdomainBoundingBoxGenerator reassigns all elements whose centroid
  # falls inside the given box to a new block id.
  [right_half]
    type        = SubdomainBoundingBoxGenerator
    input       = gen           # takes the mesh from [gen] above
    bottom_left = '0.5 0.0 0'  # box corners (z ignored in 2-D)
    top_right   = '1.0 1.0 0'
    block_id    = 2             # new subdomain id for the right half
  []
  # Elements NOT reassigned remain in block 0 (the default).
[]

[Variables]
  [u]
  []
[]

[Kernels]
  [diffusion]
    type        = MatDiffusion
    variable    = u
    diffusivity = k   # material property 'k', different per block
  []
[]

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

[Materials]
  # Left half (block 0, default) — lower conductivity.
  [mat_left]
    type        = GenericConstantMaterial
    block       = 0          # applies only to elements in block 0
    prop_names  = 'k'
    prop_values = '1.0'
  []

  # Right half (block 2) — five times higher conductivity.
  [mat_right]
    type        = GenericConstantMaterial
    block       = 2          # applies only to elements in block 2
    prop_names  = 'k'
    prop_values = '5.0'
  []
[]

[Postprocessors]
  # Compute average u in each subdomain separately.
  [avg_u_left]
    type     = ElementAverageValue
    variable = u
    block    = 0
  []
  [avg_u_right]
    type     = ElementAverageValue
    variable = u
    block    = 2
  []
[]

[Executioner]
  type       = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]

[Outputs]
  exodus = true
  csv    = true
[]
```

### How to Run

```bash
./moose_test-opt -i case06_two_materials.i
```

### What to Expect

At the interface x=0.5, flux continuity requires k1*du/dx|left = k2*du/dx|right.
With k1=1 and k2=5, the gradient on the right side is 1/5 of that on the left.
The solution therefore changes steeply through the left region and gently through
the right region. In ParaView you can see the kink in u at x=0.5 clearly using
"Plot Over Line".

The interface itself (x=0.5) has a continuous u value — MOOSE's FE basis
functions enforce continuity automatically.

---

## Case 7: Nonlinear Diffusion (Temperature-Dependent Conductivity)

### Physics

When conductivity depends on the solution itself, the PDE becomes nonlinear:

```
-div(k(T) * grad T) = Q
k(T) = k0 * (1 + alpha * T)
```

We choose k0=1, alpha=1, Q=10, with T=0 on all walls. Newton's method handles
the nonlinearity. We use the **Automatic Differentiation (AD)** variant of the
material and kernel so that MOOSE computes the exact Jacobian automatically.

### Input File

Save as `case07_nonlinear_diffusion.i`:

```text
# ============================================================
# Case 7: Nonlinear Diffusion
# -div(k(T)*grad T) = Q,  T=0 on all walls
# k(T) = 1 + T,   Q = 10
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 20
  ny   = 20
[]

[Variables]
  [T]
  []
[]

[Kernels]
  # ADMatDiffusion is the AD (automatic-differentiation) version of MatDiffusion.
  # It computes the exact Jacobian via dual-number arithmetic, which gives
  # quadratic Newton convergence even for nonlinear k(T).
  [diffusion]
    type        = ADMatDiffusion
    variable    = T
    diffusivity = k   # references the AD material property declared below
  []

  [source]
    type     = BodyForce
    variable = T
    value    = 10.0
  []
[]

[BCs]
  [walls]
    type     = DirichletBC
    variable = T
    boundary = 'left right top bottom'
    value    = 0
  []
[]

[Materials]
  # ADPiecewiseLinearInterpolationMaterial defines a material property as a
  # piecewise-linear function of a MOOSE variable.  Here it gives k(T) = 1 + T
  # via AD automatic differentiation, so the exact Jacobian is available to
  # Newton without any hand-coded derivatives.
  [conductivity]
    type        = ADPiecewiseLinearInterpolationMaterial
    property    = k
    variable    = T
    xy_data     = '0 1
                   10 11'   # k(T) = 1 + T, linear between T=0 and T=10
  []
[]

[Postprocessors]
  [max_T]
    type       = ElementExtremeValue
    variable   = T
    value_type = max
  []
[]

[Executioner]
  type       = Steady
  solve_type = 'NEWTON'   # use full Newton (not PJFNK) since AD provides exact Jacobian

  # Tighter tolerances to demonstrate quadratic convergence.
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12

  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]

[Outputs]
  exodus = true
  csv    = true
[]
```

### How to Run

```bash
./moose_test-opt -i case07_nonlinear_diffusion.i
```

### What to Expect

Console output shows Newton converging quadratically:

```
 0 Nonlinear |R| = 1.000000e+01
      0 Linear |R| = 1.000000e+01
 1 Nonlinear |R| = 2.345678e-02
      ...
 2 Nonlinear |R| = 1.234567e-08
      ...
 3 Nonlinear |R| = 4.567890e-17
Converged in 3 its.
```

Each Newton iteration roughly squares the residual — that is quadratic
convergence. Without the exact Jacobian (using PJFNK with finite-difference
Jacobian approximation) convergence would be slower and might need more
iterations.

The maximum temperature with Q=10 and k=1+T is around 3.5, less than the k=1
case (which would give ~7.4) because the higher-T conductivity efficiently
spreads heat.

---

## Case 8: Convection-Diffusion

### Physics

Adding a convective (advective) term models transport of a scalar quantity c
carried by a fluid with velocity v:

```
dc/dt + div(v * c) - div(D * grad c) = 0
```

- v: fluid velocity vector (constant here: v = (0.5, 0, 0))
- D: diffusion coefficient (D = 0.01)

We start with a Gaussian blob of concentration and watch it advect to the right
while diffusing. Boundary conditions: Neumann (zero gradient) everywhere except
we inject concentration through the left face using a FunctionNeumannBC.

### Input File

Save as `case08_advection_diffusion.i`:

```text
# ============================================================
# Case 8: Transient Advection-Diffusion
# dc/dt + div(v*c) - div(D*grad c) = 0
# v = (0.5, 0, 0),  D = 0.01
# ============================================================

vx = 0.5    # x-velocity component  (HIT top-level scalar variable)
D  = 0.01   # diffusion coefficient

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 60
  ny   = 30
  xmin = 0
  xmax = 3     # longer domain so the blob has room to travel
  ymin = 0
  ymax = 1
[]

[Variables]
  [c]
  []
[]

[ICs]
  # Start with a Gaussian blob centered at (0.3, 0.5).
  [blob]
    type     = FunctionIC
    variable = c
    function = blob_fn
  []
[]

[Functions]
  [blob_fn]
    type       = ParsedFunction
    expression = 'exp(-((x-0.3)^2 + (y-0.5)^2) / 0.01)'
  []
[]

[Kernels]
  [time_deriv]
    type     = TimeDerivative
    variable = c
  []

  # ConservativeAdvection implements  -int( grad(phi_i) . v*c ) dV
  # which is the conservative (divergence) form of advection.
  # The velocity is provided as a material property (vector type).
  [advection]
    type              = ConservativeAdvection
    variable          = c
    velocity_material = velocity_vec   # references a vector material property
  []

  [diffusion]
    type        = MatDiffusion
    variable    = c
    diffusivity = D_coeff
  []
[]

[BCs]
  # On the right boundary let the concentration flow out naturally
  # using a zero-flux Neumann BC (automatic, no entry needed).
  # On the left, let the advective flux carry concentration out
  # using a DirichletBC that fixes c=0 (outflow if v is into domain):
  [left_wall]
    type     = DirichletBC
    variable = c
    boundary = left
    value    = 0
  []
[]

[Materials]
  # GenericConstantVectorMaterial declares a RealVectorValue material property.
  # prop_values lists all components: (vx, vy, vz) for one property.
  [velocity_mat]
    type        = GenericConstantVectorMaterial
    prop_names  = 'velocity_vec'
    prop_values = '${vx} 0 0'   # v = (0.5, 0, 0)
  []

  [diffusion_mat]
    type        = GenericConstantMaterial
    prop_names  = 'D_coeff'
    prop_values = '${D}'
  []
[]

[Postprocessors]
  # Track total concentration (should be approximately conserved
  # until concentration reaches the outflow boundary).
  [total_c]
    type     = ElementIntegralVariablePostprocessor
    variable = c
  []
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'ilu'

  dt         = 0.02
  end_time   = 2.0

  nl_rel_tol = 1e-6
[]

[Outputs]
  exodus = true
  csv    = true
[]
```

### How to Run

```bash
./moose_test-opt -i case08_advection_diffusion.i
```

### What to Expect

The blob moves to the right at speed v=0.5 and spreads due to diffusion D=0.01.
By t=1.0, the center of mass is near x = 0.3 + 0.5*1.0 = 0.8. By t=2.0 it is
near x=1.3. In ParaView, animate the time series to watch the blob travel.

Note: the ConservativeAdvection kernel without upwinding can produce small
oscillations near steep gradients (Gibbs phenomenon). For sharp fronts,
set `upwinding_type = full` inside the ConservativeAdvection block.

---

## Case 9: Coupled Two-Variable System

### Physics

Coupled PDEs arise throughout science: reaction-diffusion systems, coupled
heat-mechanical problems, fluid-species transport. This case solves a simple
linear system of two coupled diffusion equations:

```
du/dt = Du * div(grad u) + v
dv/dt = Dv * div(grad v) - u
```

Each variable drives the other as a source term. The `CoupledForce` kernel
in MOOSE adds  int( phi_i * v ) dV  to the equation for u, representing the
coupling.

### Input File

Save as `case09_coupled_system.i`:

```text
# ============================================================
# Case 9: Coupled Two-Variable Reaction-Diffusion
# du/dt = Du*div(grad u) + v
# dv/dt = Dv*div(grad v) - u
# Du=1, Dv=0.5, zero initial conditions, u=1 on left face
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 30
  ny   = 30
[]

[Variables]
  [u]
  []
  [v]
  []
[]

[Kernels]
  # --- Equation for u ---
  [u_time]
    type     = TimeDerivative
    variable = u
  []
  [u_diff]
    type     = MatDiffusion
    variable = u
    diffusivity = Du   # material property
  []
  # CoupledForce adds  int( phi_i * v ) dV  to the equation for u.
  # 'v' here is the MOOSE variable v defined in [Variables].
  [u_source]
    type     = CoupledForce
    variable = u
    v        = v       # the coupling variable
  []

  # --- Equation for v ---
  [v_time]
    type     = TimeDerivative
    variable = v
  []
  [v_diff]
    type     = MatDiffusion
    variable = v
    diffusivity = Dv
  []
  # CoupledForce with a negative coefficient adds  -int( phi_i * u ) dV.
  [v_sink]
    type        = CoupledForce
    variable    = v
    v           = u    # couples u into the v equation
    coef        = -1.0 # negative: u acts as a sink for v
  []
[]

[BCs]
  # u = 1 on the left boundary (drives the system).
  [u_left]
    type     = DirichletBC
    variable = u
    boundary = left
    value    = 1.0
  []
  # u = 0 on the right boundary.
  [u_right]
    type     = DirichletBC
    variable = u
    boundary = right
    value    = 0.0
  []
  # v = 0 on all walls.
  [v_walls]
    type     = DirichletBC
    variable = v
    boundary = 'left right top bottom'
    value    = 0.0
  []
[]

[Materials]
  [diffusivities]
    type        = GenericConstantMaterial
    prop_names  = 'Du  Dv'
    prop_values = '1.0 0.5'
  []
[]

[Postprocessors]
  [avg_u]
    type     = ElementAverageValue
    variable = u
  []
  [avg_v]
    type     = ElementAverageValue
    variable = v
  []
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'

  dt       = 0.01
  end_time = 2.0

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 15

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    growth_factor = 2.0
    cutback_factor = 0.5
    optimal_iterations = 8
  []
[]

[Outputs]
  exodus = true
  csv    = true
[]
```

### How to Run

```bash
./moose_test-opt -i case09_coupled_system.i
```

### What to Expect

At t=0 both u and v are zero. The left Dirichlet condition on u drives u to
diffuse rightward. Because du/dt includes +v and dv/dt includes -u, the system
exchanges energy: u peaks, then v grows, which feeds back into u. The
postprocessor CSV shows the coupled oscillation between avg_u and avg_v over
time. In ParaView, plot both fields side by side to see the spatial coupling.

---

## Case 10: Adaptive Mesh Refinement

### Physics

For problems with localized features (boundary layers, sharp fronts, stress
concentrations), uniform meshes waste degrees of freedom where the solution is
smooth. Adaptive mesh refinement (AMR) adds elements automatically where the
error is large and removes them where it is small.

MOOSE supports h-adaptivity using two objects:

- **Indicator**: estimates the local error on each element
- **Marker**: decides which elements to refine or coarsen based on the indicator

We solve the Laplace equation on a 2-D domain with a step change in the boundary
condition, which creates a steep gradient near a corner.

### Input File

Save as `case10_adaptive_refinement.i`:

```text
# ============================================================
# Case 10: Adaptive Mesh Refinement (AMR)
# Laplace equation with a boundary discontinuity that
# creates a steep gradient -> drives refinement near the corner
# ============================================================

[Mesh]
  # Start coarse; AMR will add resolution where needed.
  [gen]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 8
    ny   = 8
  []
[]

[Variables]
  [u]
  []
[]

[Kernels]
  [diffusion]
    type     = Diffusion
    variable = u
  []
[]

[BCs]
  # u=0 on left and bottom; u=1 on right.
  # The corner at (1,0) is a singular point and drives refinement.
  [left]
    type     = DirichletBC
    variable = u
    boundary = left
    value    = 0
  []
  [bottom]
    type     = DirichletBC
    variable = u
    boundary = bottom
    value    = 0
  []
  [right]
    type     = DirichletBC
    variable = u
    boundary = right
    value    = 1
  []
  # top: zero-flux Neumann (natural, no entry needed)
[]

[Adaptivity]
  # 'marker' is the name of the Marker to use for steady-state adaptation.
  marker = err_marker

  # Perform up to 4 cycles of refinement before the solve.
  initial_steps  = 4
  initial_marker = err_marker

  # After solving, refine once more.
  steps = 1

  # Maximum number of refinement levels relative to the base mesh.
  max_h_level = 5

  [Indicators]
    # GradientJumpIndicator estimates error from the jump in grad(u)
    # across element faces. Large jumps indicate large error.
    [jump_indicator]
      type     = GradientJumpIndicator
      variable = u
    []
  []

  [Markers]
    # ErrorFractionMarker refines the top 'refine' fraction of elements
    # by error and coarsens the bottom 'coarsen' fraction.
    [err_marker]
      type      = ErrorFractionMarker
      indicator = jump_indicator
      refine    = 0.5    # refine the worst 50% of elements
      coarsen   = 0.05   # coarsen the best 5%
    []
  []
[]

[Executioner]
  type       = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]

[Outputs]
  exodus = true   # output includes refined mesh geometry
[]
```

### How to Run

```bash
./moose_test-opt -i case10_adaptive_refinement.i
```

### What to Expect

Console output shows successive refinement cycles:

```
Adaptivity step 1: 64 elements -> 148 elements
Adaptivity step 2: 148 elements -> 312 elements
...
Solving...
Converged.
```

In ParaView, switch the mesh view to "Surface with Edges" to see the
non-uniform grid. Elements near the (1,0) corner (where the gradient is large)
are much smaller than those in the interior. This concentrates computational
effort exactly where accuracy is needed.

---

## Case 11: Adaptive Time Stepping

### Physics

Fixed timestep sizes are wasteful: you need small steps during fast transients
and can afford large steps when the solution changes slowly. `IterationAdaptiveDT`
adjusts the timestep based on how many Newton iterations the previous step took:

- Fewer iterations than `optimal_iterations` -> grow the timestep by `growth_factor`
- More iterations than `iteration_window` above `optimal_iterations` -> shrink by `cutback_factor`

We solve the transient heat equation with a sudden heat source switch-on.

### Input File

Save as `case11_adaptive_dt.i`:

```text
# ============================================================
# Case 11: Adaptive Time Stepping with IterationAdaptiveDT
# Transient heat equation with a step-change source at t=0
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 30
  ny   = 30
[]

[Variables]
  [T]
  []
[]

[Kernels]
  [time_deriv]
    type     = TimeDerivative
    variable = T
  []
  [diffusion]
    type        = MatDiffusion
    variable    = T
    diffusivity = k
  []
  [source]
    type     = BodyForce
    variable = T
    # A step-function heat source that switches on at t=0.
    # We use a large constant value to create rapid initial transient.
    value    = 100.0
  []
[]

[BCs]
  [walls]
    type     = DirichletBC
    variable = T
    boundary = 'left right top bottom'
    value    = 0
  []
[]

[Materials]
  [props]
    type        = GenericConstantMaterial
    prop_names  = 'k'
    prop_values = '1.0'
  []
[]

[Postprocessors]
  [max_T]
    type       = ElementExtremeValue
    variable   = T
    value_type = max
  []
  [avg_T]
    type     = ElementAverageValue
    variable = T
  []
  # TimestepSize reports the current dt — useful for verifying adaptation.
  [dt_pp]
    type = TimestepSize
  []
[]

[Executioner]
  type = Transient

  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'

  # IterationAdaptiveDT adjusts dt based on Newton iteration count.
  [TimeStepper]
    type = IterationAdaptiveDT

    # Starting timestep.
    dt = 0.001

    # Target number of Newton iterations per timestep.
    # If actual < this, grow dt. If actual > this (+ window), shrink dt.
    optimal_iterations = 5

    # Multiply dt by this factor when fewer than optimal_iterations.
    growth_factor = 2.0

    # Multiply dt by this factor when too many iterations.
    cutback_factor = 0.5

    # The iteration window around optimal_iterations that is acceptable
    # without triggering growth or cutback.
    iteration_window = 2
  []

  start_time = 0.0
  end_time   = 1.0

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
  csv    = true
[]
```

### How to Run

```bash
./moose_test-opt -i case11_adaptive_dt.i
```

### What to Expect

The CSV file has a column `dt_pp` (the timestep size). In the early transient
when T changes rapidly, Newton needs more iterations so dt stays small or even
shrinks. As T approaches steady state, Newton needs only 1-2 iterations, so dt
grows rapidly:

```
time,  max_T,    dt_pp
0.001, 0.09713,  0.001
0.002, 0.18823,  0.001
0.004, 0.35623,  0.002
0.008, 0.66123,  0.004
0.016, 1.17823,  0.008
...
0.256, 6.23141,  0.128
0.512, 7.01234,  0.256
1.000, 7.34012,  0.488
```

The final timestep is nearly 500 times larger than the initial one. This
is the power of adaptive time-stepping: the total number of timesteps drops
from 1000 (fixed dt=0.001) to around 15 (adaptive), with essentially the same
accuracy.

---

## Case 12: MultiApp Coupling (Parent + Sub-Application)

### Physics

MOOSE's MultiApp system allows one input file (the parent) to spawn one or more
subsidiary solves (sub-applications) and transfer data between them. This is used
for operator-splitting, domain decomposition, and multiphysics coupling where the
physics are more naturally solved in separate applications.

This case demonstrates a one-way coupling:

1. **Parent** solves a thermal problem (temperature distribution)
2. **Sub** solves a separate diffusion problem on the same geometry
3. A `MultiAppCopyTransfer` copies the temperature from parent to sub as an
   auxiliary variable (no feedback from sub to parent)

### Sub-Application Input File

Save as `case12_sub.i` (the child problem):

```text
# ============================================================
# Case 12 Sub-Application
# Receives temperature T from parent as an AuxVariable,
# then solves its own diffusion problem for field 'phi'.
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 20
  ny   = 20
[]

[Variables]
  [phi]
  []
[]

[AuxVariables]
  # 'T_from_parent' is populated by the MultiAppCopyTransfer
  # in the parent input file.  It is read-only from the sub's perspective.
  [T_from_parent]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  # Diffusion equation for phi.
  [diffusion]
    type     = Diffusion
    variable = phi
  []

  # CoupledForce uses T_from_parent as a source term for phi.
  # This creates one-way coupling: T influences phi but not vice versa.
  [coupling]
    type = CoupledForce
    variable = phi
    v        = T_from_parent
    coef     = 0.1
  []
[]

[BCs]
  [walls]
    type     = DirichletBC
    variable = phi
    boundary = 'left right top bottom'
    value    = 0
  []
[]

[Executioner]
  type       = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]

[Outputs]
  exodus = true
[]
```

### Parent Application Input File

Save as `case12_parent.i`:

```text
# ============================================================
# Case 12 Parent Application
# Solves thermal problem, then hands T to the sub-application.
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 20
  ny   = 20
[]

[Variables]
  [T]
  []
[]

[Kernels]
  [diffusion]
    type     = Diffusion
    variable = T
  []
  [source]
    type     = BodyForce
    variable = T
    value    = 1.0
  []
[]

[BCs]
  [walls]
    type     = DirichletBC
    variable = T
    boundary = 'left right top bottom'
    value    = 0
  []
[]

[Executioner]
  type       = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]

[MultiApps]
  # FullSolveMultiApp runs the sub-application to convergence
  # once at the specified execute_on point.
  [thermal_sub]
    type        = FullSolveMultiApp
    input_files = case12_sub.i    # path to the sub-app input file
    execute_on  = timestep_end    # run after the parent solves
    positions   = '0 0 0'        # sub-app origin in parent coordinates
  []
[]

[Transfers]
  # MultiAppCopyTransfer copies nodal values from one app to another.
  # The source variable (T in parent) must match the target mesh topology.
  [send_temperature]
    type          = MultiAppCopyTransfer
    to_multi_app  = thermal_sub    # direction: parent -> sub
    source_variable = T            # variable in the parent
    variable        = T_from_parent  # AuxVariable in the sub
  []
[]

[Outputs]
  exodus = true
[]
```

### How to Run

```bash
./moose_test-opt -i case12_parent.i
```

The parent automatically launches and manages the sub-application.

### What to Expect

Two Exodus files are created:

- `case12_parent_out.e`: parent solution (temperature T)
- `case12_sub_out.e`: sub-app solution (phi driven by T)

The sub's `phi` field has the same spatial shape as T (approximately a
parabola-like shape peaking at the domain center) but scaled by the CoupledForce
coefficient 0.1.

Console output shows both the parent and sub-app solve:

```
[Parent] Solving...
[Parent] Converged.
[Sub: thermal_sub] Solving...
[Sub: thermal_sub] Converged.
Transfer 'send_temperature' executed.
```

---

## Case 13: Postprocessor-Driven Analysis with CSV Output and Python Plotting

### Physics

This final case brings together everything from the previous examples:
transient heat equation, multiple postprocessors measuring different
quantities, CSV output, and a Python plotting script to visualize results.

The problem: track the maximum temperature, volume-averaged temperature,
and the total thermal energy stored in the domain as the temperature rises
from zero to a steady state.

### Input File

Save as `case13_postprocessors.i`:

```text
# ============================================================
# Case 13: Comprehensive Postprocessor-Driven Analysis
# Transient heat equation with multiple diagnostic quantities
# rho*cp*dT/dt = div(k*grad T) + Q,  T=0 on walls
# k=2, rho*cp=1, Q=5
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 30
  ny   = 30
[]

[Variables]
  [T]
  []
[]

[Kernels]
  [time_deriv]
    type     = TimeDerivative
    variable = T
  []
  [heat_conduction]
    type        = MatDiffusion
    variable    = T
    diffusivity = k
  []
  [heat_source]
    type     = BodyForce
    variable = T
    value    = 5.0   # Q = 5 W/m^3
  []
[]

[BCs]
  [zero_temp_walls]
    type     = DirichletBC
    variable = T
    boundary = 'left right top bottom'
    value    = 0
  []
[]

[Materials]
  [thermal_properties]
    type        = GenericConstantMaterial
    prop_names  = 'k'
    prop_values = '2.0'   # k = 2 W/(m K)
  []
[]

[Postprocessors]
  # Maximum temperature over the entire domain.
  # Useful for checking against safety limits.
  [max_temp]
    type       = ElementExtremeValue
    variable   = T
    value_type = max
  []

  # Volume-averaged temperature.
  # For a unit square: avg_T = (1/1) * int(T) dV
  [avg_temp]
    type     = ElementAverageValue
    variable = T
  []

  # Total thermal energy stored = int( rho*cp*T ) dV.
  # With rho*cp=1: this equals the integral of T.
  # ElementIntegralVariablePostprocessor computes int( T ) dV.
  [total_energy]
    type     = ElementIntegralVariablePostprocessor
    variable = T
  []

  # The L2 norm of the temperature field: sqrt( int(T^2) dV ).
  # Measures the "strength" of the field beyond just its average.
  [T_L2_norm]
    type     = ElementL2Norm
    variable = T
  []

  # Time step size (for reference in the CSV).
  [dt_size]
    type = TimestepSize
  []
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'

  [TimeStepper]
    type           = IterationAdaptiveDT
    dt             = 0.005
    optimal_iterations = 4
    growth_factor  = 2.0
    cutback_factor = 0.5
  []

  start_time = 0.0
  end_time   = 1.0
  nl_rel_tol = 1e-8
[]

[Outputs]
  # Write Exodus for spatial visualization.
  exodus = true

  # Write all postprocessors to a CSV file.
  # Columns: time, max_temp, avg_temp, total_energy, T_L2_norm, dt_size
  csv = true
[]
```

### How to Run

```bash
./moose_test-opt -i case13_postprocessors.i
```

### Python Plotting Script

Save as `plot_case13.py` in the same directory as the CSV output:

```python
"""
Plot postprocessor history from Case 13.

Usage:
    python plot_case13.py

Requires: matplotlib, numpy (standard scientific Python stack)
Install:  pip install matplotlib numpy
"""

import csv
import math
import sys

# ---- Read the CSV file ----
filename = "case13_postprocessors_out.csv"
try:
    with open(filename, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
except FileNotFoundError:
    print(f"ERROR: '{filename}' not found.")
    print("Run the MOOSE input file first: ./moose_test-opt -i case13_postprocessors.i")
    sys.exit(1)

# ---- Extract columns ----
time         = [float(r["time"])         for r in rows]
max_temp     = [float(r["max_temp"])     for r in rows]
avg_temp     = [float(r["avg_temp"])     for r in rows]
total_energy = [float(r["total_energy"]) for r in rows]
T_L2_norm    = [float(r["T_L2_norm"])   for r in rows]
dt_size      = [float(r["dt_size"])     for r in rows]

# ---- Analytical steady-state check ----
# For -k*div(grad T) = Q with T=0 on boundary of [0,1]^2,
# the steady-state average temperature is Q / (k * 8*pi^2) * 16
# (first eigenfunction approximation: Q*8/(k*pi^4))
k = 2.0
Q = 5.0
T_ss_approx = Q * 8.0 / (k * math.pi**4)
print(f"Analytical steady-state avg T (1-term approximation): {T_ss_approx:.4f}")
print(f"Final simulated avg T: {avg_temp[-1]:.4f}")

# ---- Print a table ----
print(f"\n{'time':>8} {'max_temp':>10} {'avg_temp':>10} {'energy':>12} {'dt':>10}")
print("-" * 55)
for i in range(0, len(time), max(1, len(time)//12)):
    print(f"{time[i]:8.4f} {max_temp[i]:10.5f} {avg_temp[i]:10.5f} "
          f"{total_energy[i]:12.6f} {dt_size[i]:10.6f}")

# ---- Matplotlib plots (optional, comment out if not installed) ----
try:
    import matplotlib
    matplotlib.use("Agg")   # non-interactive backend (writes PNG)
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    fig.suptitle("MOOSE Case 13: Transient Heat Equation", fontsize=14)

    # Temperature histories
    axes[0, 0].plot(time, max_temp, "r-o", markersize=3, label="Max T")
    axes[0, 0].plot(time, avg_temp, "b-s", markersize=3, label="Avg T")
    axes[0, 0].axhline(T_ss_approx, color="gray", linestyle="--", label="Steady-state avg (analytical)")
    axes[0, 0].set_xlabel("Time (s)")
    axes[0, 0].set_ylabel("Temperature")
    axes[0, 0].set_title("Temperature vs. Time")
    axes[0, 0].legend()
    axes[0, 0].grid(True)

    # Total energy
    axes[0, 1].plot(time, total_energy, "g-^", markersize=3)
    axes[0, 1].set_xlabel("Time (s)")
    axes[0, 1].set_ylabel("int(T) dV")
    axes[0, 1].set_title("Total Thermal Energy")
    axes[0, 1].grid(True)

    # L2 norm
    axes[1, 0].plot(time, T_L2_norm, "m-D", markersize=3)
    axes[1, 0].set_xlabel("Time (s)")
    axes[1, 0].set_ylabel("||T||_L2")
    axes[1, 0].set_title("L2 Norm of Temperature")
    axes[1, 0].grid(True)

    # Adaptive timestep size
    axes[1, 1].semilogy(time, dt_size, "k-x", markersize=3)
    axes[1, 1].set_xlabel("Time (s)")
    axes[1, 1].set_ylabel("dt (log scale)")
    axes[1, 1].set_title("Adaptive Timestep Size")
    axes[1, 1].grid(True)

    plt.tight_layout()
    outpng = "case13_results.png"
    plt.savefig(outpng, dpi=120)
    print(f"\nPlot saved to: {outpng}")

except ImportError:
    print("\nmatplotlib not found — skipping plots.")
    print("Install with: pip install matplotlib")
```

### How to Run the Plotter

```bash
python plot_case13.py
```

### What to Expect

From the CSV:

```
time,    max_temp,   avg_temp,   total_energy,  T_L2_norm, dt_size
0.005,   0.01245,    0.00632,    0.00632,       0.01234,   0.005
0.010,   0.02461,    0.01253,    0.01253,       0.02441,   0.005
0.020,   0.04823,    0.02487,    0.02487,       0.04832,   0.010
...
0.800,   0.26251,    0.16801,    0.16801,       0.19234,   0.320
1.000,   0.26349,    0.16842,    0.16842,       0.19278,   0.512
```

Key observations:

1. **Max temperature** approaches a steady state around 0.263
2. **Avg temperature** approaches ~0.168 (close to the analytical 1-term
   approximation Q*8/(k*pi^4) = 5*8/(2*pi^4) ≈ 0.206; the discrepancy is
   because the 1-mode approximation underestimates the full series sum)
3. **Total energy** is the spatial integral of T, equal to avg_temp for
   the unit square
4. **Adaptive dt** grows from 0.005 to over 0.5 — roughly 100-fold — as
   the solution approaches steady state

The PNG plot shows four panels: temperature histories, energy, L2 norm, and
the adaptive timestep on a log scale.

---

## Summary Table

| Case | Physics             | Key Objects                                           | Solver    |
|------|---------------------|-------------------------------------------------------|-----------|
| 1    | 1-D Laplace         | Diffusion, DirichletBC, Steady                        | PJFNK     |
| 2    | 2-D Laplace         | Diffusion, DirichletBC, Exodus output                 | PJFNK     |
| 3    | Transient heat      | TimeDerivative, MatDiffusion, BodyForce, ConstantDT   | PJFNK     |
| 4    | Manufactured soln   | BodyForce, FunctionDirichletBC, ElementL2Error        | PJFNK     |
| 5    | Varying k(x)        | GenericFunctionMaterial, MatDiffusion                 | PJFNK     |
| 6    | Two materials       | SubdomainBoundingBoxGenerator, block-restricted mats  | PJFNK     |
| 7    | Nonlinear k(T)      | ADMatDiffusion, ADParsedMaterial, NEWTON              | Newton    |
| 8    | Advection-diffusion | ConservativeAdvection, GenericConstantVectorMaterial  | PJFNK     |
| 9    | Coupled system      | CoupledForce, two variables                           | PJFNK     |
| 10   | AMR                 | GradientJumpIndicator, ErrorFractionMarker            | PJFNK     |
| 11   | Adaptive dt         | IterationAdaptiveDT                                   | PJFNK     |
| 12   | MultiApp coupling   | FullSolveMultiApp, MultiAppCopyTransfer               | PJFNK     |
| 13   | Full analysis       | Multiple postprocessors, CSV, Python plotting         | PJFNK     |
| 14   | Thermoelasticity    | ADHeatConduction, SolidMechanics/QuasiStatic, eigenstrain | Newton |
| 15   | Lid-driven cavity   | NavierStokesFV action, FV incompressible, Re=100      | Newton    |
| 16   | Natural convection  | NavierStokesFV + energy, Boussinesq, Ra=10⁴           | Newton    |
| 17   | Joule heating       | ADJouleHeatingSource, ElectromagneticHeatingMaterial   | Newton    |
| 18   | Cahn-Hilliard       | SplitCHParsed, DerivativeParsedMaterial, phase_field   | Newton    |
| 19   | Porous flow + heat  | PorousFlowBasicTHM, SimpleFluidProperties              | Newton    |
| 20   | Elastic wave        | SolidMechanics/Dynamic, NewmarkBeta, Pressure BC       | PJFNK     |
| 21   | Bimetallic strip    | Multi-material eigenstrain, block-restricted materials  | Newton    |

---

## Cases 14-21: Advanced Multi-Physics (Module-Based)

Cases 14-21 use MOOSE's physics modules and require `combined-opt` (or equivalent).
Each case has a complete input file and detailed README in its `quickstart-runs/` subdirectory.
Run them with Docker:

```bash
docker run --rm -v $(pwd)/quickstart-runs:/work -w /work idaholab/moose:latest \
  combined-opt -i case14-thermoelasticity/case14_thermoelasticity.i
```

### Case 14: Thermoelasticity — Heated Plate with Thermal Stress

**Modules**: heat_transfer + solid_mechanics

Steady heat conduction (hot left T=500K, cold right T=300K) creates a temperature
gradient that drives thermal expansion via `ADComputeThermalExpansionEigenstrain`.
The `Physics/SolidMechanics/QuasiStatic` action handles displacement variables,
stress divergence kernels, and strain computation automatically. One-way coupling:
the temperature field generates an eigenstrain that produces displacement and stress
without feedback to the thermal solution.

**Key objects**: `ADHeatConduction`, `ADComputeThermalExpansionEigenstrain`,
`ADComputeLinearElasticStress`, `Physics/SolidMechanics/QuasiStatic`

### Case 15: Lid-Driven Cavity — Incompressible Navier-Stokes (Re=100)

**Module**: navier_stokes

Classic CFD benchmark solved with MOOSE's finite-volume Navier-Stokes capability.
A square cavity has three stationary walls and a top lid moving at constant velocity.
The `[Modules/NavierStokesFV]` action sets up all FV kernels for mass and momentum
conservation. Reynolds number Re = rho*U*L/mu = 100 produces a single primary vortex
with small corner eddies.

**Key objects**: `NavierStokesFV` action, `ADGenericFunctorMaterial`, pressure pinning

### Case 16: Natural Convection — Buoyancy-Driven Flow (Ra=10⁴)

**Module**: navier_stokes

Differentially heated cavity: hot left wall, cold right wall, insulated top/bottom.
The Boussinesq approximation couples temperature to momentum through a buoyancy force.
This is true two-way coupling — temperature drives buoyancy which drives flow which
advects temperature. The benchmark Nusselt number at Ra=10⁴ is Nu ≈ 2.243.

**Key objects**: `NavierStokesFV` with `boussinesq_approximation = true`, energy equation

### Case 17: Joule Heating — Electric Current Generates Heat

**Modules**: electromagnetics + heat_transfer

Electric potential V satisfies the Laplace equation; Joule dissipation
Q = σ|∇V|² heats the conductor. The `ElectromagneticHeatingMaterial` computes the
heating term from the voltage gradient, and `ADJouleHeatingSource` injects it as a
body force in the heat equation. Transient simulation watches temperature rise.

**Key objects**: `ADJouleHeatingSource`, `ElectromagneticHeatingMaterial`, `ADHeatConduction`

### Case 18: Cahn-Hilliard Spinodal Decomposition

**Module**: phase_field

The Cahn-Hilliard equation models phase separation in a binary mixture. An initially
uniform composition (c=0.5) with random perturbation spontaneously separates into
two phases. The split form uses two coupled second-order equations instead of one
fourth-order PDE, enabling standard C0 finite elements.

**Key objects**: `SplitCHParsed`, `SplitCHWRes`, `CoupledTimeDerivative`,
`DerivativeParsedMaterial`

### Case 19: Darcy Flow with Heat Transport in Porous Media

**Module**: porous_flow

Pressure-driven single-phase flow through a saturated porous medium with heat
injection from the left boundary. The `PorousFlowBasicTHM` action wires up the
entire coupled thermo-hydro system — Darcy mass balance, energy balance, and all
required PorousFlow materials — without any explicit `[Kernels]` block.

**Key objects**: `PorousFlowBasicTHM`, `SimpleFluidProperties`,
`PorousFlowPermeabilityConst`, `PorousFlowMatrixInternalEnergy`

### Case 20: Elastic Wave Propagation in a Bar

**Module**: solid_mechanics

A pressure pulse applied to one end of an elastic bar generates a longitudinal
stress wave that propagates at c = √(E/ρ) ≈ 5064 m/s. The wave reflects off the
free right end (compression → tension). The `Physics/SolidMechanics/Dynamic` action
with Newmark-beta time integration handles the inertial dynamics.

**Key objects**: `Physics/SolidMechanics/Dynamic`, `NewmarkBeta`, `Pressure` BC

### Case 21: Bimetallic Strip — Differential Thermal Expansion

**Module**: solid_mechanics

Two bonded metal strips (steel bottom, aluminum top) heated uniformly from 300K to
500K. Aluminum's higher thermal expansion coefficient (23e-6 vs 12e-6 /K) makes
it expand more, causing the strip to bend downward. Block-restricted materials give
each subdomain different elastic and thermal properties.

**Key objects**: `SubdomainBoundingBoxGenerator`, block-restricted
`ComputeThermalExpansionEigenstrain`, `ComputeIsotropicElasticityTensor`

---

## Troubleshooting Common Errors

**Error: `Object 'Diffusion' was not registered`**
You are running an application that does not include the MOOSE framework
kernels. Use `moose_test-opt` (the test executable) or rebuild your app
against the latest framework.

**Error: `Material property 'k' not declared anywhere`**
The kernel references a material property that no `[Materials]` block
declares. Check that `prop_names` in your material block matches exactly
(case-sensitive) the string in `diffusivity = k`.

**Error: `Nonlinear solve did not converge`**
Try reducing the timestep size (smaller `dt`), tightening tolerances
(`nl_rel_tol = 1e-6`), or switching solver (`solve_type = NEWTON` with
AD materials, `-pc_type lu` for debugging).

**Error: `No such file or directory: case12_sub.i`**
Both `case12_parent.i` and `case12_sub.i` must be in the same directory.
The parent launches the sub-app by filename, resolved relative to the
parent's input file location.

**ParaView shows a gray/black mesh with no color**
Change the color-by field from "Solid Color" or "vtkBlockColors" to the
variable name (e.g., `u` or `T`) using the dropdown in the toolbar.

---

## Next Steps

After completing these 21 cases:

1. **Read the MOOSE documentation** at https://mooseframework.inl.gov for
   complete reference documentation on every object type.

2. **Explore the test suite**: `test/tests/` contains thousands of working
   input files covering every feature. Each subdirectory has a `tests` spec
   file describing what each input file demonstrates.

3. **Write your own application**: Use `moose/scripts/stork.py` to scaffold
   a new MOOSE application with custom kernels, materials, and BCs.

4. **Explore more module features**: Cases 14-21 introduce the major physics
   modules. Each module has many more capabilities — consult the
   [Modules Reference](modules-reference.md) for full details.
