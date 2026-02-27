# MOOSE Quick-Start Guide: 58 Working Examples

This guide walks a complete beginner through 58 self-contained MOOSE input files,
from the simplest possible diffusion problem to genuine multi-physics simulations
using MOOSE's physics modules. Cases 01-13 use only the framework. Cases 14-58
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

## Cases 14-29: Advanced Multi-Physics (Module-Based)

Cases 14-29 use MOOSE's physics modules and require `combined-opt` (which includes
all 25+ modules). If you are on Windows, run them with Docker:

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/path/to/quickstart-runs:/work" \
  -w /work/case14-thermoelasticity \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case14_thermoelasticity.i 2>&1 | tail -30'
```

On Linux/macOS with a local build:

```bash
combined-opt -i case14_thermoelasticity.i
```

---

## Case 14: Thermoelasticity — Heated Plate with Thermal Stress

### Physics

This case couples two physics: **heat conduction** and **solid mechanics**. A 2D
steel plate has its left edge held at 500 K and its right edge at 300 K. The
resulting temperature gradient causes the material to expand non-uniformly — the
hot side expands more than the cold side. This differential expansion generates
internal stresses and deformation even though no external mechanical load is applied.

The coupling is **one-way**: the temperature field drives the mechanical response
through a thermal eigenstrain, but the displacement does not feed back into the
temperature solution.

Governing equations:

```
Heat:       -div( k * grad(T) ) = 0           (steady-state conduction)
Mechanics:  div( sigma ) = 0                   (quasi-static equilibrium)
            epsilon_thermal = alpha * (T - T_ref) * I
            sigma = C : (epsilon_total - epsilon_thermal)
```

### Input File

Save as `case14_thermoelasticity.i`:

```text
# ============================================================
# Case 14: Thermoelasticity — Heated Plate with Thermal Stress
# Steady-state heat conduction drives thermal expansion in a
# 2D elastic solid.  Hot left (T=500K), cold right (T=300K).
# One-way coupling: T field -> eigenstrain -> displacement/stress
# Requires: combined-opt  (heat_transfer + solid_mechanics modules)
# ============================================================

# Material constants for structural steel
E     = 200e9   # Young's modulus, Pa
nu    = 0.3     # Poisson's ratio, dimensionless
alpha = 12e-6   # coefficient of thermal expansion, 1/K
k_th  = 50      # thermal conductivity, W/(m K)
cp    = 500     # specific heat, J/(kg K)
T_ref = 300     # stress-free (reference) temperature, K

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 20
  ny   = 20
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
[]

[Variables]
  [T]
    [InitialCondition]
      type  = FunctionIC
      function = '300 + 200*(1-x)'
    []
  []
[]

[Kernels]
  [heat_conduction]
    type     = ADHeatConduction
    variable = T
  []
[]

# The QuasiStatic action automatically creates disp_x and disp_y variables,
# stress divergence kernels, strain computation, and vonmises_stress output.
[Physics/SolidMechanics/QuasiStatic]
  [solid]
    strain             = SMALL
    add_variables      = true
    eigenstrain_names  = 'thermal_eigenstrain'
    generate_output    = 'vonmises_stress'
    use_automatic_differentiation = true
  []
[]

[BCs]
  [T_hot]
    type     = DirichletBC
    variable = T
    boundary = left
    value    = 500
  []
  [T_cold]
    type     = DirichletBC
    variable = T
    boundary = right
    value    = 300
  []
  [pin_bottom_x]
    type     = DirichletBC
    variable = disp_x
    boundary = bottom
    value    = 0
  []
  [pin_bottom_y]
    type     = DirichletBC
    variable = disp_y
    boundary = bottom
    value    = 0
  []
[]

[Materials]
  [thermal_props]
    type                = ADHeatConductionMaterial
    thermal_conductivity = ${k_th}
    specific_heat        = ${cp}
  []
  [elasticity_tensor]
    type          = ADComputeIsotropicElasticityTensor
    youngs_modulus = ${E}
    poissons_ratio = ${nu}
  []
  [thermal_eigenstrain]
    type                    = ADComputeThermalExpansionEigenstrain
    temperature             = T
    thermal_expansion_coeff = ${alpha}
    stress_free_temperature  = ${T_ref}
    eigenstrain_name        = thermal_eigenstrain
  []
  [stress]
    type = ADComputeLinearElasticStress
  []
[]

[Postprocessors]
  [max_temperature]
    type       = ElementExtremeValue
    variable   = T
    value_type = max
  []
  [max_disp_x]
    type       = ElementExtremeValue
    variable   = disp_x
    value_type = max
  []
  [max_vonmises]
    type       = ElementExtremeValue
    variable   = vonmises_stress
    value_type = max
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       mumps'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Outputs]
  exodus = true
  csv = true
[]
```

### How to Run

```bash
combined-opt -i case14_thermoelasticity.i
```

### What to Expect

The steady solver converges in a small number of Newton iterations. The temperature
field is a linear ramp from 500 K (left) to 300 K (right). The plate bows upward
because the hot side expands while the bottom is pinned. Von Mises stress is highest
near the constrained bottom edge where the thermal expansion is resisted.

Output files: `case14_thermoelasticity_out.e` (spatial fields: T, disp_x, disp_y,
vonmises_stress), `case14_thermoelasticity_out.csv` (postprocessor values).

---

## Case 15: Lid-Driven Cavity — Incompressible Navier-Stokes (Re=100)

### Physics

The lid-driven cavity is the most widely used benchmark problem in computational
fluid dynamics. A square cavity filled with viscous fluid has three stationary walls
and a top lid that moves horizontally at constant velocity U=1. The moving lid drags
fluid along, creating a recirculating flow pattern.

This case introduces MOOSE's **finite volume (FV)** Navier-Stokes solver, which is
fundamentally different from the finite element kernels used in Cases 1-14. The FV
formulation uses cell-centered unknowns and flux-based discretization.

```
Continuity:   div(v) = 0
Momentum:     rho * (v . grad) v = -grad(p) + mu * laplacian(v)
Re = rho * U * L / mu = 1 * 1 * 1 / 0.01 = 100
```

### Input File

Save as `case15_lid_driven_cavity.i`:

```text
# ============================================================
# Case 15: Lid-Driven Cavity — Incompressible Navier-Stokes
# Classic CFD benchmark: Re = rho*U*L/mu = 1*1*1/0.01 = 100
# Steady state, 2D square cavity, top wall moves at U=1
# ============================================================

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

[Modules]
  [NavierStokesFV]
    compressibility          = 'incompressible'
    porous_medium_treatment  = false
    add_energy_equation      = false

    density          = 'rho'
    dynamic_viscosity = 'mu'

    initial_velocity = '0 0 0'
    initial_pressure = 0.0

    # No-slip on left, right, and bottom walls
    wall_boundaries   = 'left right bottom'
    momentum_wall_types = 'noslip noslip noslip'

    # Moving lid (top) treated as a fixed-velocity inlet
    inlet_boundaries        = 'top'
    momentum_inlet_types    = 'fixed-velocity'
    momentum_inlet_functors = '1 0'

    mass_advection_interpolation     = 'average'
    momentum_advection_interpolation = 'average'

    # Pin pressure to remove the null space
    pin_pressure       = true
    pinned_pressure_type  = average
    pinned_pressure_value = 0
  []
[]

[FunctorMaterials]
  [fluid_properties]
    type        = ADGenericFunctorMaterial
    prop_names  = 'rho mu'
    prop_values = '1   0.01'
  []
[]

[Postprocessors]
  [max_vel_x]
    type       = ElementExtremeValue
    variable   = vel_x
    value_type = max
  []
  [max_vel_y]
    type       = ElementExtremeValue
    variable   = vel_y
    value_type = max
  []
  [avg_pressure]
    type     = ElementAverageValue
    variable = pressure
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 50
  l_tol     = 1e-6
  l_max_its = 200
[]

[Outputs]
  exodus = true
  csv    = true
[]
```

### How to Run

```bash
combined-opt -i case15_lid_driven_cavity.i
```

### What to Expect

Newton converges in roughly 10-20 iterations. The velocity field shows a large
primary vortex centered slightly above and to the right of the cavity center, with
small counter-rotating eddies in the bottom corners. This matches the published
Ghia et al. (1982) benchmark results at Re=100.

Key observations:
- All variables (vel_x, vel_y, pressure) are **element variables** — FV does not
  produce nodal data. Use element-centered plotting.
- The pressure is determined only up to an additive constant; the `pin_pressure`
  parameter fixes the average to zero.

Output files: `case15_lid_driven_cavity_out.e`, `case15_lid_driven_cavity_out.csv`.

---

## Case 16: Natural Convection — Buoyancy-Driven Flow (Ra=10⁴)

### Physics

A square cavity has a hot left wall (T=1) and a cold right wall (T=0), with
insulated top and bottom. The temperature difference drives fluid motion through
the Boussinesq approximation: density varies linearly with temperature, creating
buoyancy forces that drive a recirculating flow.

This is **true two-way coupling** — temperature drives buoyancy → buoyancy drives
flow → flow advects temperature. The problem is parameterized by two dimensionless
numbers:

```
Rayleigh number:  Ra = g * alpha * dT * L^3 / (nu * kappa) = 10000
Prandtl number:   Pr = nu / kappa = 0.71  (air)
```

The benchmark result by de Vahl Davis (1983) gives an average Nusselt number
Nu ≈ 2.243 at Ra=10⁴.

### Input File

Save as `case16_natural_convection.i`:

```text
# ============================================================
# Case 16: Natural Convection — Buoyancy-Driven Flow
# Incompressible Navier-Stokes + energy, Boussinesq approximation
# Ra = 10000, Pr = 0.71 (air)
# ============================================================

nu    = 0.008426   # = mu  since rho = 1
kappa = 0.011867   # = k   since rho*cp = 1

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim  = 2
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
    nx   = 25
    ny   = 25
  []
[]

[Modules]
  [NavierStokesFV]
    compressibility          = 'incompressible'
    porous_medium_treatment  = false
    add_energy_equation      = true

    density              = 'rho'
    dynamic_viscosity    = 'mu'
    thermal_conductivity = 'k'
    specific_heat        = 'cp'

    initial_velocity    = '1e-15 1e-15 0'
    initial_pressure    = 0.0
    initial_temperature = 0.5

    wall_boundaries       = 'left right top bottom'
    momentum_wall_types   = 'noslip noslip noslip noslip'

    energy_wall_types    = 'fixed-temperature fixed-temperature heatflux heatflux'
    energy_wall_functors = '1 0 0 0'

    boussinesq_approximation = true
    gravity                  = '0 -1 0'
    ref_temperature          = 0.5
    thermal_expansion        = 'alpha'

    pin_pressure         = true
    pinned_pressure_type = average
    pinned_pressure_value = 0

    mass_advection_interpolation     = 'average'
    momentum_advection_interpolation = 'upwind'
    energy_advection_interpolation   = 'upwind'
  []
[]

[FunctorMaterials]
  [const]
    type        = ADGenericFunctorMaterial
    prop_names  = 'rho  mu       k         cp  alpha'
    prop_values = '1.0  ${nu}   ${kappa}  1.0  1.0'
  []
[]

[Postprocessors]
  [max_vel_x]
    type    = ADElementExtremeFunctorValue
    functor = vel_x
  []
  [max_vel_y]
    type    = ADElementExtremeFunctorValue
    functor = vel_y
  []
  [avg_T]
    type     = ElementAverageValue
    variable = T_fluid
  []
[]

[Executioner]
  type       = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 30
  automatic_scaling = true
[]

[Outputs]
  exodus = true
  csv    = true
[]
```

### How to Run

```bash
combined-opt -i case16_natural_convection.i
```

### What to Expect

The solver converges to a steady-state flow with a single clockwise recirculation
cell. Hot fluid rises along the left wall, crosses the top, descends along the cold
right wall, and returns along the bottom. The temperature field shows thermal
boundary layers near the hot and cold walls, with a nearly stratified interior.

The average temperature should remain near 0.5 (energy conservation). Like Case 15,
all variables are element-centered (FV).

Output files: `case16_natural_convection_out.e`, `case16_natural_convection_out.csv`.

---

## Case 17: Joule Heating — Electric Current Generates Heat

### Physics

When electric current flows through a conductor, electrons collide with the lattice
and dissipate kinetic energy as heat. This is Joule heating (also called resistive
heating or ohmic heating). The case solves two coupled equations:

```
Electric potential:  -div( sigma * grad(V) ) = 0          (Laplace equation)
Heat equation:       rho*cp * dT/dt = div(k*grad(T)) + sigma*|grad(V)|^2
```

The Joule heat source Q = sigma * |grad(V)|^2 couples the electric field into the
thermal equation. This is a **one-way coupling** — the temperature does not affect
the electrical conductivity in this simplified model.

A new technique appears here: `ADHeatConduction` is reused for the voltage equation
by pointing its `thermal_conductivity` parameter at `electrical_conductivity`. This
works because both equations have the same mathematical form: -div(k * grad(u)) = 0.

### Input File

Save as `case17_joule_heating.i`:

```text
# ============================================================
# Case 17: Joule Heating — Electric Current Generates Heat
# -div(sigma * grad(V)) = 0          (voltage)
# rho*cp * dT/dt = div(k*grad(T)) + sigma*|grad(V)|^2  (heat)
# Domain: 2D rectangle, x in [0,2], y in [0,1]
# V = 10 V on left, V = 0 V on right
# T = 300 K on left and right (electrodes as heat sinks)
# ============================================================

[Mesh]
  [gen]
    type  = GeneratedMeshGenerator
    dim   = 2
    nx    = 40
    ny    = 20
    xmin  = 0
    xmax  = 2
    ymin  = 0
    ymax  = 1
  []
[]

[Variables]
  [V]
    initial_condition = 0.0
  []
  [T]
    initial_condition = 300.0
  []
[]

[Kernels]
  # Reuse ADHeatConduction for voltage: same math, different material property
  [V_diff]
    type                 = ADHeatConduction
    variable             = V
    thermal_conductivity = electrical_conductivity
  []
  [T_time]
    type     = ADHeatConductionTimeDerivative
    variable = T
  []
  [T_diff]
    type     = ADHeatConduction
    variable = T
  []
  # Joule heating source: Q = sigma * |grad(V)|^2
  [T_joule]
    type         = ADJouleHeatingSource
    variable     = T
    heating_term = electric_field_heating
  []
[]

[BCs]
  [V_left]
    type     = ADDirichletBC
    variable = V
    boundary = left
    value    = 10.0
  []
  [V_right]
    type     = ADDirichletBC
    variable = V
    boundary = right
    value    = 0.0
  []
  [T_left]
    type     = ADDirichletBC
    variable = T
    boundary = left
    value    = 300.0
  []
  [T_right]
    type     = ADDirichletBC
    variable = T
    boundary = right
    value    = 300.0
  []
[]

[Materials]
  [thermal]
    type        = ADGenericConstantMaterial
    prop_names  = 'thermal_conductivity specific_heat density'
    prop_values = '50.0              500.0        8000.0'
  []
  [electrical]
    type        = ADGenericConstantMaterial
    prop_names  = 'electrical_conductivity'
    prop_values = '1e6'
  []
  [joule_material]
    type                       = ElectromagneticHeatingMaterial
    electric_field             = V
    electric_field_heating_name = electric_field_heating
    electrical_conductivity    = electrical_conductivity
    formulation                = time
    solver                     = electrostatic
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
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
  dt       = 0.25
  end_time = 5.0
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
combined-opt -i case17_joule_heating.i
```

### What to Expect

The voltage field establishes immediately (Laplace equation, no time dependence) as
a linear ramp from 10 V (left) to 0 V (right). The current density is uniform:
J = sigma * dV/dx = 1e6 * 10/2 = 5e6 A/m^2. The Joule source Q = sigma * (dV/dx)^2
= 1e6 * 25 = 25 MW/m^3 heats the conductor uniformly.

Temperature rises from 300 K toward a steady state where Joule input balances
conduction to the cold electrodes. The maximum temperature occurs at the center of
the conductor (x=1), forming a parabolic profile symmetric about the midpoint.

Output files: `case17_joule_heating_out.e`, `case17_joule_heating_out.csv`.

---

## Case 18: Cahn-Hilliard Spinodal Decomposition

### Physics

The Cahn-Hilliard equation models **phase separation** in a binary mixture (alloy,
polymer blend, or any two-component system). Starting from a nearly uniform
composition c = 0.5 with small random perturbations, the system spontaneously
separates into two distinct phases (c ≈ 0 and c ≈ 1). This happens because the
free energy F(c) = c^2 * (1-c)^2 has a double-well shape with an unstable region
between the two minima.

The equation is fourth-order, so it is split into two coupled second-order equations
to allow standard linear (C0) finite elements:

```
dc/dt = div( M * grad(w) )              (mass transport)
w     = dF/dc - kappa * laplacian(c)    (chemical potential)
```

Periodic boundary conditions eliminate edge effects. The gradient energy coefficient
kappa controls the interface width between phases.

### Input File

Save as `case18_cahn_hilliard.i`:

```text
# ============================================================
# Case 18: Cahn-Hilliard Spinodal Decomposition
# Split form: dc/dt = div(M*grad(w)), w = dF/dc - kappa*lap(c)
# Free energy: F(c) = c^2*(1-c)^2, periodic BCs
# Domain: [0,25] x [0,25], 40x40 mesh
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 40
  ny   = 40
  xmin = 0
  xmax = 25
  ymin = 0
  ymax = 25
[]

[Variables]
  [c]
    order  = FIRST
    family = LAGRANGE
    [InitialCondition]
      type = RandomIC
      min  = 0.44
      max  = 0.56
    []
  []
  [w]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  # Equation 1 (on w): dc/dt - div(M*grad(w)) = 0
  [c_dot]
    type     = CoupledTimeDerivative
    variable = w
    v        = c
  []
  [w_res]
    type     = SplitCHWRes
    variable = w
    mob_name = M
  []

  # Equation 2 (on c): w - dF/dc + kappa*laplacian(c) = 0
  [c_res]
    type       = SplitCHParsed
    variable   = c
    f_name     = F
    kappa_name = kappa_c
    w          = w
  []
[]

[BCs]
  [Periodic]
    [all]
      auto_direction = 'x y'
    []
  []
[]

[Materials]
  [free_energy]
    type             = DerivativeParsedMaterial
    property_name    = F
    coupled_variables = 'c'
    expression       = 'c^2*(1-c)^2'
    derivative_order = 2
    disable_fpoptimizer = true
    enable_jit          = false
  []
  [const]
    type        = GenericConstantMaterial
    prop_names  = 'kappa_c M'
    prop_values = '1.0     1.0'
  []
[]

[Postprocessors]
  [avg_c]
    type       = ElementAverageValue
    variable   = c
    execute_on = 'initial timestep_end'
  []
  [bulk_energy]
    type          = ElementIntegralMaterialProperty
    mat_prop      = F
    execute_on    = 'initial timestep_end'
  []
[]

[Executioner]
  type       = Transient
  solve_type = 'NEWTON'
  scheme     = bdf2
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       mumps'
  nl_max_its = 30
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-11
  end_time   = 100.0
  [TimeStepper]
    type           = IterationAdaptiveDT
    dt             = 0.1
    growth_factor  = 1.2
    cutback_factor = 0.5
    optimal_iterations = 8
  []
[]

[Outputs]
  exodus = true
  csv    = true
  [checkpoint]
    type      = Checkpoint
    num_files = 2
  []
[]
```

### How to Run

```bash
combined-opt -i case18_cahn_hilliard.i
```

### What to Expect

The simulation starts from near-uniform c ≈ 0.5 with random noise. Within the first
few time steps, the spinodal instability amplifies the perturbations. By t ≈ 10,
distinct domains of c ≈ 0 and c ≈ 1 have formed. Over time, these domains coarsen —
small regions shrink and large ones grow to reduce interface energy.

Key checks:
- `avg_c` should remain at ~0.5 throughout (mass conservation)
- `bulk_energy` should decrease monotonically (thermodynamic consistency)
- The `DerivativeParsedMaterial` requires `disable_fpoptimizer = true` and
  `enable_jit = false` when running in Docker (no `mpicxx` available)

Output files: `case18_cahn_hilliard_out.e`, `case18_cahn_hilliard_out.csv`.

---

## Case 19: Darcy Flow with Heat Transport in Porous Media

### Physics

This case models **thermo-hydro (TH) coupling** in a saturated porous medium. A
pressure gradient drives groundwater flow through porous rock (Darcy's law), and hot
fluid injected at the inlet carries heat downstream (advection-diffusion). This is
relevant to geothermal energy extraction, groundwater contamination, and CO2 storage.

```
Mass balance:  d(rho*phi)/dt + div(rho * q) = 0
               q = -k/mu * grad(P)                        (Darcy velocity)
Energy:        (rho*cp)_eff * dT/dt + rho_f*cp_f * q . grad(T) = div(lambda*grad(T))
```

The `PorousFlowBasicTHM` action automatically creates all kernels, time derivatives,
and the PorousFlow material hierarchy. This drastically reduces input file complexity
compared to specifying each kernel individually.

### Input File

Save as `case19_porous_flow.i`:

```text
# ============================================================
# Case 19: Darcy Flow with Heat Transport in Porous Media
# PorousFlowBasicTHM action: coupled thermo-hydro
# Domain: 3 m x 2 m, pressure-driven flow, hot injection
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 30
  ny   = 20
  xmin = 0
  xmax = 3
  ymin = 0
  ymax = 2
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[Variables]
  [porepressure]
    initial_condition = 1e6
  []
  [temperature]
    initial_condition = 300
    scaling = 1e-6
  []
[]

[PorousFlowBasicTHM]
  porepressure    = porepressure
  temperature     = temperature
  coupling_type   = ThermoHydro
  gravity         = '0 0 0'
  fp              = simple_fluid
  multiply_by_density = true
[]

[FluidProperties]
  [simple_fluid]
    type                = SimpleFluidProperties
    density0            = 1000
    viscosity           = 0.001
    thermal_expansion   = 0
    cp                  = 4186
    cv                  = 4186
    thermal_conductivity = 0.6
  []
[]

[Materials]
  [porosity]
    type          = PorousFlowPorosity
    porosity_zero = 0.3
    mechanical    = false
    thermal       = false
    fluid         = false
  []
  [biot_modulus]
    type                 = PorousFlowConstantBiotModulus
    biot_coefficient     = 1.0
    solid_bulk_compliance = 1e-10
    fluid_bulk_modulus   = 2e9
  []
  [thermal_expansion]
    type                = PorousFlowConstantThermalExpansionCoefficient
    biot_coefficient    = 1.0
    drained_coefficient = 0.0
    fluid_coefficient   = 0.0
  []
  [permeability]
    type        = PorousFlowPermeabilityConst
    permeability = '1e-12 0 0  0 1e-12 0  0 0 1e-12'
  []
  [rock_heat]
    type               = PorousFlowMatrixInternalEnergy
    density            = 2600
    specific_heat_capacity = 800
  []
  [thermal_conductivity]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '2.0 0 0  0 2.0 0  0 0 2.0'
  []
[]

[BCs]
  [p_inlet]
    type     = DirichletBC
    variable = porepressure
    boundary = left
    value    = 1.1e6
  []
  [p_outlet]
    type     = DirichletBC
    variable = porepressure
    boundary = right
    value    = 1.0e6
  []
  [T_inlet]
    type     = DirichletBC
    variable = temperature
    boundary = left
    value    = 350
  []
  [T_outlet]
    type     = DirichletBC
    variable = temperature
    boundary = right
    value    = 300
  []
[]

[Postprocessors]
  [avg_T]
    type     = ElementAverageValue
    variable = temperature
  []
  [max_T]
    type = ElementExtremeValue
    variable = temperature
  []
[]

[Executioner]
  type       = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       mumps'
  dt       = 100
  end_time = 5000
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
combined-opt -i case19_porous_flow.i
```

### What to Expect

The pressure field establishes quickly as a linear gradient from 1.1 MPa (left) to
1.0 MPa (right). The Darcy velocity is approximately q = k/mu * dP/dx =
1e-12 / 0.001 * (0.1e6 / 3) ≈ 3.3e-8 m/s.

Hot fluid (350 K) enters from the left and slowly displaces the cool ambient fluid
(300 K). The thermal front advances slower than the fluid because heat is also stored
in the rock matrix (thermal retardation). By t = 5000 s, a thermal plume extends
partway across the domain. The `avg_T` postprocessor rises gradually as the domain
warms.

Output files: `case19_porous_flow_out.e`, `case19_porous_flow_out.csv`.

---

## Case 20: Elastic Wave Propagation in a Bar

### Physics

This case introduces **dynamic solid mechanics** — the equation of motion includes
inertia (mass times acceleration), making the problem truly time-dependent rather
than quasi-static. A short pressure pulse applied to the left end of an elastic bar
generates a compressive stress wave that propagates at the longitudinal wave speed:

```
c = sqrt(E / rho) = sqrt(200e9 / 7800) ≈ 5064 m/s
```

The wave reaches the free right end at t ≈ L/c ≈ 2 ms, where it reflects with a
sign reversal (compression becomes tension) and the displacement doubles momentarily.
The reflected wave then travels back to the left.

Newmark-beta time integration (beta=0.25, gamma=0.5 — the trapezoidal rule) provides
unconditionally stable implicit time stepping for this dynamic problem.

### Input File

Save as `case20_elastic_wave.i`:

```text
# ============================================================
# Case 20: Elastic Wave Propagation in a Bar
# Dynamic solid mechanics with Newmark-beta time integration
# Thin 2D strip (100x5 elements) approximates a 1D bar
# Pressure pulse at left, free right end (wave reflection)
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 100
  ny   = 5
  xmin = 0.0
  xmax = 10.0
  ymin = 0.0
  ymax = 0.5
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

# The Dynamic action sets up displacement variables, velocity/acceleration
# AuxVariables, StressDivergence kernels, and InertialForce kernels
# with Newmark-beta integration.
[Physics]
  [SolidMechanics]
    [Dynamic]
      [all]
        add_variables  = true
        newmark_beta   = 0.25
        newmark_gamma  = 0.5
        strain         = SMALL
        density        = 7800
        generate_output = 'stress_xx stress_yy vonmises_stress'
      []
    []
  []
[]

[BCs]
  # Pressure pulse on the left face
  [Pressure]
    [pulse_left]
      boundary = left
      function = pressure_pulse
      factor = 1
    []
  []
  [fix_y_bottom]
    type     = DirichletBC
    variable = disp_y
    boundary = bottom
    value    = 0.0
  []
  [fix_y_top]
    type     = DirichletBC
    variable = disp_y
    boundary = top
    value    = 0.0
  []
  [fix_y_left]
    type     = DirichletBC
    variable = disp_y
    boundary = left
    value    = 0.0
  []
  # Right end is FREE — no BC on disp_x
[]

# 1 MPa trapezoidal pulse: ramp up in 0.2 ms, hold 0.2 ms, ramp down 0.2 ms
[Functions]
  [pressure_pulse]
    type = PiecewiseLinear
    x = '0.0    0.0002  0.0004  0.0006'
    y = '0.0    1.0e6   1.0e6   0.0'
  []
[]

[Materials]
  [elasticity]
    type           = ComputeIsotropicElasticityTensor
    youngs_modulus = 200.0e9
    poissons_ratio = 0.0
  []
  [stress]
    type = ComputeLinearElasticStress
  []
[]

[Postprocessors]
  [disp_x_right]
    type  = PointValue
    variable = disp_x
    point = '10.0 0.25 0'
  []
  [disp_x_left]
    type  = PointValue
    variable = disp_x
    point = '0.0 0.25 0'
  []
  [avg_stress_xx]
    type     = ElementAverageValue
    variable = stress_xx
  []
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'
  dt       = 2.0e-5
  end_time = 0.006
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-12
[]

[Outputs]
  csv = true
  [exodus]
    type     = Exodus
    time_step_interval = 5
  []
[]
```

### How to Run

```bash
combined-opt -i case20_elastic_wave.i
```

### What to Expect

The simulation runs 300 time steps (dt = 20 us, end_time = 6 ms). The pressure
pulse creates a compressive wave packet that takes about 2 ms to traverse the 10 m
bar. Key events visible in the `disp_x_right` postprocessor:

- t ≈ 0 - 0.4 ms: Pulse applied at left end
- t ≈ 2 ms: Wave arrives at free right end — displacement doubles momentarily
- t ≈ 4 ms: Reflected wave returns to left end
- t ≈ 6 ms: Second reflection from left boundary

The `[exodus]` sub-block writes every 5th step using `time_step_interval = 5` to
keep the output file manageable.

Output files: `case20_elastic_wave_exodus.e` (note: named sub-block, not `_out.e`),
`case20_elastic_wave_out.csv`.

---

## Case 21: Bimetallic Strip — Differential Thermal Expansion

### Physics

Two metal strips bonded together (steel on the bottom, aluminum on top) are heated
uniformly from 300 K to 500 K. Because aluminum has a higher coefficient of thermal
expansion (23e-6 /K vs 12e-6 /K for steel), it tries to expand more than the steel.
The bond constrains both materials to deform together, producing internal stresses
and causing the strip to bend — this is the classic bimetallic thermostat mechanism.

This case combines the multi-material technique from Case 6 with the thermoelasticity
from Case 14. `SubdomainBoundingBoxGenerator` splits the mesh into two blocks, and
block-restricted materials give each subdomain different elastic and thermal properties.

### Input File

Save as `case21_bimetallic_strip.i`:

```text
# ============================================================
# Case 21: Thermo-Mechanical Bimetallic Strip
# Steel (bottom, block 0): alpha=12e-6, E=200 GPa, nu=0.3
# Aluminum (top, block 1): alpha=23e-6, E=70 GPa, nu=0.33
# Uniform heating from T_ref=300K to T=500K, pinned at left
# ============================================================

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 40
    ny   = 8
    xmin = 0
    xmax = 10
    ymin = 0
    ymax = 1
  []
  # Reassign top half (y > 0.5) to block 1 (aluminum)
  [top_block]
    type        = SubdomainBoundingBoxGenerator
    input       = gen
    bottom_left = '0 0.5 0'
    top_right   = '10 1.0 0'
    block_id    = 1
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        strain              = SMALL
        add_variables       = true
        generate_output     = 'vonmises_stress stress_xx stress_yy'
        eigenstrain_names   = 'thermal_eigenstrain'
      []
    []
  []
[]

# Temperature is prescribed (not solved), set to 500 K everywhere
[AuxVariables]
  [T]
    initial_condition = 500
  []
[]

[Materials]
  # Steel (block 0)
  [elasticity_steel]
    type           = ComputeIsotropicElasticityTensor
    block          = 0
    youngs_modulus = 200e9
    poissons_ratio = 0.3
  []
  [thermal_expansion_steel]
    type                    = ComputeThermalExpansionEigenstrain
    block                   = 0
    eigenstrain_name        = thermal_eigenstrain
    stress_free_temperature = 300
    thermal_expansion_coeff = 12e-6
    temperature             = T
  []
  [stress_steel]
    type  = ComputeLinearElasticStress
    block = 0
  []

  # Aluminum (block 1)
  [elasticity_aluminum]
    type           = ComputeIsotropicElasticityTensor
    block          = 1
    youngs_modulus = 70e9
    poissons_ratio = 0.33
  []
  [thermal_expansion_aluminum]
    type                    = ComputeThermalExpansionEigenstrain
    block                   = 1
    eigenstrain_name        = thermal_eigenstrain
    stress_free_temperature = 300
    thermal_expansion_coeff = 23e-6
    temperature             = T
  []
  [stress_aluminum]
    type  = ComputeLinearElasticStress
    block = 1
  []
[]

[BCs]
  [fix_x_left]
    type     = DirichletBC
    variable = disp_x
    boundary = left
    value    = 0
  []
  [fix_y_left]
    type     = DirichletBC
    variable = disp_y
    boundary = left
    value    = 0
  []
[]

[Postprocessors]
  [tip_disp_y]
    type     = PointValue
    variable = disp_y
    point    = '10 0.5 0'
  []
  [max_vonmises]
    type     = ElementExtremeValue
    variable = vonmises_stress
  []
[]

[Executioner]
  type       = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Outputs]
  exodus = true
  csv    = true
[]
```

### How to Run

```bash
combined-opt -i case21_bimetallic_strip.i
```

### What to Expect

The steady solver converges in a few Newton iterations (the problem is linear).
The strip bends downward because the aluminum top layer (alpha = 23e-6 /K) expands
more than the steel bottom layer (alpha = 12e-6 /K), pushing the free right tip
downward. The `tip_disp_y` postprocessor gives the vertical deflection at the
right end midline.

Von Mises stress is highest near the bonded interface (y = 0.5) and near the fixed
left edge, where the differential expansion is most constrained. The deformed shape
resembles a cantilever beam under a distributed load.

Output files: `case21_bimetallic_strip_out.e` (contains fields for both blocks:
disp_x, disp_y, vonmises_stress, stress_xx, stress_yy), `case21_bimetallic_strip_out.csv`.

---

## Case 22: Charge Relaxation in an Ohmic Medium

### Physics

When free electric charge is deposited inside a conducting medium it does not stay
there — conduction currents sweep it to the surfaces on a characteristic time scale
called the **charge relaxation time** tau_e = epsilon / sigma (Melcher,
*Continuum Electromechanics*, MIT Press 1981, §5.9). For a good conductor this time
is nanoseconds; for a resistive dielectric liquid it can be milliseconds to seconds.

Two coupled equations govern the field problem:

```
dρ_e/dt + (σ/ε)·ρ_e = 0        (charge relaxation ODE at every point)
-div(ε·grad(φ))      = ρ_e     (Poisson's equation for the potential)
```

The first equation has the exact solution rho_e(x,y,t) = rho_e(x,y,0)·exp(−t/tau_e).
Every spatial point decays at the same rate, so the Gaussian blob initial condition
simply shrinks in amplitude without changing shape.

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `ADTimeDerivative` | drho_e/dt term in the relaxation equation |
| `ADReaction` | +(sigma/epsilon)·rho_e linear volumetric decay |
| `ADHeatConduction` | Repurposed as −div(epsilon·grad(phi)) for Poisson |
| `ADCoupledForce` | Injects rho_e as source on the RHS of Poisson's equation |
| `ADGenericConstantMaterial` | Permittivity property read by ADHeatConduction |
| `SMP full = true` | Builds the full coupled Jacobian (rho_e–phi off-diagonal blocks) |

The HIT top-level variable `sigma_over_eps = 10.0` sets the decay rate in a single
place. The natural (zero-flux Neumann) condition on rho_e is applied automatically
when no BC block is provided — the physically correct choice for charge that decays
in place by Ohmic conduction.

### What Makes This Case Interesting

This case shows that `ADHeatConduction` is a general Laplacian operator, not just a
"heat" kernel. Mapping `thermal_conductivity -> permittivity` repurposes it for
Poisson's equation with no new code. It also introduces `ADReaction` as the standard
pattern for any linear volumetric sink or source term, and demonstrates one-way
coupling via `ADCoupledForce`.

### Expected Results

Running on a 30x30 mesh from t = 0 to t = 0.5 s (five relaxation times) with
dt = 0.01 s, the `max_rho_e` postprocessor follows exp(−10·t) precisely. A log-linear
plot of max_rho_e vs. time is a perfect straight line with slope −10 s⁻¹. The
electric potential phi mirrors the charge decay at the same exponential rate.

---

## Case 23: Magnetic Diffusion into a Conducting Slab

### Physics

When a magnetic field is suddenly applied to the surface of a conductor, eddy currents
near the surface shield the interior. The field penetrates progressively deeper over
time — a process called **magnetic diffusion** (Melcher §6.2–6.3). The governing
equation is:

```
∂B/∂t = D_m · ∂²B/∂x²        D_m = 1/(μ₀·σ)   [magnetic diffusivity]
```

This is mathematically identical to the heat equation. With D_m = 1/(μ₀σ) playing
the role of thermal diffusivity, the MOOSE setup is a direct reuse of the transient
diffusion pattern from Cases 3–4. The exact solution for a step-applied surface field
on a semi-infinite slab is B(x,t) = erfc(x / (2√(D_m·t))), and the penetration
depth scales as delta ~ 2√(D_m·t).

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `ADTimeDerivative` | ∂B/∂t (magnetic flux density time rate) |
| `ADMatDiffusion` | D_m·∇²B with diffusivity read from material |
| `ADGenericConstantMaterial` | Provides `diffusivity = D_m` (AD-compatible) |
| BoomerAMG preconditioner | Efficient AMG solver for scalar diffusion |

The HIT top-level variable `D_m = 0.01` is referenced by the material block so that
changing the magnetic diffusivity requires editing one line. A 50x2 mesh on a
1 m × 0.04 m strip makes the domain quasi-1D while remaining formally 2D.

### What Makes This Case Interesting

The analogy between magnetic diffusion and heat conduction is one of the most useful
conceptual tools in continuum electromechanics. This case makes the analogy concrete:
the same MOOSE kernel (`ADMatDiffusion`), the same material (`ADGenericConstantMaterial`),
and the same boundary conditions produce the correct magnetic field penetration profile.
The skin-effect physics — exponential penetration with depth — is visible directly in
the Exodus output.

### Expected Results

The B field advances from the left boundary (B = 1) toward the right (B = 0) following
the erfc profile. The `avg_B` postprocessor increases from 0 toward the steady-state
value of 0.5 as flux fills the slab. At t = 5 the analytical average is approximately
0.108; MOOSE gives approximately 0.107, with the slight underestimate due to backward
Euler temporal smoothing.

---

## Case 24: Charge Drift-Diffusion Between Parallel Plates

### Physics

Positive ions are injected at the left electrode of a parallel-plate gap and drift
rightward under an applied electric field while also spreading by Fickian diffusion
(Melcher §5.5–5.7). The drift-dominated transport at Peclet number Pe = 100 creates a
near-step charge front that advances at the drift speed v = mu_i·E:

```
∂ρ_e/∂t + div(v·ρ_e) = D_i·∇²ρ_e        (drift-diffusion, conservative form)
-div(ε·grad(φ))       = ρ_e              (Poisson's equation)
```

As charge accumulates it distorts the electric potential away from the linear Laplace
solution, with a parabolic correction at steady state.

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `ConservativeAdvection` | Conservative div(v·rho_e) with full upwinding |
| `ADMatDiffusion` | Fickian diffusion D_i·∇²rho_e |
| `ADHeatConduction` | Poisson equation for the electric potential |
| `ADCoupledForce` | Charge density as source on Poisson RHS |
| `GenericConstantVectorMaterial` | Prescribed constant drift velocity vector |
| `SMP full = true` | Coupled rho_e–phi Jacobian |

Full upwinding (`upwinding_type = full`) is essential at Pe = 100 to prevent spurious
oscillations at the sharp charge front. The `GenericConstantVectorMaterial` (non-AD)
supplies the constant velocity vector to the non-AD `ConservativeAdvection` kernel.

### What Makes This Case Interesting

This is the prototype for all self-consistent field transport problems in
electromechanics, plasma physics, and electrochemistry. Solving Poisson and the
transport equation simultaneously in one Newton iteration is more robust than operator
splitting. The case also demonstrates the important distinction between
drift-dominated and diffusion-dominated regimes through the Peclet number, and
introduces upwinding as the stabilization strategy for hyperbolic transport.

### Expected Results

The charge front advances at approximately v = 1 m/s, reaching the cathode (x = 1)
at t ~ 1 s. The `avg_rho` postprocessor rises from 0 toward approximately 0.5 as the
gap fills. The `phi` field bends away from the linear Laplace solution as space charge
accumulates; the midpoint potential exceeds the linear-interpolated value by
approximately 0.125 V at steady state (the parabolic space-charge correction).

---

## Case 25: Induction Heating — Magnetic Diffusion Generates Heat

### Physics

An alternating magnetic field applied to the surface of a conducting slab drives eddy
currents in the skin-depth layer, and those currents dissipate Joule heat
(Melcher §6.7–6.8). The two coupled equations are:

```
∂B/∂t = D_m · ∂²B/∂x²                    (magnetic diffusion)
ρcp·∂T/∂t = k·∂²T/∂x² + Q_coeff·(∂B/∂x)²  (heat equation with eddy source)
```

The heating power density Q ~ (∂B/∂x)² is exponentially concentrated within one skin
depth delta = sqrt(2·D_m/omega) of the surface. With D_m = 0.005 and oscillation
period tau = 0.5 s, the skin depth is approximately 0.028 m, so interior regions
(x > 0.1 m) see negligible heating.

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `ADTimeDerivative` + `ADMatDiffusion` | Magnetic diffusion equation for B |
| `ADHeatConductionTimeDerivative` + `ADHeatConduction` | Transient heat equation for T |
| `VariableGradientComponent` (AuxKernel) | Extracts ∂B/∂x into a CONSTANT MONOMIAL AuxVariable |
| `ParsedAux` | Computes Q_coeff·(dBdx)² from the gradient AuxVariable |
| `CoupledForce` | Injects the eddy heat AuxVariable as a source in the T equation |
| `ADFunctionDirichletBC` | Oscillating BC: B = sin(2*pi*t/tau) at the left surface |

The AuxVariable pair (`dBdx`, `eddy_heat`) implements the coupling with a one-timestep
lag: the gradient is extracted and squared at TIMESTEP_END, then used as a forcing
term at the next step. This lagged pattern avoids a strongly nonlinear system at the
cost of first-order accuracy in the coupling, which is acceptable for slowly varying
heating.

### What Makes This Case Interesting

This case introduces the `VariableGradientComponent` + `ParsedAux` pattern for
extracting derived quantities from primary variable gradients. It also demonstrates
how two previously independent physics — magnetic diffusion (Case 23) and heat
conduction (Case 17) — are composed into a genuine multi-physics problem by bridging
them with AuxVariables, without writing any new C++ code.

### Expected Results

The B field oscillates at the surface with unit amplitude and decays exponentially
with depth. The eddy_heat AuxVariable oscillates at twice the field frequency
(it goes as sin²). Temperature at the left surface rises by approximately 1–2 K per
oscillation period; the interior stays near 300 K. After 10 periods (t = 5 s) the
peak temperature is approximately 314 K.

---

## Case 26: EHD Pumping — Coulomb Force Drives Fluid Flow

### Physics

In an EHD pump, a prescribed Coulomb body force f = rho_e·E acts on a dielectric
liquid, driving recirculating flow without any moving parts (Melcher §9.11–9.12).
This case prescribes the force analytically:

```
f_x(x, y) = A·(1 − x)·sin(π·y)        A = 500
f_y = 0
```

The factor (1−x) represents the charge density decaying from the injection electrode
(x = 0) to the collecting electrode (x = 1). The sin(pi·y) variation ensures zero
force at the no-slip walls and maximum force at mid-height. This drives a single
large recirculating loop: rightward in the interior, returning along the walls.

The governing equations are the incompressible Navier-Stokes equations plus the body
force:

```
rho·(v · grad)v = −grad p + mu·lap(v) + f_Coulomb(x, y)
div(v) = 0
```

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `[Modules/NavierStokesFV]` action | Creates all FV NS variables and kernels automatically |
| `INSFVBodyForce` | Adds the Coulomb body force to the x-momentum FV equation |
| `ADParsedFunctorMaterial` | Defines the analytical force expression f_x = A(1-x)sin(pi*y) |
| `ADGenericFunctorMaterial` | Supplies constant fluid properties rho and mu |
| `pin_pressure = true` | Removes pressure null space in the closed cavity |

The `rhie_chow_user_object` parameter of `INSFVBodyForce` must reference the Rhie-Chow
interpolation object created by the action (`ins_rhie_chow_interpolator`). Omitting
this causes pressure-velocity decoupling and incorrect results.

### What Makes This Case Interesting

EHD pumping is the electromagnetic analog of natural convection: a body force that
depends on a field quantity drives recirculating flow in a closed cavity. This case
demonstrates the pattern for adding any custom body force to an action-based NS setup
— a pattern that applies equally to Lorentz forces, dielectrophoretic forces, and
other electromechanical body forces. The prescribed-force approach separates the
fluid-mechanics coupling from the charge-transport problem.

### Expected Results

The converged steady solution shows a single clockwise recirculation loop with the
strongest rightward flow in the interior (near x = 0) and return flow along the walls.
The peak x-velocity is in the range 1–10 (dimensionless), much less than the Stokes
estimate because inertia limits the driven velocity. The pressure is higher near the
right wall (where the body force pushes fluid into the wall) and lower near the
left wall.

---

## Case 27: MHD Hartmann Flow — Magnetic Braking of Channel Flow

### Physics

An electrically conducting fluid flows through a channel under a uniform transverse
magnetic field B0. The fluid motion generates induced currents j = sigma*(v x B0),
and the resulting Lorentz force j x B0 opposes the flow (Melcher §9.9–9.10). This
braking redistributes momentum and flattens the velocity profile from the parabolic
Poiseuille shape into the characteristic flat-topped **Hartmann profile**:

```
0 = −dp/dx + mu·d²v_x/dy² − sigma·B0²·v_x
```

The Hartmann number Ha = B0·L·sqrt(sigma/mu) characterises the ratio of
electromagnetic to viscous forces. At Ha = 5 (used here, Ha² = 25), the flat core
and thin boundary layers (thickness ~ L/Ha = 0.2) are clearly visible.

The key implementation insight is the Darcy-Lorentz equivalence: the Lorentz drag
−sigma·B0²·v is mathematically identical to Darcy friction −(mu/(rho·K))·v when
K = mu/(rho·sigma·B0²) = 1/Ha². Setting the Darcy coefficient to Ha² in the
`NavierStokesFV` porous medium treatment models MHD braking with no custom kernels.

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `[Modules/NavierStokesFV]` with `porous_medium_treatment = true` | Activates Darcy friction term |
| `friction_types = 'darcy'` | Selects linear drag model (Lorentz equivalent) |
| `ADGenericVectorFunctorMaterial` | Supplies isotropic Darcy coefficient vector [Ha², Ha², Ha²] |
| `ADGenericFunctorMaterial` | Fluid properties rho and mu; porosity = 1.0 |
| Inlet/outlet BCs | Fixed velocity inlet (v_x = 1), fixed pressure outlet (p = 0) |

The HIT variable `Ha2 = 25.0` appears in a single place and propagates to the Darcy
coefficient. The positive Ha² diagonal contribution to the Jacobian improves
conditioning, so Newton typically converges in just 2 iterations.

### What Makes This Case Interesting

This case demonstrates that the porous-medium friction mechanism in MOOSE is a general
linear drag model, not restricted to porous-medium physics. By mapping Lorentz drag
onto Darcy friction, the same action handles MHD channel flow without any custom code.
The Hartmann profile is a canonical benchmark for low-magnetic-Reynolds-number MHD
solvers, and the case introduces open-channel (inlet-outlet) boundary conditions in
contrast to the closed-cavity setups of Cases 15–16.

### Expected Results

In the developed flow region (downstream from the inlet development length) the
cross-section velocity profile shows no-slip at the Hartmann walls, thin boundary
layers of thickness ~0.2, and a flat core at approximately 1.38 — higher than the
mean inlet velocity of 1.0 because mass conservation requires the core to accelerate
to compensate for the slow near-wall fluid. The `max_vel_x` postprocessor is
approximately 1.378; `avg_vel_x` is approximately 1.0 (mass conservation).

---

## Case 28: Two-Way Joule Heating — Temperature-Dependent Conductivity

### Physics

This case extends Case 17 (Joule Heating) by making the electrical conductivity
sigma a function of temperature, creating **two-way coupling**:

```
T changes → σ(T) changes → V distribution changes → Q changes → T changes → ...
```

The metallic conductivity model is sigma(T) = sigma_0 / (1 + alpha·(T − T_ref)),
where alpha = 0.004 /K. As temperature rises, conductivity decreases (phonon
scattering), which reduces Joule heating — a **negative feedback** that self-limits
the temperature rise. This contrasts with semiconductors (dσ/dT > 0), where positive
feedback can lead to thermal runaway (Melcher §10.1–10.3).

The two coupled PDEs are:

```
−div(σ(T)·grad(V)) = 0                              (current conservation)
ρ·cp·∂T/∂t = div(k·grad(T)) + σ(T)·|grad(V)|²     (heat equation with Joule source)
```

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `ADPiecewiseLinearInterpolationMaterial` | Tabulated sigma(T) — AD-compatible, no JIT dependency |
| `ADHeatConduction` | Both the V equation (with sigma as conductivity) and T diffusion |
| `ADJouleHeatingSource` | sigma(T)·|grad(V)|² heat source in the T equation |
| `ADHeatConductionTimeDerivative` | rho·cp·dT/dt |
| `SMP full = true` | Full Jacobian including dR_V/dT and dR_T/dV off-diagonal blocks |

`ADPiecewiseLinearInterpolationMaterial` is used instead of `ADParsedMaterial` because
the Docker MOOSE image lacks the JIT compilation toolchain that parsed materials
require at runtime. The tabulated interpolation is fully AD-compatible — the slope of
each linear segment provides dσ/dT automatically to the Jacobian. The only structural
change from Case 17 is the `[Materials]` block.

### What Makes This Case Interesting

Two-way coupling means all four Jacobian blocks (dR_V/dV, dR_V/dT, dR_T/dV, dR_T/dT)
are nonzero. This case shows how MOOSE handles strong multi-field coupling through the
AD chain: the derivative of sigma(T) with respect to T propagates automatically from
the material to both the V and T residuals, with no hand-coded cross-terms. The
metallic negative-feedback result is visible as a slower temperature rise and lower
peak temperature compared to Case 17.

### Expected Results

At t = 5 s the peak temperature is approximately 334 K, a 34 K rise above the 300 K
electrode BCs. At that temperature the conductivity has dropped to approximately
0.88·sigma_0 — a 12% reduction. The avg_T rises approximately linearly with time
(the overall Joule power changes little over this modest temperature range). The V
field develops a slight spatial non-uniformity as current routes around the
more-resistive hot spot in the centre.

---

## Case 29: Electroconvection — EHD-Enhanced Natural Convection

### Physics

This case extends the differentially heated square cavity of Case 16 (natural
convection) with an additional electrohydrodynamic body force. In a dielectric liquid
whose permittivity depends on temperature, a non-uniform electric field exerts a
dielectrophoretic force on the fluid (Melcher §9.12, §10.4). The leading-order model
reduces to:

```
f_EHD_y = Fe·(T − T_ref)
```

Because this has exactly the same form as the Boussinesq buoyancy force
f_y = alpha·(T − T_ref), the two forces combine into one effective coefficient:

```
f_total_y = (alpha + Fe)·(T − 0.5) = alpha_eff·(T − 0.5)
```

With the default Fe = 5 and alpha = 1, alpha_eff = 6.0. No additional kernels are
needed — the entire EHD contribution is absorbed into the `thermal_expansion`
property consumed by the built-in Boussinesq mechanism of the `NavierStokesFV`
action.

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `[Modules/NavierStokesFV]` with `boussinesq_approximation = true` | Built-in buoyancy mechanism handles both forces |
| `ADGenericFunctorMaterial` with `alpha = alpha_eff` | Combines buoyancy and EHD into one coefficient |
| `add_energy_equation = true` | Activates the coupled thermal transport equation |
| `pin_pressure = true` (average type) | Removes pressure null space in the closed cavity |

The HIT top-level variable `alpha_eff = 6.0` is the only change from Case 16. To
explore different EHD strengths, change the second term: alpha_eff = 1.0 + Fe.
Command-line HIT overrides (`-i case29... alpha_eff=0.5`) allow parameter sweeps
without editing the file.

### What Makes This Case Interesting

The mathematical equivalence between the dielectrophoretic force and Boussinesq
buoyancy means the entire new physics is captured by changing a single number. This
is the payoff of recognising structural analogies between different fields of physics.
The case also demonstrates EHD as a practical control mechanism for convective heat
transfer: increasing alpha_eff strengthens circulation and heat transfer;
decreasing it (Fe < 0) suppresses convection. The EHD forcing number Fe therefore
acts as a parameter that tunes the heat transfer coefficient (Nusselt number).

### Expected Results

With Fe = 5 (alpha_eff = 6.0) the flow is a single counter-clockwise recirculation
cell, qualitatively identical to Case 16 but driven with effectively six times the
buoyancy force. The postprocessors report max_vel_x ~ 0.409, max_vel_y ~ 0.613, and
avg_T = 0.5 (preserved by the (T − 0.5) symmetry of the forcing). Setting
alpha_eff = 1.0 reproduces the Case 16 pure natural-convection result exactly.

---

### Cases 30-36: Electromagnetic Noise and Quantum Optical Measurements

> Inspired by **Professor Herman A. Haus**'s masterful textbook
> [*Electromagnetic Noise and Quantum Optical Measurements*](https://doi.org/10.1007/978-3-642-57250-0) (Springer, 2000) —
> a unified treatment of electromagnetic theory from Maxwell's equations
> through waveguides, resonators, and optical fibers to noise, solitons,
> and quantum measurement. Haus was Institute Professor at MIT and a
> pioneer of laser physics, fiber soliton communication, and coupled
> mode theory.

Haus's book spans 13 chapters, from classical Maxwell theory (Chs 1-5)
through quantum noise and photon statistics (Chs 6-9) to solitons and
squeezing (Chs 10-13). These seven cases draw from the **classical chapters
only** (Chs 1-5 and 10). The quantum chapters — covering the quantum
theory of the electromagnetic field, photon operators, homodyne and
heterodyne detection, squeezed states of the radiation field, and the
quantum theory of solitons and squeezing — describe inherently
quantum-mechanical phenomena (operator commutation relations, vacuum
fluctuations, photon number statistics) that have no classical PDE
representation. A finite-element solver like MOOSE operates on classical
field equations; the quantum chapters are therefore outside its scope.

The classical chapters, however, map directly onto PDE problems that MOOSE
handles naturally: Helmholtz eigenvalue problems (Ch 2), driven resonant
cavities (Ch 3), wave scattering from dielectric interfaces (Ch 1),
coupled-mode beating (Ch 3), spectral relaxation of thermal noise (Ch 5),
dispersive pulse propagation (Ch 4), and nonlinear soliton dynamics (Ch 10).

---

## Case 30: Rectangular Waveguide Cutoff Frequencies

### Physics

Every hollow metallic waveguide has a set of resonant transverse patterns called
modes. Each mode propagates only above a characteristic **cutoff frequency**
set by the waveguide's cross-section geometry. Below cutoff the mode is
evanescent. The cutoff wavenumbers k_c are the eigenvalues of the 2D Helmholtz
equation on the waveguide cross-section with Dirichlet (PEC) boundary conditions
(Haus Ch 2):

```
∇²ψ + k_c² ψ = 0,   ψ = 0 on walls
```

For a rectangular waveguide of width a and height b the analytical eigenvalues
are k_c²(m,n) = (mπ/a)² + (nπ/b)² with m,n = 1,2,3,… for TM modes. A 2:1
aspect ratio (a=2, b=1) gives TM₁₁ at k_c² ≈ 12.337.

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `Diffusion` | Laplacian stiffness matrix ∫∇ψ·∇v dV |
| `CoefReaction` (coefficient=-1, extra_vector_tags='eigen') | Mass matrix for eigenvalue problem |
| `DirichletBC` + `EigenDirichletBC` | PEC walls: ψ=0 in both stiffness and mass systems |
| `PotentialToFieldAux` | Computes E_x, E_y from -∇ψ |
| `VectorPostprocessors/Eigenvalues` | Reports computed eigenvalues |
| `Eigenvalue` executioner | SLEPc/KRYLOVSCHUR eigenvalue solver |

### What Makes This Case Interesting

This is MOOSE's eigenvalue mode — instead of solving Ax=b for a known
right-hand side, it finds values of λ for which Ax=λBx has a non-trivial
solution. The eigenvalue infrastructure is a specialised executioner that
wraps SLEPc. The `extra_vector_tags = 'eigen'` label on CoefReaction tells
the framework which kernel contributions go into the B (mass) matrix versus
the A (stiffness) matrix.

### Expected Results

The first six eigenvalues match the analytical values within ~0.5%:

| Mode | Analytical k_c² | Computed (40×20 mesh) |
|------|------------------|-----------------------|
| TM₁₁ | 12.337 | ~12.36 |
| TM₂₁ | 22.207 | ~22.27 |
| TM₃₁ | 37.011 | ~37.15 |

The mode shape for TM₁₁ shows a single central peak.

---

## Case 31: Driven Resonant Cavity — Frequency Response and Q

### Physics

A resonant cavity excited by an oscillating source produces a large field when
the driving frequency matches an eigenfrequency (resonance) and a weak field
otherwise. The frequency response has a Lorentzian shape characterised by the
quality factor Q. This is the simplest model of cavity resonance from Haus Ch 3.

The time-harmonic Helmholtz equation governs the electric field:

```
∇²E + k²E = −J_source(x,y),   E = 0 on PEC walls
```

where k² = ω²με is the squared wavenumber. When k² approaches an eigenvalue
k_c² from Case 30, the system is nearly singular and the response diverges.

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `Diffusion` | Laplacian ∇²E |
| `CoefReaction` (coefficient = −k²) | Helmholtz reaction term k²E |
| `BodyForce` | Gaussian source at (0.7, 0.4) |
| `DirichletBC` | PEC walls E=0 |
| LU direct solver | Handles near-singular matrix |

### What Makes This Case Interesting

The sign of `CoefReaction` is negative (−k²) because MOOSE's weak-form
convention is: residual = ∫∇E·∇v dV + coefficient·∫E·v dV − ∫J·v dV = 0,
giving strong form −∇²E − k²E = −J, i.e., ∇²E + k²E = −J. A positive
coefficient would give the wrong sign. This is a common pitfall.

### Expected Results

At k²=12.3 (near TM₁₁ resonance at 12.337): large field amplitude matching
the TM₁₁ mode shape — a single central peak. At k²=15.0 (off-resonance):
weak, distorted field. The ratio of max_E at resonance vs off-resonance is
approximately 10-50×.

---

## Case 32: EM Wave Reflection from a Dielectric Slab

### Physics

A plane electromagnetic wave incident on a dielectric slab is partially
reflected and partially transmitted. This is the most fundamental wave
scattering problem in electromagnetics (Haus Ch 1). In the frequency domain
the fields satisfy:

```
d²E/dx² + k₀²εᵣ(x)·E = 0
```

where εᵣ(x) is the position-dependent relative permittivity (4 inside the
slab, 1 in vacuum). The complex field is split into real and imaginary
parts, each satisfying the same equation. An EMRobinBC at the vacuum
boundary imposes the incoming wave and absorbs the scattered wave.

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `Diffusion` | d²E/dx² for E_real and E_imag |
| `ADMatReaction` | −k₀²εᵣ(x)·E reaction term |
| `EMRobinBC` | Port BC — injects incident wave, absorbs scattered |
| `DirichletBC` | PEC wall at x=0 |
| `ReflectionCoefficient` PP | Computes \|R\| from boundary fields |
| `ADGenericFunctionMaterial` | Space-dependent coefficient from ParsedFunction |

### What Makes This Case Interesting

This case uses the electromagnetics module for a genuine wave problem. The
real/imaginary splitting is the frequency-domain alternative to time-domain
wave propagation. The `EMRobinBC` is a first-order absorbing boundary
condition that couples the two field components. The lossless simplification
(no imaginary permittivity) eliminates the `ADMatCoupledForce` terms,
making the two field equations uncoupled except at the port boundary.

### Expected Results

Standing wave pattern in vacuum (wavelength λ₀ = 2π/k₀ ≈ 15 m), shorter
wavelength inside the slab (λ = λ₀/√εᵣ ≈ 7.5 m). The `ReflectionCoefficient`
postprocessor reports |R|. For a single-interface half-space, the Fresnel
formula gives |R| = |(n−1)/(n+1)| = 1/3. The finite slab shows interference
fringes modifying this value.

---

## Case 33: Coupled Resonator Beating — Energy Exchange

### Physics

Two coupled optical resonators exchange energy periodically, a phenomenon
called beating. Haus's coupled mode theory (Ch 3) describes this with two
amplitude equations:

```
∂u/∂t = D·∇²u − γ·u + κ·v
∂v/∂t = D·∇²v − γ·v + κ·u
```

The coupling coefficient κ causes energy to slosh between modes u and v
at the beat frequency 2κ, while the damping rate γ causes the total energy
to decay exponentially. With κ=3.0 and γ=0.5, the beating period is
T = π/κ ≈ 1.05 s and the energy e-folding time is 1/γ = 2 s.

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `ADTimeDerivative` | ∂u/∂t and ∂v/∂t |
| `ADMatDiffusion` | D·∇²u spatial smoothing |
| `ADReaction` | −γ·u linear decay |
| `CoupledForce` | +κ·v cross-coupling |
| `SMP` preconditioning | Full Jacobian for coupled system |

### What Makes This Case Interesting

This is the simplest demonstration of energy exchange between coupled modes.
The spatial diffusion (D=0.01) is weak and primarily smooths the fields;
the dynamics are dominated by the local ODE coupling. Setting D=0 would give
pure ODE beating with analytic solution u(t) = u₀·cos(κt)·exp(−γt).

### Expected Results

The CSV postprocessors show avg_u and avg_v oscillating out of phase with
period ≈1.05 s. The sum avg_u² + avg_v² decays as exp(−2γt) = exp(−t).
By t=3 s the total energy has dropped to ~5% of its initial value.

---

## Case 34: Thermal Noise Relaxation — Fluctuation-Dissipation

### Physics

Thermal noise in a physical system can be modelled as random initial
conditions on a diffusion equation. The Nyquist theorem (Haus Ch 5) relates
the spectral density of thermal fluctuations to the dissipation (damping
rate) of each mode. In this case the "noise" is a random temperature field:

```
∂T/∂t = D·∇²T,   T = T_eq on walls
```

Each spatial Fourier mode (m,n) decays at rate λ_mn = D·π²(m² + n²).
High-frequency noise (small features) decays fast; only the fundamental
mode m=n=1 survives at late times.

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `ADTimeDerivative` | ∂T/∂t |
| `ADMatDiffusion` | D·∇²T |
| `ADDirichletBC` | T=0.5 on all walls |
| `RandomIC` | Random initial T ∈ [0.3, 0.7] |

### What Makes This Case Interesting

The `RandomIC` object sets each nodal value independently from a uniform
distribution. The simulation then acts as a spectral filter: watching the
random speckle smooth out in time is a direct visual demonstration of
how diffusion preferentially damps short wavelengths.

### Expected Results

t=0: random speckle. t=0.1: small features gone. t=0.5: only the
fundamental mode (single smooth bump) remains. t=2.0: nearly uniform
T=0.5. The late-time decay rate of (max_T − 0.5) approaches D·2π² ≈ 1.97.

---

## Case 35: Dispersive Pulse Broadening in an Optical Fiber

### Physics

An optical pulse propagating in a fiber experiences group velocity
dispersion (GVD): different frequency components travel at different
speeds, causing the pulse envelope to broaden. Haus Ch 4 derives the
envelope equation:

```
∂A/∂t + v_g·∂A/∂x = D_gvd·∂²A/∂x²
```

This is mathematically identical to advection-diffusion (Case 08), but
physically describes a pulse translating at group velocity v_g while
broadening due to GVD coefficient D_gvd. A Gaussian pulse of initial
width w₀ broadens as w(t) = w₀√(1 + (v_g·t/z_d)²) where the dispersion
length z_d = w₀²/(2·D_gvd).

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `ADTimeDerivative` | ∂A/∂t |
| `ADConservativeAdvection` | v_g·∂A/∂x group velocity transport |
| `ADMatDiffusion` | D_gvd·∂²A/∂x² dispersive broadening |
| `ADGenericConstantVectorMaterial` | Velocity vector (v_g, 0, 0) |

### What Makes This Case Interesting

This is Case 08 (advection-diffusion) reinterpreted as fiber optics.
The same mathematical structure describes a completely different physical
system. The GVD broadening is a real limitation in fiber-optic
telecommunications — it limits the maximum data rate for a given fiber
length. Case 36 shows how nonlinearity can counteract this broadening.

### Expected Results

The pulse translates to the right at v_g=1.0 while broadening. The peak
amplitude drops inversely with the width. The integral ∫A dV is conserved
(total energy constant). The dispersion length z_d = 0.04/(2·0.01) = 2.0,
so significant broadening is visible by t=2.

---

## Case 36: Soliton Pulse Propagation — Nonlinear Balance

### Physics

Adding Kerr nonlinearity to the dispersive pulse equation creates the
possibility of soliton propagation: a pulse that propagates without
changing shape because the dispersive broadening is exactly balanced by
nonlinear self-compression (Haus Ch 10):

```
∂A/∂t + v_g·∂A/∂x = D·∂²A/∂x² − α·A³
```

The soliton condition is α·A₀²·w₀² = 2D. With A₀=1, w₀=1, D=0.05:
α=0.1 gives a fundamental soliton with a sech profile that propagates
without distortion.

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `ADTimeDerivative` | ∂A/∂t |
| `ADConservativeAdvection` | v_g·∂A/∂x group velocity |
| `ADMatDiffusion` | D·∂²A/∂x² dispersion |
| `ADMatReaction` | Nonlinear −α·A³ from material rate |
| `DerivativeParsedMaterial` | rate = −α·A², with JIT disabled for Docker |
| `IterationAdaptiveDT` | Adaptive timestepping for nonlinear dynamics |

The key trick is that `ADMatReaction` contributes −rate·A·v to the
residual. Setting rate = −α·A² gives +α·A³·v, so the strong form gets
−α·A³ (self-focusing nonlinearity) — the correct sign.

### What Makes This Case Interesting

This is the most advanced case in the series. It combines advection,
diffusion, and nonlinearity in a single equation, demonstrating the
`DerivativeParsedMaterial` pattern for expressing a nonlinear reaction
rate. The soliton is one of the most beautiful results in nonlinear wave
theory — discovered in water waves by John Scott Russell in 1834 and
central to modern fiber-optic telecommunications via the work of Hasegawa
and Tappert (1973).

### Expected Results

With α=0.1 (soliton balance): peak amplitude stays at 1.0 and pulse
width stays constant as it translates across the domain. With α=0 (Case 35):
the pulse broadens. With α=0.3 (over-nonlinear): the pulse compresses
initially and then oscillates. The CSV shows max_A as a function of time —
a flat line at 1.0 for the soliton case, a decaying curve for dispersive,
and an oscillating curve for over-nonlinear.

---

## Case 37: Rayleigh-Benard Convection Onset

### Physics

Fluid heated from below between two rigid plates becomes unstable to
convective rolls above the critical Rayleigh number Ra_c = 1708. This is
the classical Benard problem (Rieutord Ch 7, Sec 7.5). The Boussinesq
incompressible Navier-Stokes equations with buoyancy coupling are:

```
div(u) = 0
rho*(du/dt + u.grad(u)) = -grad(p) + mu*lap(u) - rho*alpha*(T-T_ref)*g
dT/dt + u.grad(T) = kappa*lap(T)
```

With Ra = 2000 (just above Ra_c = 1708) and Pr = 0.71 (air), convective
rolls form and the Nusselt number Nu > 1 indicates heat transport
exceeding pure conduction.

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `[Modules/NavierStokesFV]` | Incompressible FV NS with Boussinesq buoyancy and energy equation |
| `boussinesq_approximation = true` | Activates buoyancy coupling rho_eff = rho*(1 - alpha*(T - T_ref)) |
| `FunctorMaterials/ADGenericFunctorMaterial` | Constant rho, mu, k, cp, alpha |
| `IterationAdaptiveDT` | Adaptive timestepping for convection onset |
| `ElementAverageValue` | Conservation check: avg(T) should stay at 0.5 |

### What Makes This Case Interesting

This is the first case to use Boussinesq buoyancy coupling with heating
from below. Case 16 (natural convection) heats from the side at Ra = 10000
in steady state. Here, the heating is from below, Ra is near critical, and
the simulation captures the transient onset of convective rolls — a
fundamentally different instability mechanism. The initial condition includes
a small sinusoidal perturbation to seed the instability.

### Expected Results

Average temperature remains at 0.5 (conservation). Convective velocities
develop slowly since Ra = 2000 is only ~17% above critical. max(vel_y)
grows from the initial perturbation level as the rolls form. The temperature
field shows weak convective rolls at the final time.

---

## Case 38: Kelvin-Helmholtz Instability — Shear Layer Rollup

### Physics

Two counterflowing streams separated by a shear layer roll up into
Kelvin-Helmholtz billows (Rieutord Ch 6, Sec 6.3.1). The velocity profile
is a tanh shear layer with half-thickness delta = 0.05, and a passive scalar
(temperature repurposed as dye) marks the two streams:

```
vel_x = tanh((y - 0.5) / 0.05)
T     = 0.5*(1 + tanh((y - 0.5) / 0.05))
```

A sinusoidal perturbation vel_y = 0.01*sin(2*pi*x) seeds the instability.

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `[Modules/NavierStokesFV]` | Incompressible FV NS with energy equation (passive scalar) |
| `ParsedFunction` | tanh velocity profile and scalar IC at inlet and initial condition |
| `inlet_boundaries` / `outlet_boundaries` | Inlet/outlet BCs to approximate periodic flow |
| `momentum_wall_types = 'slip slip'` | Free-slip top/bottom walls (symmetry planes) |
| `upwind` advection interpolation | Stabilizes steep shear layer gradients |

### What Makes This Case Interesting

This demonstrates a classical hydrodynamic instability with inlet/outlet
boundary conditions (rather than periodic). The energy equation is repurposed
as a passive scalar tracer — T acts as a dye that is advected by the flow
without affecting it. The initial tanh profile and sinusoidal perturbation
are set through MOOSE functions referenced directly by the NavierStokesFV
action's `initial_velocity` and `initial_temperature` parameters.

### Expected Results

The shear layer rolls up into vortex structures that grow over time.
The passive scalar contours show billow formation at t ~ 0.5-1.0 and
continued evolution through t = 2.0. Average temperature remains at
~0.5 (scalar conservation). max(vel_y) grows as the perturbation amplifies.

---

## Case 39: Blasius Boundary Layer — Flat Plate Laminar Flow

### Physics

The laminar boundary layer over a flat plate follows the Blasius similarity
solution (Rieutord Ch 4, Sec 4.3). Uniform flow U = 1 enters from the left
and encounters a no-slip wall at y = 0, creating a growing viscous shear
layer:

```
delta_99(x) = 5.0 * sqrt(mu * x / (rho * U))
```

The velocity profile u/U collapses onto the universal Blasius function
f'(eta) when plotted against eta = y*sqrt(U/(nu*x)). Key results:
f''(0) = 0.332, C_f = 0.664/sqrt(Re_x).

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `[Modules/NavierStokesFV]` | Incompressible FV NS, steady state, no energy equation |
| `GeneratedMeshGenerator` with `bias_y = 0.5` | Wall-refined mesh (more cells near y = 0) |
| `inlet_boundaries` / `outlet_boundaries` | Uniform inflow, pressure outlet |
| `wall_boundaries` (no-slip) + `slip` top | Flat plate at y = 0, free-stream at top |
| `PointValue` postprocessors | Sample velocity at specific stations along the plate |

### What Makes This Case Interesting

This is the first steady-state FV Navier-Stokes case with a non-trivial
spatial structure. The biased mesh (bias_y = 0.5) concentrates cells near
the wall where velocity gradients are steepest. Velocity profiles at
multiple x-stations should collapse onto the self-similar Blasius solution
when normalized by the local boundary layer thickness.

### Expected Results

The streamwise velocity contour shows the boundary layer growing along
the plate. Velocity profiles at x = 0.5, 1.0, and 1.5 show the
characteristic S-shaped profile transitioning from u = 0 at the wall to
u = U in the free stream. The boundary layer thickness at x = 2 is
delta_99 ~ 0.5 for Re_L = 400.

---

## Case 40: Turbulent Channel Flow — RANS k-epsilon Model

### Physics

Fully-developed turbulent flow in a 2D channel with RANS k-epsilon closure
and wall functions (Rieutord Ch 9, Sec 9.8). The mean velocity profile
follows the log-law of the wall:

```
u+ = (1/kappa) * ln(y+) + B,   kappa = 0.41, B = 5.2
```

The k-epsilon model solves two additional transport equations for turbulent
kinetic energy (TKE) and its dissipation rate (TKED), which together
determine the turbulent viscosity mu_t = rho*C_mu*k^2/epsilon.

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `MooseLinearVariableFVReal` | Linear FV variables: vel_x, vel_y, pressure, TKE, TKED |
| `SIMPLE` executioner | Segregated pressure-velocity coupling algorithm |
| `LinearWCNSFVMomentumFlux` | Momentum advection and diffusion |
| `LinearFVTKESourceSink` / `LinearFVTKEDSourceSink` | k-epsilon production and destruction |
| `kEpsilonViscosityAux` | Computes mu_t from k and epsilon fields |
| `LinearFVTurbulentViscosityWallFunctionBC` | Wall-function boundary conditions |
| `RhieChowMassFlux` | Pressure-velocity coupling in the mass equation |

### What Makes This Case Interesting

This is the only case using the SIMPLE segregated solver and linear FV
variables. Unlike Newton-coupled systems (Cases 15-16), each equation is
solved independently and iterated to convergence. The k-epsilon turbulence
model with wall functions is the industry-standard RANS approach. The
stitched mesh uses two blocks with opposing y-bias for wall refinement
on both channel walls.

### Expected Results

The SIMPLE solver converges in ~340 iterations. The centerline velocity
is approximately 1.1 * U_bulk (bulk velocity = 1), consistent with fully
developed turbulent channel flow. The velocity profile shows the
characteristic blunt shape of turbulent flow (much flatter than laminar
parabolic). TKE peaks near the walls and is lowest at the centerline.

---

## Case 41: Rayleigh-Taylor Instability — Heavy over Light

### Physics

A heavy fluid sitting atop a light fluid in a gravitational field is
inherently unstable (Rieutord Ch 6, Sec 6.3). A small perturbation at the
interface grows into mushroom-shaped fingers as the heavy fluid sinks
through the light fluid. The Boussinesq approximation uses temperature as
a density marker:

```
rho_eff = rho * (1 - alpha*(T - T_ref))
```

With alpha = 1, T_ref = 0.5: T = 0 (cold, heavy) on top, T = 1 (hot,
light) on bottom. The initial interface at y = 1 is perturbed by a
sinusoidal displacement.

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `[Modules/NavierStokesFV]` | Incompressible FV NS with Boussinesq buoyancy and energy |
| `ParsedFunction` (tanh interface IC) | Smooth initial density step at y = 1 |
| `ParsedFunction` (cosine perturbation) | Seeds the single-mode RT instability |
| `wall_boundaries` (all no-slip) | Confines the instability in a closed box |
| `IterationAdaptiveDT` | Adaptive timestepping as fingers accelerate |

### What Makes This Case Interesting

This case demonstrates a Rayleigh-Taylor instability in a tall domain
[0,1] x [0,2] with 25 x 50 elements. The initial condition combines a
tanh density step with a sinusoidal perturbation, seeding a single-mode
instability. The Boussinesq buoyancy drives the heavy fluid downward and
light fluid upward, creating the classic mushroom-finger morphology.

### Expected Results

max(vel_y) grows from ~0.03 at early times to ~0.33 at t = 3 as the RT
fingers develop. avg(T) remains at 0.5 throughout (scalar conservation).
The temperature contour snapshots at t = 0, 1, 2, 3 show progressive
development of the mushroom-shaped instability pattern.

---

## Case 42: Sod Shock Tube — 1D Riemann Problem

### Physics

A membrane at x = 0.5 separates high-pressure gas (left) from low-pressure
gas (right). At t = 0 the membrane bursts, producing three distinct waves:
a leftward rarefaction fan, a rightward contact discontinuity, and a
rightward shock wave (Rieutord Ch 5, Sec 5.5). The Rankine-Hugoniot jump
conditions govern the shock speed and post-shock state:

```
Left state:   rho = 1.0,   p = 1.0,   u = 0
Right state:  rho = 0.125, p = 0.1,   u = 0
gamma = 1.4 (ideal diatomic gas)
```

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `CNSFVMassHLLC`, `CNSFVMomentumHLLC`, `CNSFVFluidEnergyHLLC` | HLLC Riemann flux splitting for compressible Euler equations |
| `FVTimeKernel` | Time derivative of conservative variables (rho, rho*u, rho*E) |
| `ExplicitSSPRungeKutta` (order 2) | Explicit time integration for hyperbolic system |
| `IdealGasFluidProperties` | Equation of state: p = (gamma-1)*rho*e |
| `ConservedVarValuesMaterial` | Converts conservative to primitive variables |
| `AuxVariables` (pressure, vel_x) | Post-processes primitive fields from conservative variables |

### What Makes This Case Interesting

This is the only compressible flow case. Unlike the incompressible FV NS
cases (Cases 15-16, 37-41), it solves the conservative Euler equations
using HLLC flux splitting and an explicit Runge-Kutta time integrator —
there is no nonlinear Newton solve. The Sod shock tube is the benchmark
problem for validating any compressible flow solver against an exact
analytical solution.

### Expected Results

At t = 0.2: shock at x ~ 0.85, contact discontinuity at x ~ 0.69,
rarefaction fan spanning x ~ 0.26 to 0.49. Total mass is conserved
exactly (0.5625). The density profile shows all three wave structures
cleanly resolved on the 200-cell mesh.

---

## Case 43: Ekman Spiral — Rotating Boundary Layer

### Physics

Steady viscous flow near a wall in a rotating frame produces the Ekman
spiral (Rieutord Ch 8, Sec 8.4). The Coriolis force couples the x and y
velocity components, causing the velocity vector to rotate through the
boundary layer of thickness delta_E = sqrt(nu/Omega):

```
nu * d²vx/dz² + 2*Omega*vy = 0
nu * d²vy/dz² - 2*Omega*vx = -2*Omega*U_g
```

The analytical solution is:
vx(z) = U_g*(1 - exp(-z/delta_E)*cos(z/delta_E)),
vy(z) = U_g*exp(-z/delta_E)*sin(z/delta_E).

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `ADMatDiffusion` | Viscous diffusion nu*d²/dz² for both vx and vy |
| `CoupledForce` (coef = +2*Omega for vx, -2*Omega for vy) | Coriolis cross-coupling between velocity components |
| `BodyForce` (value = 2*Omega*U_g) | Geostrophic pressure gradient driving the flow |
| `DirichletBC` | No-slip at wall (z=0), geostrophic flow at z_max |
| `PointValue` at z = delta_E | Validates against analytical vx = 0.801, vy = 0.310 |

### What Makes This Case Interesting

This returns to the coupled scalar PDE pattern of Case 9, but now the
coupling has physical meaning: Coriolis acceleration in a rotating reference
frame. The careful derivation of CoupledForce signs (documented in the input
file) illustrates how MOOSE's residual conventions determine the coefficient
values. The Ekman hodograph (vy vs vx) traces a beautiful spiral from
(0, 0) at the wall to (U_g, 0) in the geostrophic interior.

### Expected Results

The solver converges in a single Newton iteration (linear problem). Computed
values match the analytical Ekman solution to 3 significant figures:
vx(delta_E) = 0.801, vy(delta_E) = 0.310, max(vy) = 0.322 at
z = pi/4 * delta_E = 0.0785.

---

## Case 44: Alfven Wave Propagation — MHD Elsasser Variables

### Physics

A transverse MHD wave propagates in a conducting fluid with a background
magnetic field B_0 (Rieutord Ch 10, Sec 10.4). Using Elsasser variables
d+ = vy + by/sqrt(mu_0*rho) and d- = vy - by/sqrt(mu_0*rho), the MHD
equations decouple into two advection-diffusion equations:

```
d(d+)/dt + v_A * d(d+)/dx = D_eff * d²(d+)/dx²
d(d-)/dt - v_A * d(d-)/dx = D_eff * d²(d-)/dx²
```

A Gaussian pulse in d+ propagates rightward at the Alfven speed v_A while
d- remains identically zero.

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `ADTimeDerivative` | d(d±)/dt |
| `ADConservativeAdvection` | Alfven wave advection ±v_A * d(d±)/dx |
| `ADMatDiffusion` | Resistive/viscous dissipation D_eff * d²(d±)/dx² |
| `ADGenericConstantVectorMaterial` | Advection velocity vectors (+v_A, 0, 0) and (-v_A, 0, 0) |
| `ElementExtremeValue` | Peak tracking for d+ decay and d- baseline |

### What Makes This Case Interesting

This is the final case in the series and the first to touch magnetohydrodynamics.
The Elsasser decomposition reduces the vector MHD system to two decoupled
advection-diffusion equations — the same pattern as Case 36 (soliton pulse)
but now with two fields advecting in opposite directions. The rightward d+
pulse propagates and diffusively decays while d- remains exactly zero,
confirming the clean decoupling of the Elsasser variables.

### Expected Results

The d+ Gaussian peak propagates rightward at v_A = 1.0, reaching x = 9 at
t = 6. Its amplitude decays as 1/sqrt(1 + 4*D*t/w^2) due to diffusion —
from 1.0 at t = 0 to ~0.36 at t = 6. d- stays at machine-zero throughout,
confirming perfect decoupling. Total integral of d+ decays monotonically
(diffusive dissipation).

---

## Case 45 — Monte Carlo UQ: Uncertainty in Thermal Conductivity

### Physics

Propagate uncertainty in thermal conductivity k ~ Uniform(0.5, 3.0) through 2D steady heat conduction with volumetric source. Uses the `stochastic_tools` module's MultiApp sampling architecture.

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `[StochasticTools]` | Creates dummy main problem (no PDE solved in main app) |
| `MonteCarlo` sampler | Generates 30 random k values from Uniform distribution |
| `SamplerFullSolveMultiApp` | Runs one full sub-app solve per sample |
| `SamplerParameterTransfer` | Sends sampled k to sub-app material property |
| `SamplerReporterTransfer` | Collects postprocessor values from sub-apps |
| `SamplerReceiver` | Control in sub-app that accepts transferred parameters |

### What Makes This Case Interesting

First case to use the `stochastic_tools` module. Introduces the MultiApp sampling architecture where the main app orchestrates sub-solves without solving a PDE itself. Shows how uncertainty in material properties propagates to output uncertainty.

### Expected Results

avg_T ranges from 1.09 to 6.26 (right-skewed distribution because T ~ 1/k). Strong linear correlation between avg_T and max_T (slope ≈ 2.32).

---

## Case 46 — Polynomial Chaos Expansion: Surrogate Modeling

### Physics

Build a polynomial chaos surrogate for a 1D diffusion-reaction problem with two uncertain parameters (D and sigma). Train on 36 deterministic quadrature points, then evaluate on 100 new random samples with zero additional PDE solves.

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `Quadrature` sampler | Generates deterministic integration points for training |
| `PolynomialChaosTrainer` | Trains order-5 PCE from quadrature data |
| `PolynomialChaos` surrogate | Loaded surrogate model for evaluation |
| `EvaluateSurrogate` reporter | Evaluates surrogate on new MC samples |
| `MatReaction` + `DerivativeParsedMaterial` | Handles sign convention for absorption term |

### What Makes This Case Interesting

Demonstrates surrogate-based UQ — 36 training solves replace 100+ MC solves. The polynomial chaos expansion captures the full input-output relationship analytically. Shows the Docker JIT workaround (`disable_fpoptimizer = true`, `enable_jit = false`) for `DerivativeParsedMaterial`.

### Expected Results

Surrogate predictions: mean = 0.177, std = 0.055, range 0.109–0.308. Training data shows smooth monotonic trend across 36 quadrature points.

---

## Case 47 — Heat Source Inversion: PDE-Constrained Optimization

### Physics

Recover an unknown volumetric heat source q from 4 sparse temperature measurements using adjoint-based gradient optimization. Three-file architecture: main (TAO optimizer), forward (heat equation), adjoint (sensitivity equation).

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `[Optimization]` | Creates dummy main problem for optimizer |
| `Optimize` executioner | Drives TAO L-BFGS optimizer |
| `GeneralOptimization` | Manages parameters, bounds, and objective |
| `ParsedOptimizationFunction` | Links optimization parameter to PDE source term |
| `OptimizationData` | Computes misfit and objective at measurement points |
| `ReporterPointSource` | Applies misfit as point loads in adjoint equation |
| `ElementOptimizationSourceFunctionInnerProduct` | Computes gradient dJ/dq |

### What Makes This Case Interesting

First case using the `optimization` module. Introduces adjoint methods — the gradient is computed by solving one additional PDE (the adjoint), not by finite differences. The optimizer recovers q = 1000 from initial guess 500 in just 1 L-BFGS step because the problem is linear.

### Expected Results

Recovered q = 1000.0 (exact match to true value). Objective drops from ~2×10⁴ to 1.6×10⁻²³ in a single iteration.

---

## Case 48 — Latin Hypercube Parameter Study: Multi-Parameter UQ

### Physics

Multi-parameter uncertainty study with 3 uncertain inputs (k, T_left, T_right) using the high-level `[ParameterStudy]` action. 50 Latin Hypercube samples. The entire main input file is just 15 lines.

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `[ParameterStudy]` action | Creates all samplers, multiapps, transfers, and output automatically |
| Latin Hypercube sampler | Space-filling sampling for efficient parameter coverage |

### What Makes This Case Interesting

Contrast with Case 45's manual setup — the `[ParameterStudy]` action replaces ~50 lines of configuration with ~15 lines. Shows that T_right dominates avg_T, T_left has moderate effect, and k has secondary effect through the volumetric source term.

### Expected Results

50 LHS samples. avg_T ranges from ~105 to ~248 (mean ~179). Strong correlation between T_right and avg_T.

---

---

## Case 49 — J2 Plasticity: Uniaxial Tension with Isotropic Hardening

### Physics

A 3D bar is pulled in uniaxial tension past the yield stress. The `solid_mechanics` module's `ComputeMultipleInelasticStress` object drives J2 (von Mises) plasticity with linear isotropic hardening: once the effective stress exceeds the initial yield stress, further straining accumulates plastic strain and hardens the material. The elastic-perfectly-plastic and hardening phases produce a bilinear stress-strain curve that is verified against the exact analytic result.

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `Physics/SolidMechanics/QuasiStatic` | Sets up displacement variables, strain, and stress computation |
| `ComputeMultipleInelasticStress` | Integrates the constitutive update with inelastic return mapping |
| `IsotropicPlasticityStressUpdate` | J2 return-mapping algorithm with isotropic hardening modulus |
| `FunctionDirichletBC` | Ramps applied displacement linearly with time |

---

## Case 50 — Finite Strain: Large Deformation Compression

### Physics

A rubber-like block is compressed to 50% of its original height. Finite-strain kinematics (using the multiplicative decomposition F = Fe Fp) are required because small-strain assumptions break down above ~5% deformation. The `solid_mechanics` module's `ComputeFiniteStrainElasticStress` path tracks the deformed configuration, producing correct geometric stiffening and accurate Cauchy stresses even at large strains.

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `ComputeFiniteStrain` | Computes the deformation gradient and Green-Lagrange strain |
| `ComputeFiniteStrainElasticStress` | Rotated Cauchy stress in the deformed configuration |
| `PresetDisplacement` | Applies a prescribed displacement history to the top surface |
| `RankTwoScalarAux` | Extracts scalar invariants (effective strain, volumetric strain) |

---

## Case 51 — Power-Law Creep: Column Under Sustained Compression

### Physics

A ceramic column carries a constant compressive load at elevated temperature. Over time, creep strain accumulates according to the power-law model: strain_rate = A * sigma^n * exp(-Q/RT). The `PowerLawCreepStressUpdate` object integrates this rate law implicitly at each timestep. Post-processing tracks the creep strain history and compares the mid-column strain evolution to the closed-form uniaxial solution.

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `PowerLawCreepStressUpdate` | Implicit return-mapping for steady-state power-law creep |
| `ComputeMultipleInelasticStress` | Manages inelastic strain decomposition alongside creep |
| `IterationAdaptiveDT` | Grows the timestep during slow creep, shrinks during rapid accumulation |
| `ElementAverageValue` | Tracks average creep strain for comparison to analytic solution |

---

## Case 52 — Phase-Field Fracture: Notched Specimen Under Tension

### Physics

A pre-notched rectangular specimen is pulled in tension until a crack nucleates and propagates. The phase-field fracture model couples a damage variable d (0 = intact, 1 = fully broken) to the displacement field: energy is stored elastically until the local fracture energy criterion is met, at which point d grows and the stiffness degrades smoothly. This avoids mesh-dependent crack paths that arise with sharp-crack models.

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `PhaseFieldFractureMechanicsOffDiag` | Off-diagonal Jacobian coupling between damage and displacement |
| `ComputeLinearElasticPFFractureStress` | Degraded stiffness tensor (1-d)^2 * C_e |
| `ADPFFracture` | Residual for the Allen-Cahn-type damage evolution equation |
| `SMP` | Full off-diagonal preconditioning for the coupled damage-mechanics system |

---

## Case 53 — Pressure Vessel: Thick-Walled Cylinder (Lame Solution)

### Physics

An internally pressurized thick-walled cylinder (inner radius a, outer radius b) under plane-strain conditions. The Lame analytic solution gives radial and hoop stresses as functions of r: sigma_r = A + B/r^2, sigma_theta = A - B/r^2, where A and B are determined by the boundary conditions. This case verifies that the `solid_mechanics` module reproduces the exact stress distribution to machine precision, making it an ideal regression test for the elastic kernel.

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `Physics/SolidMechanics/QuasiStatic` | Axisymmetric plane-strain formulation |
| `ComputeIsotropicElasticityTensor` | Isotropic material with given Young's modulus and Poisson's ratio |
| `ADComputeSmallStrain` | Small-strain assumption (valid for typical pressure vessel loading) |
| `ElementL2Error` | Quantifies deviation of FEM stress from the Lame analytic solution |

---

## Cases 54-58: Nuclear Reactor Physics

These five cases introduce neutronics and reactor thermal-hydraulics. Cases 54-55 solve
the neutron diffusion eigenvalue equation (1-group and 2-group). Case 56 couples the heat
equation to a volumetric source in axisymmetric geometry. Case 57 adds xenon-135 kinetics
to model flux depression after a power transient. Case 58 computes eigenvalue shift from a
spatially varying control-rod absorber. All use `combined-opt` via Docker.

---

## Case 54 — 1-Group Neutron Diffusion: Bare Slab Criticality

### Physics

The one-group neutron diffusion equation on a bare slab of thickness L describes the
neutron flux distribution at criticality. The governing eigenvalue problem is:

```
-D * d²phi/dx² + (Sigma_a - nu*Sigma_f/k_eff) * phi = 0
```

where D is the diffusion coefficient, Sigma_a is the macroscopic absorption cross section,
nu*Sigma_f is the neutron production cross section, and k_eff is the effective multiplication
factor. For a bare slab with zero-flux boundary conditions at the extrapolated boundaries,
the critical buckling is B² = (pi/L)² and the analytic k_eff = nu*Sigma_f / (D*B² + Sigma_a).
This case uses MOOSE's `Eigenvalue` executioner with `KRYLOVSCHUR` to recover k_eff = 1.0000
and verify the cosine flux shape phi(x) = cos(pi*x/L).

```bash
cd quickstart-runs/case54-neutron-diffusion-bare-slab
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case54-neutron-diffusion-bare-slab \
  --entrypoint /bin/bash idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case54_neutron_diffusion_bare_slab.i 2>&1 | tail -20'
```

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `Eigenvalue` executioner | Drives the generalized eigenvalue solve for k_eff |
| `DiffusionKernel` / `Diffusion` | Implements `-D * div(grad phi)` neutron leakage term |
| `CoefReaction` (absorption) | Implements `+Sigma_a * phi` loss term |
| `CoefReaction` (fission, eigen tag) | Implements `nu*Sigma_f * phi` as the B-matrix for the eigenvalue |
| `EigenDirichletBC` | Enforces zero flux on both slab faces in the eigen system |

---

## Case 55 — 2-Group Neutron Diffusion: Fast/Thermal Coupling with Fission

### Physics

The two-group diffusion model separates neutrons into a fast group (group 1, high energy)
and a thermal group (group 2, low energy). Fast neutrons are produced by fission and
slow down into the thermal group; thermal neutrons drive most of the fission. The coupled
eigenvalue system is:

```
-D1*div(grad phi1) + (Sigma_a1 + Sigma_12)*phi1 = (nu*Sigma_f1*phi1 + nu*Sigma_f2*phi2) / k_eff
-D2*div(grad phi2) + Sigma_a2*phi2               = Sigma_12*phi1
```

The down-scattering term Sigma_12 couples the two equations. With typical PWR-like cross
sections this system has k_eff = 1.342, indicating a supercritical assembly that would
require control rods or boron to reach criticality. The two-group formulation is the
simplest model that captures the physics of fast-to-thermal neutron slowing-down and is
the standard starting point for reactor core analysis.

```bash
cd quickstart-runs/case55-two-group-diffusion
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case55-two-group-diffusion \
  --entrypoint /bin/bash idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case55_two_group_diffusion.i 2>&1 | tail -20'
```

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `Eigenvalue` executioner | Solves the two-variable generalized eigenvalue problem |
| `Diffusion` (×2) | Neutron leakage for each group |
| `CoefReaction` (removal, ×2) | Group removal cross sections (absorption + scatter-out) |
| `CoupledForce` (down-scatter) | Sigma_12 * phi1 source term in group-2 equation |
| `CoupledForce` (fission source, eigen tag) | nu*Sigma_f contribution from both groups to group 1 |
| `EigenDirichletBC` | Zero flux at domain boundaries for both flux variables |

---

## Case 56 — Fuel Pin Heat Transfer: Radial Temperature Profile

### Physics

A cylindrical nuclear fuel pin generates heat uniformly at volumetric rate Q''' (W/m³).
The steady-state heat equation in axisymmetric (RZ) cylindrical coordinates is:

```
-(1/r) * d/dr (r * k_fuel * dT/dr) = Q'''   in fuel (0 < r < r_fuel)
-(1/r) * d/dr (r * k_clad * dT/dr) = 0      in cladding (r_fuel < r < r_clad)
```

with contact resistance at the fuel-clad interface and a convective boundary condition at
the outer cladding surface (h, T_coolant). The exact solution is a parabolic temperature
profile in the fuel and a logarithmic profile in the cladding. Typical PWR conditions
(Q''' = 300 MW/m³, k_fuel = 2.5 W/m·K, k_clad = 16 W/m·K) give a centerline temperature
of approximately 1200 °C with a coolant temperature of 300 °C — close to the UO₂ melting
point, illustrating why fuel-pin thermal analysis is safety-critical.

```bash
cd quickstart-runs/case56-fuel-pin-heat-transfer
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case56-fuel-pin-heat-transfer \
  --entrypoint /bin/bash idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case56_fuel_pin_heat_transfer.i 2>&1 | tail -20'
```

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `GeneratedMesh` (coord_type = RZ) | 1D radial mesh with axisymmetric (cylindrical) coordinates |
| `HeatConduction` | Radial heat conduction `-div(k*grad T)` in cylindrical geometry |
| `HeatSource` | Uniform volumetric source Q''' in the fuel subdomain |
| `SubdomainBoundingBoxGenerator` | Splits mesh into fuel and cladding subdomains |
| `ConvectiveHeatFluxBC` | Newton cooling `h*(T - T_coolant)` at the outer cladding surface |
| `SideAverageValue` | Reports peak centerline and interface temperatures |

---

## Case 57 — Xenon-135 Poisoning Transient: Reactor Flux Decay

### Physics

After a reactor is shut down (or power is reduced), the neutron flux drops but iodine-135
(which decays to xenon-135) continues to be produced from the decay of short-lived fission
products. Xenon-135 has an enormous thermal neutron absorption cross section (2.65 × 10⁶ barns)
and its buildup severely depresses the neutron flux — a phenomenon called xenon poisoning.
The coupled ODE/PDE system governing this transient is:

```
dI/dt   = gamma_I * Sigma_f * phi - lambda_I * I
dXe/dt  = gamma_Xe * Sigma_f * phi + lambda_I * I - (lambda_Xe + sigma_Xe * phi) * Xe
dphi/dt = D*div(grad phi) - (Sigma_a + sigma_Xe*Xe)*phi + nu*Sigma_f*phi/k0 - phi/tau_prompt
```

where I is the iodine-135 concentration, Xe is the xenon-135 concentration, phi is the
thermal neutron flux, gamma_I and gamma_Xe are fission yields, lambda_I and lambda_Xe are
decay constants, and sigma_Xe is the xenon absorption cross section. After a step-down in
power at t=0, xenon peaks at approximately 6-8 hours and recovers by 24 hours, producing
the classic "xenon transient" that operators must manage to avoid iodine pit startup failures.
This 24-hour transient simulation uses three coupled variables and demonstrates MOOSE's
ability to handle stiff ODE/PDE systems with multiple timescales.

```bash
cd quickstart-runs/case57-xenon-poisoning
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case57-xenon-poisoning \
  --entrypoint /bin/bash idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case57_xenon_poisoning.i 2>&1 | tail -20'
```

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `TimeDerivative` (×3) | Time-derivative terms for iodine, xenon, and flux variables |
| `CoupledForce` | Cross-coupling terms (iodine decay into xenon, xenon absorption on flux) |
| `ParsedMaterial` / `DerivativeParsedMaterial` | Spatially uniform nuclear data (cross sections, yields, decay constants) |
| `IterationAdaptiveDT` | Adaptive timestep that shortens during the rapid xenon peak |
| `ElementAverageValue` | Tracks spatially averaged I-135, Xe-135, and flux over 24 hours |

---

## Case 58 — Control Rod Worth: Eigenvalue Shift from Absorber

### Physics

A control rod is modeled as a spatially localized region of enhanced neutron absorption.
Inserting the rod raises Sigma_a in the rod region, shifting the effective multiplication
factor from the unrodded k_eff to a lower rodded k_eff. The difference

```
Delta_k = k_eff_unrodded - k_eff_rodded
```

is called the control rod worth and is the key parameter for reactor control and shutdown
margin calculations. This case runs two eigenvalue solves — one without the rod and one
with the rod inserted — and reports both k_eff values and the worth. The spatial flux
depression near the rod tip demonstrates the neutron flux "peaking" that must be accounted
for in fuel-rod power peaking factor analyses.

```bash
cd quickstart-runs/case58-control-rod-worth
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case58-control-rod-worth \
  --entrypoint /bin/bash idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case58_control_rod_worth.i 2>&1 | tail -20'
```

### Key MOOSE Objects

| Object | Role |
|--------|------|
| `SubdomainBoundingBoxGenerator` | Creates a rod-shaped subdomain with elevated absorption cross section |
| `GenericConstantMaterial` | Assigns group-wise nuclear data to fuel and rod subdomains |
| `Eigenvalue` executioner | Solves the one-group diffusion eigenvalue for each rod position |
| `CoefReaction` (absorption) | `Sigma_a(x) * phi` using block-restricted material for the rod region |
| `ElementIntegral` (flux) | Spatial flux distribution showing flux depression in the rod region |

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

After completing these 58 cases:

1. **Read the MOOSE documentation** at https://mooseframework.inl.gov for
   complete reference documentation on every object type.

2. **Explore the test suite**: `test/tests/` contains thousands of working
   input files covering every feature. Each subdirectory has a `tests` spec
   file describing what each input file demonstrates.

3. **Write your own application**: Use `moose/scripts/stork.py` to scaffold
   a new MOOSE application with custom kernels, materials, and BCs.

4. **Explore more module features**: Cases 14-58 introduce the major physics
   modules — from solid mechanics and heat transfer through Navier-Stokes,
   electrodynamics, MHD, nonlinear solid mechanics, and nuclear reactor physics.
   Each module has many more capabilities — consult the [Modules Reference](modules-reference.md)
   for full details.
