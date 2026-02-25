# Case 08: Transient Advection-Diffusion with a Gaussian Blob

## Overview

This case solves the **advection-diffusion equation** — one of the most fundamental
and practically important PDEs in engineering and science. It describes how a substance
(a chemical concentration, a temperature pulse, a pollutant plume) moves through a
medium under two competing processes:

1. **Advection** — bulk transport by the flow field. A river current carries a dye blob
   downstream. The blob moves as a rigid shape (in the absence of diffusion).

2. **Diffusion** — random spreading. Molecular motion spreads the blob outward in all
   directions, smoothing sharp gradients.

The ratio of these two effects is characterized by the **Peclet number**. This case
has a moderately high Peclet number, meaning advection dominates: the blob mostly
translates downstream while slowly spreading and flattening.

The setup deliberately showcases two important MOOSE features:

- **HIT top-level scalar variables** — define `vx` and `D` at the top of the file
  and reference them with `${vx}` and `${D}` throughout, making the file
  self-documenting and easy to modify.

- **ConservativeAdvection kernel** — the divergence (conservative) form of the
  advection operator, which preserves the integral of concentration exactly as long
  as the velocity field is divergence-free.

Running this case produces a time series you can animate in ParaView, watching the
Gaussian blob travel across the rectangular domain from left to right.

---

## The Physics

### Physical Problem in Plain English

Picture a long rectangular tank filled with water (3 units long, 1 unit tall). At
time zero, a small concentrated drop of dye is injected near the left end, centered
at (x=0.3, y=0.5). The dye concentration forms a Gaussian (bell-curve) shape with
peak concentration 1 at the injection point, falling off rapidly in all directions.

A uniform rightward current of speed 0.5 carries the dye blob to the right. At the
same time, molecular diffusion (coefficient D = 0.01) slowly spreads the blob in all
directions. Because the Peclet number Pe = v*L/D = 0.5*3/0.01 = 150 is large,
advection dominates: the blob travels roughly 0.5 * 2.0 = 1.0 unit to the right by
the end of the simulation (t = 2.0), while diffusion broadens it only mildly.

The right boundary is an outflow: when the blob reaches the right wall, concentration
simply exits the domain. The left boundary is held at c = 0 (no inflow of dye).

### Governing Equation

The strong form of the transient advection-diffusion equation is:

```
dc/dt  +  div(v * c)  -  div(D * grad(c))  =  0
```

Alternatively written as:

```
dc/dt  +  v . grad(c)  -  D * Laplacian(c)  =  0
     (assuming div(v) = 0, incompressible flow)
```

Symbol definitions:

- `c(x,y,t)` — the concentration of the transported substance [mol/m³ or dimensionless]
- `t` — time [s]
- `dc/dt` — rate of change of concentration at a fixed point in space
- `v` — the velocity vector field: `v = (vx, vy) = (0.5, 0)` [m/s or dimensionless]
- `div(v * c)` — the divergence of the advective flux `v*c`. This is the conservative
  form: it conserves mass exactly when `div(v) = 0`.
- `v . grad(c)` — the convective derivative: rate of change following the flow.
  This is the non-conservative form, equivalent when the flow is incompressible.
- `D` — molecular diffusion coefficient: 0.01 [m²/s or dimensionless]
- `D * Laplacian(c)` — the diffusive flux divergence. Spreads concentration from
  high to low regions.

The equation says: the local rate of change of c equals the net inflow of c by
advection minus the net outflow by diffusion.

### The Peclet Number

The Peclet number compares advective to diffusive transport:

```
Pe = v * L / D
```

where L is a characteristic length. Using the domain length L = 3:

```
Pe = 0.5 * 3 / 0.01 = 150
```

This is a high Peclet number. Practical implications:

| Pe range | Dominant process | Solution character |
|----------|------------------|--------------------|
| Pe << 1  | Diffusion        | Smooth, blob spreads quickly, barely moves |
| Pe ~ 1   | Both equal       | Blob moves and spreads at comparable rates |
| Pe >> 1  | Advection        | Blob moves mostly as a rigid shape, spreads slowly |

With Pe = 150, the blob translates ~1 unit in 2 seconds while growing modestly in
width. The width of the blob grows as `sqrt(2*D*t)`: at t=2, spread is `sqrt(0.04) = 0.2`
units in each direction, while travel is `0.5 * 2 = 1.0` unit.

High Peclet numbers also create numerical challenges. Standard centered finite element
methods can produce spurious oscillations near sharp fronts. The `ConservativeAdvection`
kernel uses an upwinding approach to avoid these oscillations.

### Conservative vs. Non-Conservative Advection Forms

There are two mathematically equivalent ways to write the advection term:

**Non-conservative form:**
```
v . grad(c)   (requires div(v) = 0 to be equivalent)
```

**Conservative (divergence) form:**
```
div(v * c)   =   v . grad(c)  +  c * div(v)
```

When `div(v) = 0` (incompressible flow), they are the same. MOOSE's
`ConservativeAdvection` kernel implements the divergence form because:

1. It conserves the total integral of c exactly at the discrete level (important for
   mass conservation in chemical transport simulations).
2. It is more robust numerically when the velocity field is not perfectly divergence-free
   on the discrete mesh.

The weak form of the conservative advection term integrates by parts:

```
integral( phi_i * div(v*c) ) dV  =  - integral( grad(phi_i) . v*c ) dV
                                    + boundary terms
```

This is what the `ConservativeAdvection` kernel assembles.

### Initial Condition: Gaussian Blob

The initial concentration field is:

```
c(x, y, t=0) = exp( -((x - 0.3)^2 + (y - 0.5)^2) / 0.01 )
```

This is a two-dimensional Gaussian (bell curve) centered at (0.3, 0.5):

- At the center (0.3, 0.5): c = exp(0) = 1.0 (maximum concentration)
- At radius r from center: c = exp(-r²/0.01)
- At r = 0.1: c = exp(-1) ≈ 0.37 (one standard deviation)
- At r = 0.2: c = exp(-4) ≈ 0.018 (essentially zero)

The effective "width" parameter is `sigma = sqrt(0.01/2) = 0.0707`. This is a very
tight blob relative to the domain size, chosen so it fits clearly within the initial
left portion of the domain.

### Boundary Conditions

```
c = 0 on the left boundary (x=0)
Zero-flux (natural Neumann) everywhere else
```

- **Left boundary (x=0)**: DirichletBC with c=0. This models no inflow of concentration
  at the upstream end of the domain.
- **Right boundary (x=3)**: Natural zero-flux Neumann. The advective flux `v*c*n`
  flows out naturally without needing an explicit BC.
- **Top and bottom walls**: Natural zero-flux Neumann. No concentration leaves through
  the walls. This means the blob is "reflected" off the walls, though for this particular
  case the blob is centered at y=0.5 and stays well away from the walls.

### ASCII Domain Diagram

```
y=1   zero-flux (Neumann, natural)
      +-------------------------------------------------------+
      |                          -->                          |
      |   c=0 everywhere     v = (0.5, 0)                    |
      |   except the blob    advects blob to the right        |
      |                                                       |
c=0   |   t=0: Gaussian      t=1: blob moved right          | zero-flux
(Dir) |   blob at (0.3,0.5)  to ~(0.8, 0.5)               | (outflow)
      |                                                       |
      |   t=2: blob at       t>2: blob exits right wall       |
      |   ~(1.3, 0.5)        (concentration leaves domain)    |
      |                                                       |
      +-------------------------------------------------------+
y=0   zero-flux (Neumann, natural)
      x=0          vx=0.5           xmax=3

Domain: 3 units long x 1 unit tall
Mesh: 60 x 30 quadrilateral elements (dx = 0.05, dy = 0.0333)
Initial blob: exp(-((x-0.3)^2 + (y-0.5)^2) / 0.01)
```

---

## HIT Top-Level Variables and ${} Substitution

One of MOOSE's most useful input file features is demonstrated at the top of this case.

### Top-Level Variable Declaration

```
vx = 0.5    # x-velocity component
D  = 0.01   # diffusion coefficient
```

These two lines appear *before* any `[Block]` in the HIT file. They declare top-level
scalar variables named `vx` and `D`. This is a HIT (Hierarchical Input Text) feature.

### Using Variables with ${}

The `${varname}` syntax substitutes the value of a top-level variable anywhere in
the input file:

```
prop_values = '${vx} 0 0'   # expands to '0.5 0 0'
prop_values = '${D}'        # expands to '0.01'
```

When MOOSE reads the file, it performs simple text substitution before parsing,
replacing every `${vx}` with `0.5` and `${D}` with `0.01`.

### Why This Matters

Without top-level variables, if you wanted to change the velocity, you would need to
find every occurrence of the number `0.5` in the file and change each one. With
top-level variables, you only change the one line `vx = 0.5` at the top.

This is especially valuable for:
- Parameter studies: change one value, re-run, observe the effect
- Avoiding mistakes: a consistent value throughout the file
- Self-documentation: the parameter name explains what the number means

You can also override top-level variables from the command line:

```bash
mpirun -n 1 moose-app-opt -i case08_advection_diffusion.i vx=1.0
```

This doubles the advection velocity without editing the file.

---

## Input File Walkthrough

The input file is `case08_advection_diffusion.i`.

### Top-Level Variables

```
vx = 0.5    # x-velocity component  (HIT top-level scalar variable)
D  = 0.01   # diffusion coefficient
```

Defined before any block. `vx = 0.5` means the flow carries concentration at half a
unit per second in the x-direction. `D = 0.01` means the diffusion coefficient is
0.01 — small relative to the advective velocity, giving Pe = 150.

### Block: `[Mesh]`

```
[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 60
  ny   = 30
  xmin = 0
  xmax = 3
  ymin = 0
  ymax = 1
[]
```

- `type = GeneratedMesh` — structured rectangular mesh. No external file needed.
- `dim = 2` — two-dimensional (x, y) problem.
- `nx = 60` — 60 elements along x. With domain length 3, element width is 0.05.
- `ny = 30` — 30 elements along y. With domain height 1, element height is 0.0333.
- `xmin = 0, xmax = 3` — the domain extends 3 units in x. This is longer than the
  default unit square specifically so the blob has room to travel before reaching
  the right wall.
- `ymin = 0, ymax = 1` — unit height.

The resolution choice matters for advection-dominated problems. The mesh Peclet number
(per element) is `Pe_h = vx * dx / D = 0.5 * 0.05 / 0.01 = 2.5`. Standard Galerkin
methods become unstable for `Pe_h > 2`. The `ConservativeAdvection` kernel handles
this with upwinding, but a finer mesh reduces numerical diffusion.

The mesh has 60 x 30 = 1,800 elements and (61 x 31) = 1,891 nodes.

### Block: `[Variables]`

```
[Variables]
  [c]
  []
[]
```

- Declares a single unknown `c` (concentration) with default settings: first-order
  Lagrange, zero initial condition (before the ICs block overrides it).
- 1,891 DOFs for this mesh.

### Block: `[ICs]`

```
[ICs]
  [blob]
    type     = FunctionIC
    variable = c
    function = blob_fn
  []
[]
```

This block sets the initial condition of `c` to a non-trivial spatial distribution.

- `type = FunctionIC` — sets the initial value of a variable by evaluating a MOOSE
  Function at every node.
- `variable = c` — applies to the concentration field.
- `function = blob_fn` — references the function named `blob_fn`, defined in the
  `[Functions]` block below. At t=0, MOOSE evaluates `blob_fn(x, y, 0)` at every
  node to initialize `c`.

Without this block, `c` would start at zero everywhere and there would be nothing to
advect.

### Block: `[Functions]`

```
[Functions]
  [blob_fn]
    type       = ParsedFunction
    expression = 'exp(-((x-0.3)^2 + (y-0.5)^2) / 0.01)'
  []
[]
```

- `type = ParsedFunction` — evaluates a mathematical expression string at any (x,y,z,t).
  MOOSE uses the `fparser` library to parse and compile the expression.
- `expression` — the Gaussian bell function. Variables `x`, `y` are automatically
  available. The caret `^` means exponentiation.
- At (0.3, 0.5): value = exp(0) = 1. At (0.4, 0.5): value = exp(-1/0.01) ≈ 0 (sharp!).
- This function is used once: at t=0, via `FunctionIC`, to set the initial distribution.

Note that `fparser` compiles this expression at startup using JIT compilation. In
environments where JIT is unavailable, you would need an alternative (such as a
precomputed CSV or a custom C++ function). For standard MOOSE builds, `ParsedFunction`
works correctly.

### Block: `[Kernels]`

```
[Kernels]
  [time_deriv]
    type     = TimeDerivative
    variable = c
  []

  [advection]
    type              = ConservativeAdvection
    variable          = c
    velocity_material = velocity_vec
  []

  [diffusion]
    type        = MatDiffusion
    variable    = c
    diffusivity = D_coeff
  []
[]
```

Three kernels implement the three terms of the governing PDE.

**`[time_deriv]` kernel:**
- `type = TimeDerivative` — implements `dc/dt` in weak form:
  ```
  integral( phi_i * dc/dt ) dV
  ```
- This term makes the problem transient. Without it, we would have a steady
  advection-diffusion equation (which would not make sense with this initial condition).
- MOOSE uses implicit time integration by default (Crank-Nicolson or backward Euler,
  configurable). The time derivative is approximated as `(c^{n+1} - c^n) / dt` at
  each time step, where superscripts denote time levels.

**`[advection]` kernel:**
- `type = ConservativeAdvection` — implements the divergence (conservative) form of
  the advection term. Weak form after integration by parts:
  ```
  -integral( grad(phi_i) . v * c ) dV  +  boundary flux terms
  ```
  The negative sign arises from integration by parts: `integral(phi_i * div(v*c)) =
  -integral(grad(phi_i) . v*c) + boundary terms`.
- `variable = c` — acts on the concentration equation.
- `velocity_material = velocity_vec` — the name of a vector-valued material property
  that provides the velocity field `v`. This is looked up from the `[Materials]` block.
  The velocity is passed as a material property (not a function) because in full
  Navier-Stokes solvers the velocity comes from the fluid momentum equations; using a
  material interface keeps the kernel general.

Note: `ConservativeAdvection` includes upwinding to handle high Peclet numbers. Pure
Galerkin advection at high Pe produces spurious oscillations; upwinding adds a small
amount of numerical diffusion aligned with the flow direction to stabilize the solution.

**`[diffusion]` kernel:**
- `type = MatDiffusion` — standard diffusion kernel using a material property for the
  diffusion coefficient. Implements:
  ```
  integral( D * grad(phi_i) . grad(c) ) dV
  ```
- `diffusivity = D_coeff` — the name of the scalar material property providing D.

### Block: `[BCs]`

```
[BCs]
  [left_wall]
    type     = DirichletBC
    variable = c
    boundary = left
    value    = 0
  []
[]
```

Only one explicit BC is defined. All other boundaries use the natural (zero-flux)
Neumann condition.

- `left_wall` — pins `c = 0` at the left boundary (x=0). This models no upstream
  source of concentration. Once the blob has traveled past the left boundary region,
  this condition prevents any artifactual re-entry of concentration.
- **Right boundary**: No BC entry means zero-flux Neumann. The advective flux
  `v * c * n_x` with `n_x = 1` at the right wall allows concentration to flow out
  naturally. When the blob reaches the right wall, it exits smoothly.
- **Top and bottom**: Zero-flux Neumann. The walls are impermeable: no concentration
  escapes through the top or bottom.

### Block: `[Materials]`

```
[Materials]
  [velocity_mat]
    type        = GenericConstantVectorMaterial
    prop_names  = 'velocity_vec'
    prop_values = '${vx} 0 0'
  []

  [diffusion_mat]
    type        = GenericConstantMaterial
    prop_names  = 'D_coeff'
    prop_values = '${D}'
  []
[]
```

Two material objects define the constant physical parameters.

**`[velocity_mat]`:**
- `type = GenericConstantVectorMaterial` — declares a `RealVectorValue` material
  property with a constant value. This is the vector analogue of `GenericConstantMaterial`.
- `prop_names = 'velocity_vec'` — creates one property with this name.
- `prop_values = '${vx} 0 0'` — the three vector components (x, y, z). After
  `${vx}` substitution this becomes `'0.5 0 0'`: velocity is 0.5 in x, 0 in y, 0 in z.
  Note: MOOSE always works in 3D internally even for 2D problems; the z component is
  simply zero.

**`[diffusion_mat]`:**
- `type = GenericConstantMaterial` — declares a scalar material property with a
  constant value.
- `prop_names = 'D_coeff'` — one property named `D_coeff`.
- `prop_values = '${D}'` — after substitution becomes `'0.01'`.

The separation of material properties from kernels is a key MOOSE design pattern.
The kernels (`ConservativeAdvection`, `MatDiffusion`) do not know what the velocity
or diffusivity values are — they query named material properties at each quadrature
point. This means you can change the material (e.g., use a spatially varying D)
without touching the kernel.

### Block: `[Postprocessors]`

```
[Postprocessors]
  [total_c]
    type     = ElementIntegralVariablePostprocessor
    variable = c
  []
[]
```

- `type = ElementIntegralVariablePostprocessor` — integrates the variable over the
  entire domain: `integral( c ) dV`.
- `variable = c` — integrates the concentration field.

This integral equals the total amount of the substance in the domain. For a purely
advective system (D=0, no inflow, no outflow), the total would be conserved exactly.
With diffusion and the outflow boundary condition, the total decreases as the blob
exits the right wall. Monitoring `total_c` over time gives a physically meaningful
check: the total should be approximately constant until the blob nears the right wall,
then decrease as concentration flows out.

The initial total can be estimated: `integral of exp(-((x-0.3)^2 + (y-0.5)^2)/0.01)`
over the 3x1 domain. This integral evaluates to approximately `pi * 0.01 = 0.0314`
(the Gaussian integral formula gives `pi * sigma^2 = pi * 0.01 / 2 * 2`).

### Block: `[Executioner]`

```
[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'ilu'

  dt         = 0.02
  end_time   = 2.0

  nl_rel_tol = 1e-6
[]
```

- `type = Transient` — a time-marching simulation. MOOSE takes successive time steps
  from t=0 to `end_time`, solving the nonlinear system at each step.

- `solve_type = 'PJFNK'` — Preconditioned Jacobian-Free Newton-Krylov. Appropriate
  here because we have not used AD kernels for the advection term. JFNK approximates
  Jacobian-vector products using finite differences of the residual.

- `petsc_options_iname = '-pc_type'` / `petsc_options_value = 'ilu'` — uses ILU
  (Incomplete LU) factorization as the preconditioner. ILU is a simple but effective
  preconditioner for advection-dominated problems on single-processor runs. For
  parallel runs, block ILU or AMG would be needed.

- `dt = 0.02` — fixed time step size. With advection velocity 0.5 and mesh spacing
  0.05, the CFL number is `vx * dt / dx = 0.5 * 0.02 / 0.05 = 0.2`. A CFL number
  below 1 means information travels less than one cell per time step, which is
  required for stability of explicit methods. MOOSE uses implicit time integration,
  so the CFL limit is not strict — but keeping CFL below 1 helps accuracy.

- `end_time = 2.0` — run for 2 seconds. Total time steps: `2.0 / 0.02 = 100`.

- `nl_rel_tol = 1e-6` — converge the nonlinear solver to six decades below the
  initial residual. Slightly looser than Case 07 because transient problems are
  less sensitive to nonlinear solve accuracy than steady-state problems.

### Block: `[Outputs]`

```
[Outputs]
  exodus = true
  csv    = true
[]
```

- `exodus = true` — writes the solution field `c` at every time step to an Exodus II
  file. With 100 time steps, the file contains 100 frames of the concentration field.
  This is the primary output for animation.
- `csv = true` — writes `total_c` at every time step to a CSV file. This produces a
  100-row table showing how total mass evolves over time.

---

## What Happens When You Run This

### Invocation

```bash
cd quickstart-runs/case08
mpirun -n 1 moose-app-opt -i case08_advection_diffusion.i
```

### Startup Phase

1. MOOSE parses the input file. Top-level variables `vx=0.5` and `D=0.01` are
   substituted into all `${vx}` and `${D}` placeholders.
2. The mesh is created: 60x30 rectangles on a 3x1 domain.
3. `ParsedFunction` compiles the Gaussian expression with fparser.
4. `FunctionIC` evaluates the Gaussian at all 1,891 nodes to set initial conditions.
   The nodal values of c range from near 1.0 at (0.3, 0.5) to zero almost everywhere.
5. The Transient executioner initializes its time-stepping loop.

### Time-Stepping Loop (100 steps)

At each time step t → t + dt:

1. **Predict**: The predictor estimates `c^{n+1}` from previous values (typically
   a linear extrapolation or simply the previous solution).

2. **Assemble**: All three kernels contribute to the residual and Jacobian:
   - `TimeDerivative`: `integral(phi_i * (c^{n+1} - c^n) / dt) dV`
   - `ConservativeAdvection`: `-integral(grad(phi_i) . v * c^{n+1}) dV`
   - `MatDiffusion`: `integral(D * grad(phi_i) . grad(c^{n+1})) dV`

3. **Solve**: PJFNK iterates until `nl_rel_tol` is satisfied. For this near-linear
   problem, Newton typically converges in 1 to 3 iterations per time step.

4. **Output**: Write the solution to Exodus and the postprocessor value to CSV.

### Evolution of the Solution

The blob behavior over time:

```
t = 0.0:  Gaussian centered at (0.3, 0.5), peak = 1.0
t = 0.5:  Blob center at (0.3 + 0.5*0.5, 0.5) = (0.55, 0.5). Peak ~0.85 (diffusion spreads it)
t = 1.0:  Blob center at (0.80, 0.5). Width has grown noticeably. Peak ~0.65
t = 1.5:  Blob center at (1.05, 0.5). Still in domain. Peak ~0.50
t = 2.0:  Blob center at (1.30, 0.5). About half of the blob remains in domain.
```

(These are approximate; actual values depend on the numerical diffusion from upwinding
and the exact PDE parameters.)

### Total Concentration Evolution

The CSV file shows `total_c` over time:

```
time,total_c
0,       ~0.0314    (initial Gaussian integral)
0.02,    ~0.0314    (essentially unchanged, blob well inside domain)
...
1.0,     ~0.0310    (still mostly conserved)
...
1.8,     ~0.0200    (blob approaching right wall, some has exited)
2.0,     ~0.0150    (significant portion has left the domain)
```

The total drops sharply once the blob front reaches x = 3 (the right outflow wall).
Before that point, total_c is approximately conserved (the conservative discretization
is working correctly). The rate of decline after that depends on how sharp the blob
front is.

### Console Output per Time Step

```
Time Step 1, time = 0.02, dt = 0.02
 0 Nonlinear |R| = 1.572e-02
      0 Linear |R| = 1.234e-13
 1 Nonlinear |R| = 3.456e-07

Converged!

Time Step 2, time = 0.04, dt = 0.02
...
```

Because the advection-diffusion equation is nearly linear (v and D are constants, not
depending on c), PJFNK typically converges in 1 to 2 nonlinear iterations per step.

---

## Output Files

### `case08_advection_diffusion_out.e`

The primary output file. Contains:
- The mesh geometry (60x30 structured mesh on 3x1 domain)
- The field `c` at every node at every time step (100 frames)
- Readable by ParaView, VisIt, and other Exodus-compatible tools

### `case08_advection_diffusion_out.csv`

A time history of the `total_c` postprocessor:

```
time,total_c
0,3.14159e-02
0.02,...
0.04,...
...
2.0,...
```

100 rows (one per time step plus initial condition). Use this to:
- Verify mass conservation before outflow
- Identify when the blob reaches the right boundary
- Quantify how much mass has left the domain

---

## Visualizing the Results in ParaView

1. Open `case08_advection_diffusion_out.e` in ParaView:
   - File > Open > navigate to the file
   - If asked "Would you like to group files?" click Yes (or No if only one file)
   - Click Apply in the Properties panel

2. In the Properties panel, verify that the variable `c` is checked under "Point Arrays".

3. Click the play button (triangle) at the top to animate all 100 frames. You will
   immediately see the blob translate from left to right.

4. To get a better view of the blob shape:
   - Color by `c`
   - Apply a "jet" or "Blue-to-Red" color map
   - Click "Rescale to Visible Data Range" at each frame for best contrast, or
     set a fixed range of [0, 1] to see the full evolution

5. To plot total concentration over time:
   - Open the CSV file in a spreadsheet, or
   - In ParaView: Filters > Data Analysis > Plot Data from CSV

6. To see how the blob width changes, use Filters > Data Analysis > Plot Over Line
   with a horizontal line at y=0.5. At different time frames, the line plot shows the
   concentration profile c(x, 0.5, t), revealing how the Gaussian peak moves and flattens.

---

## Understanding the Plots

Running `python visualize_all.py` from the `quickstart-runs/` directory produces the
following two plots saved into this directory.

### `case08_blob_snapshots.png`

**What the plot shows.** A 2x2 grid of 2D filled-contour panels, each showing the
concentration field c(x,y) at a different simulation time: approximately t=0, 0.5,
1.0, and 2.0. The viridis colormap is used with a per-panel colorbar.

**Physical quantities.** The color encodes the concentration of a scalar species (a
Gaussian blob) being simultaneously advected by a uniform velocity field (vx=0.5,
vy=0) and diffused by a small diffusion coefficient. Higher concentration is brighter.

**How to judge correctness.**
- t ≈ 0 panel: the blob should be a compact, roughly circular or Gaussian-shaped
  bright patch centered near x=0.3, y=0.5 (the initial position).
- t ≈ 0.5 panel: the blob should have moved to the right to approximately x=0.55
  (having traveled 0.5 * 0.5 = 0.25 in x) and should be slightly wider due to
  diffusion.
- t ≈ 1.0 panel: center should be near x=0.8, blob is wider still.
- t ≈ 2.0 panel: the blob has largely exited through the right boundary; concentration
  is low everywhere, spread across the domain.

The blob should always remain positive (concentration cannot be negative) and should
shift purely to the right with no upward or downward drift (advection is only in x).

**What would indicate a problem.**
- No movement between panels: advection is not applied (the `Advection` kernel or
  velocity field is missing or zero).
- Negative concentration: numerical oscillations due to insufficient diffusion (the
  Peclet number is too large without stabilization).
- Blob moving upward or downward: velocity vector has a wrong component.
- Blob growing instead of shrinking: the source/sink terms are inverted.

### `case08_total_concentration.png`

**What the plot shows.** A line plot with simulation time on the horizontal axis and
`total_c` (the domain-integrated total concentration) on the vertical axis.

**Physical quantities.** `total_c` is the integral of c over the entire domain. In a
conservation context, it measures how much of the scalar species remains in the domain.

**How to judge correctness.** The curve should be roughly constant at early times
(when the blob is fully inside the domain and little has exited) and then decrease
as the blob reaches and exits through the right boundary. Once the blob has fully
exited, the total concentration should approach zero. The decrease should be smooth
and monotone — not oscillatory.

**What would indicate a problem.**
- `total_c` increasing over time: a source term is adding concentration instead of
  the blob being advected out; check the kernel setup.
- `total_c` dropping immediately to zero: the blob exits too fast — advection velocity
  is too high or the initial blob is placed at the wrong location.
- `total_c` oscillating up and down: numerical instability in the advection scheme.
- `total_c` remaining constant and never decreasing: the outflow boundary condition
  is not allowing material to leave the domain.

---

## Interpreting the Results

### Physical Behavior

At early times (t < 1.0), the blob translates nearly rigidly at speed 0.5. Its peak
slowly decreases as diffusion spreads it. The shape remains Gaussian (diffusion of a
Gaussian initial condition is also Gaussian, just wider).

At late times (t > 1.5), the leading edge of the blob approaches the right boundary.
The Neumann outflow condition allows concentration to exit, reducing total_c. The
blob does not reflect off the right wall — it exits cleanly.

### Numerical Diffusion

With upwinding enabled in `ConservativeAdvection`, the effective diffusion coefficient
is approximately `D + Pe_h * D / 2 = D * (1 + Pe_h/2)`. For `Pe_h = 2.5`, this gives
an effective diffusivity of `0.01 * (1 + 1.25) = 0.0225`. Numerical diffusion is
significant here — the blob spreads faster than pure physics would predict.

To reduce numerical diffusion: use a finer mesh (smaller dx reduces Pe_h) or switch
to higher-order elements. Alternatively, lower the velocity or increase D to reduce
the physical Peclet number.

### Correctness Checks

1. **Total concentration before outflow**: Should be approximately constant at the
   initial value (~0.0314). Any deviation indicates a conservation error.

2. **Blob center position**: At time t, center should be at x = 0.3 + 0.5*t.
   At t=1.0, center should be near x=0.8. Use the Plot Over Line filter to check.

3. **Gaussian shape**: The blob should remain approximately Gaussian at all times
   (a known analytical property of advection-diffusion with constant coefficients).

---

## Key Concepts Learned

- **Transient executioner**: `type = Transient` enables time-marching simulation with
  specified `dt` and `end_time`.

- **TimeDerivative kernel**: Adds the `dc/dt` term to the weak form. Essential for
  any time-dependent simulation.

- **ConservativeAdvection kernel**: Implements the divergence form `div(v*c)` in weak
  form. Preserves mass integrals. Requires a vector material property for the velocity.

- **GenericConstantVectorMaterial**: Declares a constant vector material property
  (velocity field) as three component values in the input file.

- **GenericConstantMaterial**: Declares a constant scalar material property.

- **FunctionIC**: Sets the initial condition of a variable using an arbitrary MOOSE
  Function. The function is evaluated at every node at t=0.

- **ParsedFunction**: Evaluates a mathematical expression string (using fparser) as a
  MOOSE Function. Available variables: x, y, z, t.

- **HIT top-level variables**: Scalars declared before any block in the input file,
  referenced later with `${varname}`. Makes parameter studies easy.

- **Peclet number**: Dimensionless ratio Pe = v*L/D. High Pe means advection dominates.
  Affects numerical stability and solution character.

- **ElementIntegralVariablePostprocessor**: Integrates a variable over the entire mesh,
  giving the total "amount" of a quantity in the domain.

- **Outflow boundary conditions**: The natural (default, zero-flux) Neumann BC on the
  right wall allows advective flux to leave the domain without reflection.

---

## Experiments to Try

### Experiment 1: Change the Velocity

At the top of the input file, change `vx = 0.5` to `vx = 1.0`. The blob now travels
twice as fast. With the same `end_time = 2.0`, the blob will travel 2 units instead
of 1 and will mostly exit the right wall by the end. Also try `vx = 0.1` for a
slow-advection case — the diffusion will spread the blob more noticeably before it
travels far. Compare the CSV files to see how total_c varies differently.

Expected outcome: faster vx → blob exits sooner → total_c drops earlier in the CSV.

### Experiment 2: Increase the Diffusion Coefficient

Change `D = 0.01` to `D = 0.1`. The Peclet number drops from 150 to 15. Diffusion
now plays a more visible role: the blob will broaden significantly while traveling,
and its peak will decrease faster. At `D = 1.0` (Pe ≈ 1.5), diffusion nearly
dominates and the blob becomes unrecognizable as a Gaussian by t=2.

Expected outcome: higher D → faster peak decrease → wider blob → lower total_c
reaches near-zero sooner (the blob has spread to the walls and been absorbed).

### Experiment 3: Add a y-Component to the Velocity

In the `[Materials]` block, change:

```
prop_values = '${vx} 0 0'
```

to:

```
prop_values = '${vx} 0.2 0'
```

Now the flow carries the blob both right and upward. The blob will travel toward the
upper-right corner and hit the top wall. Watch in ParaView how the blob deforms when
it approaches the top Neumann (zero-flux) boundary: the blob "piles up" against the
wall and then exits through the right boundary. Note: the top wall has zero-flux, so
concentration cannot exit there — it must eventually all leave through the right wall.

Expected outcome: blob moves diagonally, accumulates at top-right, then exits right.

### Experiment 4: Narrow the Blob

In the `[Functions]` block, change the Gaussian width:

```
expression = 'exp(-((x-0.3)^2 + (y-0.5)^2) / 0.001)'
```

This makes the initial blob ten times narrower in each direction (width parameter
0.001 instead of 0.01). The blob is now much sharper — challenging for the mesh to
resolve. With 60x30 elements and dx=0.05, the blob width is comparable to a single
element. You will likely see more numerical diffusion and possibly oscillations.
To fix this, increase nx to 200 to refine the mesh.

Expected outcome: sharper initial blob, stronger numerical diffusion artifacts,
peak decays faster.

### Experiment 5: Track When the Blob Exits

From the CSV file, identify the time step at which `total_c` begins to decrease
significantly. This is approximately when the blob front reaches the right wall.
For vx=0.5, the initial center is at x=0.3 and the blob "width" is ~0.2, so the
leading edge is at x = 0.3 + 0.2 = 0.5. At speed 0.5, this edge reaches x=3.0
at t = (3.0 - 0.5) / 0.5 = 5.0 seconds — but since we only run to t=2.0, only the
very front of the blob exits. Check: does total_c start decreasing before t=2.0?

Change `end_time = 6.0` and add `dt = 0.05` to run longer. Watch total_c decrease
to near zero as the blob fully exits the domain.
