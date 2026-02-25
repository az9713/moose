# Case 12: MultiApp Coupling (Parent-Sub Application System)

## Overview

This case introduces MOOSE's MultiApp system — one of the most architecturally
distinctive features of the framework. Instead of solving everything in a single
input file, MOOSE allows you to compose a simulation from multiple independent
applications: a parent application that manages the overall execution, and one or
more sub-applications that each solve their own physics. Data is exchanged between
them via Transfers.

The physical scenario is a one-way coupled system:

1. The **parent application** solves a steady-state Poisson/diffusion problem for
   temperature T on a 20x20 mesh. This represents a thermal analysis.

2. The **sub-application** receives T from the parent (via a Transfer), then uses
   it as a forcing function to drive a diffusion equation for a scalar field phi.
   This represents a downstream physics calculation — perhaps neutron flux, species
   concentration, or structural displacement — that is influenced by temperature but
   does not feed back into the thermal calculation.

The coupling is one-way: T influences phi, but phi does not change T. This is the
simplest and most common form of multi-physics coupling and is a natural starting
point before tackling two-way (tightly coupled) systems.

---

## The Physics

### The Parent Problem (Thermal)

The parent solves a steady Poisson equation for temperature T:

    -div( grad(T) ) = Q      in Omega = [0,1] x [0,1]

where Q = 1.0 (uniform source) and T = 0 on all boundaries. This is the same
problem type seen in earlier cases. The solution T is a smooth hill shape: T = 0
at all walls, rising to a maximum at the center.

### The Sub Problem (Coupled Diffusion)

The sub solves another steady Poisson equation, this time for phi:

    -div( grad(phi) ) = 0.1 * T_from_parent      in Omega

with phi = 0 on all boundaries. The source term `0.1 * T` means the sub problem's
forcing is proportional to the temperature field received from the parent. Because
T is largest at the center, the source for phi is also largest at the center, so phi
also develops a hill-shaped profile — but it is driven by the parent's solution, not
by a user-specified constant.

### Why Multi-Application Coupling?

The MultiApp system solves a fundamental problem in large-scale simulation: different
physics often have very different numerical requirements (different meshes, different
solvers, different time scales, different programming teams) and should not be
forcefully merged into one monolithic system. Instead, MOOSE allows each physics
to live in its own input file with its own solver settings, and data passes between
them through well-defined transfer interfaces.

Real-world examples of one-way coupled analyses:
- Thermal -> structural (thermal expansion from temperature field)
- Neutronics -> thermal (heat deposition from neutron flux)
- Fluid flow -> species transport (advection field drives concentration)
- Macro-scale -> micro-scale (homogenized quantities drive fine-scale sub-models)

### Domain Relationship

Both parent and sub use identical 20x20 meshes on the unit square. The sub-app's
origin in parent coordinates is (0, 0, 0), meaning the meshes overlap exactly. This
is important for the `MultiAppCopyTransfer`, which copies nodal values directly from
parent nodes to sub nodes — the meshes must be topologically compatible for this
transfer type.

### ASCII Diagram of the Two-Application System

```
PARENT APPLICATION (case12_parent.i)
+--------------------------------------+
|  Mesh: 20x20 on [0,1]x[0,1]         |
|  Variable: T                         |
|  Kernels:                            |
|    -div(grad T) = 1.0                |
|  BCs: T=0 on all walls               |
|                                      |
|  Solve --> T field (smooth hill)     |
+-------------------+------------------+
                    |
                    | [Transfers]
                    | MultiAppCopyTransfer
                    | Parent T --> Sub T_from_parent
                    |
                    v
SUB-APPLICATION (case12_sub.i)
+--------------------------------------+
|  Mesh: 20x20 on [0,1]x[0,1]         |
|  Variable: phi                       |
|  AuxVariable: T_from_parent (read-only, set by Transfer)|
|  Kernels:                            |
|    -div(grad phi) = 0.1 * T_from_parent|
|  BCs: phi=0 on all walls             |
|                                      |
|  Solve --> phi field                 |
|  (hill shape driven by T)            |
+--------------------------------------+

Execution order:
  Step 1: Parent solves for T
  Step 2: Transfer copies T -> T_from_parent
  Step 3: Sub solves for phi using T_from_parent as source
```

---

## Input File Walkthrough: Parent (case12_parent.i)

### `[Mesh]` Block

```
[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 20
  ny   = 20
[]
```

A 20x20 mesh on the unit square. This is the parent's computational domain. The
boundaries are named left, right, top, bottom automatically.

### `[Variables]` Block

```
[Variables]
  [T]
  []
[]
```

The parent declares one degree of freedom: T (temperature). This is a standard
first-order Lagrange nodal variable.

### `[Kernels]` Block

```
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
```

Two kernels:

**`Diffusion`**: Assembles `integral( grad(T) . grad(v) ) dOmega`. Represents
`-div(grad T)` in the equation.

**`BodyForce`**: Assembles `integral( 1.0 * v ) dOmega`. Represents the constant
source Q = 1.0 on the right-hand side.

Together: solve `-div(grad T) = 1.0` with T=0 BCs.

### `[BCs]` Block

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

All four walls held at T = 0. The solution T will be a smooth dome peaking at the
center.

### `[Executioner]` Block

```
[Executioner]
  type       = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]
```

Steady-state solve. This runs once and produces the T field that will be transferred
to the sub-application.

### `[MultiApps]` Block — Launching the Sub-Application

```
[MultiApps]
  [thermal_sub]
    type        = FullSolveMultiApp
    input_files = case12_sub.i
    execute_on  = timestep_end
    positions   = '0 0 0'
  []
[]
```

This is where the sub-application is configured. Let's examine each parameter:

**`type = FullSolveMultiApp`**: Runs the sub-application to full convergence (complete
solve) at the specified execution point. This is appropriate for steady-state sub-apps.
For transient sub-apps, `TransientMultiApp` would be used instead, allowing the parent
and sub to march forward in time together.

**`input_files = case12_sub.i`**: The path to the sub-application's input file. MOOSE
reads this file, creates a separate `MooseApp` instance, and runs it as a child
process. From the sub-app's perspective, it is a completely normal MOOSE simulation —
it does not know it is being run as a sub-app.

**`execute_on = timestep_end`**: Controls when the sub-app runs. For a `Steady`
executioner, `timestep_end` means "after the parent's solve is complete." This ensures
the parent has a converged T field before the sub-app starts. Other options include
`initial` (before any parent solve) or `timestep_begin` (before each parent solve).

**`positions = '0 0 0'`**: The position of the sub-app's origin in the parent's
coordinate system. With `0 0 0`, the sub-app's mesh is placed with its lower-left
corner at the parent's origin — the two meshes overlap exactly. If you specified
`positions = '2 0 0'`, the sub-app's coordinate system would be shifted 2 units
in x relative to the parent.

Multiple sub-apps can be run by providing multiple positions:
`positions = '0 0 0    1 0 0    2 0 0'` would run three copies of the sub-app at
three different positions. Each gets its own independent mesh and solve.

### `[Transfers]` Block — Moving Data Between Apps

```
[Transfers]
  [send_temperature]
    type            = MultiAppCopyTransfer
    to_multi_app    = thermal_sub
    source_variable = T
    variable        = T_from_parent
  []
[]
```

The Transfer block defines how data moves from the parent to the sub (or back).

**`type = MultiAppCopyTransfer`**: Performs a direct node-to-node copy of variable
values. This requires that the parent and sub meshes are topologically identical (same
number of nodes, same connectivity). It is the fastest and simplest transfer type.
Other transfer types (like `MultiAppInterpolationTransfer`) handle non-matching meshes
by interpolating.

**`to_multi_app = thermal_sub`**: Specifies the direction of transfer. `to_multi_app`
means data flows FROM the parent INTO the sub-app named `thermal_sub`. The counterpart
direction would be `from_multi_app`, for sub-to-parent transfer.

**`source_variable = T`**: The variable in the parent that provides the data. This
is the T field that the parent just solved for.

**`variable = T_from_parent`**: The variable in the sub-application that receives the
data. This must be declared in the sub-app's input file — see the sub-app's
`[AuxVariables]` block.

The transfer executes automatically at `timestep_end` (because that's when the
sub-app executes), before the sub-app begins its solve.

### `[Outputs]` Block

```
[Outputs]
  exodus = true
[]
```

The parent writes its T field to an Exodus file: `case12_parent_out.e`.

---

## Input File Walkthrough: Sub-Application (case12_sub.i)

### `[Mesh]` Block

```
[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 20
  ny   = 20
[]
```

Identical mesh to the parent. This is required for `MultiAppCopyTransfer` to work
correctly.

### `[Variables]` Block

```
[Variables]
  [phi]
  []
[]
```

The sub-app has its own independent degree of freedom `phi`. This is the field the
sub-app solves for. It has nothing to do with the parent's T; it is the sub-app's
primary unknown.

### `[AuxVariables]` Block — Receiving Data from the Parent

```
[AuxVariables]
  [T_from_parent]
    order  = FIRST
    family = LAGRANGE
  []
[]
```

This is the critical connection between parent and sub. `AuxVariables` are variables
that are NOT solved for — they are either computed by `AuxKernels` or, in this case,
populated by the parent's Transfer.

`T_from_parent` must match the basis (FIRST LAGRANGE) of the parent's T variable so
that the node-to-node copy transfer works. From the sub-app's perspective, `T_from_parent`
is just a read-only field that magically contains the parent's T solution after the
transfer runs.

The `order = FIRST` and `family = LAGRANGE` declarations are explicit here (they are
also the default) to make the coupling requirement visible and clear.

### `[Kernels]` Block

```
[Kernels]
  [diffusion]
    type     = Diffusion
    variable = phi
  []
  [coupling]
    type     = CoupledForce
    variable = phi
    v        = T_from_parent
    coef     = 0.1
  []
[]
```

Two kernels:

**`Diffusion`**: The standard `-div(grad phi)` term.

**`CoupledForce`**: Adds `coef * v` as a source term, where `v` is another variable
(or AuxVariable) in the same system. Here it adds `0.1 * T_from_parent` as the right-
hand side source for phi. This is the coupling: the thermal field T drives the phi
field.

The weak form is: `integral( grad(phi) . grad(v) ) = integral( 0.1 * T_from_parent * v )`

This is a linear system: for a given T_from_parent, solving for phi is a single linear
solve (one Newton iteration). The "nonlinearity" would only arise if phi fed back into T,
which it does not in this one-way coupling.

### `[BCs]` Block

```
[BCs]
  [walls]
    type     = DirichletBC
    variable = phi
    boundary = 'left right top bottom'
    value    = 0
  []
[]
```

phi = 0 on all walls. Same structure as the parent's BCs for T.

### `[Executioner]` Block

```
[Executioner]
  type       = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]
```

Identical to the parent's executioner. Each sub-app has its own independent executioner.
In a more complex setup, the sub-app could be Transient while the parent is Steady.

### `[Outputs]` Block

```
[Outputs]
  exodus = true
[]
```

The sub-app writes its phi field (and T_from_parent as an aux variable) to a separate
Exodus file.

---

## Data Flow: Step by Step

Understanding exactly when things happen is essential for MultiApp problems. Here is
the complete execution sequence:

```
1. Parent reads case12_parent.i and initializes its mesh, variables, and objects.

2. Parent runs [Executioner] type=Steady:
   - Assembles stiffness matrix and load vector for T
   - Calls PJFNK solver
   - Converges: T field now contains the dome-shaped solution

3. Parent reaches "timestep_end":
   - [MultiApps] block: FullSolveMultiApp launches case12_sub.i
     * Sub-app creates its own mesh (20x20)
     * Sub-app initializes phi = 0, T_from_parent = 0 everywhere
   - [Transfers] block: MultiAppCopyTransfer executes
     * Reads T values at all 441 nodes (21x21 = 441) in the parent
     * Writes those exact values into T_from_parent in the sub-app
     * T_from_parent now contains the parent's T solution

4. Sub-app runs [Executioner] type=Steady:
   - Assembles stiffness matrix for phi with source = 0.1 * T_from_parent
   - Calls PJFNK solver
   - Converges: phi field contains the solution driven by T

5. Sub-app writes case12_parent_out_thermal_sub0.e
   (contains both phi and T_from_parent fields)

6. Parent writes case12_parent_out.e
   (contains the T field)

7. Simulation complete.
```

---

## What Happens When You Run This

Run only the parent input file:

```bash
./moose_test-opt -i case12_parent.i
```

MOOSE automatically handles the sub-app. You only specify the parent; the parent
launches the sub-app internally via the `[MultiApps]` block.

The console output shows both the parent and sub-app solves:

```
Framework Information:
  MOOSE version: ...

Mesh Information:
  ...

Solving...
  0 Nonlinear |R| = 2.5e+00
  1 Nonlinear |R| = 8.1e-09
  Nonlinear solve converged.

Executing MultiApps...
  Executing thermal_sub (FullSolveMultiApp) on timestep_end

  Sub-app: thermal_sub (App 0)
    Mesh Information: ...

    Transferring T -> T_from_parent

    Solving sub-app...
      0 Nonlinear |R| = 1.8e-01
      1 Nonlinear |R| = 5.2e-10
      Nonlinear solve converged.

Simulation Complete.
```

The parent solve appears first (it completes in ~1 Newton iteration since this is
a linear problem). Then MOOSE launches the sub-app, executes the transfer, and runs
the sub-app's solve.

---

## Output Files

| File | Description |
|------|-------------|
| `case12_parent_out.e` | Parent Exodus: T field on the 20x20 mesh |
| `case12_parent_out_thermal_sub0.e` | Sub-app Exodus: phi and T_from_parent fields |

### Output File Naming Convention

The sub-app output file is named using the pattern:
`{parent_output_base}_{multiapp_name}{app_index}.e`

- Parent output base: `case12_parent_out`
- MultiApp name: `thermal_sub`
- App index: `0` (first/only instance)

If three sub-apps were running (e.g., `positions = '0 0 0   1 0 0   2 0 0'`), the
outputs would be `..._thermal_sub0.e`, `..._thermal_sub1.e`, `..._thermal_sub2.e`.

### Visualizing in ParaView

To compare parent and sub solutions side by side:

1. Open both `case12_parent_out.e` and `case12_parent_out_thermal_sub0.e` in ParaView.
2. Click Apply for each.
3. For the parent file: color by T. You will see the dome-shaped temperature field.
4. For the sub file: color by phi. You will see a similar dome shape, scaled by 0.1
   (because the coupling coefficient is 0.1) and slightly different in shape due to
   the diffusion smoothing.
5. Also color the sub file by T_from_parent to verify it matches the parent's T.

---

## Understanding the Plots

Running `python visualize_all.py` from the `quickstart-runs/` directory produces the
following two plots saved into this directory.

### `case12_parent_temperature.png`

**What the plot shows.** A 2D filled-contour of the parent application's temperature
field T(x,y) at the final timestep, using the coolwarm colormap.

**Physical quantities.** The color encodes the temperature T in the parent simulation.
The parent solves a Poisson equation with a volumetric heat source and zero-Dirichlet
walls. The parent T field is also transferred to the sub-application as boundary data.

**How to judge correctness.** The plot should show a smooth dome peaked at the domain
center (0.5, 0.5) with zero temperature on all four walls. This is the same shape
as Case 03/11. The peak value depends on the source strength used in the parent input
file.

**What would indicate a problem.**
- Flat field (T=0 everywhere): the parent solve did not run, or the output file was
  not written.
- Non-symmetric dome: the boundary conditions in the parent are incorrect on one side.

### `case12_sub_phi.png`

**What the plot shows.** A 2D filled-contour of the sub-application's field phi(x,y)
at the final timestep, using the viridis colormap.

**Physical quantities.** The sub-application receives the parent temperature T as
input (via a MultiAppInterpolationTransfer) and uses it to drive its own field phi.
The physics of the sub-app determines how phi responds to T — typically phi is
proportional to T, for example phi = 0.1 * T or a similarly scaled version.

**How to judge correctness.** The phi field should have the same shape as the parent
T field — a dome peaked at the center, zero on the walls. The magnitude should be
scaled relative to T by whatever relationship the sub-app physics implements. If the
sub-app sets phi proportional to T with a coefficient of 0.1, the peak of phi should
be 0.1 times the peak of T.

Crucially: if the MultiApp transfer is working correctly, the spatial pattern of phi
should clearly mirror the pattern of T. A phi field that is completely flat or uniform
(all the same value) means the transfer failed and the sub-app is not receiving T
from the parent.

**What would indicate a problem.**
- phi is uniformly zero or uniformly constant: the transfer from parent to sub-app
  is not working. Check that the sub-app output file was created (`case12_parent_out_thermal_sub0.e`)
  and that the `MultiAppInterpolationTransfer` block is correctly configured.
- phi has the wrong shape relative to T: the transfer is delivering values to the
  wrong variable, or the sub-app is computing phi from something other than T.
- phi file is missing entirely: the sub-app did not run or its output is in a
  different file than expected.

---

## Interpreting the Results

### Parent Solution (T field)

The parent solves `-div(grad T) = 1` with T=0 on all walls. The solution is the
classical "tent function" or "drum" solution: T is zero at all four edges and rises
to a maximum at the center (0.5, 0.5). For a unit square with Q=1 and unit
conductivity, the maximum temperature is approximately T_max ~ 0.073.

The solution is symmetric about both the x=0.5 and y=0.5 lines.

### Sub Solution (phi field)

The sub solves `-div(grad phi) = 0.1 * T` with phi=0 on all walls. Since T peaks
at the center, the source for phi also peaks at the center. The resulting phi field
is another dome shape, but:

- Peak value of phi is smaller than peak of T (the coupling coefficient 0.1 and the
  additional diffusion both reduce the amplitude)
- The phi dome is "smoother" than T because it results from diffusing an already-
  smooth source (T itself is smooth)
- phi is also symmetric about x=0.5 and y=0.5

### The T_from_parent Field

The sub-app's Exodus file contains T_from_parent as an auxiliary variable. It should
be visually identical to the parent's T field. Verifying this confirms that the
transfer worked correctly.

### Correctness Verification

Three checks:
1. `T_from_parent` in sub-app should match `T` in parent (within floating-point precision)
2. Both T and phi should be symmetric about x=0.5 and y=0.5
3. Both fields should be zero at all boundaries
4. The maximum value of phi should be approximately 0.1 * max(T) / 8 ~ 0.001 (rough
   estimate, actual value depends on the Green's function of the diffusion operator)

---

## Key Concepts Learned

| Concept | What It Is | Where |
|---------|-----------|-------|
| MultiApp system | Multiple coupled MOOSE applications | `[MultiApps]` block in parent |
| FullSolveMultiApp | Runs sub-app to convergence once | `type = FullSolveMultiApp` |
| execute_on | Controls when sub-app runs | `execute_on = timestep_end` |
| positions | Sub-app origin in parent coordinates | `positions = '0 0 0'` |
| Transfers | Data exchange between parent and sub | `[Transfers]` block in parent |
| MultiAppCopyTransfer | Node-to-node copy on matching meshes | `type = MultiAppCopyTransfer` |
| to_multi_app | Parent-to-sub transfer direction | `to_multi_app = thermal_sub` |
| AuxVariables | Variables populated externally (not solved for) | `[AuxVariables]` in sub |
| CoupledForce | Source term driven by another variable | `[Kernels]` in sub |
| One-way coupling | T influences phi; phi does not affect T | Physical design choice |

---

## Experiments to Try

**Experiment 1: Change the coupling coefficient**
In `case12_sub.i`, change `coef = 0.1` in the `CoupledForce` kernel to `coef = 1.0`.
The phi solution should be about 10x larger. Does the shape change? (It should not —
scaling the source scales the solution linearly.)

**Experiment 2: Add a from_multi_app Transfer (two-way coupling)**
Add a second Transfer to `case12_parent.i`:
```
[send_phi_back]
  type           = MultiAppCopyTransfer
  from_multi_app = thermal_sub
  source_variable = phi
  variable        = phi_from_sub
[]
```
And add `[AuxVariables] [phi_from_sub] []` to the parent. After running, the parent's
Exodus file will contain the sub's phi solution. This is the first step toward two-way
coupling.

**Experiment 3: Use a different mesh in the sub-app**
Change the sub-app's mesh from `nx=20 ny=20` to `nx=10 ny=10`. You will need to
change the transfer type to `MultiAppInterpolationTransfer` (which handles non-matching
meshes) instead of `MultiAppCopyTransfer`. This demonstrates the power of the transfer
abstraction.

**Experiment 4: Run multiple sub-app instances**
Change the `positions` parameter in `case12_parent.i` to:
`positions = '0 0 0'`
and add a second sub-app input file that solves a different problem (e.g., with a
Neumann BC on one wall). This demonstrates running multiple parallel sub-apps.

**Experiment 5: Change execute_on**
Change `execute_on = timestep_end` to `execute_on = initial` in `case12_parent.i`.
Now the sub-app runs before the parent solves. The transfer will send T=0 (the initial
condition) to the sub. The sub will solve with T_from_parent = 0 everywhere, so its
solution will be trivially zero. This illustrates how the execute_on timing controls
the coupling direction.
