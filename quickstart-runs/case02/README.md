# Case 02: Steady-State 2D Diffusion

## Overview

Case 02 extends the 1D diffusion problem from Case 01 into two spatial dimensions.
The governing equation and solver configuration are identical — only the mesh changes
from a line segment to a square, and the boundary conditions are adjusted to match the
2D domain.

Despite that small change, a great deal happens internally: the mesh is now a grid of
quadrilateral elements instead of line segments, the stiffness matrix is larger and
sparser, and the 2D visualization requires different tools than a simple line plot.

What you will learn:
- How MOOSE generalizes seamlessly from 1D to 2D (and beyond) without changing the kernel
- The concept of **natural boundary conditions** (zero-flux Neumann) and how they arise
  automatically when you omit a `[BCs]` entry for a boundary
- The named boundary convention for 2D `GeneratedMesh` (`left`, `right`, `bottom`, `top`)
- How to read and interpret a 2D color-map visualization in ParaView
- Why the exact solution `u(x,y) = x` has no y-dependence even though the domain is 2D

The exact solution to this case is again `u(x,y) = x` — a planar ramp that varies
only in the x-direction. The top and bottom boundaries impose no constraint because
the solution naturally has zero slope in the y-direction. This provides another
opportunity to verify the result by inspection.

---

## The Physics

### The Physical Problem in Plain English

Imagine a flat square plate, one meter on each side, lying in the (x, y) plane.
The left edge of the plate (x = 0) is held at 0 degrees — perhaps cooled by a cold
water pipe running along that edge. The right edge (x = 1) is held at 1 degree — heated
by a warm pipe. The top edge and the bottom edge are **insulated**: no heat can flow
through them perpendicular to the plate edge.

No heat is generated inside the plate. After enough time, the temperature reaches a
steady state: every point in the plate stops changing temperature.

What is the temperature distribution at steady state?

Because the top and bottom edges are insulated (no heat escapes from them), heat flows
purely in the x-direction — from the hot right edge toward the cold left edge. There is
no reason for the temperature to vary in the y-direction, because nothing breaks the
symmetry in y. The temperature profile is the same at every height y: a linear ramp
from 0 to 1 as x goes from 0 to 1.

The solution is `u(x, y) = x` — identical to the 1D case, just extended across the
two-dimensional plate. The insulated top and bottom walls are fully consistent with
this solution: the heat flux normal to these walls is `du/dy = 0`, which is exactly
what "insulated" means.

### The Governing Equation

The PDE is the 2D Laplace equation:

```
-div(grad u) = 0    on the domain [0,1] x [0,1]
```

Written out in components:

```
  d²u     d²u
- ---- - ----- = 0
  dx²     dy²
```

Every symbol:
- `u(x, y)` — the scalar field, a function of both x and y
- `grad u` — the gradient vector: `(du/dx, du/dy)`. Points in the direction of steepest
  increase of u.
- `div(grad u)` — the Laplacian: sum of second partial derivatives. Measures how much u
  at a point differs from the average of u in a small neighborhood.
- `-div(grad u) = 0` — the field u is harmonic: it has no local maxima or minima in the
  interior (maximum principle), and satisfies the averaging property.

The 2D operator reduces to the 1D operator when u does not depend on y:
- `d²u/dy² = 0` (since u = x, which has no y-dependence)
- `d²u/dx² = 0` (since u = x is linear in x, its second derivative is zero)
- The equation is satisfied at every interior point

### Boundary Conditions

The domain has four boundaries. Only two have explicit Dirichlet conditions:

```
u = 0    on x=0 (left boundary)
u = 1    on x=1 (right boundary)
```

The other two boundaries have no explicit condition:

```
du/dn = 0    on y=0 (bottom boundary)   [natural / Neumann]
du/dn = 0    on y=1 (top boundary)      [natural / Neumann]
```

**What is a natural boundary condition?**

When you omit a boundary from the `[BCs]` block entirely, the finite element method
automatically imposes a **zero-flux Neumann condition** on that boundary. This is not
an accident — it falls out of the mathematics of the weak form.

The weak form of `-div(grad u) = 0` is derived by multiplying by a test function `phi_i`
and integrating by parts. Integration by parts introduces a boundary term:

```
integral( grad(phi_i) . grad(u)  dV ) = integral( phi_i * (grad u . n)  dS )
```

The right side is the boundary integral of the normal flux `grad u . n`. When no BC
is applied on a boundary, this term is simply left as zero in the assembly. That
implicitly enforces `grad u . n = 0` — zero flux normal to the boundary.

"Natural" means it arises naturally from the weak form without any special treatment.
Neumann (zero-flux) is the default, and Dirichlet is the exception that requires active
enforcement.

**Physical interpretation of zero flux:**

`grad u . n = 0` on the top and bottom means no heat (or concentration, or anything)
flows through those walls. They are perfectly insulated. This is the appropriate
condition for a problem where the driving force (temperature difference) is purely
horizontal.

### Exact Solution

The exact solution is:

```
u(x, y) = x
```

This is a planar surface (not a 2D curve) that:
- Equals 0 at x=0 for all y — matches left Dirichlet BC
- Equals 1 at x=1 for all y — matches right Dirichlet BC
- Has `du/dy = 0` everywhere — satisfies the natural Neumann condition on top and bottom
- Satisfies `-div(grad u) = -(0 + 0) = 0` — satisfies the PDE

Because `u = x` is linear in x and constant in y, it is exactly representable by
first-order Lagrange elements on any mesh. The numerical solution should match the exact
solution to within machine precision.

### ASCII Domain Diagram

```
  y=1  +--------------------+  Insulated (du/dn = 0)
       |                    |
       |                    |
 u=0   |                    |  u=1
(left) |   u(x,y) = x       | (right)
       |                    |
       |                    |
       |                    |
  y=0  +--------------------+  Insulated (du/dn = 0)
      x=0                  x=1

Color map (top view):
  Blue (u=0)   ----> Cyan ----> Yellow ----> Red (u=1)
  x=0                                        x=1
  Perfectly vertical stripes: u depends only on x

Flux arrows (heat flow direction):
  <---- <---- <---- <---- <---- <---- <---- <----
  Heat flows from right (hot) to left (cold)
  Arrows are horizontal, uniform, everywhere the same
```

---

## Input File Walkthrough

The input file is `case02_diffusion_2d.i`. The structure is nearly identical to Case 01,
but dimensionality is increased from 1D to 2D.

### Header Comments

```
# ============================================================
# Case 2: Steady-State 2-D Diffusion
# Solves -div(grad u) = 0 on the unit square
# Exact solution: u(x,y) = x
# ============================================================
```

The comments document the domain shape ("unit square"), the operator, and the exact
solution. Note the move from the 1D notation `-d²u/dx²` to the coordinate-free
`-div(grad u)`, which is valid in any number of dimensions. This is the notation MOOSE
uses internally.

---

### Block: `[Mesh]`

```
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
```

**`dim = 2`**

This is the only required change from Case 01 to make a 2D mesh. Setting `dim = 2` tells
`GeneratedMesh` to build a 2D quadrilateral mesh instead of a 1D line mesh.

What changes internally:
- Elements are 4-node quadrilaterals (QUAD4) instead of 2-node lines (EDGE2)
- Each element has a 2D Jacobian (area, not length)
- Shape functions are bilinear products of 1D hat functions
- Quadrature uses a 2x2 tensor product of Gauss points per element
- The stiffness matrix grows: each element contributes a 4x4 block instead of 2x2

**`nx = 20` and `ny = 20`**

The mesh has 20 elements in each direction, totaling:
- 20 x 20 = 400 quadrilateral elements
- 21 x 21 = 441 nodes
- 441 degrees of freedom (one DOF per node for variable `u`)

Each element is a square of side 1/20 = 0.05 units. The element aspect ratio is 1:1
(perfectly square), which is ideal for numerical accuracy — distorted elements
degrade accuracy.

**`xmin`, `xmax`, `ymin`, `ymax`**

These four parameters explicitly set the domain bounds. In Case 01, `xmin` and `xmax`
were already present. Adding `ymin = 0` and `ymax = 1` completes the 2D domain
specification. The domain is the unit square [0,1] x [0,1].

**Named boundaries in 2D**

`GeneratedMesh` in 2D automatically assigns four named boundary sets:

| Name | Location | Boundary type |
|------|----------|---------------|
| `left` | x = xmin = 0 | Left edge (a line from y=0 to y=1) |
| `right` | x = xmax = 1 | Right edge (a line from y=0 to y=1) |
| `bottom` | y = ymin = 0 | Bottom edge (a line from x=0 to x=1) |
| `top` | y = ymax = 1 | Top edge (a line from x=0 to x=1) |

Each named boundary is a set of edges (in 2D, a boundary is a 1D curve — the side
of an element). MOOSE identifies which elements touch each boundary and applies BCs
to those edges.

**Scale of the system**

Comparing 1D (Case 01) to 2D (Case 02):

| Quantity | Case 01 (1D) | Case 02 (2D) |
|----------|-------------|-------------|
| Elements | 20 | 400 |
| Nodes | 21 | 441 |
| DOFs | 21 | 441 |
| Matrix size | 21 x 21 | 441 x 441 |
| Matrix nonzeros | ~60 | ~2200 |

The 2D problem has 21x more DOFs and a proportionally larger sparse matrix, but still
solves in well under a second.

---

### Block: `[Variables]`

```
[Variables]
  [u]
  []
[]
```

Identical to Case 01. The variable `u` uses default first-order Lagrange elements.
The only difference is that MOOSE now allocates 441 DOFs (one per node in the 21x21
node grid) instead of 21.

In 2D, the basis functions are bilinear:
- On each quadrilateral element, u is approximated as `u ≈ a + b*x + c*y + d*x*y`
- The four coefficients (a, b, c, d) are determined by the four nodal values
- The function is continuous across element boundaries

The `x*y` term makes bilinear elements slightly richer than tensor products of 1D linear
elements. They can represent tilted planes exactly.

---

### Block: `[Kernels]`

```
[Kernels]
  [diffusion]
    type     = Diffusion   # same kernel, works in any dimension
    variable = u
  []
[]
```

**The kernel is identical to Case 01.** This is one of MOOSE's key design features:
the `Diffusion` kernel is written to work in any spatial dimension without modification.

Internally, the kernel computes:

```
R_i = integral( grad(phi_i) . grad(u)  dV )
```

In 2D, `grad(phi_i)` and `grad(u)` are 2D vectors, and the dot product gives:

```
R_i = integral( d(phi_i)/dx * du/dx  +  d(phi_i)/dy * du/dy  dV )
```

The kernel does not contain any explicit `if (dim == 2)` logic. Instead, MOOSE provides
`_grad_phi[_i][_qp]` and `_grad_u[_qp]` as vectors whose length matches the problem
dimension. The dot product `.dot()` works correctly in any dimension.

The Jacobian contribution per element is:

```
J_ij = d(R_i)/d(u_j) = integral( grad(phi_i) . grad(phi_j)  dV )
```

This is the element stiffness matrix — a 4x4 matrix for each quadrilateral element.
Assembling all 400 elements gives the global 441x441 sparse stiffness matrix.

---

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
  # No BCs on 'top' and 'bottom' => zero-flux Neumann (the natural condition).
[]
```

Two Dirichlet conditions, exactly as in Case 01, but now each boundary is a line (an
edge set) rather than a single point.

**`boundary = left`**

The `left` boundary contains 21 edges (the left sides of the 20 elements in the leftmost
column, plus the corner nodes). The Dirichlet BC applies to the 21 nodes on this edge:
all have their DOF equation replaced by `u = 0`. These nodes form the column x=0,
y = 0, 0.05, 0.10, ..., 1.00.

**`boundary = right`**

Similarly, 21 nodes at x=1 are fixed to u=1.

**`top` and `bottom` are absent**

The comment makes this explicit: not listing `top` and `bottom` in `[BCs]` imposes the
natural zero-flux Neumann condition on those edges. No action is needed in the input
file to get this behavior — it is automatic.

The 21 nodes on the bottom edge (y=0) and 21 nodes on the top edge (y=1) are free
interior DOFs (from the solver's perspective). Their values are determined by the PDE,
not pinned by BCs.

---

### Block: `[Executioner]`

```
[Executioner]
  type       = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]
```

Identical to Case 01. `Steady` + `PJFNK` + `BoomerAMG` is the standard configuration.

The system is still linear (Laplace equation), so Newton converges in one iteration.
The only difference is that the linear system is now 441x441 instead of 21x21.
BoomerAMG handles this without any change to the configuration — it scales to millions
of DOFs with the same settings.

**Why does BoomerAMG work so well for this problem?**

Algebraic multigrid builds a hierarchy of increasingly coarser representations of the
matrix. At each level, errors at the corresponding spatial scale are efficiently damped.
By combining corrections from all levels, AMG eliminates error at all scales in each
iteration. For the Laplace equation, which is a prototypical elliptic PDE, AMG is
provably optimal: the solve time grows linearly with the number of DOFs (O(N)
complexity).

---

### Block: `[Outputs]`

```
[Outputs]
  exodus = true
[]
```

Only the Exodus output is enabled this time (no CSV). Without postprocessors, there
are no scalar quantities to write to CSV anyway.

The output file `case02_diffusion_2d_out.e` stores the full 2D field `u(x, y)` over
all 441 nodes. In 2D, visualization is much more informative than in 1D — you can see
color maps, vector fields, contour lines, and surface plots.

---

## What Happens When You Run This

### Invocation

```bash
cd quickstart-runs/case02
../../test/moose_test-opt -i case02_diffusion_2d.i
```

### Step-by-Step Internally

**1. Mesh Generation**

`GeneratedMesh` creates a 20x20 structured quad mesh. The 441 nodes are at
positions `(i/20, j/20)` for `i, j = 0, 1, ..., 20`. The 400 elements are ordered
row by row. Boundary sidesets `left`, `right`, `bottom`, and `top` are tagged.

**2. DOF Numbering**

MOOSE's DofMap assigns global equation numbers 0 through 440 to the 441 nodes.
The boundary conditions will later constrain equations for the 21 left-edge nodes
and 21 right-edge nodes.

**3. Assembly Loop**

For each of the 400 elements:
- 4 Gauss quadrature points are used (2x2 tensor product)
- At each quadrature point, 4 shape function gradients are computed
- The `Diffusion` kernel contributes a 4x4 element stiffness matrix
- These are accumulated into the 441x441 global sparse matrix

Total assembly: 400 elements x 4 quadrature points x 4x4 entries = 6400 operations.
The global matrix has approximately 2200 nonzero entries (each node connects to its
neighbors, typically 9 for interior nodes in a structured quad mesh).

**4. Boundary Condition Application**

For the 21 left-edge nodes: row `i` of the system is replaced by `1*u_i = 0`.
For the 21 right-edge nodes: row `j` of the system is replaced by `1*u_j = 1`.
The 399 interior and top/bottom nodes have their equations from the kernel assembly.

**5. Linear Solve**

GMRES + BoomerAMG solves the 441x441 system. For this small problem, convergence
is essentially instantaneous. Typical iteration count: 5-10 Krylov iterations.

**6. Solution Check**

At every node (i, j), the computed value should be `u = i/20` (the x-coordinate).
The y-coordinate plays no role.

### Typical Console Output

```
Mesh Information:
  Spatial dimension:      2
  Mesh:                   400 elements, 441 nodes

Nonlinear System:
  Num DOFs:               441

 0 Nonlinear |R| = 1.000000e+00
      0 Linear |R| = 1.000000e+00
      1 Linear |R| = 3.824958e-16
 1 Nonlinear |R| = 5.551115e-16

Solve Converged!
```

The pattern is identical to Case 01: one Newton step, linear residual drops to machine
precision. The nonlinear residual at step 1 is different from Case 01 because the
system is larger (different matrix norm), but convergence is equally fast.

---

## Output Files

### `case02_diffusion_2d_out.e`

The Exodus II file contains the 2D mesh and the solution field `u(x,y)` at all 441 nodes.

**Visualizing in ParaView:**

1. File > Open, select `case02_diffusion_2d_out.e`, click OK, click Apply
2. The mesh appears as a colored square. The default coloring may not show u.
3. In the Properties panel (or pipeline), set the Color field to `u`
4. Click the "Rescale to Data Range" button (the magnifying glass icon)
5. You should see a color gradient: blue on the left (u=0) fading to red on the right (u=1)
6. The color bands are perfectly vertical — u depends only on x, not y

**Additional visualizations in ParaView:**

- **Contour lines**: Filters > Common > Contour, set field to `u`, add values 0.2, 0.4,
  0.6, 0.8. The contours are perfectly vertical straight lines.
- **Plot Over Line (x-direction)**: Filters > Data Analysis > Plot Over Line,
  from (0, 0.5, 0) to (1, 0.5, 0). The plot shows u vs x: a straight line from 0 to 1.
- **Plot Over Line (y-direction)**: from (0.5, 0, 0) to (0.5, 1, 0). The plot shows
  u vs y: a flat horizontal line at u=0.5. Confirms no y-dependence.
- **Warp By Scalar**: Filters > Common > Warp By Scalar, field `u`, scale factor 0.5.
  Visualizes the solution as a 3D surface height — a flat inclined plane.

**Reading the data programmatically:**

```python
import numpy as np
# Read nodal values from the Exodus file
# Using the NetCDF library or meshio:
# u_values.shape = (441,)  -- one value per node
# All values u_i should satisfy: u_i == x_coordinate_of_node_i
```

---

## Interpreting the Results

### What the Solution Looks Like

The solution `u(x,y) = x` is a planar surface that:
- Is 0 along the entire left edge (x=0)
- Is 1 along the entire right edge (x=1)
- Varies linearly from 0 to 1 across the plate in the x-direction
- Is perfectly constant in the y-direction

Color map appearance: perfectly vertical color bands (stripes running from bottom to top
of the plate). The color at any point depends only on how far right that point is.

### Verifying Correctness

Four verification checks:

1. **Visual**: vertical color stripes in ParaView, with left edge uniformly blue and right
   edge uniformly red.

2. **Plot over line (y-direction)**: at any x-location, u should be constant in y.
   A horizontal plot at y=const should show a flat line.

3. **Nodal value check**: for the node at position (x_i, y_j), the value should be
   exactly x_i. The y-coordinate y_j has no effect.

4. **Convergence**: solve completes in 1 Newton iteration, confirming linearity and
   correct setup.

### Comparing to Case 01

The two cases have the same exact solution `u = x`. The only difference is:
- Case 01: 1D problem on a line segment, 21 DOFs
- Case 02: 2D problem on a square, 441 DOFs

The extra DOFs in Case 02 correspond to DOFs at different y-values, all of which
satisfy the same linear ramp. Moving from 1D to 2D added no new physics — the
solution is the same. This confirms that the 2D `Diffusion` kernel is a correct
extension of the 1D version.

### Physical Insight

This result has a deep physical meaning: **in a uniform medium with no internal sources,
steady-state temperature profiles have no local extrema in the interior** (maximum
principle). The temperature can only be as high as the maximum boundary temperature
and as low as the minimum.

Furthermore, **symmetry determines structure**: because the problem has no y-dependence
in the BCs, forcing, or material properties, the solution cannot vary in y. Any
y-variation would break a symmetry that the equations do not break. This argument —
symmetry of the solution follows from symmetry of the problem — is a powerful tool for
understanding FEM results before running the solver.

---

## Key Concepts Learned

- **Dimension independence**: MOOSE kernels like `Diffusion` are dimension-agnostic;
  changing `dim` in the mesh block is all that is needed to move from 1D to 2D
- **2D GeneratedMesh**: creates structured quadrilateral meshes with four named
  boundaries: `left`, `right`, `bottom`, `top`
- **Natural boundary conditions**: boundaries absent from `[BCs]` automatically get
  zero-flux Neumann conditions — the physically meaningful default for insulated walls
- **Bilinear elements**: 2D first-order Lagrange elements are bilinear (linear in x and y,
  plus the x*y cross-term) and can represent affine surfaces exactly
- **Sparse matrix structure**: the 2D stiffness matrix connects each interior node to its
  8 neighbors (in a structured quad mesh), giving a banded sparse structure that AMG
  handles efficiently
- **Symmetry principle**: when the BCs, geometry, and material properties are symmetric
  in a direction, the solution is symmetric (or constant) in that direction
- **2D visualization**: ParaView color maps, contour plots, and plot-over-line are the
  primary tools for interpreting 2D FEM results
- **Scaling**: moving from 1D to 2D multiplies DOF count by approximately N (the number
  of elements per edge), but the solver's O(N) AMG complexity keeps run time manageable

---

## Experiments to Try

### Experiment 1: Mesh Refinement and Verification

Double and halve the mesh resolution:

```
nx = 10, ny = 10    # 100 elements, 121 nodes
nx = 20, ny = 20    # 400 elements, 441 nodes (baseline)
nx = 40, ny = 40    # 1600 elements, 1681 nodes
```

For this linear-solution problem, all three produce the exact same answer (u=x at each
node) to machine precision. Observe that:
- Solve time increases slightly with refinement
- Newton still converges in 1 iteration
- The solution field looks smoother in ParaView as the mesh is refined (more colored
  squares visible)

For a practical convergence study, you would need a problem with a non-polynomial exact
solution (see Case 04).

### Experiment 2: Non-Square Domain

Change the domain to a 2:1 rectangle:

```
[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 40      # twice as many x-elements to keep aspect ratio 1:1 per element
  ny   = 20
  xmin = 0
  xmax = 2       # domain is now [0,2] x [0,1]
  ymin = 0
  ymax = 1
[]
```

Keep the same BCs (u=0 on left, u=1 on right). The exact solution becomes:

```
u(x, y) = x / 2    (ranging from 0 at x=0 to 1 at x=2)
```

The temperature ramp is shallower because the same temperature difference (1 degree)
is spread over a longer rod (2 meters). Observe this in ParaView: the color gradient
transitions more slowly from left to right.

### Experiment 3: Add BCs on Top and Bottom

Convert the natural Neumann conditions on top and bottom to Dirichlet conditions:

```
[BCs]
  [left]
    type = DirichletBC; variable = u; boundary = left; value = 0
  []
  [right]
    type = DirichletBC; variable = u; boundary = right; value = 1
  []
  [bottom]
    type = DirichletBC; variable = u; boundary = bottom; value = 0.5
  []
  [top]
    type = DirichletBC; variable = u; boundary = top; value = 0.5
  []
[]
```

Now the corners are over-constrained (e.g., the bottom-left corner must simultaneously
be 0, 1, and 0.5 — impossible). MOOSE will resolve corner conflicts using priority rules.
The interior solution will now be a blend of all four boundary values and will no longer
be the simple ramp u=x. The solution becomes a 2D Laplace problem with four non-trivial
boundaries, which has no simple closed-form expression (it involves an infinite series
of sinusoids). Observe the curved contours in ParaView.

### Experiment 4: Anisotropic Mesh

Try a very coarse mesh in one direction:

```
nx = 20, ny = 2
```

This creates only 2 elements in the y-direction. Since the exact solution has no
y-dependence, the result should still be exactly `u = x`. This experiment confirms
that y-resolution is irrelevant for a problem with no y-variation — you need elements
where the solution changes, not where it is constant.

### Experiment 5: Add a Postprocessor to Quantify the 2D Solution

Add a `[Postprocessors]` block:

```
[Postprocessors]
  [avg_u]
    type     = ElementAverageValue
    variable = u
  []
  [max_u]
    type          = ElementExtremeValue
    variable      = u
    value_type    = max
  []
[]
```

Also add `csv = true` to `[Outputs]`. For the solution `u = x` on the unit square:
- `avg_u` should be exactly 0.5 (the spatial average of x over [0,1]x[0,1])
- `max_u` should be exactly 1.0 (the maximum of x on [0,1]x[0,1], at x=1)

These postprocessors appear in the CSV output file and in the console. They provide
scalar summary metrics that are useful for quick automated checking of results.
