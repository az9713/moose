# Case 06: Two-Region Domain with Different Conductivities

## Overview

Cases 01 through 05 used a single material throughout the domain. Most engineering
problems involve **multiple materials** — a wall made of insulation and concrete, a
nuclear fuel pellet surrounded by cladding, a composite panel with a metal face sheet
and foam core. Each material occupies a distinct region with its own properties.

This case introduces MOOSE's approach to multi-material domains:
1. **Mesh subdomains** (also called blocks): the mesh is tagged so different regions
   carry different identifiers.
2. **Block-restricted materials**: each material definition applies only to elements
   in specified blocks.
3. **MeshGenerator pipeline**: a sequence of mesh operations is chained together to
   build the final mesh from simpler pieces.

The problem solved is steady-state heat conduction through a two-layer wall:

```
-div( k * grad(u) ) = 0
k = 1.0  in the left half  (x in [0, 0.5])
k = 5.0  in the right half (x in [0.5, 1.0])
u = 0    at x=0  (left wall)
u = 1    at x=1  (right wall)
```

The right half is five times more conductive than the left. This creates a kink in the
temperature profile at the interface x=0.5, and the majority of the temperature drop
occurs in the less conductive left half.

This case demonstrates how MOOSE handles multi-material problems without requiring
separate meshes — a single conforming mesh is subdivided into labeled regions, and
the material system applies different properties to each region automatically.

---

## The Physics

### The Physical Problem

Consider a composite wall between two temperature-controlled surfaces. The wall has
two layers of equal thickness (each 0.5 m wide). The left layer is a poor conductor
(like insulation foam, k=1 W/m/K). The right layer is a good conductor (like aluminum,
k=5 W/m/K).

The left wall is held at u=0 (cold side), the right wall at u=1 (hot side). Heat flows
from hot to cold: right to left. At steady state, the same heat flux passes through
both layers (energy conservation). Because the materials have different conductivities,
they develop different temperature gradients to carry that common flux.

### The Governing Equation

```
-div( k * grad(u) ) = 0    on [0,1] x [0,1]
```

where `k` is now piecewise constant:

```
k(x) = 1.0    for x in [0.0, 0.5)   (left half, block 0)
k(x) = 5.0    for x in [0.5, 1.0]   (right half, block 2)
```

The equation is the same diffusion equation as previous cases. The only difference is
that `k` is discontinuous at x=0.5. This discontinuity is the material interface.

### Interface Conditions

At the interface between the two materials (x=0.5), two physical conditions hold:

1. **Temperature continuity**: the temperature is the same on both sides of the interface.
   (No temperature jump — no thermal contact resistance in this model.)
   ```
   u_left(0.5) = u_right(0.5)
   ```

2. **Flux continuity**: the heat flux through the interface is the same on both sides.
   (Energy conservation — what flows into the interface from the right must flow out
   to the left.)
   ```
   k_left * du/dx |_{x=0.5^-} = k_right * du/dx |_{x=0.5^+}
   ```
   ```
   1.0 * du/dx |_{left side}  = 5.0 * du/dx |_{right side}
   ```

These conditions are automatically satisfied by the finite element method on a
conforming mesh (no special treatment needed), provided the mesh has nodes exactly
at the interface. `SubdomainBoundingBoxGenerator` places the block boundary at x=0.5,
which aligns with element boundaries since nx=40 gives element edges at multiples of
0.025, including 0.5.

### Analytical Solution

For the 1D version, the two-layer problem has a simple exact solution. Let the heat flux
be `q` (constant, same in both layers). In each layer:

```
Left layer (0 <= x <= 0.5):
  q = -k_left * du/dx = -1.0 * du/dx
  du/dx = -q
  u(x) = u_0 - q*x = -q*x   (since u(0) = 0)

Right layer (0.5 <= x <= 1.0):
  q = -k_right * du/dx = -5.0 * du/dx
  du/dx = -q/5
  u(x) = u_0.5 - (q/5)*(x - 0.5)
```

From temperature continuity at x=0.5: `u_0.5 = -q*0.5 = -0.5*q`

From the right BC u(1) = 1:
```
1 = u_0.5 - (q/5)*0.5
1 = -0.5*q - 0.1*q
1 = -0.6*q
q = -5/3
```

Therefore:

```
Left half (0 <= x <= 0.5):
  u(x) = (5/3)*x

Right half (0.5 <= x <= 1.0):
  u(x) = 5/6 + (1/3)*(x - 0.5)
```

Key values:
- u(0.0) = 0.000  (left BC)
- u(0.5) = 5/6 * 0.5 * (1/0.3?) ... let us compute directly:
  u(0.5) = (5/3)*0.5 = 5/6 ≈ 0.833
- u(1.0) = 1.000  (right BC)

The interface temperature is 5/6 ≈ 0.833, much closer to the hot side than to the cold
side. This makes physical sense: the highly conductive right layer quickly equilibrates
to near-uniform temperature, so most of the temperature drop (0.833 out of 1.0 total)
occurs in the poorly conducting left layer.

Temperature gradient in each layer:
- Left: du/dx = 5/3 ≈ 1.667
- Right: du/dx = 1/3 ≈ 0.333

The left gradient is five times steeper. This is precisely the flux continuity condition:
`1.0 * 1.667 = 5.0 * 0.333 = 1.667` (same flux in both layers).

### ASCII Domain Diagram

```
  y=1  (insulated)
       +--------------------+--------------------+
       |                    |                    |
       |  Block 0           |  Block 2           |
       |  k = 1.0           |  k = 5.0           |
  u=0  |  (low conduct.)    | (high conduct.)    |  u=1
(cold) |                    |                    | (hot)
       |  steep gradient    |  shallow gradient  |
       |                    |                    |
       +--------------------+--------------------+
  y=0  (insulated)
       x=0               x=0.5                x=1

  Temperature profile along x:
  1.0 |                              *   *
  0.8 |                    ** <- u(0.5) = 0.833
  0.6 |               *
  0.4 |          *
  0.2 |     *
  0.0 +-------------------+--------------------
      x=0               x=0.5                x=1
      Block 0              |  Block 2
      slope = 5/3          |  slope = 1/3

  The kink at x=0.5 is where conductivity jumps from 1 to 5.
  Steeper slope in block 0 (harder to conduct) -> larger gradient.
  Shallower slope in block 2 (easier to conduct) -> smaller gradient.
```

### What Multi-Material Interfaces Mean Physically

A material interface is not a special thing mathematically — it is just a discontinuity
in the coefficient k. The finite element method handles this naturally when:
1. The mesh is **conforming** at the interface: element boundaries align with the
   material boundary (nodes lie on the interface, not cutting through elements)
2. No elements span the interface (each element is entirely inside one material)

When both conditions hold, the FEM solution automatically satisfies both interface
conditions (temperature and flux continuity) without any special treatment. This is
one of the strengths of the finite element method.

If the mesh is not conforming at the interface (element cut by material boundary), the
flux continuity condition is violated and the solution is inaccurate. For this reason,
MOOSE uses mesh generators that place element edges at material boundaries.

---

## Input File Walkthrough

The input file is `case06_two_materials.i`. This case has the most complex mesh setup
of the series.

### Block: `[Mesh]`

```
[Mesh]
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

  [right_half]
    type        = SubdomainBoundingBoxGenerator
    input       = gen
    bottom_left = '0.5 0.0 0'
    top_right   = '1.0 1.0 0'
    block_id    = 2
  []
[]
```

This case uses the **MeshGenerator pipeline** rather than the legacy `type = GeneratedMesh`
syntax. The mesh is built in stages, each stage transforming or extending the previous one.

#### Sub-block: `[gen]`

```
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
```

- `type = GeneratedMeshGenerator` — the MeshGenerator version of `GeneratedMesh`.
  The difference from previous cases: this is a MeshGenerator object that can be
  chained with other generators.
- `dim = 2` — 2D mesh.
- `nx = 40` — 40 elements in x. Choosing an even number ensures element edges fall
  exactly at x=0.5 (element edges are at 0, 0.025, 0.05, ..., 0.5, ..., 1.0).
  This is critical: the material boundary must coincide with element boundaries.
- `ny = 20` — 20 elements in y. Fewer than x because the y-direction is less important
  for this effectively 1D problem (chosen for visual aspect ratio).
- `xmin/xmax = 0/1`, `ymin/ymax = 0/1` — explicit domain bounds (same as default).

At this stage, all 800 elements belong to **block 0** (the default block ID).

#### Sub-block: `[right_half]`

```
[right_half]
  type        = SubdomainBoundingBoxGenerator
  input       = gen
  bottom_left = '0.5 0.0 0'
  top_right   = '1.0 1.0 0'
  block_id    = 2
[]
```

- `type = SubdomainBoundingBoxGenerator` — reassigns elements to a new subdomain
  based on whether their centroid falls inside a specified bounding box.
- `input = gen` — takes the mesh produced by the `[gen]` generator as input.
  This is how the pipeline is connected: each generator names its input.
- `bottom_left = '0.5 0.0 0'` — the lower-left corner of the bounding box.
  Format: 'x y z' (z is ignored for 2D meshes).
- `top_right = '1.0 1.0 0'` — the upper-right corner of the bounding box.
  Together these define the right half of the domain.
- `block_id = 2` — all elements whose centroid is inside the box [0.5,1.0] x [0.0,1.0]
  are reassigned from block 0 to block 2.

After this generator runs:
- **Block 0**: 400 elements in the left half (x in [0, 0.5])
- **Block 2**: 400 elements in the right half (x in [0.5, 1.0])

The **final generator** in the `[Mesh]` block is used as the active mesh. MOOSE
determines the pipeline automatically from the `input = ...` connections.

#### Why block_id = 2 and not 1?

Block IDs are user-chosen integers. The choice of 2 is intentional to avoid confusion:
block 1 is never used. In more complex problems, the gap leaves room to insert
additional intermediate blocks without renumbering. Any positive integer can be used.
The block ID is just a label — it has no physical meaning by itself.

#### The MeshGenerator Pipeline Concept

Think of MeshGenerators as a Unix pipe:

```
GeneratedMeshGenerator -----> SubdomainBoundingBoxGenerator -----> [Final Mesh]
(creates raw mesh)             (labels the right half)
       "gen"                         "right_half"
```

Additional generators could be chained:
- `RenameBlockGenerator` — rename block IDs or names
- `SideSetsFromBoundingBoxGenerator` — create new boundary sidesets
- `TransformGenerator` — rotate/scale/translate
- `MeshRepairGenerator` — fix mesh quality issues
- `RefineBlockGenerator` — selectively refine one block

The pipeline approach keeps each operation simple and composable.

### Block: `[Variables]`

```
[Variables]
  [u]
  []
[]
```

- Single unknown `u`, first-order Lagrange, default initial condition of zero.
- The variable spans the entire mesh (all blocks). MOOSE does not restrict variables
  to blocks by default. The same `u` is solved across both materials.

Total DOFs: (41*21) = 861 nodes, so 861 degrees of freedom.

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

- `type = MatDiffusion` — variable-coefficient diffusion, reads `k` from materials.
- `diffusivity = k` — material property name. MOOSE will look up `k` in the material
  system for each element. The material system knows which block each element belongs
  to, so it automatically provides k=1 for block-0 elements and k=5 for block-2 elements.

Note: only **one kernel** is needed, despite having two materials. The kernel is
general — it asks the material system for `k`, and the material system dispatches to
the appropriate material object based on the element's block ID. This is the elegance
of MOOSE's architecture: PDE operators (kernels) and material properties are completely
decoupled.

The kernel is not block-restricted. It applies to the entire domain (both blocks). If
you wanted different kernels in different blocks (e.g., a different PDE in each region),
you would add `block = 0` or `block = 2` parameters to each kernel.

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

- Left wall: u=0 (cold).
- Right wall: u=1 (hot).
- Same as Case 05. No BCs on top and bottom (zero-flux, insulated).
- The named boundaries `left`, `right`, `top`, `bottom` are assigned by
  `GeneratedMeshGenerator` in the same way as `GeneratedMesh`.

### Block: `[Materials]`

```
[Materials]
  [mat_left]
    type        = GenericConstantMaterial
    block       = 0
    prop_names  = 'k'
    prop_values = '1.0'
  []

  [mat_right]
    type        = GenericConstantMaterial
    block       = 2
    prop_names  = 'k'
    prop_values = '5.0'
  []
[]
```

This is the core of the multi-material approach. Two material objects are defined,
each restricted to a different block.

**`[mat_left]`:**
- `type = GenericConstantMaterial` — the simplest material type: sets a material
  property to a constant value everywhere in its applicable region.
- `block = 0` — this material applies **only** to elements in block 0 (the left half).
  MOOSE will never call this material for block-2 elements.
- `prop_names = 'k'` — declares the property named `k`.
- `prop_values = '1.0'` — k = 1.0 (constant) for all block-0 elements.

**`[mat_right]`:**
- `type = GenericConstantMaterial` — same type, different value.
- `block = 2` — applies only to block-2 elements (the right half).
- `prop_names = 'k'` — same property name `k`.
- `prop_values = '5.0'` — k = 5.0 for all block-2 elements.

Both materials define the same property name `k`. MOOSE resolves which one to use
based on the current element's block ID. This works because the materials are
block-restricted: there is no ambiguity about which k to use for any given element.

If you forgot the `block` parameter, MOOSE would detect that two materials claim to
provide the same property `k` on the same elements and throw an error. Block-restriction
is what makes them compatible.

**GenericConstantMaterial vs GenericFunctionMaterial:**
- `GenericConstantMaterial`: k = constant number. Simple and efficient (value computed
  once, not per quadrature point).
- `GenericFunctionMaterial`: k = evaluated function at each quadrature point (Case 05).
  For piecewise-constant k (as here), the constant material is both simpler and faster.

### Block: `[Postprocessors]`

```
[Postprocessors]
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
```

Two postprocessors, one for each subdomain:

- `[avg_u_left]` — average of `u` over block 0 (left half) only.
- `[avg_u_right]` — average of `u` over block 2 (right half) only.
- `block = 0` / `block = 2` — restricts the averaging to elements in the specified block.

Postprocessors, like kernels and materials, can be block-restricted.

Expected values from the analytical solution:
- Average over left half: average of `(5/3)*x` for x in [0, 0.5]
  = `(5/3) * integral_0^{0.5} x dx / 0.5 = (5/3) * (0.25/2) / 0.5 = (5/3) * 0.25 = 5/12 ≈ 0.4167`

  Wait, let us be more careful. The average of u over [0, 0.5] is:
  ```
  (1/0.5) * integral_0^{0.5} (5/3)*x dx = 2 * (5/3) * [x^2/2]_0^{0.5}
  = 2 * (5/3) * (0.125) = 2 * 0.2083 = 0.4167
  ```

- Average over right half: average of `u(0.5) + (1/3)*(x-0.5)` for x in [0.5, 1.0]
  = `(1/0.5) * integral_{0.5}^{1.0} [5/6 + (1/3)*(x-0.5)] dx`

  Let s = x - 0.5, s in [0, 0.5]:
  = `2 * integral_0^{0.5} [5/6 + (1/3)*s] ds`
  = `2 * [5/6 * 0.5 + (1/3) * 0.125]`
  = `2 * [0.4167 + 0.0417]`
  = `2 * 0.4583 = 0.9167`

So expected: avg_u_left ≈ 0.417, avg_u_right ≈ 0.917.

These large values confirm the asymmetry: the right (conductive) half is at nearly the
hot-side temperature, while the left (resistive) half spans most of the temperature range.

### Block: `[Executioner]`

```
[Executioner]
  type       = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]
```

- Identical to previous cases. Steady-state, PJFNK with AMG preconditioning.
- The problem is linear, so Newton converges in one iteration.

### Block: `[Outputs]`

```
[Outputs]
  exodus = true
  csv    = true
[]
```

- Exodus output: `case06_two_materials_out.e`
- CSV output: `case06_two_materials_out.csv` (contains both postprocessor values)

---

## What Happens When You Run This

### Invocation

```bash
cd quickstart-runs/case06
mpirun -n 1 moose-app-opt -i case06_two_materials.i
```

### Mesh Generation Pipeline

MOOSE processes the mesh in order:
1. Creates the initial 40x20 mesh with all elements in block 0.
2. Runs `SubdomainBoundingBoxGenerator`: for each element, computes its centroid,
   checks if it falls inside [0.5,1.0] x [0.0,1.0], and if so, reassigns it to block 2.
   - Element centroids for the 20 rightmost columns are at x = 0.5125, 0.5375, ..., 0.9875
   - All fall inside the bounding box, so all 400 right-half elements become block 2.
3. The final mesh is ready: 400 elements in block 0, 400 in block 2.

### Material System Setup

At initialization, MOOSE validates that:
- Every element in the mesh has a material that provides property `k`.
- Block 0 is covered by `mat_left`.
- Block 2 is covered by `mat_right`.
- No element is uncovered (would cause a "no material defined" error).
- No element is covered by two materials claiming the same property (would be ambiguous).

### Assembly with Two Materials

During assembly, for each element:
- If element is in block 0: `mat_left.computeQpProperties()` runs, setting `k = 1.0`
  at all quadrature points.
- If element is in block 2: `mat_right.computeQpProperties()` runs, setting `k = 5.0`
  at all quadrature points.
- The `MatDiffusion` kernel then uses whichever value was set.

The kink in the temperature profile at x=0.5 arises naturally because the assembled
stiffness matrix encodes the discontinuous conductivity. No special interface handling
is required.

### Console Output

```
Mesh Information:
  Spatial dimension:    2
  Mesh:                 800 elements, 861 nodes
  Subdomains:           [0, 2]

Nonlinear System:
  Num DOFs:             861

 0 Nonlinear |R| = 1.000000e+00
      0 Linear |R| = 7.32e-15
 1 Nonlinear |R| = 7.32e-15

Converged!

Postprocessor Values:
  avg_u_left  = 4.1667e-01
  avg_u_right = 9.1667e-01
```

- MOOSE lists the two subdomains: [0, 2].
- One Newton iteration (linear problem).
- avg_u_left ≈ 0.417 and avg_u_right ≈ 0.917 match the analytical values exactly.

---

## Output Files

### `case06_two_materials_out.e`

The Exodus file contains:
- The mesh with block information (blocks 0 and 2 are stored separately)
- The solution field `u` at all nodes

In ParaView:
1. Open the file, click Apply.
2. Color by `u`. The kink at x=0.5 should be visible: a sharp change in the color
   gradient where the material changes.
3. Use "Plot Over Line" from (0,0.5,0) to (1,0.5,0) to extract the 1D profile.
   Verify: straight line from 0 to 0.833 in the left half, straight line from 0.833
   to 1 in the right half. The two slopes differ by a factor of 5.
4. Use Threshold or Extract Block to isolate each subdomain and view them separately.

### `case06_two_materials_out.csv`

```
time,avg_u_left,avg_u_right
0,4.1667e-01,9.1667e-01
```

Both postprocessor values are in a single CSV file, one column per postprocessor.
Compare these to the analytical predictions to verify correctness.

---

## Understanding the Plots

Running `python visualize_all.py` from the `quickstart-runs/` directory produces the
following two plots saved into this directory.

### `case06_contour_2d.png`

**What the plot shows.** A 2D filled-contour of the solution `u(x,y)` on the unit
square with a material interface at x=0.5. The viridis colormap is used. A white
dashed vertical line marks the interface at x=0.5.

**Physical quantities.** The color encodes temperature in a two-material wall: left
half has k=1, right half has k=5. The white dashed line marks the material interface.

**How to judge correctness.** The contour bands should be vertical (no y-variation)
but with a clear difference in band spacing on either side of the interface. The left
half (k=1, poor conductor) needs a steeper gradient to carry the same flux and
therefore has closely spaced bands. The right half (k=5, good conductor) has a
shallower gradient and widely spaced bands. The transition between band spacings
should be sharp at x=0.5.

**What would indicate a problem.**
- Uniformly spaced bands across the entire domain: the block-restricted materials are
  not working — both halves are using the same conductivity.
- The kink in band spacing appears at a location other than x=0.5: the subdomain
  boundaries in the mesh are incorrectly set.
- Curved bands: incorrect top/bottom boundary conditions.

### `case06_line_interface.png`

**What the plot shows.** A line plot of u along the midline y ≈ 0.5, with x on the
horizontal axis (0 to 1) and u on the vertical axis (0 to 1). Blue circles are MOOSE
data; a vertical dashed red line marks the material interface at x=0.5.

**Physical quantities.** The curve shows the temperature profile through the two-
material wall. Two distinct linear segments are visible, joined at the interface.

**How to judge correctness.** The curve should consist of two straight-line segments
that meet at a kink at x=0.5. The left segment (k=1) should be steeper, and the
right segment (k=5) should be shallower. The exact interface temperature is
approximately u(0.5) ≈ 0.833 (the majority of the temperature drop happens in the
less conductive left material). The endpoints must be exactly u=0 at x=0 and u=1
at x=1.

A useful rule: the temperature drop in each layer is inversely proportional to its
conductivity. With k1=1 and k2=5, the left layer carries 5/6 of the total temperature
drop and the right carries 1/6. So u at the interface is 1 - 5/6 = 1/6 from the
hot side, which gives u(0.5) = 5/6 ≈ 0.833.

**What would indicate a problem.**
- A single straight line from 0 to 1 (no kink): both materials have the same k, the
  block-restricted materials are not working.
- The kink at the wrong location: the interface is not at x=0.5 in the mesh.
- The interface temperature above 0.9 or below 0.7: the conductivity ratio is wrong.
- Non-linear segments (curved instead of straight on each side): the solver did not
  converge, or a nonlinear material property was accidentally used.

---

## Interpreting the Results

### The Solution Field

The temperature profile is **piecewise linear with a kink at x=0.5**:

```
Left half (block 0, k=1):    u(x) = (5/3) * x
Right half (block 2, k=5):   u(x) = 5/6 + (1/3) * (x - 0.5)
```

Key values:
- u(0.0)  = 0.000  (left BC, exactly enforced)
- u(0.5-) = (5/3)*0.5 = 0.833  (approached from the left)
- u(0.5+) = 5/6 = 0.833  (same from the right — temperature is continuous)
- u(1.0)  = 1.000  (right BC, exactly enforced)

The kink is at x=0.5 where the slope changes from 5/3 to 1/3. This is physically the
material interface. The temperature is continuous (no jump) but the gradient is not.

### Verifying Correctness

1. **Interface temperature**: approximately 0.833. From the exact solution, exactly 5/6.
   Verify in ParaView with a "Plot Over Line" filter at y=0.5.
2. **Average values**: avg_u_left ≈ 0.417, avg_u_right ≈ 0.917. These match the
   analytical computation (both should match to near-machine precision).
3. **Slope ratio**: Extract the gradient from the left and right halves. The ratio of
   gradients should be exactly 5 (= k_right/k_left). MOOSE captures this correctly
   because the mesh is conforming at the interface.
4. **Flux conservation**: at any x-plane, `k * du/dx` should be the same value.
   Left: `1.0 * (5/3) = 5/3`. Right: `5.0 * (1/3) = 5/3`. Both match.

### Physical Insight

The two-material result reveals a fundamental principle: **thermal resistance** controls
temperature distribution. The thermal resistance of each layer is `R = L / (k * A)`
where L is thickness and A is area. For unit area:

```
R_left  = 0.5 / 1.0 = 0.5
R_right = 0.5 / 5.0 = 0.1
R_total = 0.6
```

The temperature drop across each layer is proportional to its resistance:

```
Delta_T_left  = 1.0 * (R_left / R_total)  = 1.0 * (0.5/0.6) = 0.833
Delta_T_right = 1.0 * (R_right / R_total) = 1.0 * (0.1/0.6) = 0.167
```

Interface temperature = u(0) + Delta_T_left = 0 + 0.833 = 0.833. This matches.

This is exactly the electrical resistor analogy for heat conduction: resistors in series
divide the voltage (temperature difference) in proportion to their resistance.

### Comparing to Case 05

| Property | Case 05 (k = 1+x) | Case 06 (k = 1 then k = 5) |
|----------|-------------------|----------------------------|
| Profile  | Smooth logarithm  | Piecewise linear with kink |
| Midpoint | 0.585             | 0.833                      |
| avg_u    | 0.557             | (0.417+0.917)/2 = 0.667    |
| Interface| No interface      | Sharp kink at x=0.5        |

Case 05's smoothly varying k gives a smooth curved profile. Case 06's step-change in k
gives a profile with a sharp kink. Both are above the 0.5 of the uniform-k case because
both concentrate resistance on the left side.

---

## Key Concepts Learned

- **Mesh subdomains (blocks)**: regions of the mesh tagged with integer IDs that allow
  different properties and physics to be applied to different geometric regions
- **MeshGenerator pipeline**: chain mesh generators together, each taking the previous
  as `input`, to build complex meshes from simple operations in a readable sequence
- **GeneratedMeshGenerator**: the MeshGenerator version of GeneratedMesh, usable in
  a pipeline
- **SubdomainBoundingBoxGenerator**: reassign elements within a bounding box to a new
  subdomain ID, the simplest way to create a two-region mesh for a rectangular domain
- **Block-restricted materials**: use the `block` parameter on material objects to
  apply different properties in different regions of the mesh, enabling multi-material simulations
- **GenericConstantMaterial**: assign a constant value to a named material property
  over a specified block
- **Block-restricted postprocessors**: compute averages and other quantities
  independently for each subdomain using the `block` parameter
- **Interface conditions**: temperature continuity and flux continuity are automatically
  enforced by the FEM on a conforming mesh; no special interface treatment is needed
- **Thermal resistance analogy**: series resistors divide temperature differences in
  proportion to resistance, explaining the highly asymmetric temperature distribution
- **Conforming mesh requirement**: the material interface must coincide with element
  boundaries; choosing nx=40 ensures element edges fall exactly at x=0.5

---

## Experiments to Try

### Experiment 1: Vary the Conductivity Ratio

Change `prop_values = '5.0'` for `mat_right` to different values: try 2, 10, 100.
Predict the interface temperature before running:

```
u_interface = R_left / (R_left + R_right) = (0.5/1) / (0.5/1 + 0.5/k_right)
            = 1 / (1 + 1/k_right)
```

For k_right = 2: u_interface = 1/(1+0.5) = 0.667
For k_right = 10: u_interface = 1/(1+0.1) = 0.909
For k_right = 100: u_interface = 1/(1+0.01) = 0.990

As k_right increases, the right half becomes nearly isothermal (at u=1), and essentially
all the temperature drop occurs in the left half.

### Experiment 2: Three-Layer Domain

Add a third material region (e.g., a thin middle strip from x=0.45 to x=0.55 with a
very low or high conductivity). This requires:
1. A second `SubdomainBoundingBoxGenerator` to carve out the middle strip (block 3)
2. Careful ordering: first carve the right half as block 2, then carve the middle as block 3
   (or vice versa; order matters for overlapping generators)
3. A third material definition `block = 3`

Alternatively, change `block_id = 2` to use block 1 and add a thin strip as block 2,
with the right bulk as block 3. Experiment to understand how the pipeline ordering works.

### Experiment 3: Asymmetric Domain

Move the interface by changing the bounding box from `0.5` to a different x-coordinate,
say `0.3` (thinner left layer, thicker right layer). Predict the new interface temperature
using the resistance formula. Verify that the simulation matches.

For interface at x=x_0:
```
R_left = x_0 / k_left = x_0 / 1.0
R_right = (1 - x_0) / k_right = (1 - x_0) / 5.0
u_interface = R_left / (R_left + R_right)
```

For x_0 = 0.3: R_left = 0.3, R_right = 0.14, u_interface = 0.3/0.44 ≈ 0.682

### Experiment 4: Add Block-Restricted Kernels

Suppose you want a source term only in the left half (simulating heat generation in
one layer, like nuclear fuel). Add:

```
[Kernels]
  [source_left]
    type    = BodyForce
    variable = u
    value   = 2.0
    block   = 0
  []
[]
```

The `block = 0` parameter restricts this kernel to block-0 elements. The right half
has no source. This breaks the linear profile in the left half (it becomes parabolic)
and changes the interface temperature. Predict the new profile analytically if you can.

### Experiment 5: Non-rectangular Interface

For a non-rectangular interface (e.g., a circle or diagonal), `SubdomainBoundingBoxGenerator`
is insufficient. Use `SubdomainBoundingBoxGenerator` with a rotated bounding box or
explore `ParsedSubdomainMeshGenerator`, which uses a mathematical expression to classify
elements. For example, to put all elements with centroid satisfying `x + y < 1` into
block 2, use:

```
[inside]
  type              = ParsedSubdomainMeshGenerator
  input             = gen
  combinatorial_geometry = 'x + y < 1'
  block_id          = 2
[]
```

This creates a triangular subdomain and demonstrates MOOSE's flexibility for complex geometries.
