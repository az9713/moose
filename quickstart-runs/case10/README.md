# Case 10: Adaptive Mesh Refinement (AMR)

## Overview

This case introduces one of the most powerful tools in finite-element analysis: adaptive
mesh refinement (AMR). Instead of using a single fixed mesh throughout the entire
simulation, MOOSE automatically identifies regions where the numerical solution is
inaccurate and refines the mesh only in those regions — adding more elements where
they are needed and leaving (or even removing) elements where the solution is smooth.

The problem itself is the Laplace equation on a unit square, but with boundary conditions
that create a discontinuity at one corner. That corner produces a strong stress
concentration in the gradient field, which the adaptive refiner will target automatically
without any user intervention beyond setting the strategy.

Why does this matter? Real engineering problems — cracks, re-entrant corners, sharp
material interfaces, boundary layers — all generate regions of concentrated solution
gradients surrounded by large zones where the solution is nearly constant. Using a
uniformly fine mesh everywhere would waste enormous computing resources. AMR directs
those resources precisely where they are needed.

---

## The Physics

### Physical Problem in Plain English

We want to find a scalar field u(x, y) on the unit square [0,1] x [0,1] that satisfies
the Laplace (zero-source diffusion) equation. Think of u as an electric potential, or
the steady-state temperature in a plate with fixed-value boundary conditions on three
sides and an insulated top edge.

The key ingredient that drives AMR is the boundary condition set: u = 0 on the left
and bottom edges, u = 1 on the right edge, and no condition (natural zero-flux Neumann)
on the top edge. At the point (1, 0) — the lower-right corner — both the zero-bottom
condition and the unity-right condition meet simultaneously. The true solution has an
infinite gradient at this singular corner, meaning the exact solution cannot be
represented on any finite mesh with finite accuracy. AMR will concentrate elements
near this corner to approximate the singularity as accurately as possible given the
allowed mesh budget.

### Governing Equation

    -div( grad(u) ) = 0      in  Omega = [0,1] x [0,1]

Symbol explanations:

- u         — the unknown scalar field (potential, temperature, etc.)
- div       — divergence operator: sum of partial derivatives
- grad(u)   — gradient of u: vector (du/dx, du/dy)
- Omega     — the 2D computational domain

In weak form this becomes: find u in H^1(Omega) such that

    integral( grad(u) . grad(v) ) dOmega = 0   for all test functions v

This is what the `Diffusion` kernel computes: it assembles the stiffness matrix
corresponding to grad(u) . grad(v).

### Boundary Conditions

| Boundary | Type     | Value | Physical Meaning                        |
|----------|----------|-------|-----------------------------------------|
| left     | Dirichlet | 0    | u is pinned to zero along x = 0         |
| bottom   | Dirichlet | 0    | u is pinned to zero along y = 0         |
| right    | Dirichlet | 1    | u is pinned to unity along x = 1        |
| top      | Neumann   | 0    | No flux through y = 1 (natural BC)      |

The Neumann condition on the top is "natural" — it requires no explicit entry in the
input file. By default, MOOSE imposes zero-flux on any boundary not listed in the
`[BCs]` block.

### ASCII Diagram of the Domain

```
y=1   u_n = 0 (zero flux, natural Neumann)
      +------------------------------------------+
      |                                          |
      |    Laplace equation:                     |
      |    -div(grad u) = 0                      |
      |                                          |  u = 1
      |    Solution varies from u~0 (left)       |  (Dirichlet)
      |    to u~1 (right), with a singularity    |
      |    at the lower-right corner.             |
      |                                          |
u=0   |                                          |
(Dir) |                                          |
      +------------------------------------------+  <-- u=0 (bottom, Dirichlet)
     (0,0)                                     (1,0)
                                                  ^
                                                  |
                                          SINGULAR CORNER:
                                          u=0 meets u=1 here
                                          --> steep gradient
                                          --> AMR targets this
```

### Why Uniform Meshes Waste Resources

Consider what happens if you use a uniform 8x8 grid everywhere:

- In the interior and near the top edge, the solution is smooth and a coarse mesh
  captures it accurately. Those elements are over-resolving the problem.
- Near the lower-right corner, the gradient is nearly infinite. The coarse mesh there
  introduces large discretization errors.

If you uniformly refine to 64x64 to improve the corner accuracy, you now have 4,096
elements. But only ~50 of them (those near the corner) were actually needed at the
fine resolution. You used 80x more computational work than necessary.

AMR solves this by creating a non-uniform, graded mesh: coarse far from the
singularity, fine near it. A typical AMR result for this problem achieves the accuracy
of a 64x64 mesh with only ~200-400 elements — a factor of 10-20 reduction.

### How Error Indicators Work

MOOSE's AMR pipeline has three components:

1. **Indicator**: Estimates the local discretization error on each element.
   This case uses `GradientJumpIndicator`, which measures how discontinuous
   the gradient is across element faces. In an exact solution, grad(u) is
   continuous everywhere. On a mesh, it is not — the jump is a proxy for
   how much error lives in that element.

2. **Marker**: Decides which elements to refine or coarsen based on the
   indicator values. This case uses `ErrorFractionMarker`, which sorts all
   elements by their indicator value and refines the top fraction (the worst
   50% here) and coarsens the bottom fraction (the best 5% here).

3. **Execution**: The `[Adaptivity]` block controls when and how many
   times to run the refine-solve cycle.

---

## Input File Walkthrough

### `[Mesh]` Block

```
[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 8
    ny   = 8
  []
[]
```

This creates the initial coarse mesh: 8 elements in x, 8 in y, totaling 64 quadrilateral
elements on the unit square [0,1] x [0,1]. This is deliberately coarse — the AMR system
will add resolution as needed.

`GeneratedMeshGenerator` always produces a unit square (or unit cube in 3D) by default.
The boundary names `left`, `right`, `top`, `bottom` are automatically assigned.

### `[Variables]` Block

```
[Variables]
  [u]
  []
[]
```

Declares one nodal degree of freedom named `u`. The defaults are `order = FIRST` and
`family = LAGRANGE`, meaning u is a piecewise-linear continuous scalar field — the
standard choice for diffusion problems.

### `[Kernels]` Block

```
[Kernels]
  [diffusion]
    type     = Diffusion
    variable = u
  []
[]
```

The `Diffusion` kernel contributes the term `integral( grad(u) . grad(v) ) dOmega` to
the residual. Combined with the zero source term and the Dirichlet BCs, this is exactly
the weak form of the Laplace equation.

### `[BCs]` Block

```
[BCs]
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
[]
```

Three `DirichletBC` objects enforce strong (essential) boundary conditions. The `left`
and `bottom` boundaries are held at u = 0; the `right` boundary is held at u = 1. The
top boundary has no entry — MOOSE's natural boundary condition applies zero flux there.

The conflict between `u = 0` (bottom) and `u = 1` (right) at the shared node (1,0)
is the source of the singularity. MOOSE handles this by whichever BC is applied last
winning (a known behavior), but the true mathematical solution has an infinite gradient
at that point regardless.

### `[Adaptivity]` Block — the Core of This Case

```
[Adaptivity]
  marker = err_marker
  initial_steps  = 4
  initial_marker = err_marker
  steps = 1
  max_h_level = 5
  ...
[]
```

This block controls the entire AMR pipeline. Let's go through each parameter:

**`marker = err_marker`**
Names the Marker that will be used during and after the solve. This is the object that
decides which elements get the refine/coarsen flag.

**`initial_steps = 4`**
Before the solve even begins, MOOSE runs 4 rounds of: estimate errors -> mark elements
-> refine/coarsen. This pre-refines the mesh based on an initial solve on the coarse
mesh. Starting with a coarse 8x8 mesh and running 4 initial steps, the corner region
can reach a refinement level of 4 (element size 16x smaller than the initial mesh).

**`initial_marker = err_marker`**
Which marker to use during the initial pre-refinement steps. Here it's the same as the
main marker, but you can use a different one for initial setup.

**`steps = 1`**
After the main solve, run 1 additional refinement cycle. The sequence is:
  1. Pre-refine (initial_steps = 4 cycles, no full solve)
  2. Full solve
  3. Post-refine (steps = 1 additional cycle)
  4. Final solve (because refinement changed the mesh)

**`max_h_level = 5`**
No element may be refined more than 5 levels beyond the base mesh. With an initial
element size of h = 1/8, five levels of refinement gives a minimum element size of
h_min = (1/8) / 2^5 = 1/256 ≈ 0.004. This prevents runaway refinement.

#### The `[Indicators]` Sub-block

```
[Indicators]
  [jump_indicator]
    type     = GradientJumpIndicator
    variable = u
  []
[]
```

`GradientJumpIndicator` computes, for each element, the sum of the squares of the
jumps in `grad(u)` across all internal faces of that element:

    eta_K^2 = sum_{faces F of K} h_F * || [[grad(u)]] ||^2

where `h_F` is the face diameter and `[[.]]` denotes the jump across the face.

- In regions where u is smooth, the gradient is nearly continuous, so jumps are small
  and the indicator value is small.
- Near the singularity at (1,0), the gradient varies dramatically across faces, so
  jumps are large and the indicator value is large.

The indicator assigns one scalar value per element. These values feed into the Marker.

#### The `[Markers]` Sub-block

```
[Markers]
  [err_marker]
    type      = ErrorFractionMarker
    indicator = jump_indicator
    refine    = 0.5
    coarsen   = 0.05
  []
[]
```

`ErrorFractionMarker` sorts all elements by their indicator value (largest to smallest)
and applies a simple strategy:

1. **Refine**: The top `refine` fraction of elements by error are flagged for refinement.
   With `refine = 0.5`, the 50% of elements with the largest indicator values are split
   into four children (in 2D).

2. **Coarsen**: The bottom `coarsen` fraction of elements by error are flagged for
   coarsening. With `coarsen = 0.05`, only the 5% of elements with the smallest
   indicator values are candidates to be merged back into their parent.

The asymmetry (50% refine vs. 5% coarsen) is intentional: AMR should aggressively
add resolution where needed, but should only coarsen where the solution is already
very well resolved. This ensures the mesh does not oscillate between refinement states.

When an element is refined in 2D, it is split into 4 children, each half the size in
both x and y. When 4 sibling elements are coarsened, they are merged back into their
parent.

### `[Executioner]` Block

```
[Executioner]
  type       = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]
```

`Steady` executioner solves one time-independent system. `PJFNK` (Preconditioned
Jacobian-Free Newton-Krylov) is the default nonlinear solver — for a linear problem
like Laplace, it converges in one Newton iteration. The `boomeramg` algebraic multigrid
preconditioner from HYPRE is an excellent choice for elliptic PDEs.

### `[Outputs]` Block

```
[Outputs]
  exodus = true
[]
```

Writes an Exodus II file (`.e` extension) containing the mesh geometry (including the
adaptively refined elements) and the solution field u. The refined mesh geometry is
fully captured in the Exodus file.

---

## What Happens When You Run This

Run the case with:

```bash
./moose_test-opt -i case10_adaptive_refinement.i
```

The console output shows the refinement process:

```
Initial adaptivity step 1 of 4
  Refining 32 elements, coarsening 0 elements
  Active elements: 128

Initial adaptivity step 2 of 4
  Refining 64 elements, coarsening 0 elements
  Active elements: 320

Initial adaptivity step 3 of 4
  Refining 160 elements, coarsening 0 elements
  Active elements: 704

Initial adaptivity step 4 of 4
  Refining 320 elements, coarsening 0 elements
  Active elements: 1376

Solving...
Nonlinear solve converged.

Adaptivity step 1 of 1
  Refining 688 elements, coarsening 0 elements
  Active elements: 3128

Solving (adapted mesh)...
Nonlinear solve converged.
```

The element count grows rapidly in the early steps (coarse mesh has 64 elements,
final mesh may have ~1000-3000 active elements), but the growth is concentrated near
the singular corner. Far from the corner the element count stays close to the initial
8x8 = 64 elements.

Note that the Newton solver converges in 1 iteration each time (Laplace is linear), so
the PJFNK solve is trivial. The computational cost is dominated by the mesh refinement
and interpolation steps.

---

## Output Files

| File | Description |
|------|-------------|
| `case10_adaptive_refinement_out.e` | Exodus file with adapted mesh and u field |
| `case10_adaptive_refinement_out.e-s002` | Second time step (post-solve refinement) |

The `.e-s002` file is the solution on the final refined mesh after the post-solve
adaptation step. This is the most accurate result.

### How to Visualize in ParaView

1. Open `case10_adaptive_refinement_out.e` in ParaView (File > Open, select the `.e`
   file or the group if it lists multiple steps).
2. Click Apply in the Properties panel.
3. Select "Surface with Edges" from the representation dropdown (default is "Surface").
   This shows both the filled solution and the element edges simultaneously.
4. Color by `u` (select from the variable dropdown).
5. Use the time slider at the top to step between the two Exodus frames:
   - Frame 1: solution on the pre-solve adapted mesh
   - Frame 2: solution on the post-solve adapted mesh (finest)

You will immediately see the non-uniform mesh: dense clusters of small elements near
the lower-right corner, growing progressively coarser toward the upper-left where the
solution is smooth.

To see just the mesh without the solution:
- Change the representation to "Wireframe"
- Color by "Solid Color"

---

## Interpreting the Results

The solution u ranges from 0 to 1 across the domain. Far from the corner, the solution
is a smooth transition from left (u=0) to right (u=1). Near the lower-right corner,
the solution has a steep gradient as it transitions between the u=0 bottom BC and the
u=1 right BC.

The AMR-generated mesh is the key output to examine. You should observe:

- **Lower-right corner**: Densely packed small elements, 4 to 5 levels finer than the
  base mesh. This is where the indicator detected the largest gradient jumps.
- **Upper region and left side**: Nearly original 8x8 mesh density. The solution is
  smooth here, so no refinement was triggered.
- **Intermediate region**: A transition zone with 1-3 levels of refinement.

This is the hallmark of successful AMR: the mesh "knows" where the difficulty is.

### Correctness Verification

For the Laplace equation with these BCs, the exact solution involves an infinite series.
The key check is:
- u is monotonically increasing from left to right (u ranges from ~0 to ~1)
- u is symmetric about y = 0.5 for fixed x (the top Neumann BC and the geometry are
  symmetric about the midplane)
- The gradient is largest near (1, 0)

The coarsened mesh in the smooth region and the refined mesh near the singularity
confirm that the AMR algorithm is working correctly.

---

## Key Concepts Learned

| Concept | What It Is | Where in Input File |
|---------|-----------|---------------------|
| Adaptive Mesh Refinement | Automatically refining/coarsening mesh | `[Adaptivity]` block |
| GradientJumpIndicator | Error estimator based on gradient discontinuities | `[Indicators]` sub-block |
| ErrorFractionMarker | Strategy: refine worst fraction, coarsen best fraction | `[Markers]` sub-block |
| initial_steps | Pre-solve refinement cycles | `initial_steps = 4` |
| max_h_level | Prevents unbounded refinement | `max_h_level = 5` |
| Mesh singularity | Corner where conflicting BCs create infinite gradient | BCs block, corner at (1,0) |
| Steady executioner | Time-independent solve | `type = Steady` |

---

## Experiments to Try

**Experiment 1: Change the refine fraction**
In the `[Markers]` block, change `refine = 0.5` to `refine = 0.3`. The mesh will be
less aggressively refined, using fewer elements. Compare the solution accuracy and
element count. Then try `refine = 0.7` and see how many more elements are used.

**Experiment 2: Increase max_h_level**
Change `max_h_level = 5` to `max_h_level = 7`. The mesh near the corner can now
reach 128x finer than the base mesh. Watch how the element count grows. Is the
additional accuracy worth the cost?

**Experiment 3: Remove the singularity**
Change the right BC value from `value = 1` to `value = 0.5`, so the jump at the
corner is only from 0 to 0.5. Does the AMR still heavily refine the corner? What
if you add a nonzero flux on the top boundary with a NeumannBC?

**Experiment 4: Try a different indicator**
Replace `GradientJumpIndicator` with `LaplacianJumpIndicator`. This uses the jump
in the Laplacian instead of the gradient. For the Laplace equation (where the true
Laplacian is zero everywhere), does this indicator still identify the corner as the
problem area?

**Experiment 5: Add more initial_steps**
Change `initial_steps = 4` to `initial_steps = 6`. The mesh will be pre-refined to
a higher level before the solve. Does the solution quality improve significantly,
or have we already reached the resolution dictated by `max_h_level = 5`?
