# Case 99: Conducting Cylinder in Uniform Electric Field — Polar BVP with Penalty Method

## Overview

A perfectly conducting circular cylinder of radius R = 0.5 m is placed in a uniform background electric field E₀ = 1 V/m directed along the y-axis. Outside the cylinder the potential satisfies Laplace's equation, and the conducting surface enforces an equipotential boundary condition. This is the canonical electrostatic scattering problem from MIT 6.641 Lecture 11, covering boundary value problems in polar coordinates.

The analytic solution in polar coordinates is:

```
Phi(r, theta) = -E0 (r - R²/r) sin(theta)
```

Expressed in Cartesian coordinates (using y = r sin θ):

```
Phi(x, y) = -E0 * y * (1 - R² / (x² + y²))
```

The interior of the conductor (r < R) sits at Phi = 0, and the induced surface charge density on the cylinder is σ_s(θ) = 2 ε₀ E₀ sin θ. Because MOOSE's `GeneratedMesh` cannot punch a hole in a structured mesh, the simulation models the full square domain [-3, 3]² and uses a smooth penalty term to drive the potential to zero inside the cylinder, avoiding the need for a conforming mesh with a circular void.

This case demonstrates the **penalty method** for encoding conductor constraints without mesh surgery, **ParsedFunction** for analytic boundary conditions and penalty fields, and quantitative verification via `ElementL2Error` against the known exact solution.

---

## The Physics

**Governing equation** (outside the cylinder):

```
nabla^2 Phi = 0
```

Inside the cylinder the conductor is enforced by adding a large penalty reaction:

```
nabla^2 Phi + p(x,y) * Phi = 0
```

where p(x,y) is a smooth step function that is 1000 inside r < R and 0 outside.

**Boundary conditions:**

| Boundary | Condition | Expression |
|----------|-----------|-----------|
| All outer walls (x or y = ±3) | Analytic Dirichlet | Phi = -E₀ y (1 - R²/(x²+y²)) |
| Cylinder surface (r = R, implicit) | Conductor, Phi = 0 | Enforced by penalty |

**Material properties:**

| Parameter | Value | Units |
|-----------|-------|-------|
| E₀ (background field) | 1.0 | V/m |
| R (cylinder radius) | 0.5 | m |
| L (half-domain size) | 3.0 | m |
| Penalty strength | 1000 | 1/m² |

**Domain geometry:**

- Square domain: [-3, 3] × [-3, 3], 40 × 40 QUAD4 elements
- Cylinder at origin, fully contained in domain; L/R = 6, so the outer boundary is in the far field where the analytic solution is accurate

---

## Input File Walkthrough

### Variables

```
[Variables]
  [phi]
    order  = FIRST
    family = LAGRANGE
  []
[]
```

A single scalar variable `phi` representing the electric potential Phi [V]. FIRST-order Lagrange elements are appropriate for Laplace's equation.

### Kernels

```
[Kernels]
  [laplacian]
    type     = Diffusion
    variable = phi
  []
  [conductor_penalty]
    type          = MatReaction
    variable      = phi
    reaction_rate = conductor_penalty
  []
[]
```

`Diffusion` contributes the weak form of nabla^2 Phi. `MatReaction` adds the term `conductor_penalty * Phi` to the residual. When `conductor_penalty` is large inside the cylinder, this forces Phi to zero there.

### Materials

```
[Materials]
  [conductor_mat]
    type        = GenericFunctionMaterial
    prop_names  = 'conductor_penalty'
    prop_values = 'penalty_func'
  []
[]
```

`GenericFunctionMaterial` evaluates the `penalty_func` ParsedFunction at each quadrature point to produce the material property. This avoids subdomains and works directly on the unmodified structured mesh.

### Functions

```
[Functions]
  [analytic_bc]
    type       = ParsedFunction
    expression = '-${E0} * y * (1.0 - ${R}*${R} / (x*x + y*y + 1e-20))'
  []
  [penalty_func]
    type       = ParsedFunction
    expression = '1000.0 * 0.5 * (1.0 - tanh(20.0 * (sqrt(x*x + y*y) - ${R})))'
  []
  [analytic_full]
    type       = ParsedFunction
    expression = 'if(x*x+y*y < ${R}*${R}, 0, -${E0} * y * (1.0 - ${R}*${R}/(x*x+y*y+1e-20)))'
  []
[]
```

`analytic_bc` gives the far-field potential for the outer boundary DirichletBC. `penalty_func` uses a tanh transition of width 1/20 m centred at r = R; the factor 0.5(1 - tanh) is 1 for r < R and 0 for r > R, multiplied by the penalty magnitude 1000. `analytic_full` encodes the exact solution in both regions for L2 error computation.

### Boundary Conditions

```
[BCs]
  [outer_analytic]
    type     = FunctionDirichletBC
    variable = phi
    boundary = 'left right top bottom'
    function = analytic_bc
  []
[]
```

All four outer walls receive the exact analytic potential as a Dirichlet condition. At distance L = 3 m (with R = 0.5 m), the correction term R²/r² is at most 0.028, so the far-field approximation is very accurate.

### Postprocessors

`ElementL2Error` computes the global L2 norm of (phi - phi_exact). `PointValue` at three radii along the y-axis provides spot checks against the analytic formula. `ElementExtremeValue` records the maximum and minimum potentials.

### Executioner

```
[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]
```

Steady-state Newton solve with BoomerAMG algebraic multigrid preconditioning. Tight tolerances ensure the L2 error comparison is meaningful.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case99-cylinder-uniform-field \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case99_cylinder_uniform_field.i 2>&1 | tail -30'
```

---

## Expected Results

**Analytic values at postprocessor points:**

| Location | Exact Phi [V] | Notes |
|----------|--------------|-------|
| (0, 0.5) — cylinder surface | 0 | conductor equipotential |
| (0, 1.0) | -0.75 | -E₀ · 1 · (1 - 0.25) |
| (0, 2.0) | -1.875 | -E₀ · 2 · (1 - 0.0625) |

**Global extrema:**

The potential ranges from approximately -E₀ · L to +E₀ · L (roughly ±3 V), with the asymmetric distortion of the field lines around the cylinder producing slightly enhanced extrema near the top and bottom of the cylinder where the field is strongest.

**L2 error:**

With 40 × 40 elements and the tanh penalty transition, the L2 error is expected in the range 0.01–0.05 V, dominated by the smeared transition region near r = R. Refining the mesh or sharpening the tanh (increasing the coefficient from 20 to, say, 50) reduces this error at the cost of a stiffer linear system.

**Physical interpretation:**

The potential contours outside the cylinder are dipole-distorted circles. The electric field lines — perpendicular to equipotentials — are deflected around the cylinder and terminate normally on the conductor surface. The induced charge σ_s = 2ε₀E₀ sinθ is positive on the top half (y > 0) and negative on the bottom half (y < 0), consistent with a dipole moment p = 4πε₀ R² E₀ pointing in the +y direction.

---

## Key Takeaways

- **Penalty method for conductors**: Rather than carving a hole in the mesh or enforcing Dirichlet BCs on a curved boundary, a large `MatReaction` coefficient drives the interior potential to zero — a simple and robust approach for immersed-boundary problems.
- **GenericFunctionMaterial**: Spatial variation of a material property can be encoded directly from a ParsedFunction, removing the need for separate subdomains when the geometry is simple.
- **tanh regularisation**: The smooth tanh transition avoids the numerical difficulty of a discontinuous penalty coefficient, allowing Newton's method to converge without special treatment of the interface.
- **Far-field Dirichlet BCs**: Applying the analytic far-field solution as a Dirichlet condition replaces the need for absorbing or mapped-infinity boundaries, provided the domain is large enough (here L/R = 6).
- **L2 error verification**: `ElementL2Error` against a ParsedFunction analytic solution provides a single scalar metric to quantify accuracy, useful for convergence studies.
- **Connection to MIT 6.641 Lecture 11**: This is the textbook BVP in polar coordinates (Zahn, Section 2.5), solved here by a Cartesian FEM approach rather than the separation-of-variables method.
