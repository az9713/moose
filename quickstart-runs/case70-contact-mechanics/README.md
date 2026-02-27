# Case 70: Contact Mechanics — Mortar Frictionless Contact

## Overview

Contact mechanics — the study of stresses and deformations when two bodies come into physical contact — is one of the most challenging problems in computational solid mechanics. The difficulty is threefold: the contact region is not known in advance (it depends on the solution), the contact constraint is inherently nonlinear (surfaces can only push, not pull), and the stress field is singular at the contact edge. MOOSE's `contact` module implements mortar-based contact formulation, which uses Lagrange multipliers defined on an auxiliary mortar mesh to enforce the non-penetration constraint in a mathematically rigorous, variationally consistent manner.

This case models two elastic blocks compressed together by a prescribed displacement. The left block ([-1,0] x [-0.5,0.5]) is pushed rightward by a linearly increasing displacement, closing an initial gap of 0.01 units and then compressing the contact interface. The right block ([0,1] x [-0.6,0.6]) is held fixed at its right boundary. Both blocks share identical isotropic elastic material properties (E = 1e6 Pa, nu = 0.3). The mortar formulation handles the geometrically nonlinear finite-strain contact between their opposing faces.

Key concepts demonstrated:

- Multi-block mesh construction using `MeshCollectionGenerator`, `SubdomainIDGenerator`, and `RenameBlockGenerator`
- The `Contact` action block for mortar frictionless contact — a high-level interface that automatically creates the mortar mesh, Lagrange multiplier variables, and constraint objects
- `Physics/SolidMechanics/QuasiStatic` action for automatic setup of finite-strain solid mechanics variables and kernels
- `ComputeFiniteStrainElasticStress` for geometrically nonlinear large-deformation elasticity
- Block-restricted postprocessors to avoid the extra mortar subdomains created by the contact action
- `FunctionDirichletBC` with `preset = true` for prescribed displacement loading

---

## The Physics

### Governing Equations — Finite-Strain Elasticity

Each block satisfies the equilibrium equations in the reference configuration (quasi-static, no inertia):

```
div(P) = 0
```

where P is the first Piola-Kirchhoff stress tensor. For finite strains, the Green-Lagrange strain tensor is:

```
E = (1/2) * (F^T * F - I)
```

where F = I + grad(u) is the deformation gradient. The constitutive law (St. Venant-Kirchhoff) gives the second Piola-Kirchhoff stress S in terms of the Lame parameters:

```
S = lambda * tr(E) * I + 2*mu * E
```

with lambda = nu*E/(1+nu)/(1-2*nu) = 576923 Pa and mu = E/(2*(1+nu)) = 384615 Pa for E = 1e6 Pa, nu = 0.3.

### Mortar Contact Formulation

The frictionless contact constraint requires that:
1. The normal gap between the surfaces is non-negative: g_N >= 0
2. The contact pressure (normal Lagrange multiplier) is non-positive: lambda_N <= 0 (compression only)
3. The complementarity condition: lambda_N * g_N = 0 (either gap is open and pressure is zero, or surfaces are in contact and gap is zero)

The mortar formulation enforces these constraints weakly on a mortar mesh created at the interface. The secondary surface (`lb_right`, the right face of the left block) has nodes constrained relative to the primary surface (`rb_left`, the left face of the right block). Lagrange multipliers represent the contact pressure and are solved simultaneously with the displacements.

### Loading Protocol

The left block starts with an initial displacement of -0.01 (shifted left by the gap amount). The left boundary then moves right at 0.1 units per second:

```
disp_x(left boundary) = 0.1 * t
```

Timeline:
- t = 0 to 0.1 s: left block moves right, gap closes (0.01 / 0.1 m/s = 0.1 s to close)
- t = 0.1 s: surfaces touch, contact initiates
- t = 0.1 to 0.5 s: contact pressure builds as left block compresses right block
- t = 0.5 s: total applied displacement = 0.05 m; interface compression = 0.05 - 0.01 = 0.04 m

### Domain and Mesh

The two blocks have different heights (0.5 vs. 0.6 half-heights) by design — this ensures the meshes are non-conforming at the contact interface, which exercises the mortar interpolation machinery. The left block uses 10x10 elements; the right block uses 10x12 elements.

---

## Input File Walkthrough

The input file is `case70_contact_mechanics.i`.

### `[Mesh]` — Multi-Block Construction

The mesh is assembled from two separately-generated blocks:

```
[left_block]
  type = GeneratedMeshGenerator
  xmin = -1.0  xmax = 0.0  ymin = -0.5  ymax = 0.5
  nx = 10  ny = 10
  boundary_name_prefix = lb
[]
[left_block_id]
  type = SubdomainIDGenerator
  input = left_block
  subdomain_id = 1
[]
[right_block]
  type = GeneratedMeshGenerator
  xmin = 0.0  xmax = 1.0  ymin = -0.6  ymax = 0.6
  nx = 10  ny = 12
  boundary_name_prefix = rb
  boundary_id_offset = 10
[]
...
[combined]
  type = MeshCollectionGenerator
  inputs = 'left_block_id right_block_id'
[]
[block_rename]
  type = RenameBlockGenerator
  input = combined
  old_block = '1 2'
  new_block = 'left_block right_block'
[]
```

`boundary_name_prefix` prefixes all boundary names of each generated mesh (`lb_left`, `lb_right`, `rb_left`, `rb_right`, etc.). `boundary_id_offset = 10` ensures the right block's boundary IDs do not collide with the left block's. `MeshCollectionGenerator` merges the two meshes into a single combined mesh. `RenameBlockGenerator` gives the subdomain IDs human-readable names.

The top-level variable `gap = 0.01` sets the initial separation between the blocks. It is used in the `ICs` block to pre-displace the left block.

### `[Physics/SolidMechanics/QuasiStatic]`

```
[all]
  strain = FINITE
  incremental = true
  add_variables = true
  block = 'left_block right_block'
[]
```

This action automatically creates `disp_x` and `disp_y` variables and the corresponding `StressDivergenceTensors` kernels for quasi-static mechanics. `strain = FINITE` selects the finite-strain (geometrically nonlinear) formulation. `incremental = true` uses an incremental update consistent with `ComputeFiniteStrainElasticStress`. Block restriction to the two physical blocks prevents the action from operating on the mortar subdomains created by the `Contact` block.

### `[AuxVariables]` and `[AuxKernels]`

```
[vonmises]
  order = CONSTANT
  family = MONOMIAL
  block = 'left_block right_block'
[]
[vonmises_kernel]
  type = RankTwoScalarAux
  rank_two_tensor = stress
  scalar_type = VonMisesStress
  block = 'left_block right_block'
[]
```

`RankTwoScalarAux` with `scalar_type = VonMisesStress` extracts the von Mises equivalent stress from the full stress tensor. CONSTANT MONOMIAL shape functions store one value per element — appropriate for element-averaged stress. Block restriction is required here because the `Contact` action adds mortar subdomains that do not have a `stress` material property.

### `[Contact]`

```
[leftright]
  secondary = lb_right
  primary = rb_left
  model = frictionless
  formulation = mortar
[]
```

This single block replaces hundreds of lines of manual mortar constraint setup. The `Contact` action creates the mortar mesh at the interface, declares the Lagrange multiplier variable for contact pressure, and adds the constraint objects. `secondary = lb_right` is the contact surface that will be constrained (the face of the left block that contacts the right block). `primary = rb_left` is the reference surface. `model = frictionless` means only normal forces are transmitted (no tangential friction). `formulation = mortar` selects the mortar method over the simpler node-on-segment approach.

### `[ICs]`

```
[disp_x_left]
  type = ConstantIC
  block = left_block
  variable = disp_x
  value = -${gap}
[]
```

The left block starts shifted 0.01 units to the left (negative x). This creates the initial gap between the contact surfaces. The `${gap}` syntax references the top-level variable defined as `gap = 0.01`. The right block starts with zero displacement — it is already at its natural position.

### `[Materials]`

```
[elasticity_left]
  type = ComputeIsotropicElasticityTensor
  block = left_block
  youngs_modulus = 1.0e6
  poissons_ratio = 0.3
[]
[stress_left]
  type = ComputeFiniteStrainElasticStress
  block = left_block
[]
```

`ComputeIsotropicElasticityTensor` computes the 4th-order elasticity tensor C from E and nu. `ComputeFiniteStrainElasticStress` computes the rotated Cauchy stress using the Jaumann rate, consistent with the `FINITE` strain formulation. Both materials are defined separately for each block to allow future differentiation of material properties.

### `[Executioner]`

```
petsc_options_iname = '-pc_type -mat_mffd_err -pc_factor_shift_type -pc_factor_shift_amount'
petsc_options_value  = 'lu       1e-5          NONZERO               1e-15'
dt = 0.05
end_time = 0.5
nl_rel_tol = 1e-8
nl_abs_tol = 1e-6
nl_max_its = 50
```

The contact problem requires a direct solver (LU factorization) for robustness — iterative methods often struggle with the saddle-point system created by Lagrange multiplier contact. The `NONZERO` shift prevents zero pivots that arise when nodes are not yet in contact. Ten time steps of dt = 0.05 s span the full loading history.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case70-contact-mechanics \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case70_contact_mechanics.i 2>&1 | tail -30'
```

Output files:
- `case70_contact_mechanics_out.e` — Exodus with displacement and von Mises stress fields over time
- `case70_contact_mechanics_out.csv` — Time series of max_disp_x, max_vonmises, avg_vonmises

---

## Expected Results

### Stress vs. Time

| Time (s) | max_disp_x (m) | max_vonmises (Pa) | avg_vonmises (Pa) |
|----------|----------------|-------------------|-------------------|
| 0.0      | 0.000          | 0                 | 0                 |
| 0.05     | 0.005          | 3,382             | 2,339             |
| 0.10     | 0.010          | 6,756             | 4,684             |
| 0.15     | 0.015          | 10,122            | 7,036             |
| 0.20     | 0.020          | 13,481            | 9,394             |
| 0.30     | 0.030          | 20,182            | 14,128            |
| 0.40     | 0.040          | 26,863            | 18,889            |
| 0.50     | 0.050          | 33,528            | 23,675            |

The von Mises stress grows approximately linearly with the applied displacement, consistent with linear elastic behavior at moderate strains. The ratio max_vonmises / avg_vonmises ~ 1.42 reflects the stress concentration at the contact edge (corners of the contact patch experience higher stress than the average).

### Contact Initiation

At t = 0, the left block is displaced -0.01 m (initial gap), so disp_x ranges from -0.01 (initial IC) to 0 (the gap amount). As the left boundary moves rightward at 0.1 m/s, disp_x increases. Contact initiates at t ~ 0.1 s when the gap closes. The stress remains near zero for t < 0.1 s (gap phase) and then rises steeply once contact is established. The CSV shows that all steps have nonzero stress because the initial IC already puts the left block close to the right block, and the first dt = 0.05 s step partially closes the gap.

### Spatial Distribution

The von Mises stress concentrates at the contact interface between the two blocks. Within each block, the stress field is relatively uniform in the bulk and peaks near the contact boundary. The stress singularity at the contact edge (where the contact region meets the free surface) causes the max_vonmises to be higher than the values deep in the interior.

---

## Key Takeaways

- The `Contact` action is the recommended high-level interface for mortar contact in MOOSE. A single block with `secondary`, `primary`, `model`, and `formulation` parameters replaces the need to manually create mortar meshes, Lagrange multiplier variables, and constraint objects. The mortar formulation provides variationally consistent contact force distributions, superior to older node-on-segment methods.
- Multi-block meshes require `boundary_name_prefix` to prevent name collisions between mesh generators, and `boundary_id_offset` to prevent ID collisions. When using `MeshCollectionGenerator`, the input blocks retain their original subdomain IDs; `RenameBlockGenerator` provides human-readable names.
- All postprocessors and AuxKernels that operate on element or node data must be block-restricted to `'left_block right_block'` when using mortar contact. The `Contact` action creates additional mortar subdomains, and objects like `RankTwoScalarAux` that reference material properties will fail on these subdomains because no stress material is defined there.
- `ComputeFiniteStrainElasticStress` requires `strain = FINITE` and `incremental = true` in the SolidMechanics action. These three settings together specify a geometrically nonlinear formulation using an incremental multiplicative decomposition of the deformation gradient — appropriate for large displacements even when the material remains elastic.
- Direct solvers (LU via PETSc's `-pc_type lu`) are typically required for contact problems. The saddle-point structure of the mortar system (displacements and Lagrange multiplier unknowns together) can cause iterative preconditioners to fail or converge slowly. The `NONZERO` pivot shift prevents factorization failures when nodes are not in contact and have zero diagonal entries.
- `preset = true` in `FunctionDirichletBC` applies the prescribed displacement value at the beginning of each time step, before the Newton solve, rather than enforcing it as a constraint during the solve. This is the preferred approach for displacement-controlled loading because it eliminates the Dirichlet DOFs from the system and simplifies the Jacobian.
- `volumetric_locking_correction = true` in `[GlobalParams]` applies the selective reduced integration correction to QUAD4 elements, preventing spurious stiffness (locking) in near-incompressible materials. With nu = 0.3, locking is not severe, but the correction is good practice for general contact simulations.
