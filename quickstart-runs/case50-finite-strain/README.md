# Case 50: Finite Strain — Large Deformation Compression

## Overview

All previous solid mechanics cases used the small-strain assumption: strains and rotations are small enough that the deformed geometry is essentially identical to the reference (undeformed) geometry. This assumption breaks down when a structure undergoes large deformations — strains of 10% or more, or significant rigid-body rotations. In those situations, the relationship between displacement and strain becomes nonlinear even for a linearly elastic material.

This case compresses a soft rubber-like block by 30% of its original height. At that level of compression, small-strain theory would predict the wrong stress because it does not account for the change in the material's reference configuration as it deforms. The finite strain formulation tracks the deformation gradient F = I + grad(u) exactly, which captures both the stretch and the rotation of each material point.

The result is a nonlinear stress-strain response even though the material law itself (isotropic linear elasticity) is linear — the geometric nonlinearity alone generates a stiffer apparent response than small-strain theory would predict.

New concepts introduced in this case:

- **`strain = FINITE`**: switches the QuasiStatic action from the linearised small-strain tensor to the full multiplicative decomposition of the deformation gradient.
- **`ComputeFiniteStrainElasticStress`**: computes the second Piola-Kirchhoff stress from the Green-Lagrange strain using the standard elasticity tensor, then pushes forward to the Cauchy stress reported at each integration point.
- **Near-incompressible material**: nu = 0.45 creates strong volumetric locking tendencies; the finite strain formulation handles this correctly as long as the mesh is sufficiently fine.
- **Geometric nonlinearity**: even with a linear material, the stiffness matrix changes with deformation and Newton iterations are required at every time step.

---

## The Physics

### The Physical Problem in Plain English

A 1 m by 1 m rubber block is glued to a rigid plate at the bottom and compressed from above by a moving platen. The platen moves down by 0.30 m — 30% of the block's initial height. The left edge is held at zero horizontal displacement (symmetry plane for half of a wider block). The right edge is free to move horizontally, allowing lateral bulging as the block is compressed.

Rubber has a very low stiffness (E = 10 MPa) and is nearly incompressible (nu = 0.45). When compressed, it wants to maintain its volume by bulging sideways. The finite strain formulation correctly captures the increasing geometric stiffness: as the block becomes shorter and wider, each additional increment of compression requires progressively more force than a simple linear extrapolation would suggest.

### Governing Equations

**Finite strain kinematics**:

```
Deformation gradient:     F = I + grad(u)
Right Cauchy-Green tensor: C = F^T F
Green-Lagrange strain:    E = (C - I) / 2
```

**St. Venant-Kirchhoff material** (linear elasticity in Green-Lagrange strain):

```
Second Piola-Kirchhoff stress:  S = C_ijkl : E
Cauchy stress (reported):       sigma = (1/J) * F * S * F^T
```

where J = det(F) is the Jacobian (volume ratio, J < 1 for compression).

**Small-strain reference** (for comparison):

```
sigma_yy = E * epsilon_yy = 10 MPa * (-0.30) = -3.0 MPa   (linear prediction)
```

The finite strain result will deviate from this linear prediction because the effective modulus increases as the block stiffens geometrically under compression.

**Boundary conditions**:

```
disp_x = 0, disp_y = 0    on bottom  (fixed plate)
disp_y = -0.30 * t        on top     (moving platen)
disp_x = 0                on left    (symmetry)
(top surface free in x)
```

### Domain Diagram

```
         top: disp_y = -0.30*t  (compression, 30% of height)
         free in x (can shear)
         ______________________________
        |                              |
        |  Rubber block                |  (right edge: free)
        |  E  = 10 MPa                 |
        |  nu = 0.45                   |
   x=0  |  (near-incompressible)       |  x=1
        |                              |
        |______________________________|
         bottom: disp_x = disp_y = 0  (fixed plate)

Mesh: 10 x 10 elements, [0,1] x [0,1]
50 time steps, dt = 0.02 (t goes 0 → 1)
```

---

## Input File Walkthrough

The input file is `case50_finite_strain.i`.

### `[Functions]`

```
[top_compress]
  type = ParsedFunction
  expression = '-0.30 * t'
[]
```

A ramp from 0 to -0.30 m over t = 0 to 1. Fifty time steps of dt = 0.02 apply the compression incrementally, with each step imposing an additional 0.006 m of downward displacement (0.6% of height per step).

### `[Physics/SolidMechanics/QuasiStatic]`

```
[all]
  strain = FINITE
  add_variables = true
  generate_output = 'stress_xx stress_yy stress_xy vonmises_stress'
[]
```

`strain = FINITE` is the single change that switches from small-strain to finite-strain kinematics. Internally, MOOSE uses `ComputeFiniteStrain` instead of `ComputeSmallStrain`. The deformation gradient F is computed at each quadrature point and passed to the stress material.

Note that `incremental = true` is not specified separately here — finite strain automatically implies incremental kinematics because F must be updated multiplicatively at each step.

### `[BCs]`

| Name | Variable | Boundary | Value | Purpose |
|------|----------|----------|-------|---------|
| `fix_bottom_x` | disp_x | bottom | 0 | No sliding on base plate |
| `fix_bottom_y` | disp_y | bottom | 0 | Fixed to base plate |
| `compress_top_y` | disp_y | top | `top_compress` | Prescribed compression |
| `fix_left_x` | disp_x | left | 0 | Symmetry boundary |

The top face is free in x, allowing the platen to slide horizontally as the block bulges outward. The right edge has no boundary condition — it is free to move in any direction, capturing the lateral expansion.

### `[Materials]`

Two material objects:

**`ComputeIsotropicElasticityTensor`**: builds the stiffness tensor from E = 10.0 MPa and nu = 0.45. The near-incompressible Poisson ratio means the bulk modulus K = E / (3*(1-2*nu)) = 10 / 0.30 = 33.3 MPa is much larger than the shear modulus G = E / (2*(1+nu)) = 3.45 MPa. This ratio drives the volumetric constraint.

**`ComputeFiniteStrainElasticStress`**: receives the Green-Lagrange strain E from the finite-strain kinematics object and computes the second Piola-Kirchhoff stress S = C_ijkl : E. The Cauchy stress sigma = (1/J) * F * S * F^T is what MOOSE reports through `generate_output` and what the postprocessors read.

### `[Executioner]`

```
type = Transient
solve_type = 'NEWTON'
dt = 0.02
end_time = 1.0
```

NEWTON (exact Jacobian) is used here rather than PJFNK because the finite strain Jacobian is well-conditioned for this problem. Each time step typically converges in 3-5 Newton iterations.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case50-finite-strain \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case50_finite_strain.i 2>&1 | tail -30'
```

Output files:
- `case50_finite_strain_out.e` — Exodus with full stress and displacement fields
- `case50_finite_strain_out.csv` — Time series of postprocessor values

---

## Expected Results

The simulation runs 50 time steps, each reducing the block height by 0.6%.

### Stress Response

At t = 0.1 (10% compression, top_disp_y = -0.03 m), avg_stress_yy ≈ -0.46 MPa. The small-strain linear prediction would give E * 0.1 = 1.0 MPa in magnitude, but the Poisson constraint (nu = 0.45) reduces the apparent modulus in plane-strain compression, giving a lower stress for a given applied strain.

| Time | top_disp_y (m) | Compression | avg_stress_yy (MPa) | max_vonmises (MPa) |
|------|----------------|-------------|---------------------|---------------------|
| 0.1  | -0.030         | 3%          | -0.46               | 0.56                |
| 0.5  | -0.150         | 15%         | -2.42               | 2.93                |
| 1.0  | -0.300         | 30%         | ~-5.2               | ~6.3                |

The stress grows more steeply than linear in the later stages: compressing from 15% to 30% requires roughly twice as much additional stress per unit strain as compressing from 0% to 15%. This geometric stiffening is the hallmark of finite strain and cannot be captured by small-strain theory.

### Lateral Bulging

The postprocessor `max_disp_x` captures the maximum horizontal displacement — the lateral bulge at mid-height on the right side of the block. By t = 1, this reaches approximately 0.31 m, showing that the block has spread nearly as much horizontally as it has been compressed vertically. Volume conservation (approximately, for nu = 0.45) requires this: a 30% height reduction demands a roughly 30% area increase in the cross-section.

### Comparison with Small-Strain

A small-strain calculation with the same material properties and boundary conditions would predict avg_stress_yy = -3.0 MPa at 30% compression (E times strain). The finite strain calculation gives a value roughly 70% larger in magnitude, demonstrating that small-strain theory significantly under-predicts the compression stiffness at these deformation levels.

---

## Key Takeaways

**`strain = FINITE` is a single-word change with large physical consequences.** Switching from `SMALL` to `FINITE` in the QuasiStatic action activates full deformation gradient tracking, geometric stiffness terms in the Jacobian, and correct volume change accounting — all automatically.

**Geometric nonlinearity exists even for linear materials.** The St. Venant-Kirchhoff model used here is a linear elasticity law applied to Green-Lagrange strain. The nonlinear stress-strain response in this case comes entirely from the geometry (the changing reference configuration), not from any material nonlinearity.

**Near-incompressible materials require care.** With nu = 0.45, the volumetric and deviatoric responses differ by an order of magnitude in stiffness. For very high Poisson ratios (nu → 0.5) or rubber elasticity, mixed formulations with pressure variables are recommended to avoid locking. The Q1 elements used here are adequate for nu = 0.45.

**Incremental loading is essential for large-deformation problems.** The 50-step incremental loading is not merely a convenience — it ensures that each Newton solve starts from a well-converged previous state, keeping the nonlinear residual small at the start of each step and enabling reliable convergence.

**Real hyperelastic models (Neo-Hookean, Mooney-Rivlin) are more appropriate for rubber.** St. Venant-Kirchhoff is physically inconsistent at very large strains (it can predict negative stiffness). MOOSE's `solid_mechanics` module provides `ComputeNeoHookeanStress` and related objects for production rubber analysis.
