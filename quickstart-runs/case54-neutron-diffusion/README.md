# Case 54: 1-Group Neutron Diffusion — Bare Slab Criticality

## Overview

This case introduces the most fundamental problem in reactor physics: finding the critical size of a bare slab reactor. A reactor is "critical" when the neutron chain reaction is exactly self-sustaining — each fission event produces exactly one neutron that causes another fission. The effective multiplication factor k_eff is the ratio of neutrons produced to neutrons lost, so criticality means k_eff = 1.

The governing equation is the one-speed neutron diffusion equation, an eigenvalue problem where the fission source term is scaled by 1/k. MOOSE's `Eigenvalue` executioner finds k_eff and the corresponding fundamental mode (the flux shape). For a bare slab with vacuum boundaries, the analytical solution is a half-cosine, and the critical thickness is exactly determined by the material properties.

This case demonstrates how MOOSE handles generalized eigenvalue problems in reactor physics. The key technique is tagging the fission kernel with `extra_vector_tags = 'eigen'`, which places it in the B-matrix of the generalized eigenproblem Ax = (1/k)Bx rather than in the standard residual. The `EigenDirichletBC` objects work alongside ordinary `DirichletBC` to enforce zero-flux boundary conditions in both the A and B systems.

## The Physics

### Governing Equation

The one-group steady-state neutron diffusion equation is:

```
-D * laplacian(phi) + Sigma_a * phi = (1/k) * nu_Sigma_f * phi
```

where:
- `phi` — scalar neutron flux (neutrons/cm^2/s, arbitrary normalization)
- `D = 1.0 cm` — diffusion coefficient (related to mean free path)
- `Sigma_a = 0.10 /cm` — macroscopic absorption cross section
- `nu_Sigma_f = 0.15 /cm` — nu (neutrons per fission) times fission cross section
- `k` — effective multiplication factor (the eigenvalue)

This is a generalized eigenvalue problem: find k (and the flux shape phi) such that the equation is satisfied with phi > 0 everywhere inside the domain.

### Analytical Criticality Condition

Rearranging for a uniform infinite medium, the geometric buckling B^2 must satisfy:

```
B^2 = (nu_Sigma_f - Sigma_a) / D = (0.15 - 0.10) / 1.0 = 0.05 /cm^2
```

For a bare slab with vacuum boundaries (phi = 0 at x = 0 and x = L), the fundamental mode is:

```
phi(x) = A * sin(pi * x / L)
```

which requires B = pi/L, giving the critical half-thickness:

```
L_crit = pi / sqrt(B^2) = pi / sqrt(0.05) = 14.049 cm
```

This case uses L = 14.05 cm, so k_eff should come out to exactly 1.0000.

### Boundary Conditions

Both ends use vacuum boundary conditions, which model a reactor-to-air interface. In the diffusion approximation, the flux extrapolates to zero at a small distance beyond the physical boundary; for this educational case we place the boundary directly at the extrapolated zero:

```
phi(0) = 0       (left vacuum BC)
phi(L) = 0       (right vacuum BC)
```

The top and bottom surfaces (y = 0, y = 0.5) use the natural (zero-flux-gradient) boundary condition automatically — this is correct for the quasi-1D geometry where nothing varies in y.

### Material Properties

| Parameter | Value | Units | Physical meaning |
|-----------|-------|-------|-----------------|
| D | 1.0 | cm | Diffusion coefficient |
| Sigma_a | 0.10 | /cm | Absorption cross section |
| nu_Sigma_f | 0.15 | /cm | Fission production cross section |
| k_eff (expected) | 1.0000 | dimensionless | Exactly critical |

### Domain and Mesh

- Geometry: 1D slab, x in [0, 14.05] cm, modeled as 2D quasi-1D (100 x 2 elements)
- The y-dimension (height 0.5 cm, 2 elements) exists only because MOOSE requires 2D geometry here; all physics is in x
- Element size: ~0.14 cm in x

## Input File Walkthrough

### `[Mesh]`

```
[gmg]
  type = GeneratedMeshGenerator
  dim = 2
  nx = 100
  ny = 2
  xmin = 0
  xmax = ${L}    # 14.05 cm
  ymin = 0
  ymax = 0.5
[]
```

The top-level variable `L = 14.05` is referenced with `${L}`. The mesh is purely structured quads. Using `ny = 2` gives exactly two elements in y, enough to avoid degenerate geometry while keeping the problem essentially 1D.

### `[Variables]`

```
[phi]   # scalar neutron flux
[]
```

A single variable for the scalar flux, using the default first-order Lagrange basis (linear elements). No initial condition is needed — the Eigenvalue solver initializes from a uniform guess.

### `[Kernels]`

Three kernels form the weak form of the equation:

**Diffusion** (`MatDiffusion`): Implements the term integral(D * grad(phi) . grad(test) dV). The diffusivity `D` is a material property defined in `[Materials]`.

**Absorption** (`Reaction`): Implements integral(Sigma_a * phi * test dV) = integral(0.10 * phi * test dV). This is a loss term — neutrons removed by absorption. Note the positive sign; `Reaction` with a positive rate is a loss.

**Fission** (`Reaction` with `extra_vector_tags = 'eigen'`): Implements the fission source as integral(nu_Sigma_f * phi * test dV). The rate is -0.15 (negative) because `Reaction` adds rate * phi * test to the residual; we want this on the right-hand side (source), so the negative sign moves it there. The `extra_vector_tags = 'eigen'` tag tells MOOSE to place this contribution in the B-matrix rather than the A-matrix, which is how the generalized eigenproblem Ax = (1/k)Bx is constructed.

### `[Materials]`

```
[diffusion_coeff]
  type = GenericConstantMaterial
  prop_names = 'D'
  prop_values = '1.0'
[]
```

`GenericConstantMaterial` is the simplest material in MOOSE — it just declares a constant-valued material property. Here D = 1.0 cm is registered so that `MatDiffusion` can look it up by name.

### `[BCs]`

Four boundary conditions are required for eigenvalue problems with Dirichlet constraints:

- `DirichletBC` at left and right: enforces phi = 0 in the standard (A-matrix) system
- `EigenDirichletBC` at left and right: enforces phi = 0 in the eigen (B-matrix) system

Both sets are required. Without `EigenDirichletBC`, the B-matrix does not have the Dirichlet constraints applied correctly, and the eigenvalue solve produces spurious modes.

### `[Executioner]`

```
type = Eigenvalue
solve_type = PJFNK
petsc_options_iname = '-pc_type -pc_hypre_type'
petsc_options_value = 'hypre    boomeramg'
```

`Eigenvalue` uses SLEPc (the Scalable Library for Eigenvalue Problem Computations, built on PETSc) to find the dominant eigenvalue. `PJFNK` (Preconditioned Jacobian-Free Newton Krylov) is the default solve type for the linear systems within each SLEPc iteration. The `boomeramg` algebraic multigrid preconditioner from Hypre is highly effective for diffusion-dominated problems.

### `[VectorPostprocessors]`

```
[flux_profile]
  type = LineValueSampler
  variable = phi
  start_point = '0 0.25 0'
  end_point = '${L} 0.25 0'
  num_points = 101
  sort_by = x
[]
```

`LineValueSampler` extracts values along a line through the mesh centerline (y = 0.25). The 101 sample points give a smooth profile for plotting. This writes a CSV file with one row per sample point.

### `[Outputs]`

Both `exodus = true` (for ParaView visualization) and `csv = true` (for the postprocessor scalars and vector sampler data) are enabled.

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case54-neutron-diffusion \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case54_neutron_diffusion.i 2>&1 | tail -30'
```

## Expected Results

The solver output will show the SLEPc eigenvalue iteration converging. The final k_eff appears in the console as the dominant eigenvalue. For L = 14.05 cm with the given cross sections:

```
Eigenvalue: 1.0000  (k_eff)
```

The flux profile sampled in `flux_profile_0001.csv` will show a perfect half-sine shape:

```
x = 0:        phi ~ 0.000   (vacuum BC)
x = 7.025:    phi ~ 1.000   (peak at center, normalized by SLEPc)
x = 14.05:    phi ~ 0.000   (vacuum BC)
```

The midpoint flux value reported by the `k_eff` postprocessor (confusingly named — it samples phi at the midpoint, not k_eff directly) should be close to the SLEPc-normalized peak.

To verify the result analytically: the eigenvalue should match the ratio nu_Sigma_f / (Sigma_a + D * (pi/L)^2) = 0.15 / (0.10 + 1.0 * 0.05) = 0.15 / 0.15 = 1.000.

## Key Takeaways

- The one-group neutron diffusion equation is a generalized eigenvalue problem: Ax = (1/k)Bx, where A contains diffusion and absorption, and B contains the fission source.
- MOOSE handles this by tagging fission kernels with `extra_vector_tags = 'eigen'`, separating them into the B-matrix automatically.
- `EigenDirichletBC` must be added alongside `DirichletBC` whenever Dirichlet conditions appear in an eigenvalue problem — one for each matrix.
- The critical slab thickness follows directly from the material buckling: L_crit = pi / sqrt((nu_Sigma_f - Sigma_a) / D).
- The fundamental mode (lowest-k eigenfunction) is always positive everywhere inside the domain — this is a physical requirement (flux cannot be negative) and is guaranteed by the Perron-Frobenius theorem for this class of problems.
- The `Eigenvalue` executioner requires SLEPc, which is included in the standard MOOSE Docker image.
