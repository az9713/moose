# Case 102: Membrane Levitation Stability — Electrostatic Pull-In Eigenvalue Problem

## Overview

A thin elastic membrane (string) of tension S is stretched between two clamped supports at x = 0 and x = l = 1 m. The membrane is suspended a distance d above a grounded electrode, and a voltage V is applied across the gap. The electrostatic attraction creates a force density that is proportional to the local displacement — a destabilising effect that reduces the effective stiffness of each membrane mode.

The linearised equation of motion for small displacement ξ about the equilibrium position is:

```
rho * d²xi/dt² = S * d²xi/dx² + (eps0 * V² / d³) * xi
```

The electrostatic term (ε₀V²/d³)ξ acts as a negative spring constant: unlike mechanical spring-back, the electric force increases with displacement, making the system conditionally stable. This is the **pull-in instability** (also called snap-through or electrostatic collapse) that limits the operating range of MEMS actuators, electrostatic speakers, and ink-jet printer nozzles.

Eigenvalue analysis (with ρ = 1, S = 1, l = 1, d = 1, ε₀ = 1 normalised) gives natural frequencies:

```
omega_n² = S * (n*pi/l)² - Fe,    Fe = eps0 * V² / d³
```

The n = 1 mode goes unstable (ω₁² < 0) when Fe > π², i.e. V > V_crit = π ≈ 3.14159. This case solves at V = 2 (so F_e = 4), which is stable, and uses the SLEPc Krylov-Schur eigensolver to compute the first four mode frequencies and shapes.

This is MIT 6.641 Lecture 18 material on stability and instabilities in electromagnetic systems. MOOSE models the problem as a generalised eigenvalue problem A·x = ω²·B·x, using the standard `Eigenvalue` executioner with SLEPc as the backend.

---

## The Physics

**Linearised equation of motion:**

```
rho * d²xi/dt² = S * d²xi/dx² + Fe * xi
```

Substituting xi ~ X(x) exp(j omega t) and separating:

```
-S * d²X/dx² - Fe * X = omega² * rho * X
```

This is the generalised eigenvalue problem (A - omega² B) X = 0 with:
- **A** = stiffness operator = -S nabla² - Fe (tension minus electric softening)
- **B** = mass matrix = rho (identity)

**Boundary conditions:**

```
xi(0) = 0,    xi(1) = 0    (clamped supports)
```

**Predicted eigenvalues** (normalised S = 1, rho = 1, F_e = 4, l = 1):

| Mode n | Formula omega_n² = (n pi)² - F_e | Value |
|--------|-----------------------------------|-------|
| 1 | pi² - 4 | 5.870 |
| 2 | 4pi² - 4 | 35.478 |
| 3 | 9pi² - 4 | 84.783 |
| 4 | 16pi² - 4 | 153.914 |

**Pull-in condition:**

The mode n = 1 becomes unstable when omega_1² < 0, i.e.:

```
Fe > pi²    =>    V > V_crit = pi * sqrt(S / eps0) * d^{3/2} / l
```

With normalised parameters V_crit = π ≈ 3.14159. At V = 2 the system is approximately 36% below the pull-in voltage.

**Domain:** 1D line segment [0, 1], 200 elements (fine enough to resolve mode shapes accurately).

---

## Input File Walkthrough

### Mesh

```
[Mesh]
  [line]
    type = GeneratedMeshGenerator
    dim  = 1
    nx   = 200
    xmin = 0
    xmax = 1.0
  []
  [rename]
    type = RenameBoundaryGenerator
    input = line
    old_boundary = 'left right'
    new_boundary = 'clamp_left clamp_right'
  []
[]
```

A 1D mesh with 200 EDGE2 elements. The `RenameBoundaryGenerator` gives semantically meaningful boundary names, improving input file readability. A 1D mesh is appropriate because the membrane displacement varies only along x.

### Variables

```
[Variables]
  [xi]
    order  = FIRST
    family = LAGRANGE
  []
[]
```

The membrane displacement ξ [m] is the eigenfunction. The eigensolver finds the normalised shapes sin(nπx/l) and the associated eigenvalues ω²_n.

### Kernels — Stiffness Matrix A

```
[Kernels]
  [tension_stiffness]
    type     = Diffusion
    variable = xi
  []
  [electric_destabilise]
    type        = CoefReaction
    variable    = xi
    coefficient = -${F_e}
  []
```

`Diffusion` contributes +S∫(dxi/dx)(dv/dx)dx to the weak form — i.e., the -d²xi/dx² strong-form stiffness. `CoefReaction` with coefficient -F_e contributes -F_e∫xi·v dx, implementing the destabilising electric term. Together these form the A matrix of the generalised eigenproblem.

### Kernels — Mass Matrix B

```
  [mass]
    type              = MatReaction
    variable          = xi
    reaction_rate     = mass_rho
    extra_vector_tags = 'eigen'
  []
[]
```

The mass kernel is tagged `'eigen'` so MOOSE routes it to the B system rather than the A system. The B matrix here is simply the identity (mass density ρ = 1), so the generalised problem A x = ω² B x reduces to A x = ω² x. The `extra_vector_tags` mechanism is the standard MOOSE way to separate contributions to the two matrices in an eigenvalue problem.

### Boundary Conditions

```
[BCs]
  [clamp_left]
    type     = DirichletBC
    variable = xi
    boundary = clamp_left
    value    = 0
  []
  [eigen_clamp_left]
    type     = EigenDirichletBC
    variable = xi
    boundary = clamp_left
  []
  ...same pair for clamp_right...
[]
```

Each clamped end requires two BC objects: a standard `DirichletBC` for the A matrix, and an `EigenDirichletBC` for the B matrix. Without `EigenDirichletBC`, the boundary nodes contribute non-zero entries to B at the constrained DOFs, which generates spurious near-zero eigenvalues (the "spectral pollution" problem). The `EigenDirichletBC` zeroes out those entries in the B system.

### Executioner

```
[Executioner]
  type = Eigenvalue
  solve_type = KRYLOVSCHUR
  n_eigen_pairs     = 4
  which_eigen_pairs = TARGET_MAGNITUDE
  petsc_options_iname = '-eps_target -st_type -st_ksp_type -st_pc_type -st_pc_factor_mat_solver_type'
  petsc_options_value = '5.0         sinvert  preonly       lu           mumps'
[]
```

The `Eigenvalue` executioner calls SLEPc's Krylov-Schur solver to find 4 eigenpairs. `TARGET_MAGNITUDE` with shift 5.0 and shift-invert (`st_type = sinvert`) finds eigenvalues closest to 5.0 first — i.e., mode 1 at ω² ≈ 5.87. MUMPS provides the direct LU factorisation required for the shift-invert spectral transformation.

### Outputs

```
[Outputs]
  exodus     = true
  csv        = true
  execute_on = FINAL
[]
```

`execute_on = FINAL` ensures the eigenvector fields are written only once (at convergence), not at every iteration. The VectorPostprocessor `Eigenvalues` writes a CSV table of all four complex eigenvalues (imaginary parts should be zero for this real symmetric problem).

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case102-membrane-levitation \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case102_membrane_levitation.i 2>&1 | tail -30'
```

---

## Expected Results

**Computed eigenvalues ω²_n (should match analytic within ~0.01%):**

| Mode | Analytic omega² | Expected MOOSE | Mode shape |
|------|----------------|----------------|-----------|
| 1 | 5.870 | ~5.870 | sin(pi x) — one half-arch |
| 2 | 35.478 | ~35.478 | sin(2 pi x) — full arch |
| 3 | 84.783 | ~84.783 | sin(3 pi x) — 3 half-arches |
| 4 | 153.914 | ~153.914 | sin(4 pi x) — 2 full arches |

All eigenvalues are positive at V = 2 < V_crit = π, confirming the system is stable.

**Mode shape diagnostics:**

The postprocessor `xi_sq_mid` measures ξ at x = 0.5. For mode 1 (antinode at centre) this is the maximum displacement; for mode 2 (node at centre) it is zero. The postprocessor `xi_sq_quarter` at x = 0.25 is an antinode for mode 2 and a node for mode 4.

**Physical interpretation of the destabilised frequencies:**

Without the electric field (F_e = 0), the natural frequencies would be ω²_n = (nπ)² = {9.87, 39.48, 88.83, 157.91}. The electric destabilisation reduces each frequency by F_e = 4, uniformly across all modes. This is because the destabilising force is proportional to ξ (spatially uniform coefficient), so it adds the same constant F_e to the eigenvalue of every mode — a rigid shift of the entire spectrum downward.

**Pull-in instability at V_crit = π:**

If the simulation were rerun with F_e = π² ≈ 9.87 (V = π), mode 1 would give ω²_1 = π² - π² = 0 — marginal stability. With F_e > π² the first eigenvalue goes negative, meaning the mode is unstable and the membrane collapses to the electrode. This is the pull-in instability that sets the maximum operating voltage for electrostatic MEMS devices.

---

## Key Takeaways

- **Generalised eigenvalue problem A x = omega² B x**: MOOSE's `Eigenvalue` executioner with `extra_vector_tags = 'eigen'` provides a clean way to partition kernel contributions between the stiffness (A) and mass (B) matrices.
- **EigenDirichletBC is mandatory**: Omitting it produces spurious near-zero eigenvalues from unconstrained B-matrix entries at boundary nodes. Both `DirichletBC` (for A) and `EigenDirichletBC` (for B) are required at every constrained DOF.
- **Shift-invert spectral transformation**: The SLEPc shift-invert approach (`st_type = sinvert`) with target near the expected smallest eigenvalue is essential for finding interior eigenvalues efficiently on a fine mesh.
- **CoefReaction vs. MatReaction**: `CoefReaction` with a negative coefficient implements a uniform softening term (-F_e ξ) that reduces stiffness. `MatReaction` with the `eigen` tag provides the mass B-matrix contribution.
- **Electrostatic pull-in**: The uniform downward shift of all modal frequencies by F_e = ε₀V²/d³ means the lowest mode goes unstable first, at the critical voltage V_crit = π sqrt(S/ε₀) d^{3/2}/l. This is the physical pull-in condition governing MEMS actuators.
- **Connection to MIT 6.641 Lecture 18**: This is Zahn's membrane levitation stability problem, canonical in the study of electromechanical system stability and pull-in phenomena.
