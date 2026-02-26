# Case 30: Rectangular Waveguide Cutoff Frequencies

## Overview

A rectangular metallic waveguide supports electromagnetic modes that can only
propagate above a characteristic **cutoff frequency**. Below this frequency the
mode is evanescent — it decays exponentially along the guide rather than
propagating. Designing filters, antennas, and microwave circuits requires knowing
these cutoff wavenumbers precisely.

This case finds the TM-mode cutoff wavenumbers of a 2:1 rectangular waveguide
(a = 2, b = 1) by solving the 2D Helmholtz eigenvalue problem on the cross-section.
The approach follows Haus, *Electromagnetic Noise and Quantum Optical Measurements*
(Springer, 2000), Ch. 2. The MOOSE `Eigenvalue` executioner drives SLEPc to find
the lowest eigenvalues; no time-marching is involved.

---

## The Physics

### Governing Equations

Inside the waveguide cross-section Ω the transverse field satisfies:

```
∇²ψ + k_c² ψ = 0    in Ω

ψ = 0                on ∂Ω  (PEC walls, TM modes)
```

where ψ is the longitudinal electric field component and k_c is the cutoff
wavenumber. This is a standard Helmholtz eigenvalue problem. The analytical
eigenvalues for a rectangle of width a and height b are:

```
k_c²(m,n) = (mπ/a)² + (nπ/b)²    m, n = 1, 2, 3, ...
```

For a = 2, b = 1, the first three TM eigenvalues are:

| Mode   | k_c² (exact) |
|--------|--------------|
| TM₁₁  | 12.337       |
| TM₂₁  | 22.207       |
| TM₁₂  | 42.088       |

### Weak Form

Written as a generalised eigenvalue problem Aψ = λBψ:

```
A = ∫ ∇ψ · ∇v dV          (stiffness — from Diffusion kernel)
B = ∫ ψ · v dV             (mass — from CoefReaction with eigen tag)
λ = k_c²                   (eigenvalue sought by SLEPc)
```

---

## MOOSE Implementation

| Object | Role |
|--------|------|
| `GeneratedMesh` 40×20 QUAD4 | 2:1 rectangular cross-section, a=2, b=1 |
| `Diffusion` kernel | Builds stiffness matrix ∫∇ψ·∇v dV |
| `CoefReaction` (coeff=-1, `eigen_kernel=true`) | Builds mass matrix -∫ψ·v dV with eigenvalue tag |
| `DirichletBC` + `EigenDirichletBC` | Enforces ψ=0 on all four PEC walls |
| `SlepcEigenvalue` executioner | Calls SLEPc to find smallest k_c² values |
| `VectorPostprocessor` on eigenvalues | Extracts computed eigenvalue list to CSV |

The `eigen_kernel = true` flag on `CoefReaction` marks it as the mass-matrix
contribution; MOOSE separates it from the stiffness matrix automatically when
assembling the generalised eigenvalue problem.

---

## Expected Results

The first three computed eigenvalues should match the analytical values within
approximately 0.5%:

| Mode  | Analytical k_c² | MOOSE (expected) |
|-------|-----------------|------------------|
| TM₁₁ | 12.337          | ~12.34           |
| TM₂₁ | 22.207          | ~22.23           |
| TM₁₂ | 42.088          | ~42.13           |

The spatial mode shape for TM₁₁ is a single arch spanning the full cross-section.
TM₂₁ shows two arches side-by-side in the long (x) direction. TM₁₂ shows one arch
in x and two in y, confirming the nodal pattern predicted by the analytical formula.

---

## How to Run

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case30-waveguide-cutoff \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case30_waveguide_cutoff.i 2>&1 | tail -30'
```

Output files produced:
- `case30_waveguide_cutoff_out.e` — Exodus file with mode-shape fields
- `case30_waveguide_cutoff_out.csv` — computed eigenvalues k_c²
