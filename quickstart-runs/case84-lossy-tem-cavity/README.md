# Case 84: Lossy TEM Cavity Resonator — Q Factor

## Overview

This case solves the 1D Helmholtz eigenvalue problem in a PEC-walled cavity filled with a **lossy dielectric** material, computing the complex resonant eigenvalues and deriving the **Q factor** of each cavity mode.

The cavity has:
- PEC walls at x = 0 and x = D = 1.0 m
- Lossy dielectric filling with complex permittivity ε_r = 4.0 − j 0.2
- Loss tangent: tan δ = ε''/ε' = 0.2/4.0 = 0.05 (moderate loss)

## Physics

### Governing Equation

Inside the cavity the electric field satisfies the 1D Helmholtz eigenvalue equation:

```
d²E/dx² + λ ε_r E = 0,    E(0) = E(D) = 0
```

where λ = (ω/c)² is the eigenvalue (squared normalised wavenumber) and ε_r = ε' − j ε'' is the complex relative permittivity.

### Complex Eigenvalues and Q Factor

With lossy material the eigenvalues become **complex**: λ = λ_real + j λ_imag.

Perturbation theory for small loss (ε'' << ε') gives:

```
λ_m ≈ (mπ/D)²/ε'  ×  [1 − j ε''/ε']
```

The **Q factor** of the m-th resonant mode is:

```
Q_m = |Re(λ_m)| / (2 |Im(λ_m)|)  ≈  ε' / (2 ε'')
```

For ε' = 4.0, ε'' = 0.2: **Q ≈ 10.0** (same for all modes in the uniform filling approximation).

### Coupled Real/Imaginary Formulation

Writing E = E_real + j E_imag and splitting the equation into real and imaginary parts yields the 2×2 coupled eigenvalue system:

```
d²E_real/dx² + λ (ε' E_real + ε'' E_imag) = 0
d²E_imag/dx² + λ (ε' E_imag − ε'' E_real) = 0
```

This is the generalised eigenvalue problem **A x = λ B x** where:
- **A**: stiffness matrix (from Diffusion kernels)
- **B**: complex ε-weighted mass matrix (from MatReaction + CoupledForce)
- **x** = [E_real; E_imag] stacked DOF vector

## MOOSE Formulation

### Variables

| Variable | Type | Description |
|----------|------|-------------|
| `E_real` | FIRST LAGRANGE | Real part of electric field phasor |
| `E_imag` | FIRST LAGRANGE | Imaginary part of electric field phasor |

### Kernels

| Kernel | Variable | Tag | Role |
|--------|----------|-----|------|
| `Diffusion` | E_real | (none) | Stiffness A: −d²E_real/dx² |
| `Diffusion` | E_imag | (none) | Stiffness A: −d²E_imag/dx² |
| `MatReaction(rate=ε')` | E_real | eigen | Mass B: −ε' ∫ E_real φ dx |
| `CoupledForce(v=E_imag, coef=+ε'')` | E_real | eigen | Mass B: −ε'' ∫ E_imag φ dx |
| `MatReaction(rate=ε')` | E_imag | eigen | Mass B: −ε' ∫ E_imag φ dx |
| `CoupledForce(v=E_real, coef=−ε'')` | E_imag | eigen | Mass B: +ε'' ∫ E_real φ dx |

**Note**: `MatReaction` (non-AD) is used because `GenericConstantMaterial` provides non-AD properties. Using `ADMatReaction` with non-AD material properties causes a runtime error.

### Boundary Conditions

Both PEC walls require two BC objects each:

| BC | Variable | Boundary | System |
|----|----------|----------|--------|
| `DirichletBC` (value=0) | E_real, E_imag | pec_left, pec_right | Stiffness A |
| `EigenDirichletBC` | E_real, E_imag | pec_left, pec_right | Mass B |

`EigenDirichletBC` is essential — without it the boundary nodes contribute spurious entries to the B matrix, producing near-zero eigenvalues that contaminate the spectrum.

## Expected Results

### Analytic Eigenvalues (lossless limit, ε'' = 0)

```
λ_m = (mπ/D)²/ε' = m² × π²/4 ≈ m² × 2.4674
```

| Mode m | λ_m (lossless) | Mode shape |
|--------|----------------|------------|
| 1 | 2.467 | sin(πx/D) |
| 2 | 9.870 | sin(2πx/D) |
| 3 | 22.21 | sin(3πx/D) |
| 4 | 39.48 | sin(4πx/D) |
| 5 | 61.68 | sin(5πx/D) |
| 6 | 88.83 | sin(6πx/D) |

### With Loss (ε'' = 0.2, perturbation theory)

```
Re(λ_m) ≈ m² × 2.4674    (resonant wavenumber squared, barely shifted by loss)
Im(λ_m) ≈ −m² × 0.1234   (negative: energy extracted from field by absorption)
Q_m ≈ ε'/(2ε'') = 10.0   (same for all modes in perturbation limit)
```

## Running the Simulation

```bash
# Inside Docker with combined-opt:
cd /path/to/case84-lossy-tem-cavity
/opt/moose/bin/combined-opt -i case84_lossy_tem_cavity.i

# Expected outputs:
#   case84_lossy_tem_cavity.e    — Exodus field file (E_real, E_imag, |E|² vs x)
#   case84_lossy_tem_cavity.csv  — Postprocessor values
#   case84_lossy_tem_cavity_eigenvalues_*.csv — Complex eigenvalues table
```

## Verification

Post-process the eigenvalues CSV with Python:

```python
import numpy as np
import pandas as pd

# Load eigenvalues (MOOSE reports real and imaginary parts)
df = pd.read_csv('case84_lossy_tem_cavity_eigenvalues_0001.csv')

# Compute Q factor for each mode
lambda_real = df['eigen_values_real'].values
lambda_imag = df['eigen_values_imag'].values
Q = np.abs(lambda_real) / (2.0 * np.abs(lambda_imag))

print("Mode  Re(lambda)   Im(lambda)   Q factor")
for m, (lr, li, q) in enumerate(zip(lambda_real, lambda_imag, Q), 1):
    print(f"  {m}   {lr:10.4f}   {li:10.4f}   {q:.3f}")

# Expected: Q ≈ 10.0 for all modes
# Expected: Re(lambda_m) ≈ m^2 * 2.4674
```

## Key Concepts

1. **Lossy eigenvalue problem**: Complex ε_r → complex eigenvalues λ. The imaginary part of λ directly encodes the energy loss rate per oscillation cycle.

2. **Q factor from eigenvalues**: Q = |Re(λ)| / (2|Im(λ)|) is derived from the complex eigenvalue without any driven-source simulation. This is the analytic continuation of the resonator quality factor.

3. **Real/imaginary splitting**: The 2N×2N real eigenproblem is mathematically equivalent to the N×N complex eigenproblem but avoids complex arithmetic in MOOSE's standard FEM assembly.

4. **Off-diagonal ε'' coupling**: The imaginary part of ε_r creates off-diagonal blocks in the B matrix coupling E_real ↔ E_imag. This is the FEM encoding of material loss in the eigenvalue problem.

## Relation to Other Cases

- **Case 30** (waveguide cutoff): 2D eigenvalue with real ε_r; diagonal B matrix.
- **Case 81** (photonic crystal): 2D eigenvalue with spatially-varying real ε_r.
- **Case 82** (3D cavity): 3D eigenvalue with uniform real ε (vacuum).
- **Case 83** (Veselago lens): Driven source problem with complex ε_r; real cross-coupling between E_real and E_imag via ADMatCoupledForce in the bulk.
- **Case 84** (this case): 1D eigenvalue with uniform complex ε_r; demonstrates Q-factor extraction from complex eigenspectrum.

## References

- Pozar, D.M., *Microwave Engineering*, 4th ed., Wiley (2012), Sections 6.7–6.8
- Collin, R.E., *Foundations for Microwave Engineering*, 2nd ed., McGraw-Hill (1992), Chapter 7
- Kong, J.A., *Electromagnetic Wave Theory*, 2nd ed., EMW Publishing (1990), Chapter 1
- MIT 6.635 Advanced Electromagnetism, Spring 2003 — Prof. Jin Au Kong, Lecture 1
  https://ocw.mit.edu/courses/6-635-advanced-electromagnetism-spring-2003/
