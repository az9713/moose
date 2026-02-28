# Case 89: Dielectric Slab Waveguide — TE Guided Mode Eigenvalues

**MIT 6.635 Lecture(s):** 3-5 (Guided waves, dielectric waveguides, dispersion relations,
eigenmode analysis)

## Overview

This case computes the propagation constants (eigenvalues beta^2) of the guided TE modes in a
symmetric planar dielectric slab waveguide. The core has refractive index n1 = 1.5 and thickness
d = 1.0 m; the cladding has n2 = 1.0. With free-space wavelength lambda_0 = 1 m (k0 = 2 pi
rad/m) the V-number is V = k0 * (d/2) * sqrt(n1^2 - n2^2) = pi * sqrt(1.25) ~ 3.51, which places
the waveguide in the three-guided-mode regime (V > pi, supporting TE0, TE1, and TE2).

The transverse field profile psi(x) satisfies a 1D Sturm-Liouville problem. MOOSE reformulates
it as the generalised eigenvalue problem A psi = lambda B psi with lambda = beta^2, using the
same kernel pattern as Case 81 (photonic crystal): Diffusion builds the stiffness part of A,
MatReaction with rate = n^2(x) k0^2 (non-AD, from GenericFunctionMaterial) adds the weighted
mass contribution to A, and CoefReaction with coefficient = -1 tagged extra_vector_tags = 'eigen'
builds the negative identity B matrix. SLEPc's Krylov-Schur solver with shift-invert near the
cladding cut-off (eps_target = 39.5 ~ n2^2 k0^2) finds all three guided modes efficiently.

This case extends Case 82 (3D cavity resonance eigenvalue) to a continuous-medium waveguide
problem with spatially varying material coefficients, and demonstrates how a non-uniform
refractive index creates a guidance mechanism through the eigenvalue spectrum.

## The Physics

Governing transverse equation for TE modes (E_y component, propagation in z):

    d^2 psi/dx^2 + (n^2(x) k0^2 - beta^2) psi = 0

where n(x) = n1 for |x| < d/2 (core) and n(x) = n2 for |x| > d/2 (cladding).

Rearranged as a generalised eigenvalue problem with eigenvalue lambda = beta^2:

    d^2 psi/dx^2 + n^2(x) k0^2 psi = lambda psi

Weak form (multiply by test phi, integrate by parts, psi=0 at domain boundaries):

    A psi = lambda B psi
    A_ij = integral( dphi_j/dx * dphi_i/dx dx ) - integral( n^2(x) k0^2 phi_j phi_i dx )
    B_ij = -integral( phi_j phi_i dx )

KERNEL MAPPING:
- Diffusion (standard):          contributes +integral(dpsi/dx dphi/dx) to A
- MatReaction(n2k0sq, standard): contributes -n^2 k0^2 integral(psi phi) to A (non-eigen-tagged)
- CoefReaction(-1, eigen-tagged): contributes -integral(psi phi) to B

Refractive index profile (step-index waveguide):

    n(x) = 1.5  for |x| < 0.5 m  (core)
    n(x) = 1.0  for |x| >= 0.5 m (cladding)

Material parameters:
- n1^2 k0^2 = 2.25 * (2pi)^2 = 88.826 m^-2  (core coefficient)
- n2^2 k0^2 = 1.00 * (2pi)^2 = 39.478 m^-2  (cladding coefficient)

Guidance condition: n2^2 k0^2 < beta^2 < n1^2 k0^2, i.e., 39.478 < lambda < 88.826

V-number analysis (V = k0 * (d/2) * sqrt(n1^2 - n2^2) = pi * sqrt(1.25) ~ 3.512):
- TE0 exists for all V > 0: V = 3.512 > 0 -> guided
- TE1 cut-off at V = pi/2 ~ 1.571: V = 3.512 > 1.571 -> guided
- TE2 cut-off at V = pi ~ 3.141: V = 3.512 > 3.141 -> guided
- TE3 cut-off at V = 3pi/2 ~ 4.712: V = 3.512 < 4.712 -> not guided

Expected eigenvalues (approximate, from graphical/numerical solutions of characteristic equations):
- TE0 (fundamental, even, no interior nodes):   lambda ~ 80  (n_eff = beta/k0 ~ 1.42)
- TE1 (first order, odd, one interior node):    lambda ~ 60  (n_eff ~ 1.23)
- TE2 (second order, even, one interior node):  lambda ~ 43  (n_eff ~ 1.04, near cut-off)

Domain: 1D, [-5, 5] m (evanescent field decays to < 2e-4 of peak at x = +/-5 m).
Mesh: nx = 400, dx = 0.025 m. Boundary conditions: DirichletBC + EigenDirichletBC at x = +/-5.

## MOOSE Objects Used

- **Kernels**: `Diffusion` (stiffness, A matrix); `MatReaction` with reaction_rate = n2k0sq
  (n^2 k0^2 weighted mass, A matrix, NOT eigen-tagged); `CoefReaction` with coefficient = -1
  and extra_vector_tags = 'eigen' (negative identity B matrix for the eigenproblem)
- **Materials**: `GenericFunctionMaterial` (NOT ADGenericFunctionMaterial) wrapping n2k0sq_fn
  as non-AD property 'n2k0sq' — required because MatReaction (non-AD) reads non-AD properties.
  The same AD/non-AD pairing rule as Case 81: GenericFunctionMaterial + MatReaction (non-AD)
- **AuxVariables / AuxKernels**: `FunctionAux` populates n_sq_field (CONSTANT MONOMIAL) from
  eps_r_vis_fn for visualisation of the refractive-index step profile; `PotentialToFieldAux`
  computes dpsi_dx = d psi/dx (transverse derivative, proportional to H_z for TE modes)
- **BCs**: `DirichletBC` psi=0 on 'left right' (A matrix BC) and `EigenDirichletBC` psi=0
  on 'left right' (B matrix BC) — both required for the eigenvalue problem, as in Case 82
- **VectorPostprocessors**: `Eigenvalues` reports all six computed lambda = beta^2 values
- **Postprocessors**: `PointValue` at core center (x=0), interface (x=0.5), core mid (x=0.25),
  near cladding (x=1), far cladding (x=3); `ElementExtremeValue` for psi_max and psi_min
  (mode symmetry check); `ElementAverageValue` for average n^2 (sanity check)
- **Executioner**: Eigenvalue type, KRYLOVSCHUR, n_eigen_pairs=6, eps_target=39.5 (just above
  n2^2 k0^2 to target guided modes), st_type=sinvert with MUMPS LU inner solver

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case89-dielectric-waveguide \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case89_dielectric_waveguide.i 2>&1 | tail -30'
```

## Expected Results

The VectorPostprocessor CSV file lists 6 eigenvalues. Three should fall in the guidance range
39.478 < lambda < 88.826, corresponding to TE0, TE1, and TE2 guided modes. The remaining
eigenvalues fall below 39.478 (radiation modes, physically leaky but captured as pseudo-modes
on the finite computational domain).

Approximate expected lambda values:
- lambda_TE2 ~ 43  (closest to shift target 39.5, converges first in shift-invert)
- lambda_TE1 ~ 60
- lambda_TE0 ~ 80  (fundamental mode, farthest from cut-off)

Postprocessor interpretation for the fundamental mode TE0 (even symmetry):
- psi_at_center (x=0): maximum amplitude (cosine peak in core)
- psi_at_interface (x=0.5): amplitude reduced from center value (continuous, but slope reverses)
- psi_at_cladding_near (x=1): much smaller, exponential decay well underway
- psi_at_cladding_far (x=3): near zero, confirming guided (bound) mode character
- psi_min / psi_max: near zero ratio -> confirms even-symmetry mode

In ParaView, the Exodus output shows the 1D mode profile psi(x). Using the "Plot Over Line"
filter along x produces the characteristic shape:
- TE0: smooth bell-shaped curve concentrated in the core, exponential tails in the cladding
- TE1: antisymmetric with a zero crossing at x=0 (after cycling to mode 2 in SLEPc output)
- TE2: symmetric with two extrema in the core and a node inside the core

The n_sq_field overlay shows the step from n^2=1.0 (cladding) to n^2=2.25 (core) at x=+/-0.5,
visually confirming the core location and width.

## Key Takeaways

- The 1D dielectric waveguide eigenvalue problem maps directly to MOOSE's A psi = lambda B psi
  framework: Diffusion (stiffness, A), MatReaction with n^2 k0^2 (weighted mass, A), and
  CoefReaction(-1, eigen-tagged) for the negative identity B matrix.
- The V-number determines how many guided modes exist. For V > m*pi/2 (m = 0,1,2,...), mode
  TE_m is guided. This waveguide (V ~ 3.51) supports exactly three guided modes: TE0, TE1, TE2.
- Use GenericFunctionMaterial (non-AD) + MatReaction (non-AD) for eigenvalue problems, not
  ADMatReaction. This is the same pattern as Case 81. Mixing AD and non-AD types causes a
  runtime error.
- EigenDirichletBC is essential alongside DirichletBC to enforce zero boundary conditions in
  BOTH the A and B matrices. Without it, spurious eigenvalues appear near the domain boundaries.
- Shift-invert (eps_target just above n2^2 k0^2) is the key to efficiently finding only the
  guided modes, avoiding the dense continuum of radiation modes below the cut-off.
- The effective refractive index n_eff = sqrt(lambda)/k0 = sqrt(lambda)/(2pi) always satisfies
  n2 < n_eff < n1 for guided modes — this is the guidance condition translated to eigenvalue form.
- This 1D eigenvalue problem is the transverse counterpart of Case 30 (2D waveguide cutoff): both
  use the same Diffusion + CoefReaction(eigen) kernel structure, but Case 89 adds the n^2(x)
  weighted MatReaction term to incorporate the spatially varying refractive index.
