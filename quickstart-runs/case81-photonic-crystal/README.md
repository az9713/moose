# Case 81: 2D Photonic Crystal Band Gap — TM Eigenvalue at Gamma Point

**MIT 6.635 Lecture(s):** 8-9 (Photonic crystals, Bloch modes, band gaps)

## Overview

This case computes the TM electromagnetic eigenfrequencies of a 2D photonic crystal at the Gamma
point (k = 0) of the first Brillouin zone. The structure is a square lattice of circular alumina
rods (epsilon_r = 8.9, radius r = 0.2a) in vacuum (epsilon_r = 1), with normalised lattice
constant a = 1. Just as electrons in a crystal develop energy band gaps from the periodic ionic
potential, electromagnetic modes in a periodic dielectric structure develop photonic band gaps —
frequency ranges where no propagating modes exist regardless of direction.

The problem is cast as a generalised eigenvalue problem A x = lambda B x, where
lambda = (omega/c)^2 is the squared normalised frequency. MOOSE's SLEPc interface with the
Krylov-Schur solver and shift-invert spectral transformation computes the first several TM
eigenvalues simultaneously. At the Gamma point, Bloch's theorem reduces to ordinary periodicity,
and Neumann (PMC) boundary conditions on the unit cell faces select the symmetric modes without
requiring an explicit BCs block.

This case extends Case 30 (2D waveguide cutoff eigenvalue) to a periodic unit cell and
introduces the photonic crystal eigenvalue formulation, the eps_r-weighted mass matrix, and
the GenericFunctionMaterial / MatReaction (non-AD) combination required for eigenvalue problems.

## The Physics

Governing equation (TM polarisation, E_z = psi, 2D unit cell):

    nabla^2(psi) + (omega/c)^2 * eps_r(x,y) * psi = 0

Weak form — generalised eigenvalue problem A psi = lambda B psi:

    A_ij = integral( nabla phi_j . nabla phi_i dV )     (stiffness, Diffusion kernel)
    B_ij = -integral( eps_r * phi_j * phi_i dV )        (weighted mass, MatReaction with 'eigen' tag)
    lambda = (omega/c)^2                                 (eigenvalue = squared normalised frequency)

Dielectric function on the unit cell [0,1] x [0,1]:
    eps_r(x,y) = 8.9  if (x-0.5)^2 + (y-0.5)^2 < 0.04  (alumina rod, r = 0.2)
               = 1.0  otherwise                           (vacuum host)
    Filling fraction f = pi r^2/a^2 = pi * 0.04 ~ 0.1257

Expected eigenvalues — normalised frequency f_a = sqrt(lambda)/(2 pi), a = 1:
- Mode 1 (trivial DC): lambda = 0
- Modes 2-3 (degenerate pair): lambda ~ 5.1, f_a ~ 0.36
- Mode 4: lambda ~ 11.1, f_a ~ 0.53

Boundary conditions: Neumann natural BC (d psi/dn = 0) on all four unit cell faces — no
explicit [BCs] block needed. Selects symmetric Gamma-point modes (PMC equivalent).

Domain: unit cell [0,1] x [0,1], a = 1 (normalised). Mesh: 40 x 40 QUAD4 elements.

## MOOSE Objects Used

- **Kernels**: `Diffusion` (stiffness matrix A); `MatReaction` with reaction_rate = eps_r and
  extra_vector_tags = 'eigen' (eps_r-weighted mass matrix B routed to SLEPc B matrix)
- **Materials**: `GenericFunctionMaterial` (NOT ADGenericFunctionMaterial) wrapping eps_r_fn as
  non-AD property eps_r — required because MatReaction (non-AD) reads non-AD properties
- **AuxVariables / AuxKernels**: `FunctionAux` populates eps_r_field (CONSTANT MONOMIAL) for
  visualisation; `PotentialToFieldAux` computes grad_x = -d(psi)/dx and grad_y = -d(psi)/dy
- **VectorPostprocessors**: `Eigenvalues` reports all computed lambda values as a CSV column
- **Postprocessors**: `ElementExtremeValue` tracks psi_max and psi_min; `ElementAverageValue`
  verifies avg(eps_r) ~ 2.00 (filling-fraction sanity check)
- **Executioner**: Eigenvalue type, KRYLOVSCHUR solve_type, n_eigen_pairs = 8,
  eps_target = 4.0 (just below the first non-trivial mode), st_type = sinvert with LU inner solver

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case81-photonic-crystal \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case81_photonic_crystal.i 2>&1 | tail -30'
```

## Expected Results

The VectorPostprocessor CSV file lists 8 eigenvalues. The first entry near zero is the trivial
constant mode (psi = const, lambda = 0). The next eigenvalues should appear near lambda = 5.1
(f_a ~ 0.36, a doubly degenerate pair corresponding to the two lowest non-trivial TM bands at
Gamma), then 11.1 (f_a ~ 0.53), consistent with published MPB band structure results for
alumina rods in vacuum.

The Exodus output shows psi, eps_r_field, grad_x, and grad_y. In ParaView, colouring by
eps_r_field shows the rod geometry; colouring by psi reveals the mode shape. The lowest
non-trivial mode concentrates inside the high-eps rod (dielectric band), while higher modes
concentrate in the surrounding air (air band) — the physical origin of the photonic band gap.

## Key Takeaways

- Photonic crystal modes are generalised eigenvalue problems: Diffusion builds A, MatReaction
  with extra_vector_tags = 'eigen' builds the eps_r-weighted B matrix.
- Use GenericFunctionMaterial (non-AD) paired with MatReaction (non-AD) — not ADMatReaction.
  Mixing AD properties with non-AD kernels or vice versa causes a runtime error.
- Neumann BCs (the FEM natural BC) require no explicit block and select the symmetric Gamma-point
  modes; the trivial lambda = 0 constant mode is handled automatically by shift-invert.
- Shift-invert with eps_target just below the first physical mode converges efficiently to the
  desired band frequencies without finding the trivial DC mode first.
- Convert eigenvalue to normalised frequency: f_a = sqrt(lambda) / (2 pi) with a = 1 (normalised).
- Modes below a band gap concentrate in high-eps regions (dielectric bands); modes above
  concentrate in low-eps regions (air bands) — this spatial separation creates the gap.
