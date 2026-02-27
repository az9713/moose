# Case 82: 3D Rectangular Cavity Resonance — TM Mode Eigenvalues

**MIT 6.635 Lecture(s):** 2 (Rectangular waveguides and cavities; resonant mode structure)

## Overview

This case computes the resonant eigenfrequencies and mode shapes of a 3D perfectly conducting
(PEC) rectangular cavity with dimensions a = 3 m, b = 2 m, d = 1 m. Inside the cavity each
Cartesian field component satisfies the scalar Helmholtz eigenvalue equation with homogeneous
Dirichlet (PEC) boundary conditions on all six faces. The exact resonant wavenumbers are given
analytically by k^2(m,n,p) = (m pi/a)^2 + (n pi/b)^2 + (p pi/d)^2 with mode indices m, n, p = 1, 2, 3, ...

The case extends Case 30 (2D rectangular waveguide cutoff modes) to a fully 3D closed cavity.
The 2D QUAD4 mesh becomes a 3D HEX8 mesh, PotentialToFieldAux now computes all three gradient
components (Ex, Ey, Ez), and the HEX8-based FEM discretisation introduces a third mode family
in the z-direction. The mathematical formulation is otherwise identical: Diffusion + CoefReaction
with extra_vector_tags = 'eigen', EigenDirichletBC on all cavity walls, and Krylov-Schur with
shift-invert.

The case relates to Pozar, Microwave Engineering, 4th ed., Section 6.7, and Collin, Foundations
for Microwave Engineering, Chapter 7.

## The Physics

Governing equation (scalar Helmholtz eigenvalue, 3D PEC cavity):

    nabla^2(psi) + k^2 * psi = 0   in Omega = [0,3] x [0,2] x [0,1]
    psi = 0                         on all six PEC walls

Separable solution:

    psi(x,y,z) = sin(m pi x/a) * sin(n pi y/b) * sin(p pi z/d)

Exact resonant eigenvalues for a = 3, b = 2, d = 1:

    k^2(m,n,p) = (m pi/a)^2 + (n pi/b)^2 + (p pi/d)^2

    TM_111: k^2 = 1.0966 + 2.4674 +  9.8696 = 13.434   (fundamental)
    TM_211: k^2 = 4.3864 + 2.4674 +  9.8696 = 16.723
    TM_121: k^2 = 1.0966 + 9.8696 +  9.8696 = 20.836
    TM_311: k^2 = 9.8696 + 2.4674 +  9.8696 = 22.207
    TM_221: k^2 = 4.3864 + 9.8696 +  9.8696 = 24.126

Generalised eigenvalue problem A x = lambda B x, lambda = k^2:
- A_ij = integral( nabla phi_j . nabla phi_i dV )  — stiffness from Diffusion
- B_ij = integral( phi_j * phi_i dV )              — mass from CoefReaction(coefficient=-1, eigen)

Boundary conditions: psi = 0 (Dirichlet) on all six faces (left, right, bottom, top, back, front).

Domain: [0,3] x [0,2] x [0,1] m. Mesh: 15 x 10 x 5 = 750 HEX8 elements.

## MOOSE Objects Used

- **Kernels**: `Diffusion` (stiffness matrix A = integral nabla phi . nabla phi); `CoefReaction`
  with coefficient = -1 and extra_vector_tags = 'eigen' (mass matrix B = integral phi phi,
  the -1 coefficient ensures B is positive-definite as required by SLEPc)
- **BCs**: `DirichletBC` on all six cavity walls enforces psi = 0 in the stiffness system;
  `EigenDirichletBC` on the same boundaries enforces psi = 0 in the eigen (mass) system —
  both are required to prevent spurious near-zero modes from boundary DOFs
- **AuxVariables / AuxKernels**: Ex, Ey, Ez (CONSTANT MONOMIAL) computed by `PotentialToFieldAux`
  with sign = negative, giving E = -nabla(psi) for 3D mode field visualisation
- **VectorPostprocessors**: `Eigenvalues` reporter outputs all 6 computed k^2 values for direct
  comparison with the analytic table
- **Executioner**: Eigenvalue type, KRYLOVSCHUR solve_type, n_eigen_pairs = 6,
  eps_target = 13.0 (just below TM_111 = 13.434), st_type = sinvert with MUMPS inner LU

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case82-3d-cavity-resonance \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case82_3d_cavity_resonance.i 2>&1 | tail -30'
```

## Expected Results

The VectorPostprocessors CSV lists 6 eigenvalues. The first five should approximate the
analytical values 13.434, 16.723, 20.836, 22.207, and 24.126 with FEM discretisation error
O(h^2). On the 15 x 10 x 5 mesh, errors of a few percent are expected; refining to 30 x 20 x 10
reduces errors by roughly a factor of four.

The Exodus output shows the 3D mode shapes in ParaView. For TM_111: psi = sin(pi x/3) sin(pi y/2)
sin(pi z) has a single maximum near the cavity centre (1.5, 1.0, 0.5) and vanishes on all walls.
Using isosurface at psi = 0.5 clearly visualises the half-sinusoid nodal structure.

## Key Takeaways

- 3D cavity resonance is a direct extension of 2D waveguide cutoff: the same Diffusion +
  CoefReaction(eigen) + EigenDirichletBC pattern applies with only a dim = 3 mesh change.
- GeneratedMeshGenerator with dim = 3 produces HEX8 elements; boundary labels include
  `back` and `front` for the z = 0 and z = d faces.
- PotentialToFieldAux works in 3D, extracting all three components of E = -nabla(psi) from
  a single scalar Lagrange variable.
- Both DirichletBC (for matrix A) and EigenDirichletBC (for matrix B) must be applied; omitting
  EigenDirichletBC introduces spurious near-zero eigenvalues from unconstrained boundary DOFs.
- The shift target eps_target = 13.0 must be placed just below the fundamental mode k^2 = 13.434
  to focus the Krylov-Schur solver on physically meaningful cavity modes.
- Analytical resonant frequencies k^2(m,n,p) provide exact verification of the FEM solution
  — a textbook benchmark for any 3D EM eigenvalue solver.
