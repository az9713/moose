# Case 74: EM Wave in a Left-Handed Material (LHM) Slab

## Overview

This case models a plane electromagnetic wave at f = 20 MHz propagating through a slab of
left-handed material (LHM) — a medium in which both the relative permittivity (epsilon_r = -1.5)
and relative permeability (mu_r = -1.5) are simultaneously negative. Such double-negative media
were first predicted theoretically by Veselago (1968) and later realised experimentally as
metamaterials using split-ring resonator arrays (Smith et al., Science 292, 2001).

The key physics that sets an LHM apart from an ordinary dielectric is the sign of the refractive
index: n = -sqrt(epsilon_r * mu_r) = -1.5. While the energy (group) velocity remains forward
(into the slab), the phase velocity is reversed. This reversal is not visible in the bulk
Helmholtz equation — n^2 = 2.25 is identical for n = +1.5 and n = -1.5 — but it manifests
through the interface condition at the vacuum-LHM boundary, where the sign of d(E)/dx flips.

This case introduces the 1/mu_r weighted diffusion formulation (ADMatDiffusion with a negative
diffusivity), which is the only FEM approach capable of distinguishing negative-index from
positive-index media. The case extends Case 32 (single dielectric slab) by replacing the
positive-epsilon material with a double-negative LHM, and by using ADMatDiffusion in place of
plain Diffusion. It relates to MIT 6.635 Lectures 1-2 on left-handed materials and negative
refraction.

## The Physics

Governing equation (frequency domain, 1D, normal incidence):

    d/dx [ (1/mu_r) dE/dx ] + k0^2 * eps_r(x) * E = 0

where k0 = 2*pi*f/c = 0.41888 rad/m. This is NOT the same as the simple scalar Helmholtz
d^2E/dx^2 + k0^2*n^2*E = 0 because the interface condition carries the physics of mu_r:

    Continuity: E continuous, (1/mu_r) dE/dx continuous

At the vacuum-to-LHM interface this gives dE/dx|LHM = -1.5 * dE/dx|vac (sign reversal).

Boundary conditions:
- x = 0: PEC wall (DirichletBC, E_real = E_imag = 0)
- x = 120 m: Robin port (EMRobinBC) injecting E_inc = exp(j*k0*x) and absorbing reflections

Material properties:
- Vacuum (x < 40 m, x > 80 m): 1/mu_r = +1.0, k0^2*eps_r = +0.17546 m^-2
- LHM slab (40 <= x <= 80 m): 1/mu_r = -0.66667, k0^2*eps_r = -0.26320 m^-2

Domain: 1D, x in [0, 120] m. Domain is 8 free-space wavelengths (lambda_0 = 15 m at 20 MHz).
The LHM slab thickness is 40 m = 2.67 lambda_0. Mesh: 600 first-order Lagrange elements.

## Input File Walkthrough

**HIT variables**: k0 = 0.41888 rad/m, L = 120 m, E0 = 1 V/m, theta = 0 (normal incidence).

**[Mesh]**: 1D GeneratedMesh from x=0 to x=120 with 600 elements. Boundaries renamed
to `metal` (x=0, PEC) and `port` (x=120, injection port).

**[Variables]**: Two real-valued scalars `E_real` and `E_imag` representing the complex
phasor E(x) = E_real + j*E_imag. Both use FIRST-order LAGRANGE elements.

**[Functions]**: Three ParsedFunctions. `inv_mu_fn` returns 1/mu_r(x): -0.66667 in the LHM
slab, +1.0 in vacuum. `k0sq_eps_fn` returns k0^2*eps_r(x): -0.26320 in the slab, +0.17546
in vacuum. `cosTheta` returns cos(0) = 1.

**[Materials]**: ADGenericFunctionMaterial wraps each ParsedFunction as an AD material
property, enabling automatic differentiation through the spatially varying coefficients.

**[Kernels]**: Two identical sets of kernels — one for E_real, one for E_imag. Each set
uses ADMatDiffusion (diffusivity = inv_mu_r) and ADMatReaction (reaction_rate = k0sq_eps_r).
No cross-coupling kernels appear because the LHM is lossless (Im(eps_r) = Im(mu_r) = 0).

**[BCs]**: DirichletBC enforces E_real = E_imag = 0 at x=0 (PEC). EMRobinBC at x=120
simultaneously injects the incident plane wave and absorbs the outgoing reflected wave.
The two field components couple only at this boundary through the j*k0 factor.

**[Postprocessors]**: ReflectionCoefficient computes |R|^2 at the port. PointValue
postprocessors sample the field at x = 20, 60, 100, 120 m for diagnostics.

**[Executioner]**: Steady Newton solver with direct LU factorisation. The non-coercive
system (negative diffusivity in LHM) makes iterative solvers unreliable; LU is robust.

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case74-left-handed-material \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case74_left_handed_material.i 2>&1 | tail -30'
```

## Expected Results

The solver converges in 1-2 Newton iterations (the system is linear). The CSV output
reports the reflection coefficient |R|^2 and field values at diagnostic points. Because
the domain ends with a PEC wall and the slab has n = -1.5, the Fabry-Perot pattern
inside the slab differs markedly from a positive-index slab of the same n^2 = 2.25.
The field at x = 20 m shows a standing wave with a phase structure that reveals the
sign reversal at the vacuum-LHM interfaces. The reflection coefficient depends on the
specific slab thickness and differs from the RHM result despite having the same |n|.

Visualise `case74_left_handed_material_out.e` in ParaView: plot E_real and E_imag along
the x-axis and compare with an equivalent positive-index slab to see the phase reversal.

## Key Takeaways

- A medium with eps_r < 0 and mu_r < 0 simultaneously supports propagation with n < 0.
- The bulk Helmholtz equation is insensitive to the sign of n (n^2 is the same); only the
  interface condition (1/mu_r) dE/dn = continuous distinguishes LHM from RHM.
- ADMatDiffusion with a negative diffusivity (1/mu_r = -0.667 in the slab) naturally
  enforces the sign-reversed interface condition through the weak form.
- Direct LU is required when the diffusion coefficient is negative (non-coercive system).
- The ReflectionCoefficient postprocessor extracts |R|^2 from the complex total field.
- Relates to MIT 6.635 Lectures 1-2: left-handed materials and backward-wave propagation.
