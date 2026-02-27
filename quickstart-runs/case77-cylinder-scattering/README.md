# Case 77: EM Scattering from a Dielectric Cylinder (2D, TE Polarisation)

## Overview

This case solves the 2D frequency-domain electromagnetic scattering problem for a plane wave
striking a dielectric cylinder with epsilon_r = 4 (n = 2) and radius R = 1 m. The problem
is formulated in TE polarisation (E-field along the cylinder axis, z-direction), which reduces
the full 3D Maxwell problem to a scalar 2D Helmholtz equation for E_z in the cross-sectional
plane. This is the first 2D scattering problem in the series, extending the 1D slab cases (32,
74, 75) to two spatial dimensions.

The case demonstrates Mie scattering physics: for an electrical size k0*R = pi/2 ~ 1.57 (Mie
regime, neither Rayleigh nor geometric optics), the scattered field shows pronounced diffraction,
field enhancement inside the cylinder due to focusing by the high-index material, and a shadow
region behind the cylinder. The FEM approach here is compared conceptually with the Electric
Field Integral Equation (EFIE) and Method of Moments (MoM) taught in MIT 6.635 Lectures 5-6.

The case introduces the 2D Helmholtz formulation, the use of four EMRobinBC absorbing boundary
conditions on all domain faces (one in port mode to inject the incident wave from the west, three
in absorbing mode), and the ParsedAux E_magnitude_sq = E_real^2 + E_imag^2 for intensity output.

## The Physics

Governing equation (2D scalar Helmholtz, TE polarisation):

    nabla^2(E_z) + k0^2 * eps_r(x,y) * E_z = 0

where:

    eps_r(x,y) = 4.0  if x^2 + y^2 < 1  (dielectric cylinder)
               = 1.0  otherwise           (free space)

    k0 = 2*pi/lambda_0 = pi/2 = 1.5708 rad/m  (lambda_0 = 4 m)
    k0^2 = 2.4674 m^-2  (free space),  k0^2 * eps_r = 9.8696 m^-2  (inside cylinder)

Electrical size: k0 * R = pi/2 ~ 1.57 (Mie regime, ka ~ O(1)).

Boundary conditions:
- West boundary (x = -4): EMRobinBC port mode, injects E_inc = exp(j*k0*x)
- East, north, south boundaries: EMRobinBC absorbing mode, implements first-order Sommerfeld ABC

Domain: 2D, [-4, 4] x [-4, 4] m (2 free-space wavelengths wide).
Mesh: 80 x 80 = 6400 quad elements, 40 elements per lambda_0 (well-resolved).

## Input File Walkthrough

**HIT variables**: k = pi/2 = 1.5708 rad/m, E0 = 1 V/m, theta = 0.

**[Mesh]**: 2D GeneratedMesh from [-4,4] x [-4,4] with 80x80 elements. Boundaries renamed
to `west`, `east`, `south`, `north`.

**[Variables]**: E_real and E_imag, both FIRST LAGRANGE. Lagrange elements enforce C^0
continuity of E_z across the cylinder boundary, implementing tangential E continuity.

**[Functions]**: `coeff_fn` is a ParsedFunction returning k0^2*eps_r(x,y): 9.8696 inside
the cylinder (x^2+y^2 < 1), 2.4674 outside. `cosTheta` = cos(0) = 1.

**[Materials]**: ADGenericFunctionMaterial wraps coeff_fn as AD property `coeff_material`.

**[Kernels]**: Standard lossless Helmholtz: Diffusion + ADMatReaction(coeff_material) for
each of E_real and E_imag. No cross-coupling kernels (eps_r is real, lossless cylinder).
The two field components are decoupled in the bulk and couple only through the port BC.

**[BCs]**: Eight EMRobinBC blocks (two per boundary face, real and imaginary components).
The west boundary uses `mode = port` with `profile_func_real = E0 = 1` to inject the
incident wave. The east, north, south boundaries use `mode = absorbing` (profile_func_real
omitted) to absorb outgoing scattered field.

**[AuxVariables / AuxKernels]**: `E_magnitude_sq` stores |E_z|^2 = E_real^2 + E_imag^2
computed via ParsedAux, proportional to time-averaged intensity.

**[Postprocessors]**: PointValue samples at cylinder centre (0,0), shadow region (2.5,0),
forward scatter (3.5,0), incident side (-2.5,0), broadside (0,2), and far corner (3.5,3.5).
`total_field_energy` integrates |E_z|^2 over the domain.

**[Executioner]**: Steady Newton with LU. The 2D Helmholtz system is indefinite; LU is
robust at this scale (12,800 DOFs).

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case77-cylinder-scattering \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case77_cylinder_scattering.i 2>&1 | tail -30'
```

## Expected Results

Newton converges in 1-2 iterations. The 2D Exodus output `case77_cylinder_scattering_out.e`
shows the E_magnitude_sq field in ParaView. Key physical features to observe:
- Inside the cylinder (|r| < 1): field enhancement from dielectric focusing (|E|^2 > 1).
- Shadow region (x > 1 on axis): reduced field amplitude below 1.
- Forward scatter direction: constructive interference produces a bright lobe.
- Standing-wave pattern on the incident side (x < -1): interference between incident and
  backscattered waves creates alternating bright/dark bands with period lambda_0/2 = 2 m.

The FEM result can be compared to the exact Mie series for the TE scattering cross-section.

## Key Takeaways

- 2D Helmholtz scattering is solved with the same Diffusion + ADMatReaction kernel pair as 1D.
- The ParsedFunction `if(x*x+y*y < 1, ...)` captures the circular cylinder geometry without
  sub-domain meshing; C^0 Lagrange elements enforce E_z continuity at the interface.
- EMRobinBC in port mode injects a plane wave; absorbing mode acts as a first-order ABC.
- The FEM (volume, sparse) and MoM/EFIE (surface, dense) approaches give the same physics
  but differ in efficiency: MoM excels for homogeneous scatterers, FEM for inhomogeneous ones.
- Field intensity |E_z|^2 = E_real^2 + E_imag^2 is efficiently computed via ParsedAux.
- Relates to MIT 6.635 Lectures 5-6: EFIE, MoM, and Mie scattering from dielectric cylinders.
