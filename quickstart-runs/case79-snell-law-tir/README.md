# Case 79: Oblique Incidence — Snell's Law and Total Internal Reflection

**MIT 6.635 Lecture(s):** 3-4 (Green's functions for layered media, TE/TM decomposition, Fresnel coefficients)

## Overview

This case demonstrates Snell's law and total internal reflection (TIR) for a TE-polarised plane
wave striking a planar dielectric interface at oblique incidence. The wave propagates in medium 1
(epsilon_r = 4, n1 = 2) at theta = 20 degrees and refracts into medium 2 (epsilon_r = 1, n2 = 1,
free space) at theta2 = arcsin(2 sin 20 deg) = 43.16 degrees. The critical angle for TIR is
theta_c = arcsin(n2/n1) = 30 degrees, so this case sits below TIR. The Fresnel power reflectance
is R = 0.194 and transmittance is T = 0.806, with R + T = 1 confirming energy conservation.

The key new technique is the **scattered-field formulation**: instead of solving for the total
field directly, the problem is decomposed as E_total = E_inc + E_scat and only E_scat is solved
for. Because E_inc satisfies the Helmholtz equation in medium 1 exactly, substituting the
decomposition produces an interior source term -k0^2 * (eps_r - 4) * E_inc that is nonzero only
in medium 2. All four domain boundaries then use purely absorbing BCs with no wave injection.

This case introduces oblique incidence, the scattered-field decomposition, phase-matching at
dielectric interfaces, and the use of ParsedAux with use_xyzt = true to reconstruct the total
field analytically.

## The Physics

Governing equation for the scattered field E_scat (TE polarisation, 2D):

    nabla^2(E_scat) + k0^2 * eps_r(x) * E_scat = -3 k0^2 * E_inc(x,y)   (x > 0, medium 2)
                                                 = 0                        (x < 0, medium 1)

where E_inc = exp(j*(kx*x + ky*y)), with kx = k0 n1 cos(20 deg) = 2.9521 rad/m and
ky = k0 n1 sin(20 deg) = 1.0745 rad/m. Free-space wavenumber k0 = pi/2 rad/m (lambda_0 = 4 m).

Fresnel coefficients for TE at theta1 = 20 degrees:
- r_s = (n1 cos theta1 - n2 cos theta2)/(n1 cos theta1 + n2 cos theta2) = 0.4408
- t_s = 2 n1 cos theta1 / (n1 cos theta1 + n2 cos theta2) = 1.4408
- Power reflectance R = |r_s|^2 = 0.194, transmittance T = 0.806

Boundary conditions:
- All four domain boundaries: EMRobinBC in absorbing mode (profile_func_real = 0), first-order
  Sommerfeld radiation condition absorbing outgoing scattered waves

Material properties:
- Medium 1 (x < 0): epsilon_r = 4, n1 = 2
- Medium 2 (x > 0): epsilon_r = 1, n2 = 1

Domain: 2D, [-5, 5] x [-5, 5] m. Dielectric interface at x = 0. Mesh: 100 x 100 QUAD4 elements
(20 elements per wavelength in medium 1, 40 per wavelength in medium 2).

## MOOSE Objects Used

- **Kernels**: `Diffusion` (Laplacian term), `ADMatReaction` (k0^2 eps_r(x) coefficient from
  ADGenericFunctionMaterial), `BodyForce` (scattered-field source term, nonzero only in medium 2)
- **BCs**: `EMRobinBC` — eight absorbing ABC blocks (two per boundary, real and imaginary
  components), all with profile_func_real = 0
- **Materials**: `ADGenericFunctionMaterial` wrapping the piecewise-constant coeff_fn as AD
  property `coeff_material` for use with ADMatReaction
- **AuxVariables / AuxKernels**: `ParsedAux` with use_xyzt = true reconstructs E_total = E_scat
  + E_inc by adding the analytical incident field back; further ParsedAux computes |E_total|^2
  and |E_scat|^2
- **Postprocessors**: `PointValue` at the interface, in medium 1, in medium 2, and at three
  y-offsets along x = 2 to verify phase matching (Snell's law)
- **Executioner**: Steady Newton with LU direct solver; SMP full = true for the off-diagonal
  coupling between E_scat_real and E_scat_imag introduced by EMRobinBC

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case79-snell-law-tir \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case79_snell_law_tir.i 2>&1 | tail -30'
```

## Expected Results

Newton converges in 1-2 iterations. Key results for Fresnel verification:

- Total field amplitude at x = 0 (the interface) should be approximately t_s = 1.44 for the
  transmitted wave component.
- In medium 1 (x < 0): total field is the superposition of incident and reflected waves, forming
  a standing-wave pattern with intensity oscillating between (1 - |r_s|)^2 = 0.31 and
  (1 + |r_s|)^2 = 2.07.
- In medium 2 (x > 0): single traveling wave propagating at theta2 = 43.16 degrees, amplitude ~1.44.
- Phase-matching check: probing E_total_real at (2,0), (2,1), (2,-1) verifies the transverse
  wavenumber ky = 1.0745 rad/m is preserved across the interface — the FEM statement of Snell's
  law.

The Exodus output shows E_total_sq in ParaView, revealing the oblique incident beam in medium 1
and the more steeply refracted beam in medium 2, with a clear phase-front rotation at x = 0.

## Key Takeaways

- The scattered-field formulation converts oblique incidence into an interior BodyForce source,
  enabling purely absorbing BCs on all domain boundaries.
- The source term -k0^2 * (eps_r - eps_inc) * E_inc is nonzero only where the background
  permittivity differs from the incident medium — here in medium 2 only.
- Snell's law is enforced implicitly by the C^0 Lagrange continuity at x = 0: the tangential
  wavenumber ky is conserved because the y-dependence of E_inc = exp(j ky y) propagates through.
- The MOOSE BodyForce residual is -integral(f v), so the source coefficient must carry a
  negative sign: source_coeff = -3 k0^2 (not +3 k0^2).
- ParsedAux with use_xyzt = true allows spatially varying analytical expressions to be evaluated
  at each node for total-field reconstruction without an extra FEM solve.
