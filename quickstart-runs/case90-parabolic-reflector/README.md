# Case 90: Parabolic Reflector — Paraxial Focusing of a Plane Wave

**MIT 6.635 Lecture themes:** Scattering, aperture antennas, focusing systems

## Overview

This case simulates a 2D parabolic reflector antenna using the frequency-domain Helmholtz
equation (TE polarisation). A plane wave travelling in the -x direction illuminates a parabolic
conducting surface defined by x = y²/(4f) with focal length f = 2.0. The reflected field
converges to a bright focal spot at the focal point (2, 0).

The parabola has the fundamental geometric property that every ray arriving parallel to the
axis of symmetry reflects through the focal point. Equivalently, all optical path lengths from
a planar wavefront to the focus are identical, so all reflected contributions arrive in phase
at (f, 0): this is the origin of the intense focal spot.

The reflector is implemented as a high-loss dielectric region (epsilon_r'' = 50) embedded in
the rectangular GeneratedMesh. This avoids any need for curved boundary meshing while accurately
capturing the PEC-like reflective behaviour: the skin depth at epsilon_r'' = 50 is only
delta_s = 1/(k0 sqrt(epsilon_r''/2)) ≈ 0.032 lambda_0, so the field is completely absorbed
within the first few mesh elements inside the reflector layer.

The case introduces the lossy-reflector approach for curved conductors on structured meshes —
the same technique used in Cases 83 (Veselago lens) and 90-93 for embedding curved material
boundaries in rectangular domains.

## The Physics

Governing equation (2D scalar Helmholtz, TE polarisation, E = E_z z-hat):

    nabla^2(E_z) + k0^2 * eps_r(x,y) * E_z = 0

With complex eps_r = eps_r' + j eps_r'' and E_z = E_real + j E_imag, this splits into:

    E_real:  nabla^2(E_r) + k0^2 eps_r' E_r - k0^2 eps_r'' E_i = 0
    E_imag:  nabla^2(E_i) + k0^2 eps_r' E_i + k0^2 eps_r'' E_r = 0

Material parameters (lambda_0 = 1.0 unit, k0 = 2*pi = 6.2832 rad/unit):

    Free space (outside reflector):
        eps_r'  = 1.0,  eps_r'' = 0.0
        k0^2 eps_r'  = 4*pi^2 = 39.478
        k0^2 eps_r'' = 0  (no cross-coupling)

    Parabolic reflector (x < y^2/8 AND |y| < 3 AND x < 1.5 AND x > -1.5):
        eps_r'  = 1.0,  eps_r'' = 50.0
        k0^2 eps_r'  = 39.478  (same real part)
        k0^2 eps_r'' = 1973.92  (large imaginary: mimics PEC)

Parabola geometry:
- Equation: x = y^2 / (4f) = y^2 / 8.0  (with f = 2.0)
- Aperture: |y| < 3.0  (total aperture = 6 lambda_0)
- Parabola tip at (0, 0); parabola reaches x = 1.125 at |y| = 3
- Focal point at (f, 0) = (2, 0)

Boundary conditions:
- East boundary (x = 10): EMRobinBC port mode, injects plane wave E_inc = exp(-j k0 x)
  propagating in the -x direction from the right
- West boundary (x = -2): EMRobinBC absorbing mode, absorbs the focused reflected beam
- North, south boundaries: EMRobinBC absorbing mode

Domain: 2D, [-2, 10] x [-5, 5] (12 x 10 lambda_0). Mesh: 120 x 100 = 12,000 QUAD4.

## MOOSE Objects Used

- **Kernels**: For each field component — `Diffusion` (standard Laplacian with unit
  diffusivity), `ADMatReaction` (k0^2 eps_r' = 39.478, spatially uniform reaction),
  `ADMatCoupledForce` (cross-coupling +/- k0^2 eps_r'' between E_real and E_imag; nonzero
  only inside the lossy reflector region)
- **BCs**: `EMRobinBC` — eight blocks total (real and imaginary for each of four boundaries).
  East boundary in `mode = port` with `profile_func_real = 1` to inject the incident wave.
  West, north, south in `mode = absorbing` with `profile_func_real = 0`
- **Functions**: Three `ParsedFunction` blocks encoding k0^2 eps_r'(x,y), +k0^2 eps_r''(x,y),
  -k0^2 eps_r''(x,y). The reflector region is defined by the if() expression
  `x < y*y/8.0 & abs(y) < 3.0 & x < 1.5 & x > -1.5`
- **Materials**: Three `ADGenericFunctionMaterial` blocks wrapping the functions as AD
  material properties for use by ADMatReaction and ADMatCoupledForce
- **AuxVariables / AuxKernels**: `ParsedAux` computes E_intensity = E_real^2 + E_imag^2
- **Postprocessors**: `PointValue` at focal point (2,0), incident field check (8,0),
  pre-focus (1,0), post-focus (-1,0), off-axis near-focus (2,1), inside reflector (0.5,0),
  above aperture edge (5,3.5); `ElementIntegralVariablePostprocessor` for total field energy
- **Executioner**: Steady Newton with MUMPS LU; `SMP full = true` required for the strongly
  coupled off-diagonal blocks from EMRobinBC and the large ADMatCoupledForce cross terms

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case90-parabolic-reflector \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case90_parabolic_reflector.i 2>&1 | tail -30'
```

## Expected Results

Newton converges in 1-2 iterations (linear system). The focal spot at (2, 0) should
have significantly higher intensity than the incident wave (|E_intensity_incident| ≈ 1).

Expected postprocessor behaviour:
- `E_intensity_incident` ≈ 1.0 (unperturbed plane wave far from reflector)
- `E_intensity_focus` > 1.0 (concentrated reflected power; 2D focusing gain ≈ 4-6)
- `E_intensity_prefocus` > 1 but < `E_intensity_focus` (converging beam)
- `E_intensity_postfocus` decreasing (diverging beam past the focus)
- `E_intensity_focus_offaxis` (at (2,1)) < `E_intensity_focus` (spot width < lambda_0)
- `E_intensity_inside_reflector` ≈ 0 (field absorbed in lossy layer)

In ParaView, colour by E_intensity:
- Uniform incident field approaching from x = 10
- Sharp drop in field where the parabolic absorbing layer blocks the wave
- Bright focal spot at (2, 0) — constructive interference of all reflected rays
- Standing-wave pattern to the left of the reflector where incident and reflected
  beams superpose
- The parabolic boundary is visible as the sharp contour of near-zero intensity

## Key Takeaways

- A parabola focuses a plane wave to a point because all path lengths from the wavefront to
  the focus are equal (definition of a parabola from focus-directrix property) — all reflected
  contributions arrive in phase → constructive interference → intensity maximum.
- Embedding a curved conductor as a high-loss dielectric (epsilon_r'' ≫ 1) allows using
  structured rectangular meshes without curved boundaries; accuracy requires the skin depth
  delta_s ≪ lambda_0, which is satisfied here with delta_s ≈ 0.032 lambda_0.
- The two-component (E_real, E_imag) Helmholtz formulation with cross-coupling
  ADMatCoupledForce kernels is the standard FEM approach for lossy media in the
  frequency domain; the same pattern appears in Cases 75, 83, and 90-93.
- A full SMP preconditioner (full = true) is essential when ADMatCoupledForce creates
  large off-diagonal Jacobian blocks; without it Newton diverges in the high-loss region.
- The 2D focusing gain scales as aperture / lambda_0 (≈ 6 for this case), in contrast to
  the 3D case where gain scales as 4*pi*A / lambda_0^2 for aperture area A.
