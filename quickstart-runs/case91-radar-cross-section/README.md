# Case 91: Radar Cross Section of a Finite Flat Conducting Plate

**MIT 6.635 Lecture themes:** Scattering, EFIE, aperture diffraction, RCS

## Overview

This case computes the 2D Radar Cross Section (RCS, or scattering width W) of a finite
flat conducting strip (plate width L = 3 lambda_0, normal incidence). The 2D frequency-domain
scalar Helmholtz equation is solved for the total electric field (TE polarisation, E_z out of
plane). The scattered field is then extracted analytically by subtracting the known incident
plane wave, and the bistatic RCS pattern is sampled at seven observation angles at R = 5
from the plate.

The case demonstrates the key physics of flat-plate scattering: a large specular main lobe
at forward scatter (phi = 0°) and backscatter (phi = 180°), and weaker edge-diffraction
lobes at broadside (phi = 90°) consistent with the Geometrical Theory of Diffraction (GTD).

The conducting plate is modelled as a thin high-loss dielectric strip (epsilon_r'' = 100,
skin depth delta_s ≈ 0.016 lambda_0 ≪ plate thickness 0.2 lambda_0) embedded in the
rectangular domain — the same absorbing-material approach used in Case 90. This avoids
internal boundary conditions on the plate surface while accurately reproducing the
specular reflection and edge-diffraction physics.

## The Physics

Governing equation (2D scalar Helmholtz, TE polarisation):

    nabla^2(E_z) + k0^2 * eps_r(x,y) * E_z = 0

Material definition:

    Free space (|x| > 0.1 OR |y| > 1.5):
        eps_r = 1.0  (real, lossless)
        k0^2 eps_r' = 39.478,  k0^2 eps_r'' = 0

    Conducting plate (|x| < 0.1 AND |y| < 1.5):
        eps_r = 1 + j * 100  (high loss, PEC model)
        k0^2 eps_r'' = 3947.84  (skin depth delta_s ≈ 0.016 lambda_0)

Numerical parameters (lambda_0 = 1.0, k0 = 2*pi):

    k0 = 6.2832 rad/unit,  k0^2 = 39.478 unit^-2

Scattered field extraction:

    E_inc(x,y) = exp(+j k0 x) = cos(k0 x) + j sin(k0 x)
    E_scat = E_total - E_inc
    |E_scat|^2 = (E_real - cos(k0 x))^2 + (E_imag - sin(k0 x))^2

2D scattering width (RCS measure) at observation point (x_obs, y_obs), R = 5:

    W(phi) = 2*pi*R * |E_scat(R, phi)|^2 / |E_inc|^2
           = 2*pi*5 * Es_sq(phi)   (since |E_inc| = 1)

Physical optics prediction (L = 3, lambda_0 = 1, normal incidence):

    Monostatic (phi = 180°):  W_back ≈ k0 L^2 / (4*pi) ≈ 1.43 lambda_0
    Forward (phi = 0°):       W_fwd  ≈ k0 L^2 / pi ≈ 5.73 lambda_0  (stronger than back)
    First null near phi = 90° ± arcsin(lambda_0/L): null at phi ≈ 90° ± 19.5°

Boundary conditions:
- West (x = -8): EMRobinBC port mode, injects E_inc = exp(+j k0 x), absorbs backscatter
- East (x = +8): EMRobinBC absorbing mode
- North, south (y = +/-8): EMRobinBC absorbing mode

Domain: 2D, [-8, 8] x [-8, 8] (16 x 16 lambda_0). Mesh: 160 x 160 = 25,600 QUAD4.

## MOOSE Objects Used

- **Kernels**: `Diffusion` + `ADMatReaction(k0sq_eps_pr)` + `ADMatCoupledForce` cross terms for
  both E_real and E_imag. The same three-kernel pattern as Case 90. The large cross-coupling
  coefficient (k0^2 eps_r'' = 3947.84 inside the plate) forces the field inside the plate
  toward zero, implementing the PEC condition
- **BCs**: Eight `EMRobinBC` blocks (real + imaginary for each of four boundaries). West
  boundary uses `mode = port` with `profile_func_real = 1` to inject the incident wave;
  east, north, south use `mode = absorbing`
- **Functions**: Three `ParsedFunction` blocks for k0^2 eps_r', +k0^2 eps_r'', -k0^2 eps_r''.
  Plate region defined by `abs(x) < 0.1 & abs(y) < 1.5`
- **Materials**: Three `ADGenericFunctionMaterial` blocks wrapping the functions as AD properties
- **AuxVariables / AuxKernels**: Four ParsedAux variables:
  - `E_intensity` = E_real^2 + E_imag^2 (total field intensity)
  - `Es_real` = E_real - cos(k0 x) (scattered field, real part)
  - `Es_imag` = E_imag - sin(k0 x) (scattered field, imaginary part)
  - `Es_sq` = Es_real^2 + Es_imag^2 (scattered field intensity, RCS proxy)
- **Postprocessors**: `PointValue` at plate centre (field suppression check), incident check
  at (5,5), and nine RCS observation points on a circle R=5 at angles 0°, 30°, 60°, 90°,
  120°, 150°, 180°, 270°, 330°; `ElementIntegralVariablePostprocessor` for total scattered
  and total field energy
- **Executioner**: Steady Newton with MUMPS LU; `SMP full = true` required for the large
  off-diagonal cross terms (k0^2 eps_r'' = 3948) and EMRobinBC coupling

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case91-radar-cross-section \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case91_radar_cross_section.i 2>&1 | tail -30'
```

## Expected Results

Newton converges in 1-2 iterations. Key postprocessor observations:

- `E_intensity_plate_centre` ≈ 0 (PEC condition: tangential E → 0 inside conductor)
- `E_intensity_check` ≈ 1.0 (unperturbed incident wave far from plate)
- `Es_sq_phi_000` (forward scatter, phi=0°): large value — forward main lobe
- `Es_sq_phi_180` (backscatter, phi=180°): large value — specular main lobe
- `Es_sq_phi_090` and `Es_sq_phi_270`: should be approximately equal (symmetry check)
- `Es_sq_phi_090` < `Es_sq_phi_180` (edge diffraction weaker than specular reflection)

Physical optics ratio estimate at R = 5:
    Es_sq_phi_180 ≈ W_back / (2*pi*R) ≈ 1.43 / (31.4) ≈ 0.046

In ParaView, visualise E_intensity and Es_sq:
- `E_intensity`: standing-wave fringes on the left (incident + backscattered), shadow behind
  the plate (reduced intensity at x > 0.1, |y| < 1.5), near-zero inside the plate
- `Es_sq`: bistatic RCS pattern as a 2D field — bright lobes in forward and backward
  directions, weaker edge-diffraction radiation at broadside; the two plate edge tips at
  (0, +/-1.5) appear as bright point sources of edge-diffracted scattered field
- To plot the RCS polar pattern: extract Es_sq on a circle of radius R = 5 using ParaView's
  Probe Line filter, then convert to polar coordinates

## Key Takeaways

- The high-loss dielectric model for a PEC flat plate (epsilon_r'' = 100, skin depth 0.016
  lambda_0) accurately reproduces specular reflection and edge diffraction on a structured
  rectangular mesh, avoiding the need for curved internal boundaries.
- The total-field approach (solve for E_total = E_inc + E_scat, then subtract E_inc
  analytically via ParsedAux) is simpler to implement than the scattered-field formulation
  for normal-incidence problems where the source can be injected via an EMRobinBC port.
- The bistatic RCS pattern of a flat strip is dominated by two mechanisms:
  (1) Specular reflection: main lobes at phi = 0° and phi = 180° with width ~ lambda_0/L
  (2) Edge diffraction (GTD): contributes at all angles, most visible at phi = 90°
- Problem symmetry (y -> -y): Es_sq at phi = +90° must equal Es_sq at phi = -90° (= 270°).
  Verifying this symmetry validates the solver and mesh.
- Increasing plate width L increases the peak RCS as L^2 (physical optics) while narrowing
  the main lobe as lambda_0/L. The side-lobe level remains approximately -13 dB below the
  main lobe (sinc-squared pattern) independent of L.
- Relates to MIT 6.635 topics: EFIE, physical optics approximation, GTD, Method of Moments
  applied to flat plate scattering.
