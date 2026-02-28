# Case 88: Single-Slit Diffraction — Fraunhofer Limit and Resolution

**MIT 6.635 Lecture(s):** 3-4 (Diffraction theory, Huygens principle, aperture fields)

## Overview

This case solves the 2D frequency-domain Helmholtz equation to compute the far-field intensity
pattern produced by a plane wave passing through a single slit of width a = 3 lambda_0 in a
perfectly conducting screen. The result is the classic Fraunhofer (far-field) diffraction pattern
with a sinc-squared intensity envelope, nulls at sin(theta) = m lambda_0/a, and secondary lobes
whose intensity is much weaker than the central maximum.

The simulation uses the Kirchhoff boundary condition approach: the left domain boundary (x=0)
represents the screen/slit plane. A FunctionDirichletBC prescribes E = 1 at the open slit
(|y| < a/2 = 1.5 m) and E = 0 on the conducting screen (|y| >= 1.5 m), injecting the aperture
field directly. The right, top, and bottom boundaries use EMRobinBC as absorbing boundary
conditions, allowing the diffracted wave to exit the 15-wavelength-long domain without significant
spurious reflections. The diffraction pattern is sampled at x = 10 m (the observation screen),
which is well into the Fraunhofer regime (x >> a^2/lambda_0 = 9 m).

## The Physics

Governing equation (TE polarisation, scalar Helmholtz, uniform free-space medium):

    nabla^2 E_z + k0^2 E_z = 0

with k0 = 2 pi / lambda_0 = 6.2832 rad/m and lambda_0 = 1.0 m.

The complex phasor E_z = E_real + j E_imag satisfies two decoupled bulk equations:

    nabla^2 E_real + k0^2 E_real = 0
    nabla^2 E_imag + k0^2 E_imag = 0

E_real and E_imag are coupled only through the EMRobinBC terms (the j k0 factor in the
Sommerfeld absorbing boundary condition mixes the two components at the boundary).

Left boundary (x=0) — Kirchhoff aperture condition:

    E_real(0, y) = 1.0  if |y| < 1.5 m   (slit aperture: incoming plane wave)
    E_real(0, y) = 0.0  if |y| >= 1.5 m  (PEC screen: zero tangential E)
    E_imag(0, y) = 0.0  everywhere        (plane wave arrives with real amplitude at x=0)

Analytical far-field result (Fraunhofer diffraction, observation at distance z >> a^2/lambda_0):

    I(theta) = I_0 * sinc^2(beta)
    where: beta = pi * a / lambda_0 * sin(theta) = 3 pi * sin(theta)
           sinc(beta) = sin(beta) / beta

Null positions at x = 10 m observation plane (sin(theta_m) = m * lambda_0 / a = m/3):
- First null (m=1):    sin(theta) = 1/3,   theta = 19.47 deg,   y ~ 3.54 m
- Second null (m=2):   sin(theta) = 2/3,   theta = 41.81 deg,   y ~ 8.94 m (near domain edge)

Secondary maximum positions (sin(theta) ~ (2m+1)/(2a)):
- m=1 secondary max:   sin(theta) ~ 1/2,   theta ~ 30 deg,   y ~ 5.77 m,   I/I_0 ~ 0.045

Domain: [0, 15] x [-8, 8] m. Slit: a = 3 m (a = 3 lambda_0). Observation: x = 10 m.
Mesh: 150 x 160 QUAD4 elements, dx = dy = 0.1 m = lambda_0/10.

## MOOSE Objects Used

- **Kernels**: `Diffusion` (builds -nabla^2 E contribution in A); `ADMatReaction` with
  reaction_rate = k0^2 (builds -k0^2 E contribution, giving Helmholtz in strong form)
- **Materials**: `ADGenericConstantMaterial` wrapping k0^2 = 39.478 m^-2 as AD property k0sq
  (uniform free-space medium — no spatial variation, so ConstantMaterial is more efficient
  than FunctionMaterial)
- **BCs**: `FunctionDirichletBC` on the left boundary using slit_fn = if(abs(y)<1.5, 1, 0)
  for E_real (the key diffraction BC); plain `DirichletBC` with value=0 for E_imag on the
  left; `EMRobinBC` absorbing (profile_func_real=0) on right, top, and bottom boundaries —
  eight BC blocks total (real and imaginary component for each of the three absorbing edges)
- **AuxVariables / AuxKernels**: `ParsedAux` computes E_intensity = E_real^2 + E_imag^2
  for field intensity visualisation (proportional to time-averaged Poynting flux)
- **Postprocessors**: `PointValue` at 13 positions along x=10 m spanning y=0 to y=7 m
  (central maximum through second lobe region); near-field probes at x=1 m and x=5 m;
  `ElementIntegralVariablePostprocessor` for total field energy
- **Executioner**: Steady Newton with LU direct solver; SMP full=true for the off-diagonal
  Robin BC coupling between E_real and E_imag

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case88-single-slit-diffraction \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case88_single_slit_diffraction.i 2>&1 | tail -30'
```

## Expected Results

Newton converges in 1-2 iterations (the system is linear). The postprocessor CSV file shows:

- E_intensity_y0 (central maximum, y=0): highest intensity value (~1 after normalisation)
- E_intensity_y3p5 (near first null, y=3.5 m): near zero (null of sinc^2 at y_null ~ 3.54 m)
- E_intensity_y5p8 (near secondary maximum, y=5.8 m): approximately 4.5% of the central peak
- E_intensity_y7p0: decreasing, approaching the second null region

The intensity ratio E_intensity_y5p8 / E_intensity_y0 ~ 0.045 confirms the sinc^2 envelope.

In ParaView, colouring by E_intensity shows:
- A bright central horizontal stripe (the central diffraction maximum) propagating in the +x direction
- The diffraction pattern expanding laterally as x increases
- At x=10-15 m, distinct vertical bands (maxima) separated by dark bands (nulls) are visible
- Use the "Plot Over Line" filter along the vertical line at x=10 to extract I(y) and compare
  with the analytic sinc^2 curve

The E_real and E_imag fields show the complex phase structure: the diffracted wave accumulates
phase as it propagates, so E_real and E_imag alternate in dominance with period lambda_0.

## Key Takeaways

- The Kirchhoff aperture BC (FunctionDirichletBC with a step function) is the simplest MOOSE
  implementation of single-slit diffraction: prescribe E=1 at the open slit, E=0 on the screen.
- The Fraunhofer (far-field) regime requires the observation distance x >> a^2/lambda_0.
  With a = 3 and lambda_0 = 1, the Fraunhofer condition is x >> 9 m; at x=10 the pattern
  is already well-developed.
- First null angle: sin(theta_1) = lambda_0 / a = 1/3. Wider slits (larger a) produce
  narrower central maxima — the trade-off between aperture size and angular resolution.
- FunctionDirichletBC with parsed if() expressions naturally implements piecewise boundary
  conditions on a single named boundary set without requiring separate mesh sidesets.
- ADGenericConstantMaterial is the correct choice for spatially uniform coefficients (avoids
  overhead of evaluating a ParsedFunction at every quadrature point).
- The EMRobinBC absorbing BC works well for near-normally-incident diffracted waves (the central
  lobe) but has increasing reflection error at large angles (the secondary lobes near y ~ 6 m
  hit the right boundary at ~30 degrees from normal, where the first-order ABC reflection
  coefficient is |cos(30)-1|/|cos(30)+1| ~ 8%). A PML or higher-order ABC would improve accuracy
  for the secondary lobe measurements.
