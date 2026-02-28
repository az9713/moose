# Case 87: Phased Array Beamforming — Two-Element Array

**MIT 6.635 Advanced Electromagnetism, Spring 2003 — Professor Jin Au Kong**
**Reference:** Balanis, "Antenna Theory: Analysis and Design", 4th Ed., Ch. 6; Van Trees, "Optimum Array Processing", Ch. 1

## Overview

This case simulates electronic beam steering using a two-element phased array antenna. Two Gaussian
point sources are placed along the y-axis, separated by d = lambda/2 (half-wavelength spacing). A
progressive phase offset delta = pi/4 is applied to the second element, steering the main radiation
lobe away from the broadside direction.

The array factor (pattern due to the array geometry alone) for two elements is:

    |AF(theta)|^2 = |1 + exp(j*(k*d*sin(theta) + delta))|^2 = 4*cos^2((k*d*sin(theta) + delta)/2)

where theta is measured from the array broadside direction (the +y axis, perpendicular to the
array axis). The maximum occurs when k*d*sin(theta) + delta = 0:

    sin(theta_max) = -delta / (k*d) = -(pi/4) / pi = -0.25
    theta_max = arcsin(-0.25) = -14.48 degrees

The beam is steered 14.5 degrees off broadside toward negative x (toward the first element at
y = -0.25). This case demonstrates that a phase shift alone — without moving the antenna physically —
is sufficient to redirect the main lobe, the fundamental principle of phased array radar and 5G
beamforming systems.

## The Physics

Governing equation (2D scalar Helmholtz, Cartesian):

    nabla^2 E + k0^2 E = -source(x, y)

where:
- k0 = 2*pi rad/unit (free-space wavenumber, lambda_0 = 1 unit)
- k0^2 = 4*pi^2 = 39.4784176 unit^-2

Two Gaussian sources model the antenna elements (sigma = 0.05, 2*sigma^2 = 0.005):

    Source 1 at (0, -0.25): amplitude A = 50, phase offset delta_1 = 0
      Drives E_imag only (purely imaginary phasor)

    Source 2 at (0, +0.25): amplitude A = 50, phase offset delta_2 = pi/4
      Drives E_imag with A*cos(pi/4) = 50*0.70711 (imaginary part of phased phasor)
      Drives E_real with A*sin(pi/4) = 50*0.70711 (real part of phased phasor)

The phase offset delta = pi/4 is implemented by distributing the excitation between E_real and
E_imag according to sin(delta) and cos(delta). This is equivalent to multiplying the source
phasor by e^{j*delta} in the complex field representation.

Splitting E = E_real + j*E_imag:

    E_real:  nabla^2 E_r + k0^2 E_r = source2_real (only Source 2's real part)
    E_imag:  nabla^2 E_i + k0^2 E_i = source1_imag + source2_imag

Boundary conditions: EMRobinBC absorbing on all four sides (profile_func_real = 0)

Domain: 2D Cartesian, x in [-8, 8], y in [-8, 8] (16 lambda_0 square)
Mesh: 80 x 80 = 6,400 elements

## Expected Array Pattern

Array factor |AF(theta)|^2 for d = lambda/2, delta = pi/4 (theta from +y axis):

| theta | x = 5*sin(theta) | y = 5*cos(theta) | |AF|^2 |
|-------|-----------------|-----------------|--------|
| 0 deg | 0.000 | 5.000 | 3.41 |
| +15 deg | 1.294 | 4.830 | 1.94 |
| +30 deg | 2.500 | 4.330 | 0.59 |
| +45 deg | 3.536 | 3.536 | 0.02 (near null) |
| +60 deg | 4.330 | 2.500 | 0.13 |
| +75 deg | 4.830 | 1.294 | 0.44 |
| +90 deg | 5.000 | 0.000 | 0.59 |
| -15 deg | -1.294 | 4.830 | 4.00 (BEAM MAX) |
| -30 deg | -2.500 | 4.330 | 3.41 |
| -45 deg | -3.536 | 3.536 | 2.27 |

Key diagnostic: E_int_m15 > E_int_0 > E_int_p15 confirms beam steering.
The ratio E_int_m15 / E_int_p45 should be approximately 4/0.02 = 200, showing a deep null
at theta = +45 degrees and the maximum at theta = -15 degrees.

## MOOSE Objects Used

- **Kernels**: `Diffusion` + `ADMatReaction` for each field component; `BodyForce` kernels for
  three source functions: `src1_imag` (Source 1, E_imag only), `src2_imag` (Source 2 imaginary
  part, amplitude A*cos(delta), E_imag), and `src2_real` (Source 2 real part, amplitude A*sin(delta), E_real)
- **Materials**: `ADGenericConstantMaterial` with `k0sq = 39.4784176` (free space throughout)
- **BCs**: `EMRobinBC` absorbing on all four boundaries (left, right, top, bottom), with
  `profile_func_real = 0` — pure absorbing, no external incident wave
- **AuxVariables / AuxKernels**: `ParsedAux` computes E_intensity = E_real^2 + E_imag^2
- **Postprocessors**: `PointValue` at 10 angular positions on the R = 5 measurement circle,
  plus source locations and total energy integral
- **Executioner**: Steady Newton with LU direct solver; SMP full = true for EMRobinBC coupling

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case87-phased-array \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case87_phased_array.i 2>&1 | tail -30'
```

## Key Takeaways

- Phased arrays steer beams electronically by applying progressive phase shifts between elements,
  without mechanical motion. The steering formula sin(theta_max) = -delta/(k*d) is the
  fundamental equation of phased array beamforming.
- Half-wavelength spacing (d = lambda/2) prevents grating lobes: with kd = pi, the maximum
  possible |sin(theta_max)| = |delta/pi| < 1 for all -pi < delta < pi, so the beam always
  steers within the physical half-space (no aliased grating lobe).
- A phase offset delta = pi/4 with d = lambda/2 gives modest steering (14.5 degrees). For
  maximum steering (endfire, theta = 90 degrees), set delta = -pi (= -kd): the beam points
  along the array axis.
- The phase offset is encoded in the MOOSE model by splitting the complex phasor A*exp(j*delta)
  into real part A*sin(delta) acting on E_real and imaginary part A*cos(delta) acting on E_imag.
- The asymmetry in the E_intensity field — higher on the theta = -15 deg side than the theta = +15
  deg side — is the direct FEM demonstration of electronic beam steering.
- For a symmetric two-element in-phase array (delta = 0), both sources drive only E_imag and the
  pattern would be exactly symmetric about the y-axis (pure broadside beam). Breaking this symmetry
  with delta = pi/4 is visible as E_int_m15 > E_int_p15 in the postprocessor output.
- Pattern multiplication theorem: the total pattern of an array of isotropic elements equals the
  element pattern times the array factor. Since Gaussian sources are approximately isotropic in 2D,
  the measured E_intensity pattern closely follows the array factor |AF(theta)|^2.
