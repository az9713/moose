# Case 86: Half-Wave Dipole Antenna — Radiation Pattern

**MIT 6.635 Advanced Electromagnetism, Spring 2003 — Professor Jin Au Kong**
**Reference:** Balanis, "Antenna Theory: Analysis and Design", 4th Ed., Ch. 4; Kong, "Electromagnetic Wave Theory", Ch. 7

## Overview

This case simulates the radiation pattern of a half-wave dipole antenna using a 2D axisymmetric
(RZ) frequency-domain Helmholtz formulation. The half-wave dipole (total length L = lambda/2) is
the canonical resonant antenna element, with an input impedance of approximately 73 Ohm and a
directivity of 1.64 dBi — slightly higher than the Hertzian (infinitesimal) dipole of Case 85.

The defining feature of the half-wave dipole is its cosine current distribution:

    I(z) = I_0 * cos(pi*z/L)    for |z| <= L/2

which differs from the uniform distribution of the Hertzian dipole and produces the characteristic
gain pattern G(theta) = cos^2((pi/2)*cos(theta)) / sin^2(theta) with maximum at theta = 90
degrees (broadside) and nulls at theta = 0 and 180 degrees (along the dipole axis).

The axisymmetric formulation (coord_type = RZ) reduces the full 3D problem to 2D, exploiting
the rotational symmetry of the dipole about the z-axis. The MOOSE RZ coordinate system maps x to
the radial coordinate r and y to the axial coordinate z.

## The Physics

Governing equation (frequency-domain Helmholtz, axisymmetric):

    (1/r) d(r dE/dr)/dr + d^2E/dz^2 + k0^2 E = -source(r, z)

where:
- k0 = 2*pi rad/unit (free-space wavenumber, lambda_0 = 1 unit)
- k0^2 = 4*pi^2 = 39.4784176 unit^-2
- L = 0.5 unit (dipole length = lambda_0/2)

The distributed current source:

    source(r, z) = 100 * cos(pi*z/L)    if r < 0.05 AND |z| < 0.25
                   0                     otherwise

This confines the cosine current to a thin strip near r = 0 (approximating the antenna wire)
with the cosine distribution I_0*cos(pi*z/0.5) along the z-axis.

Splitting E = E_real + j*E_imag with real k0^2:

    E_real:  nabla^2 E_r + k0^2 E_r = 0           (homogeneous — no source)
    E_imag:  nabla^2 E_i + k0^2 E_i = -source(r,z) (driven by distributed current)

The source acts only on E_imag because the phasor of a real current J_z is -j*omega*mu_0*J_z,
which is purely imaginary. E_real is excited by off-diagonal coupling through the EMRobinBC.

Boundary conditions:
- r = 0 (left): natural Neumann (symmetry axis, automatic in RZ formulation)
- r = 6, z = +/-6 (right, top, bottom): EMRobinBC absorbing (profile_func_real = 0)

Domain: RZ, r in [0, 6], z in [-6, 6] (6 lambda_0 in each direction)
Mesh: 60 x 120 = 7,200 elements, Δr = Δz = 0.1 (10 elements per lambda_0)

## Expected Radiation Pattern

The normalised power gain of the half-wave dipole (Balanis, Eq. 4-58a):

    G(theta) = cos^2((pi/2)*cos(theta)) / sin^2(theta)

Numerical values at the measurement arc r_meas = 4 (theta from z-axis):

| theta | r = 4*sin(theta) | z = 4*cos(theta) | G(theta) normalised |
|-------|-----------------|-----------------|---------------------|
| 10 deg | 0.6946 | 3.9392 | 0.019 (near null) |
| 30 deg | 2.0000 | 3.4641 | 0.175 |
| 45 deg | 2.8284 | 2.8284 | 0.394 |
| 60 deg | 3.4641 | 2.0000 | 0.667 |
| 75 deg | 3.8637 | 1.0353 | 0.904 |
| 90 deg | 4.0000 | 0.0000 | 1.000 (maximum — broadside) |

The simulation should produce E_intensity ratios (E_int_thetaN / E_int_theta90) that
approximately match these normalised gain values. Note that this is a near-field FEM simulation
(measurement at r = 4 = 4 lambda_0), so there will be some deviation from the pure far-field
pattern formula, particularly at small theta near the axis null.

## MOOSE Objects Used

- **Problem**: `coord_type = RZ` — activates 2D axisymmetric (cylindrical) coordinate
  formulation; MOOSE automatically applies the r Jacobian factor to all weak form integrals
- **Kernels**: `Diffusion` + `ADMatReaction` for each field component; `BodyForce` on E_imag
  with `source_fn` (cosine current distribution in thin strip near r=0)
- **Materials**: `ADGenericConstantMaterial` with `k0sq = 39.4784176` (free space throughout)
- **BCs**: `EMRobinBC` absorbing on right, top, bottom boundaries (profile_func_real = 0);
  no explicit BC on left (r=0 axis — natural Neumann is automatic in RZ)
- **AuxVariables / AuxKernels**: `ParsedAux` computes E_intensity = E_real^2 + E_imag^2
- **Postprocessors**: `PointValue` at six angles on the r = 4 measurement arc, plus on-axis
  null check and total radiated energy integral
- **Executioner**: Steady Newton with LU direct solver; SMP full = true for EMRobinBC
  off-diagonal coupling between E_real and E_imag

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case86-half-wave-dipole \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case86_half_wave_dipole.i 2>&1 | tail -30'
```

## Key Takeaways

- The RZ axisymmetric formulation (coord_type = RZ) reduces a 3D rotationally-symmetric
  problem to 2D, with MOOSE automatically handling the cylindrical Jacobian factor r.
- The r = 0 axis requires no explicit BC: the cylindrical Jacobian r = 0 makes the
  boundary flux integral vanish identically, automatically enforcing the symmetry condition.
- The cosine current distribution I(z) = I_0*cos(pi*z/L) is physically required — the current
  must vanish at the open endpoints z = +/-L/2 = +/-0.25, which cos achieves naturally.
- The BodyForce source acts on E_imag only, because the phasor of a real current source
  is purely imaginary (-j*omega*mu_0*J). E_real is indirectly driven by the Robin BCs.
- The half-wave dipole pattern G(theta) = cos^2((pi/2)*cos theta)/sin^2(theta) is narrower
  than the Hertzian dipole sin^2(theta) pattern, with greater directivity (1.64 vs 1.5 dBi).
- The Revolve filter in ParaView can rotate the 2D RZ solution into a full 3D donut-shaped
  radiation pattern, confirming the azimuthal symmetry of the antenna.
- For better far-field accuracy, increase the domain size or add PML (perfectly matched layer)
  absorbing boundaries; the first-order Robin ABC has residual reflections at oblique incidence.
