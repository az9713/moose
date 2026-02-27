# Case 75: Lossy Drude Slab — Skin Depth and Attenuation

## Overview

This case models a plane EM wave at f = 20 MHz passing through a 30 m thick metallic slab
described by the Drude free-electron model. The Drude model gives the complex permittivity
of a metal as a function of frequency; at frequencies below the plasma frequency the real part
of epsilon_r becomes negative, causing evanescent (exponentially attenuating) rather than
propagating fields inside the material. The attenuation length is the skin depth.

The case introduces complex permittivity and the resulting cross-coupling between the real
and imaginary parts of the electric field. For a lossy medium eps_r = eps_r' + j*eps_r'', the
scalar Helmholtz equation splits into two coupled real equations. The coupling terms are
proportional to eps_r'' and are implemented using the ADMatCoupledForce kernel — a pattern
absent in the lossless cases (32, 74). This case relates to Griffiths, Introduction to
Electrodynamics, Chapter 9, and Cheng, Field and Wave Electromagnetics, Chapter 8.

For the parameters chosen (plasma frequency = 2*omega, collision rate = 0.3*omega), the slab
is approximately 21 skin depths thick. The field is essentially completely extinguished inside
the slab, giving |R|^2 approximately 1 (near-total reflection) at the port.

## The Physics

Drude model permittivity at omega = 2*pi*20 MHz, omega_p = 2*omega, gamma = 0.3*omega:

    eps_r = 1 - omega_p^2 / (omega^2 + j*gamma*omega) = -2.66972 + j*1.10092

The skin depth (characteristic attenuation length):

    delta = 1 / Im(k_z)  where k_z = k0 * sqrt(eps_r)

Numerically: k_z = 0.1383 + j*0.6983 m^-1, giving delta = 1/0.6983 = 1.432 m.
The 30 m slab is ~21 skin depths thick: field amplitude attenuates by exp(-21) ~ 10^-9.

Governing equations (frequency domain, 1D, splitting into real/imaginary parts):

    d^2(E_real)/dx^2 + k0^2 * eps_r'(x) * E_real + k0^2 * eps_r''(x) * E_imag = 0
    d^2(E_imag)/dx^2 + k0^2 * eps_r'(x) * E_imag - k0^2 * eps_r''(x) * E_real = 0

Boundary conditions:
- x = 0: PEC wall (DirichletBC, E_real = E_imag = 0)
- x = 100 m: Robin port (EMRobinBC) injecting E_inc and absorbing reflections

Material properties:
- Vacuum (x <= 30 m, x > 60 m): k0^2*eps_r' = +0.17546 m^-2, eps_r'' = 0 (no coupling)
- Drude slab (30 < x <= 60 m): k0^2*eps_r' = -0.46843, k0^2*eps_r'' = +0.19317

Domain: 1D, x in [0, 100] m. Drude slab at [30, 60] m with PEC wall behind it at x=0.
Mesh: 500 first-order Lagrange elements.

## Input File Walkthrough

**HIT variables**: k = 0.41888 rad/m, L = 100 m, E0 = 1 V/m, theta = 0.

**[Mesh]**: 1D GeneratedMesh from x=0 to x=100 with 500 elements. Boundaries renamed
to `metal` (x=0) and `port` (x=100).

**[Variables]**: E_real and E_imag, both FIRST LAGRANGE.

**[Functions]**: Three ParsedFunctions. `coeff_real_fn` gives k0^2*Re(eps_r(x)).
`coeff_imag_fn` gives +k0^2*Im(eps_r(x)) for the E_real cross-coupling term.
`coeff_neg_imag_fn` gives -k0^2*Im(eps_r(x)) for the E_imag cross-coupling term.
All three are non-zero only inside the Drude slab [30, 60] m.

**[Materials]**: ADGenericFunctionMaterial wraps each of the three coefficient functions
as AD material properties for use in the AD kernels.

**[Kernels]**: Each field equation uses three kernels. `diffusion_real`/`diffusion_imag`
(plain Diffusion kernel) provide the Laplacian term. `reaction_real`/`reaction_imag`
(ADMatReaction) add the k0^2*eps_r' diagonal term. `cross_real_from_imag` (ADMatCoupledForce
with coeff_imag_material) adds the +k0^2*eps_r''*E_imag coupling into the E_real equation.
`cross_imag_from_real` (ADMatCoupledForce with coeff_neg_imag_material) adds the
-k0^2*eps_r''*E_real coupling into the E_imag equation. These cross-coupling kernels are
absent in the lossless cases and are the key new feature of this case.

**[BCs]**: DirichletBC at x=0 (PEC). EMRobinBC at x=100 (port). Both field components
are coupled at the port boundary; full SMP preconditioning is required.

**[Postprocessors]**: ReflectionCoefficient computes |R|^2. PointValue postprocessors
sample the field at x=15 (vacuum gap), 35 (near slab entrance), 45 (slab mid), 55 (deep
in slab), 80 (transmitted vacuum), and 100 (port) to trace the exponential decay.

**[Preconditioning]**: SMP with full=true is required because the ADMatCoupledForce kernels
create off-diagonal Jacobian blocks between E_real and E_imag throughout the slab.

**[Executioner]**: Steady Newton with direct LU (pc_type = lu). The Helmholtz operator is
indefinite below the plasma frequency; direct LU avoids iterative solver stagnation.

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case75-drude-slab \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case75_drude_slab.i 2>&1 | tail -30'
```

## Expected Results

Newton converges in 1-2 iterations. The key physical result is the exponential field
decay inside the Drude slab. Comparing |E|^2 at x=35, 45, 55 m should show a ratio
consistent with exp(-2*Im(k_z)*delta_x) = exp(-2*0.6983*10) = exp(-13.97) per 10 m step.
The reflection coefficient |R|^2 should be close to 1.0 (near-total reflection) because
the slab is ~21 skin depths thick and essentially opaque. Any deviation from unity
represents ohmic absorption.

Visualise `case75_drude_slab_out.e` in ParaView: plot E_real and E_imag along x to see
the abrupt transition from oscillatory to evanescent behaviour at the slab entrance (x=30).

## Key Takeaways

- The Drude model gives complex epsilon_r with Re(eps_r) < 0 below the plasma frequency.
- Complex epsilon_r couples the real and imaginary parts of E in the bulk via ADMatCoupledForce.
- The skin depth delta = 1/Im(k0*sqrt(eps_r)) characterises the exponential attenuation.
- Full SMP preconditioning is essential when cross-coupling kernels are present.
- A slab of ~21 skin depths is effectively opaque: |R|^2 ~ 1.
- Relates to Griffiths Ch. 9 (electromagnetic waves in conducting media) and Drude model.
