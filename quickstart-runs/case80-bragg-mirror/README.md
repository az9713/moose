# Case 80: Multilayer Dielectric Stack — Bragg Mirror

**MIT 6.635 Lecture(s):** 5 (Interference and multilayer stacks; distributed Bragg reflectors)

## Overview

This case models a 1D distributed Bragg reflector (DBR) — five periods of alternating high/low
permittivity dielectric layers backed by a PEC wall, operating at f0 = 20 MHz (lambda_0 = 15 m).
Each layer is designed with a quarter-wave optical thickness: n_H * d_H = n_L * d_L = lambda_0/4.
At the design frequency, partial reflections from every H-L interface arrive back at the input
with a round-trip phase shift of exactly pi, producing constructive interference that builds up
to near-perfect reflectivity. The analytic transfer-matrix method predicts |R|^2 = 1 exactly at
f0, and the FEM solution is verified against this result.

The case extends Case 32 (single dielectric slab, partial Fabry-Perot reflectance) to a
10-layer periodic stack, introducing the concept of photonic stopbands. The entire material
profile is encoded in a single nested ParsedFunction — no subdomain splitting is needed. This
demonstrates that MOOSE's pointwise material evaluation scales naturally from 1 to N layers
with no change in kernel structure.

## The Physics

Layer design — quarter-wave optical thickness at lambda_0 = 15 m:

    n_H * d_H = lambda_0/4  ->  d_H = 15/(4*2.0)         = 1.875 m  (H layer, eps_H = 4, n_H = 2)
    n_L * d_L = lambda_0/4  ->  d_L = 15/(4*sqrt(1.5))   = 3.062 m  (L layer, eps_L = 1.5, n_L = 1.225)

Governing equation (1D Helmholtz, normal incidence, lossless):

    d^2(E)/dx^2 + k0^2 * eps_r(x) * E = 0

where k0 = 2 pi f0/c = 0.41888 rad/m and eps_r(x) is piecewise constant:
- H layers (eps_r = 4):    k0^2 * 4   = 0.70184 m^-2
- L layers (eps_r = 1.5):  k0^2 * 1.5 = 0.26319 m^-2
- Vacuum  (x > 24.684 m):  k0^2 * 1   = 0.17546 m^-2

Transfer-matrix result: at f = f0 the per-period matrix is diagonal (quarter-wave condition),
and after 5 periods backed by PEC the amplitude reflection coefficient |R| = 1 exactly.

Boundary conditions:
- x = 0 (PEC wall): DirichletBC, E_real = E_imag = 0
- x = 75 m (port): EMRobinBC injects the incident wave and absorbs the reflected wave

Domain: 1D, x in [0, 75] m. Stack occupies [0, 24.684] m, vacuum gap [24.684, 75] m.
Mesh: 750 LAGRANGE/FIRST elements, dx = 0.1 m (~18 elements per H layer).

## MOOSE Objects Used

- **Kernels**: `Diffusion` (d^2E/dx^2 term), `ADMatReaction` with reaction_rate =
  coeff_real_material (k0^2 eps_r(x) coefficient); one pair per field component
- **BCs**: `DirichletBC` at the PEC wall (E = 0 for both components); `EMRobinBC` at the port
  (injects incident wave and absorbs reflected wave, couples E_real and E_imag)
- **Materials**: `ADGenericFunctionMaterial` wrapping the nested ParsedFunction coeff_real_fn
  as AD property coeff_real_material
- **Postprocessors**: `ReflectionCoefficient` at the port boundary (computes |R|^2 from the
  total field); `PointValue` at x = 0, 10.81, 50, 75 for spatial diagnostics
- **Preconditioning**: SMP full = true — required because the Robin port BC introduces
  off-diagonal Jacobian blocks coupling E_real and E_imag
- **Executioner**: Steady Newton with LU direct solver (Helmholtz is indefinite near resonance;
  LU is robust where iterative solvers stagnate)

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case80-bragg-mirror \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case80_bragg_mirror.i 2>&1 | tail -30'
```

## Expected Results

The primary result is `reflection_coefficient` ~ 1.0000 in the CSV output, confirming
near-perfect reflectivity at the design frequency. Small deviations from 1.0 arise from
FEM discretisation error (O(h^2) for first-order elements).

In the vacuum gap (x > 24.684 m), with |R| ~ 1, the field approaches a pure standing wave
with amplitude ~2 E0 = 2 V/m. Inside the Bragg stack the field is strongly attenuated between
the port-facing entrance and the PEC wall, demonstrating that no energy penetrates the stopband.

Compare with Case 32 (single slab), which shows partial transmission and Fabry-Perot oscillations.
The 5-period quarter-wave stack eliminates all transmission across the stopband by design.

## Key Takeaways

- A quarter-wave dielectric stack (n_H d_H = lambda/4) achieves constructive interference of
  Fresnel reflections, reaching |R|^2 ~ 1 at the design frequency.
- The 10-layer material profile is encoded as a single deeply nested ParsedFunction — no
  subdomain splitting, no MeshSubdomainIDs, just one if() expression evaluated pointwise.
- The transfer-matrix method gives an analytic |R| = 1 prediction that provides a direct
  quantitative test of the FEM solver.
- Changing `k` in the HIT variables shifts the operating frequency off resonance, revealing
  the finite-bandwidth stopband (the mirror is frequency selective, not broadband).
- The same Diffusion + ADMatReaction kernel pair scales from Case 32 (1 layer) to this case
  (10 layers) with no change in the MOOSE input structure.
