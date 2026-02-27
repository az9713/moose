# Case 83: Veselago Flat Lens — Point Source Focusing by a Negative-Index Slab

**MIT 6.635 Lecture(s):** 1 (Left-handed materials, perfect lens), 10 (Negative refraction in photonic crystals)

## Overview

This case demonstrates the Veselago flat lens: a slab of n = -1 material (epsilon_r = -1,
mu_r = -1) of thickness d focuses a point source at distance s from one face to an image at
distance d - s on the far side. With the slab occupying -1 < x < +1 (d = 2) and the source at
(-2, 0) (s = 1 from the left face), the predicted focus locations are at (0, 0) inside the
slab and (+2, 0) outside — the external image point.

Veselago predicted negative-index focusing in 1968; Pendry showed in 2000 that the ideal
lossless n = -1 slab also amplifies evanescent waves, enabling perfect sub-diffraction imaging.
This case focuses on the propagating-wave component of the effect, which is already visible at
the geometric-optics level. A small imaginary loss delta = 0.3 regularises the near-singular
lossless LHM system and couples E_real and E_imag through cross terms in the slab.

The key new element compared to Case 74 (1D LHM) is the 2D 1/mu_r weighted diffusion
formulation: `ADMatDiffusion` with a spatially varying diffusivity that is -1 inside the slab,
which is the FEM expression of negative permeability and negative refraction.

## The Physics

Governing equation for TE polarisation (E = E_z z-hat):

    nabla . (1/mu_r * nabla E_z) + k0^2 * eps_r * E_z = -j omega mu_0 J_z

With complex eps_r = eps_r' + j eps_r'' and splitting E_z = E_real + j E_imag:

    E_real:  nabla.(1/mu_r nabla E_r) + k0^2 eps_r' E_r + k0^2 eps_r'' E_i = 0
    E_imag:  nabla.(1/mu_r nabla E_i) + k0^2 eps_r' E_i - k0^2 eps_r'' E_r = f(x,y)

Point source Gaussian at (-2, 0):  f(x,y) = 50 exp(-((x+2)^2 + y^2) / 0.045),  sigma = 0.15

Material parameters (lambda_0 = 4, k0 = pi/2 = 1.5708 rad/unit):
- Vacuum (|x| > 1):     1/mu_r = +1.0,   k0^2 eps_r' = +2.4674,   eps_r'' = 0
- LHM slab (-1<x<+1):   1/mu_r = -1.0,   k0^2 eps_r' = -2.4674,   k0^2 eps_r'' = +0.7402

Interface condition at x = +/-1: (1/mu_r) dE/dn is continuous, meaning the normal derivative
reverses sign across the vacuum-LHM interface — the mathematical signature of negative refraction.
Impedance Z = sqrt(mu_r/eps_r) = sqrt((-1)/(-1)) = 1 matches vacuum: zero reflection in the
lossless limit.

Boundary conditions: EMRobinBC absorbing (profile_func_real = 0) on all four domain faces.
No incident wave; source is the interior BodyForce on E_imag.

Domain: 2D, [-6, 6] x [-4, 4] (3 lambda_0 x 2 lambda_0). LHM slab at -1 < x < +1 (d = 2).
Mesh: 120 x 80 = 9600 QUAD4 elements, 10 elements per unit, 40 per lambda_0.

## MOOSE Objects Used

- **Kernels**: For each field component — `ADMatDiffusion` (1/mu_r weighted diffusion, negative
  in slab), `ADMatReaction` (k0^2 eps_r' diagonal coefficient), `ADMatCoupledForce` (cross-coupling
  +/- k0^2 eps_r'' between E_real and E_imag in the lossy slab). `BodyForce` on E_imag only
  (the current source -j omega mu_0 J has only an imaginary phasor component)
- **BCs**: `EMRobinBC` — eight absorbing blocks on all four domain boundaries (real and imaginary
  components), all with profile_func_real = 0
- **Materials**: Four `ADGenericFunctionMaterial` blocks wrapping inv_mu_fn, k0sq_eps_pr_fn,
  k0sq_eps_pp_fn, and neg_k0sq_eps_pp_fn as AD properties
- **AuxVariables / AuxKernels**: `ParsedAux` computes E_intensity = E_real^2 + E_imag^2 at
  each node for field intensity visualisation
- **Postprocessors**: `PointValue` at source (-2,0), internal focus (0,0), external focus (+2,0),
  far field (+4,0), off-axis slab (0,1), and vacuum midpoints (-1.5,0) and (+1.5,0)
- **Executioner**: Steady Newton with MUMPS LU direct solver; SMP full = true is required for
  the negative-definite slab sub-blocks and the cross-coupling Robin BC off-diagonal terms

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case83-veselago-lens \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case83_veselago_lens.i 2>&1 | tail -30'
```

## Expected Results

Newton converges in 1-2 iterations. The key signature of flat-lens focusing is three intensity
peaks along y = 0 at x = {-2, 0, +2}:

- Primary peak at the source (-2, 0): highest intensity, the driving current location.
- Secondary peak at the internal focus (0, 0): inside the LHM slab at the slab centre.
- Third peak at the external image focus (+2, 0): the Veselago image point.

The intensity ratio (focus / source) is less than 1 due to loss regularisation delta = 0.3.
The far-field probe at (+4, 0) should show lower intensity than the focus at (+2, 0), confirming
a genuine local maximum at the image point rather than a monotonic decay.

In ParaView, colouring by E_intensity reveals the three-peak pattern. The phase field (E_real
or E_imag) shows reversed phase progression inside the slab, demonstrating backward phase
velocity — the hallmark of negative-index propagation.

## Key Takeaways

- The flat-lens effect requires the 1/mu_r weighted diffusion (ADMatDiffusion), not the standard
  Helmholtz. Both n = +1 and n = -1 give n^2 = 1, so the sign of refraction is only
  distinguishable through the 1/mu_r interface condition.
- A negative diffusion coefficient (inv_mu_r = -1 in the slab) is physically meaningful and
  MOOSE handles it correctly; Newton + LU converges without modification.
- Loss regularisation (Im(eps_r) = delta = 0.3) prevents a near-singular system and adds
  cross-coupling between E_real and E_imag via ADMatCoupledForce kernels.
- The BodyForce source acts on E_imag only, because the phasor of a real current J_z is -j omega mu_0 J
  which has only an imaginary part.
- Impedance matching (Z = sqrt(mu_r/eps_r) = 1 for eps_r = mu_r = -1) gives zero reflection at
  the vacuum-LHM interface in the lossless limit — all source energy enters the slab.
- The three-peak intensity pattern along the optical axis at x = {-2, 0, +2} is the visual
  signature of Veselago focusing for propagating waves, as predicted by geometric optics with n = -1.
