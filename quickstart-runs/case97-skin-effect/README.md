# Case 97: Skin Effect — Magnetic Diffusion into a Conducting Slab

## Overview

When an alternating magnetic field is applied to the surface of a conductor, it does not penetrate uniformly. Eddy currents induced by the changing flux oppose the field inside the conductor, confining most of the current — and most of the field — to a thin layer near the surface. This is the skin effect, and the thickness of that layer is the skin depth delta = sqrt(2 / (mu sigma omega)). At power-line frequencies (50-60 Hz) in copper, the skin depth is about 9 mm; at radio frequencies (1 MHz) it shrinks to 66 micrometers; at microwave frequencies it is measured in nanometers. The skin effect is the reason transformer cores are laminated, RF conductors are hollow tubes, and microwave cavities need only a micron-thick metal plating.

MIT 6.641 Lecture 9 derives the skin effect from the magnetic diffusion equation: Faraday's law and Ampere's law (in the low-frequency, displacement-current-free limit) combine to give a parabolic PDE for the magnetic field identical in form to the heat equation. The magnetic diffusivity is D_m = 1/(mu sigma), and the skin depth is related to D_m and frequency by delta = sqrt(2 D_m / omega).

This case simulates a quasi-1D conducting slab of unit length with D_m = 0.01 and omega = 2 pi, giving a skin depth of about 0.056 m. The slab is roughly 18 skin depths thick, so the field is fully attenuated before reaching the far wall. Starting from a field-free initial condition, the simulation runs for five periods (t = 5 s) until the initial transient has decayed and the oscillation has settled into sinusoidal steady state. The `FunctionDirichletBC` with a `ParsedFunction` drives the sinusoidal boundary, and `PointValue` postprocessors at x = delta, 2 delta, and 3 delta verify the exponential amplitude profile.

---

## The Physics

**Governing equation**

In the low-frequency (magneto-quasi-static) limit, displacement current is negligible and Faraday's and Ampere's laws give:

```
dH/dt = D_m * d^2 H / dx^2     where D_m = 1 / (mu sigma)
```

This is the magnetic diffusion equation — formally identical to the heat equation.

**Sinusoidal steady-state solution**

With the boundary condition H(0, t) = H0 sin(omega t), the sinusoidal steady-state solution is:

```
H(x, t) = H0 * exp(-x / delta) * sin(omega t - x / delta)
```

The field decays exponentially with depth, and also accumulates a phase lag of x/delta radians. At depth x = delta both the amplitude is 1/e ≈ 0.368 H0 and the phase lag is one radian.

**Skin depth**

```
delta = sqrt(2 D_m / omega) = sqrt(2 * 0.01 / (2 pi)) ≈ 0.0564 m
```

**Boundary conditions**

| Boundary | Condition | Physical meaning |
|----------|-----------|-----------------|
| Left (x = 0) | H = sin(omega t) | Oscillating applied field |
| Right (x = 1) | H = 0 | Far-field (18 skin depths from surface) |
| Top and bottom | Natural (zero flux) | Quasi-1D: no y-variation |

**Initial condition**

H = 0 everywhere — the conductor starts field-free. The initial transient (which has a different spatial profile than steady state) decays over roughly one magnetic diffusion time D_m / L^2 ≈ 0.01 s, much less than the period T = 1 s.

**Domain and discretization**

- 2D rectangle: x in [0, 1], y in [0, 0.04]
- 100 elements in x (fine resolution along diffusion direction), 2 elements in y (quasi-1D)

**Parameters**

| Symbol | Value | Meaning |
|--------|-------|---------|
| D_m | 0.01 m^2/s | Magnetic diffusivity 1/(mu sigma) |
| omega | 2 pi rad/s | Angular frequency (T = 1 s) |
| delta | 0.0564 m | Skin depth |
| L | 1.0 m | Slab thickness (= 17.7 delta) |
| dt | 0.02 s | Time step (50 steps per period) |
| end_time | 5.0 s | Duration (5 periods) |

---

## Input File Walkthrough

**Global parameters**

```
D_m   = 0.01
omega = 6.2831853   # 2 pi
```

The skin depth is not set explicitly — it follows from the physics as delta = sqrt(2 * D_m / omega).

**[Mesh]**

`GeneratedMesh` (the legacy shorthand) produces a 100 x 2 grid. The fine x-resolution (100 elements over 1 m, so dx = 0.01 m ≈ delta/5.6) is needed to resolve the exponential decay. The two elements in y are the minimum for a 2D mesh while keeping the problem effectively 1D.

**[Variables]**

`H` is the magnetic field in A/m (normalized with H0 = 1). It is initialized to zero.

**[Kernels]**

```
[H_time]  type = ADTimeDerivative   — dH/dt
[H_diff]  type = ADMatDiffusion     — D_m * d^2 H/dx^2
                  diffusivity = diffusivity
```

`ADMatDiffusion` reads the material property named `diffusivity` at each quadrature point.

**[BCs]**

```
[H_oscillating]  FunctionDirichletBC  boundary = left   function = H_applied
[H_far]          DirichletBC          boundary = right  value = 0
```

`H_applied` is the `ParsedFunction` `sin(omega * t)`. The top and bottom boundaries receive zero-flux Neumann by default, consistent with the quasi-1D geometry.

**[Materials]**

```
[mag_diff]  ADGenericConstantMaterial  prop_names = 'diffusivity'  prop_values = '0.01'
```

A single spatially uniform material provides D_m = 0.01 to both kernels.

**[Postprocessors]**

| Name | Point | Expected amplitude at steady state |
|------|-------|-----------------------------------|
| H_surface | (0.0, 0.02) | 1.0 (surface) |
| H_at_skin_depth | (0.056, 0.02) | exp(-1) ≈ 0.368 |
| H_at_2delta | (0.113, 0.02) | exp(-2) ≈ 0.135 |
| H_at_3delta | (0.169, 0.02) | exp(-3) ≈ 0.050 |
| avg_H | domain average | approaches 0 (symmetric oscillation) |

The postprocessors are measured at the mid-height of the slab (y = 0.02) to avoid edge effects from the top and bottom boundaries.

**[Executioner]**

Transient NEWTON with BoomerAMG algebraic multigrid. Running for 5 s (250 steps) gives the solution time to settle into sinusoidal steady state. The initial magnetic diffusion transient decays within the first period.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case97-skin-effect \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case97_skin_effect.i 2>&1 | tail -30'
```

The 250-step transient takes around 30-60 seconds depending on hardware. Exodus output captures the spatial field profile at every time step; the CSV records the five postprocessors at each step, making it easy to extract the amplitude envelope by taking the maximum over each period.

---

## Expected Results

After the initial transient (roughly the first 1-2 periods), the field enters sinusoidal steady state. Examining H_at_skin_depth in the CSV, the oscillation amplitude should stabilize near exp(-1) ≈ 0.368, with a phase lag of one radian relative to H_surface.

Snapshot of amplitudes at steady state (read from peak values in CSV after t = 3 s):

| Location | x / delta | Expected amplitude |
|----------|-----------|-------------------|
| Surface | 0 | 1.000 |
| x = delta | 1 | 0.368 |
| x = 2 delta | 2 | 0.135 |
| x = 3 delta | 3 | 0.050 |

In the Exodus output, an animation of H versus x shows the characteristic "evanescent oscillation" — a spatial sinusoid whose envelope decays exponentially from the left surface. The wavelength of the spatial oscillation is 2 pi delta ≈ 0.354 m, so roughly three spatial oscillations are visible before the amplitude drops below noise.

---

## Key Takeaways

- The magnetic diffusion equation dH/dt = D_m * nabla^2 H is the governing PDE for skin effect — it has exactly the same mathematical form as the heat equation.
- The skin depth delta = sqrt(2 D_m / omega) sets the characteristic length scale for field penetration; all fields and currents are confined within a few skin depths of the surface.
- The field decays exponentially in amplitude and accumulates a phase lag at the same rate: at depth delta the amplitude is 1/e and the phase lag is 1 radian.
- `FunctionDirichletBC` with a time-dependent `ParsedFunction` is the standard way to drive oscillating boundary conditions in MOOSE transient problems.
- The quasi-1D geometry (100 x 2 mesh) is a practical pattern: it preserves the 2D mesh requirement while concentrating all resolution in the direction of interest.
- Five periods of simulation are needed to see clean sinusoidal steady state because the initial transient decays on a timescale comparable to the period; plotting only the last two periods avoids start-up artifacts.
