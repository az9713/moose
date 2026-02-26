# Case 35: Dispersive Pulse Broadening in an Optical Fiber

## Overview

A short optical pulse launched into a fiber does not retain its shape. Different
frequency components of the pulse travel at slightly different group velocities
— a phenomenon called **group velocity dispersion (GVD)**. The result is pulse
broadening: the peak amplitude drops while the pulse duration increases, with the
total energy (pulse area squared integrated) conserved.

This case models the pulse envelope equation using advection-diffusion: the
advection term translates the pulse at group velocity v_g, while the diffusion
term (with coefficient equal to the GVD parameter) broadens it. The broadening
law is exactly Gaussian: w(z) = w₀√(1 + (z/z_d)²), where z_d is the dispersion
length. This result is derived in Haus, *Electromagnetic Noise and Quantum
Optical Measurements* (Springer, 2000), Ch. 4.

---

## The Physics

### Governing Equations

The pulse envelope amplitude A(x, t) satisfies the advection-diffusion equation:

```
∂A/∂t + v_g·∂A/∂x = D_gvd·∂²A/∂x²
```

Here time t plays the role of propagation distance z along the fiber, x is the
local time coordinate within the pulse frame, and the diffusion coefficient
D_gvd is the GVD parameter β₂/2 in standard fiber optics notation.

Parameters:

| Symbol  | Value | Meaning |
|---------|-------|---------|
| v_g     | 1.0   | Group velocity (translates pulse along x) |
| D_gvd   | 0.01  | GVD broadening coefficient |
| w₀      | 0.1   | Initial Gaussian half-width (1/e radius) |
| z_d     | w₀²/(2D_gvd) = 0.5 | Dispersion length |

### Initial Condition

A Gaussian pulse centred at x₀ = 0.5:

```
A(x, 0) = exp(−(x − 0.5)² / w₀²)      w₀ = 0.1
```

### Analytical Solution

The exact solution for the advection-diffusion equation with Gaussian IC:

```
A(x, t) = (w₀/w(t)) · exp(−(x − x₀ − v_g·t)² / w(t)²)

w(t) = w₀ · √(1 + (t/t_d)²)        t_d = w₀²/(2·D_gvd) = 0.5
```

At t = t_d = 0.5 the pulse has broadened by factor √2 and its peak dropped to
1/√2 of the initial value, while remaining Gaussian.

---

## MOOSE Implementation

| Object | Role |
|--------|------|
| `GeneratedMesh` 1D, 200 elements, [0, 3] | 1D domain long enough for pulse travel |
| `ADTimeDerivative` | ∂A/∂t |
| `ADConservativeAdvection` (v_g=1.0) | v_g·∂A/∂x group velocity transport |
| `ADMatDiffusion` (D_gvd=0.01) | GVD broadening D_gvd·∂²A/∂x² |
| `FunctionIC` (Gaussian) | Initial pulse A(x,0) = exp(-(x-0.5)²/0.01) |
| `ADDirichletBC` A=0 at x=0, x=3 | Zero-field domain boundaries |
| `ElementExtremeValue` max_A | Tracks peak amplitude decay |
| `Transient` executioner, dt=0.01, end=1.0 | Time integration |

`ADConservativeAdvection` implements the advective flux in conservative form
∂(v_g·A)/∂x, which equals v_g·∂A/∂x for constant v_g. It requires an upwind
scheme; the `upwinding_type = full` option suppresses numerical oscillations.

---

## Expected Results

| Time | Centre position | Peak amplitude | Width w(t) |
|------|----------------|----------------|------------|
| t = 0.0 | 0.50 | 1.000 | 0.100 |
| t = 0.5 | 1.00 | 0.707 | 0.141 (=w₀√2) |
| t = 1.0 | 1.50 | 0.447 | 0.224 (=w₀√5) |

Key observations from the CSV output:
- The pulse centre moves at x_centre = 0.5 + v_g·t = 0.5 + t.
- `max_A` decays as 1/√(1 + (t/t_d)²) — not exponential, but algebraic.
- The integral ∫A dx (total area) is approximately conserved (slight loss only
  from boundaries when pulse approaches the domain edge).

---

## How to Run

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case35-dispersive-pulse \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case35_dispersive_pulse.i 2>&1 | tail -20'
```

Output files produced:
- `case35_dispersive_pulse_out.e` — Exodus file with A(x,t) at all timesteps
- `case35_dispersive_pulse_out.csv` — max_A and avg_A vs. time

This case sets the stage for Case 36, which adds Kerr nonlinearity to balance
the GVD broadening and produce a soliton that propagates without changing shape.
