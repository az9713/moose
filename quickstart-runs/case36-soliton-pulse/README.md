# Case 36: Soliton Pulse Propagation — Nonlinear Balance

## Overview

In Case 35 GVD broadens every pulse regardless of amplitude. But when a Kerr
nonlinearity is present (the refractive index depends on intensity: n = n₀ + n₂|A|²),
the phase modulation it produces can exactly cancel the dispersive broadening.
The result is a **soliton**: a self-consistent pulse that propagates undistorted
because the nonlinear self-focusing precisely balances the linear GVD spreading.

This case adds a cubic nonlinear term −α·A³ to the advection-diffusion equation
of Case 35. At the **soliton balance condition** (α·A₀²·w₀² = 2D_gvd) the pulse
propagates with constant shape. Away from balance, the pulse either broadens
(too little nonlinearity) or compresses and oscillates (too much). The soliton
physics is treated in Haus, *Electromagnetic Noise and Quantum Optical
Measurements* (Springer, 2000), Ch. 10.

---

## The Physics

### Governing Equations

The nonlinear Schrödinger-type equation for the pulse envelope A(x, t):

```
∂A/∂t + v_g·∂A/∂x = D·∂²A/∂x² − α·A³
```

The nonlinear term −α·A³ models Kerr self-phase modulation in the anomalous
GVD regime. It acts as a nonlinear gain for the pulse peak and a suppression for
the wings, opposing the spreading caused by D·∂²A/∂x².

Parameters:

| Symbol | Value | Meaning |
|--------|-------|---------|
| v_g    | 1.0   | Group velocity |
| D      | 0.05  | GVD coefficient |
| α      | 0.1   | Kerr nonlinearity (soliton balance) |
| w₀     | 1.0   | Initial pulse width |
| A₀     | 1.0   | Initial peak amplitude |

### Soliton Balance Condition

Soliton condition: α·A₀²·w₀² = 2D. With the values above:

```
α·A₀²·w₀² = 0.1 × 1.0 × 1.0 = 0.1    (not yet at balance with D=0.05)
2D = 2 × 0.05 = 0.10                    ✓ balanced
```

The exact soliton initial condition is a hyperbolic secant profile:

```
A(x, 0) = A₀ · sech((x − x₀)/w₀)    x₀ = 1.5
```

which satisfies the stationary soliton equation D·A'' = α·A³ − (v_g²/(4D))·A.

### Three Regimes

| α value | Balance | Behaviour |
|---------|---------|-----------|
| 0.0     | No nonlinearity | Pure GVD broadening (Case 35 result) |
| 0.1     | Exact soliton | Constant peak/width, stable propagation |
| 0.3     | Over-nonlinear | Compression then oscillation/break-up |

---

## MOOSE Implementation

| Object | Role |
|--------|------|
| `GeneratedMesh` 1D, 300 elements, [0, 4] | 1D domain for pulse travel |
| `ADTimeDerivative` | ∂A/∂t |
| `ADConservativeAdvection` (v_g=1.0) | Group velocity transport |
| `ADMatDiffusion` (D=0.05) | GVD broadening |
| `ADMatReaction` with nonlinear material | Nonlinear term −α·A³ |
| `ADPiecewiseLinearInterpolationMaterial` | Tabulated rate = −α·A² vs. A |
| `FunctionIC` (sech profile) | A(x,0) = sech((x-1.5)/1.0) |
| `ADDirichletBC` A=0 at x=0, x=4 | Domain boundaries |
| `ElementExtremeValue` max_A | Tracks soliton peak |
| `Transient` executioner, dt=0.02, end=2.0 | Time integration |

The nonlinear term −α·A³ is implemented via `ADMatReaction` with a material
property equal to −α·A². Because `ADParsedMaterial` requires JIT compilation
(unavailable in Docker), the A-dependent reaction rate is supplied through
`ADPiecewiseLinearInterpolationMaterial` with a table of (A, −α·A²) values
sampled at sufficient resolution. This is the same Docker-compatible pattern
used in Cases 27 and 28 for nonlinear material properties.

---

## Expected Results

| α | max_A at t=0 | max_A at t=1 | max_A at t=2 | Behaviour |
|---|--------------|--------------|--------------|-----------|
| 0.0 | 1.00 | 0.71 | 0.55 | Broadening (algebraic decay) |
| 0.1 | 1.00 | ~1.00 | ~1.00 | Soliton: constant peak |
| 0.3 | 1.00 | ~1.25 | oscillating | Over-nonlinear: compression |

The soliton run (α = 0.1) should show:
- `max_A` remaining near 1.0 throughout, with only small oscillations from
  numerical dispersion.
- Pulse centre advancing at x_centre ≈ 1.5 + v_g·t = 1.5 + t.
- No visible broadening or narrowing of the pulse profile.

The dispersive run (α = 0.0) replicates Case 35 and provides a control for
comparison. Running both from the same input file using the command-line
override `alpha=0.0` versus the default `alpha=0.1` confirms the nonlinear
balance effect directly.

---

## How to Run

```bash
# Soliton run (alpha = 0.1, default)
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case36-soliton-pulse \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case36_soliton_pulse.i 2>&1 | tail -20'

# Dispersive comparison (alpha = 0.0)
  -c '/opt/moose/bin/combined-opt -i case36_soliton_pulse.i alpha=0.0 2>&1 | tail -20'

# Over-nonlinear run (alpha = 0.3)
  -c '/opt/moose/bin/combined-opt -i case36_soliton_pulse.i alpha=0.3 2>&1 | tail -20'
```

Output files produced:
- `case36_soliton_pulse_out.e` — Exodus file with A(x,t) at all timesteps
- `case36_soliton_pulse_out.csv` — max_A and avg_A vs. time

Visualise with `visualize_all.py` or load the Exodus file in ParaView and
animate to see the pulse travelling and (for α = 0.1) maintaining its shape.
