# Case 34: Thermal Noise Relaxation — Fluctuation-Dissipation

## Overview

In thermal equilibrium a system at temperature T_eq is not perfectly still —
microscopic fluctuations perturb every degree of freedom with energy ~k_B T.
The **fluctuation-dissipation theorem** links these equilibrium fluctuations to
the dissipation rate: the same mechanism that causes relaxation also sustains the
fluctuations. A concrete consequence is that high-wavenumber (short-wavelength)
fluctuations decay much faster than low-wavenumber ones — mode decay rate scales
as k².

This case demonstrates this principle by initialising a 2D heat diffusion problem
with random noise and watching the spectrum evolve. Short-wavelength speckle
disappears first; only the long-wavelength fundamental mode survives at late time.
The connection to the fluctuation-dissipation theorem is discussed in Haus,
*Electromagnetic Noise and Quantum Optical Measurements* (Springer, 2000), Ch. 5.

---

## The Physics

### Governing Equations

The 2D heat equation on the unit square with equilibrium temperature T_eq = 0.5:

```
∂T/∂t = D·∇²T        D = 0.1

T = 0.5  on all walls  (equilibrium temperature BC)
T(x,y,0) = 0.5 + noise(x,y)    (random IC from RandomIC)
```

### Modal Decay

The solution can be expanded in eigenfunctions of the Laplacian:

```
T(x,y,t) − T_eq = Σ_{m,n} A_{mn} · sin(mπx) · sin(nπy) · exp(−λ_{mn}·t)

λ_{mn} = D·π²·(m² + n²)
```

The lowest mode (m=n=1) decays at λ₁₁ = D·π²·2 ≈ 1.97 s⁻¹. Higher modes decay
at rates proportional to m² + n²:

| Mode | λ_{mn} (D=0.1) | Decay time |
|------|----------------|------------|
| (1,1) | 1.97 s⁻¹   | 0.51 s     |
| (2,1) | 4.93 s⁻¹   | 0.20 s     |
| (2,2) | 7.90 s⁻¹   | 0.13 s     |
| (3,3) | 17.8 s⁻¹   | 0.056 s    |

By t = 0.5 s only the (1,1) mode is appreciable; by t = 2.0 s the domain has
relaxed to uniform T = 0.5.

---

## MOOSE Implementation

| Object | Role |
|--------|------|
| `GeneratedMesh` 40×40 QUAD4 | Unit square with fine resolution for noise |
| `ADTimeDerivative` | ∂T/∂t |
| `ADMatDiffusion` (D=0.1) | Heat diffusion D·∇²T |
| `ADDirichletBC` T=0.5 on all walls | Equilibrium wall temperature |
| `RandomIC` (mean=0.5, std=0.2) | Random initial temperature field |
| `ADGenericConstantMaterial` | Constant diffusivity D=0.1 |
| `ElementAverageValue` avg_T | Tracks mean (should stay at ~0.5) |
| `ElementExtremeValue` max_T, min_T | Track noise amplitude decay |
| `Transient` executioner, dt=0.01, end=2.0 | Time integration |

`RandomIC` sets a spatially uncorrelated random initial condition at each node
by sampling from a specified distribution. Because all Fourier modes are excited
equally by white noise, the simulation directly demonstrates the mode-selective
relaxation predicted by the fluctuation-dissipation theorem.

---

## Expected Results

| Time | Physical observation |
|------|---------------------|
| t = 0.0 | Speckled field — all spatial scales present |
| t = 0.1 | Fine-scale (high-k) features smoothed out |
| t = 0.5 | Only the broad (1,1) arch remains |
| t = 2.0 | Domain approaches uniform T = T_eq = 0.5 |

Quantitative checks from the CSV output:
- `avg_T` remains close to 0.5 at all times (conservation in mean).
- `max_T − min_T` decreases monotonically, tracking the amplitude of the
  surviving lowest mode as exp(−λ₁₁·t) ≈ exp(−1.97·t).
- The decay of `max_T − avg_T` on a log-linear plot should approach a straight
  line with slope −1.97 s⁻¹ at late times when only the (1,1) mode remains.

---

## How to Run

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case34-thermal-noise \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case34_thermal_noise.i 2>&1 | tail -20'
```

Output files produced:
- `case34_thermal_noise_out.e` — Exodus file with T(x,y,t) at all timesteps
- `case34_thermal_noise_out.csv` — avg_T, max_T, min_T vs. time

Note: `RandomIC` uses a fixed random seed by default. To explore different
realisations, add `seed = <integer>` to the `RandomIC` block.
