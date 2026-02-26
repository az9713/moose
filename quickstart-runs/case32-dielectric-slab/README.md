# Case 32: EM Wave Reflection from a Dielectric Slab

## Overview

When an electromagnetic wave travelling in free space strikes a dielectric
interface, part of the energy is reflected and part is transmitted. For a slab
of finite thickness both interfaces contribute, and reflected waves from the
two surfaces interfere — producing a reflection coefficient that oscillates with
frequency and slab thickness (the Fabry-Pérot effect).

This case computes the 1D frequency-domain reflection from a lossless dielectric
slab (εᵣ = 4, like glass at microwave frequencies) at 20 MHz. The field is
represented as a complex-valued wave: two real MOOSE variables carry the real
and imaginary parts of E. A Robin port boundary condition injects the incident
wave and absorbs the outgoing reflected wave. The physics follows Haus,
*Electromagnetic Noise and Quantum Optical Measurements* (Springer, 2000), Ch. 1.

---

## The Physics

### Governing Equations

In 1D, the time-harmonic electric field E(x) (complex amplitude) satisfies:

```
d²E/dx² + k₀² εᵣ(x) E = 0
```

where k₀ = ω/c is the free-space wavenumber. The permittivity profile is:

```
εᵣ(x) = 4.0   for 0 ≤ x ≤ 25 m   (dielectric slab)
εᵣ(x) = 1.0   for 25 < x ≤ 75 m  (vacuum)
```

Because MOOSE solves real systems, E is split into real and imaginary parts:
E = E_r + i·E_i. This doubles the number of unknowns but keeps the system real.

### Boundary Conditions

| Boundary | Type | Physical meaning |
|----------|------|-----------------|
| x = 0   | Dirichlet E_r = 0, E_i = 0 | PEC metal wall behind slab |
| x = 75  | Robin port BC | Injects incident wave, absorbs reflected wave |

The Robin (port) BC at the right boundary couples E_r and E_i:

```
dE_r/dx + k₀ E_i = 2k₀ (incident amplitude)
dE_i/dx - k₀ E_r = 0
```

This encodes the impedance condition for a right-travelling incident wave and an
outgoing reflected wave simultaneously.

### Parameters

| Symbol | Value | Meaning |
|--------|-------|---------|
| f      | 20 MHz | Frequency |
| k₀     | 0.4189 rad/m | Free-space wavenumber |
| εᵣ (slab) | 4.0 | Relative permittivity |
| L      | 75 m | Domain length |
| slab   | 0–25 m | Dielectric region |

---

## MOOSE Implementation

| Object | Role |
|--------|------|
| `GeneratedMesh` 1D, 300 elements | 1D domain 0 to 75 m |
| `Diffusion` on E_r, E_i | d²E/dx² for each component |
| `ADMatReaction` on E_r with -k₀²εᵣ | Helmholtz term (real part) |
| `ADMatReaction` on E_i with -k₀²εᵣ | Helmholtz term (imaginary part) |
| `ADGenericFunctionMaterial` | Space-dependent εᵣ(x) via parsed function |
| `EMRobinBC` (or equivalent NeumannBC pair) | Port BC coupling E_r ↔ E_i |
| `DirichletBC` at x=0 | PEC metal wall |
| `ElementAverageValue` on E_r², E_i² | Energy density for reflection coefficient |

The reflection coefficient is computed from the complex field amplitude at the
port: |R|² = |E_reflected|² / |E_incident|².

---

## Expected Results

The Fresnel formula for reflection at a single interface between vacuum and
εᵣ = 4 gives |R| = (√εᵣ - 1)/(√εᵣ + 1) = 1/3 ≈ 0.333. For a finite slab
the reflection oscillates with slab thickness due to interference between the
two surfaces. At slab thickness = λ_slab/2 the reflections cancel (|R| → 0);
at slab thickness = λ_slab/4 they add constructively.

| Observable | Expected value |
|------------|----------------|
| Standing wave pattern in vacuum (x > 25 m) | Visible in E_r or |E|² profile |
| Shorter wavelength in slab (x < 25 m) | λ_slab = λ₀/√εᵣ = λ₀/2 |
| |R| for this slab thickness | Between 0 and 1/3 (depends on thickness/λ) |

---

## How to Run

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case32-dielectric-slab \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case32_dielectric_slab.i 2>&1 | tail -20'
```

Output files produced:
- `case32_dielectric_slab_out.e` — Exodus file with E_r(x) and E_i(x) profiles
- `case32_dielectric_slab_out.csv` — postprocessor values including |R|
