# Case 31: Driven Resonant Cavity — Frequency Response and Q

## Overview

When a cavity is driven by an internal source at a frequency near one of its
resonant modes, the field amplitude becomes very large — a phenomenon directly
analogous to mechanical resonance. Far from resonance the same source produces
only a weak, distorted field. The ratio of on-resonance to off-resonance
response characterises the cavity's quality factor Q.

This case solves the inhomogeneous Helmholtz equation for a 2D PEC cavity
driven by a Gaussian source. By comparing solutions at k² = 12.3 (near the
TM₁₁ resonance at k² ≈ 12.337) and k² = 15.0 (off-resonance), the resonant
amplification is directly observable. The physics follows Haus,
*Electromagnetic Noise and Quantum Optical Measurements* (Springer, 2000), Ch. 3.

---

## The Physics

### Governing Equations

The time-harmonic electric field in a 2D PEC cavity with a localised source:

```
∇²E + k²E = -J_source(x, y)      in Ω

E = 0                              on ∂Ω  (PEC walls)
```

where k = ω/c is the free-space wavenumber and J_source is a Gaussian source
centred at (x₀, y₀) = (0.7, 0.4):

```
J_source = exp(-((x-0.7)² + (y-0.4)²) / 0.01)
```

The source is intentionally off-centre so that it overlaps with the TM₁₁ mode
shape (which has no zero at (0.7, 0.4)), guaranteeing a nonzero coupling coefficient.

### Resonant Amplification

The field amplitude scales as 1/(k_res² - k²) near resonance. At k² = 12.3
(separation |Δk²| = 0.037 from TM₁₁) the denominator is small and the field
is large. At k² = 15.0 (separation |Δk²| = 2.663) the denominator is large and
the field is weak. The peak amplitude ratio gives a measure of the cavity Q.

### Parameters

| Symbol | Value | Meaning |
|--------|-------|---------|
| a      | 2.0   | Cavity width  |
| b      | 1.0   | Cavity height |
| k² (near resonance)  | 12.3  | Close to TM₁₁ at 12.337 |
| k² (off-resonance)   | 15.0  | Between TM₁₁ and TM₂₁  |

---

## MOOSE Implementation

| Object | Role |
|--------|------|
| `GeneratedMesh` 40×20 QUAD4 | 2:1 unit rectangle (a=2, b=1) |
| `Diffusion` kernel | Laplacian ∇²E |
| `CoefReaction` (coeff=-k²) | Helmholtz reaction term -k²E |
| `BodyForce` with Gaussian function | Localised source -J_source |
| `DirichletBC` on all walls | PEC boundary condition E=0 |
| `ElementExtremeValue` | Peak field amplitude at each k² |
| `Steady` executioner | Direct solve (no time-marching) |

The same input file is run twice with different values of the `k_sq` top-level
variable. A command-line override (`k_sq=15.0`) selects the off-resonance case
without editing the file.

---

## Expected Results

| Case | k² | max_E (expected) | Physical interpretation |
|------|----|------------------|-------------------------|
| Near resonance  | 12.3 | large (>>1) | TM₁₁ mode strongly excited |
| Off-resonance   | 15.0 | small (~0.1) | Weak, distorted response  |

At k² = 12.3 the field pattern closely resembles the TM₁₁ mode shape — a single
arch across the cavity — because the resonant mode dominates the response. At
k² = 15.0 the field is weak with no clear modal structure; the source drives a
superposition of all modes, none of which is strongly excited.

The ratio max_E(near) / max_E(off) provides a coarse estimate of Q/(2Δk/k_res).

---

## How to Run

```bash
# Near-resonance run (k² = 12.3, default)
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case31-driven-cavity \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case31_driven_cavity.i 2>&1 | tail -10'

# Off-resonance run (k² = 15.0, command-line override)
  -c '/opt/moose/bin/combined-opt -i case31_driven_cavity.i k_sq=15.0 2>&1 | tail -10'
```

Output files produced:
- `case31_driven_cavity_out.e` — Exodus file with E field
- `case31_driven_cavity_out.csv` — max_E postprocessor value
