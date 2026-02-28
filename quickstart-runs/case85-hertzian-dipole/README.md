# Case 85: Hertzian Dipole Radiation Pattern

## Overview

This case simulates a **Hertzian (infinitesimal) dipole** antenna using a 2D axisymmetric (RZ) frequency-domain Helmholtz formulation. The goal is to verify the characteristic **sin²θ radiation pattern** of a z-directed electric dipole.

The dipole is modelled as a Gaussian current source at the origin of an RZ domain, and the radiation pattern is extracted by sampling the field intensity at multiple angles on a circle of radius R = 4λ₀.

## Physics

### The Hertzian Dipole

A Hertzian dipole is a current element of infinitesimal length Δl carrying current I·exp(+jωt) along the z-axis. It is the canonical radiating antenna element. In the far field (kr >> 1), the electric field has only a θ-component:

```
E_θ(R, θ) ∝ (j k₀ I Δl) / (4π R) × sin(θ) × exp(−j k₀ R)
```

The radiated power density (Poynting magnitude) follows the famous pattern:

```
U(θ) = |E_θ|² ∝ sin²(θ)
```

where θ is the polar angle measured from the dipole (z) axis:
- **θ = 0°**: zero radiation (along the dipole axis)
- **θ = 90°**: maximum radiation (in the equatorial plane)

### RZ Axisymmetric Formulation

The problem has azimuthal symmetry (∂/∂φ = 0). In MOOSE's RZ coordinate system:
- x-axis → r (radial, r ≥ 0)
- y-axis → z (axial)

The field variable represents the z-component of the magnetic vector potential A_z, which satisfies the **scalar cylindrical Helmholtz equation**:

```
(1/r)∂/∂r(r ∂A_z/∂r) + ∂²A_z/∂z² + k₀² A_z = −μ₀ J_z(r, z)
```

MOOSE's standard `Diffusion` kernel automatically computes the cylindrical Laplacian when `coord_type = RZ` is set in the `[Problem]` block.

### Complex Phasor Splitting

Writing A_z = E_real + j·E_imag and splitting into real/imaginary parts:

```
∇²_cyl E_real + k₀² E_real = 0                  (no bulk source)
∇²_cyl E_imag + k₀² E_imag = −J(r, z)           (Gaussian source)
```

The two equations couple only through the absorbing Robin BCs on the outer boundaries.

## Setup

### Domain

| Parameter | Value | Notes |
|-----------|-------|-------|
| Radial extent R_max | 6 λ₀ = 6.0 | r ∈ [0, 6] |
| Axial half-extent Z_max | 6 λ₀ = 6.0 | z ∈ [−6, +6] |
| Wavelength λ₀ | 1.0 (normalised) | — |
| Wavenumber k₀ | 2π ≈ 6.283 rad/unit | — |
| Mesh | 60 × 120 QUAD4 | 10 elements/λ₀ |

### Source

A Gaussian approximation to the point dipole source:

```
J(r, z) = 100 × exp(−(r² + z²) / 0.005)
```

- Centre: (r=0, z=0) — the origin
- Amplitude: A = 100
- Width: σ = 0.05 λ₀ (much smaller than the wavelength)

### Boundary Conditions

| Boundary | BC Type | Notes |
|----------|---------|-------|
| Left (r = 0) | Natural Neumann (automatic) | Symmetry axis — ∂E/∂r = 0 |
| Right (r = 6) | EMRobinBC (absorbing) | First-order ABC, no injection |
| Top (z = +6) | EMRobinBC (absorbing) | First-order ABC |
| Bottom (z = −6) | EMRobinBC (absorbing) | First-order ABC |

**Critical**: Do NOT apply any BC at r = 0. The natural Neumann condition (∂E/∂r = 0) is the physically correct symmetry condition for the m = 0 azimuthal mode.

## Expected Results

### Radiation Pattern at R = 4λ₀

| θ | (r, z) | sin²(θ) | Expected E_intensity ratio |
|---|--------|---------|---------------------------|
| 0° | (0.000, 4.000) | 0.000 | ≈ 0 (small residual from ABC) |
| 30° | (2.000, 3.464) | 0.250 | ≈ 0.25 |
| 45° | (2.828, 2.828) | 0.500 | ≈ 0.50 |
| 60° | (3.464, 2.000) | 0.750 | ≈ 0.75 |
| 90° | (4.000, 0.000) | 1.000 | 1.00 (maximum) |

The pattern should be symmetric: θ and (180°−θ) give the same intensity. The symmetry probes at θ = 120° and θ = 150° should match θ = 60° and θ = 30° respectively.

### Deviations from Ideal

- **Near-field effects**: At R = 4λ₀ the far-field approximation is good but not perfect. Expect 2–5% deviation from ideal sin²θ.
- **ABC reflections**: The first-order absorbing BC has non-zero reflection for non-normal incidence. At 6λ₀ domain size, expect ~5% ABC error.
- **θ = 0 null**: The on-axis null is finite (not zero) due to the finite source width σ = 0.05λ₀. Expect a residual of ~1% of the equatorial intensity.

## Running the Simulation

```bash
# Inside Docker with combined-opt:
cd /path/to/case85-hertzian-dipole
/opt/moose/bin/combined-opt -i case85_hertzian_dipole.i

# Expected outputs:
#   case85_hertzian_dipole.e    — Exodus field file (RZ half-plane)
#   case85_hertzian_dipole.csv  — Postprocessor values (radiation pattern)
```

## Visualisation in ParaView

1. Open `case85_hertzian_dipole.e`
2. Colour by `E_intensity` — the doughnut/torus lobe shape should be visible
3. Apply **Reflect** filter (normal = [−1, 0, 0]) to show both positive and negative r half-planes
4. For 3D: apply **Rotational Extrusion** around the z-axis to generate the full torus radiation pattern
5. Use **Plot Over Arc** from (0, 4, 0) to (4, 0, 0) to extract the angular pattern

## Post-Processing in Python

```python
import numpy as np
import pandas as pd

df = pd.read_csv('case85_hertzian_dipole_0001.csv')

# Extract radiation pattern values
thetas = [0, 30, 45, 60, 90]
keys   = ['E_intensity_theta0', 'E_intensity_theta30',
          'E_intensity_theta45', 'E_intensity_theta60',
          'E_intensity_theta90']

I_max = float(df[keys[-1]])  # maximum at θ = 90°

print("θ [deg]  |E|²_norm   sin²(θ)   error [%]")
for theta, key in zip(thetas, keys):
    I_norm = float(df[key]) / I_max
    I_ideal = np.sin(np.radians(theta))**2
    err = 100 * abs(I_norm - I_ideal) / max(I_ideal, 1e-6)
    print(f"  {theta:3d}    {I_norm:.4f}     {I_ideal:.4f}    {err:.1f}")

# Check symmetry: θ_120 ≈ θ_60, θ_150 ≈ θ_30
I_60  = float(df['E_intensity_theta60'])
I_120 = float(df['E_intensity_theta120'])
print(f"\nSymmetry check: I(60°)/I(120°) = {I_60/I_120:.4f} (expected 1.0)")
```

## Coordinate System Notes

In MOOSE RZ notation:
- `x` in the input file → `r` in physics (radial coordinate, x ≥ 0)
- `y` in the input file → `z` in physics (axial coordinate)
- Point coordinates in postprocessors: `'r z 0'` → e.g., `'4.0 0.0 0'` for equatorial point

The polar angle θ from the dipole axis satisfies:
```
sin(θ) = r / sqrt(r² + z²)
cos(θ) = z / sqrt(r² + z²)
```

## Key Technical Details

### coord_type = RZ in [Problem]

This single setting makes MOOSE modify the integration measure and differential operators for cylindrical geometry:
- Volume element: dV = r · dx · dy (where x=r, y=z)
- Diffusion kernel: computes (1/r)∂/∂r(r ∂E/∂r) + ∂²E/∂z² automatically
- No special kernels are needed for the cylindrical Laplacian

### Why No BC at r = 0?

For the m = 0 (azimuthally symmetric) mode driven by a z-directed dipole, A_z(r, z) is even in r: A_z(r, z) = A_z(−r, z). Therefore ∂A_z/∂r = 0 at r = 0 — this is the natural Neumann condition, which MOOSE applies automatically. Applying a Dirichlet condition E = 0 at r = 0 would be physically wrong for this mode.

### Source Convention

The Gaussian BodyForce acts on E_imag only. This follows the phasor convention where the source −jωμ₀J has only an imaginary component (with +jωt time convention). E_real is driven indirectly through the absorbing BCs at the domain boundaries, which couple E_real ↔ E_imag via the factor j in the impedance term.

## Relation to Other Cases

- **Case 77** (cylinder scattering): 2D Cartesian, scattered field from PEC cylinder
- **Case 79** (Snell's law): 2D Cartesian, refraction at dielectric interface
- **Case 83** (Veselago lens): 2D Cartesian, point source focusing by LHM slab
- **Case 85** (this case): 2D RZ axisymmetric, radiation pattern from dipole source

The RZ formulation is unique among the Batch F cases — it provides an efficient 2D simulation of an inherently 3D radiation pattern problem.

## References

- Kong, J.A., *Electromagnetic Wave Theory*, 2nd ed., EMW Publishing (1990), Chapter 2
- Balanis, C.A., *Antenna Theory: Analysis and Design*, 3rd ed., Wiley (2005), Chapter 4
- Pozar, D.M., *Microwave Engineering*, 4th ed., Wiley (2012), Sections 9.2–9.3
- MIT 6.635 Advanced Electromagnetism, Spring 2003 — Prof. Jin Au Kong
  https://ocw.mit.edu/courses/6-635-advanced-electromagnetism-spring-2003/
