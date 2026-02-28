# Case 92: Aperture Synthesis Interferometer

## Overview

This case simulates the interference pattern created by two coherent point sources ("stars") separated by d = 2 units, observed at a baseline distance R = 12 units away. It demonstrates the fundamental Fourier relationship at the core of aperture synthesis radio astronomy:

```
Fringe spacing: Λ = λ R / d = 1 × 12 / 2 = 6 units
Angular resolution: Δθ ≈ λ / B_max
```

## Physics

**Governing equation** (2D Helmholtz, frequency domain):

```
∇²E + k₀² E = s(x,y)
```

where `s(x,y)` is the superposition of two Gaussian "star" sources:

```
s(x,y) = A·exp(−((x+6)² + (y−1)²) / 2σ²)
        + A·exp(−((x+6)² + (y+1)²) / 2σ²)
```

with A = 50, σ = 0.1, λ = 1, k₀ = 2π.

**Decomposition into real/imaginary parts** (since ε_r = 1 everywhere):

```
∇²E_real + k₀² E_real = s(x,y)    (sources drive the real part)
∇²E_imag + k₀² E_imag = 0         (imaginary part driven by Robin BCs only)
```

**Intensity (observable)**:

```
I(x,y) = E_real² + E_imag²
```

## Van Cittert–Zernike Theorem

For two equal point sources at angles ±δθ/2 from the axis (δθ = d/R = 2/12 = 0.167 rad), the mutual coherence measured at baseline spacing b is:

```
V(b) = 2 cos(π b δθ / λ)
```

The first null of the visibility occurs at `b_null = λ / (2 δθ) = λ R / (2d) = 3` units. Since the observation plane spans y ∈ [−6, 6], this null appears at y = ±3 in the fringe pattern.

## Domain and Mesh

| Parameter | Value |
|-----------|-------|
| Domain | [−8, 8] × [−6, 6] |
| Mesh | 80 × 60 elements |
| Element size | 0.2 units = λ/5 |
| Wavelength λ | 1.0 unit |
| Wavenumber k₀ | 2π = 6.283 rad/unit |

## Source Configuration

| Source | Position | Gaussian σ |
|--------|----------|-----------|
| Star 1 | (−6, +1) | 0.1 unit |
| Star 2 | (−6, −1) | 0.1 unit |
| Separation d | 2 units | — |
| Distance R (to baseline) | 12 units | — |

## Expected Fringe Pattern

At the observation baseline x = 6, the intensity varies as:

```
I(y) ∝ 1 + cos(2π y / Λ)    with Λ = λR/d = 6
```

| y (units) | Phase | Expected intensity |
|-----------|-------|-------------------|
| 0 | 0 | MAXIMUM (constructive) |
| ±1 | π/3 | ~75% of maximum |
| ±2 | 2π/3 | ~25% of maximum |
| ±3 | π | MINIMUM (destructive) |
| ±4 | 4π/3 | ~25% of maximum |
| ±5 | 5π/3 | ~75% of maximum |
| ±6 | 2π | MAXIMUM (constructive) |

## Files

- `case92_interferometer.i` — MOOSE input file (single simulation)

## Running

```bash
# From the Docker container or MOOSE installation:
cd /path/to/case92-interferometer
combined-opt -i case92_interferometer.i

# With MPI (optional, small problem):
mpirun -np 4 combined-opt -i case92_interferometer.i
```

## Output

- `case92_interferometer_out.e` — Exodus file with E_real, E_imag, E_intensity on the 2D mesh
- `case92_interferometer_out.csv` — Fringe samples at 11 points along x = 6, plus diagnostics

## Visualization

Open the Exodus file in ParaView:

1. Apply a **Slice** at x = 6 (y-z plane) to extract the baseline fringe pattern.
2. Plot `E_intensity` along the slice to see the cosine fringe with period Λ = 6.
3. The full 2D `E_intensity` field shows the interference fringes as horizontal bands.
4. Use **WarpByScalar** on E_real to visualize the cylindrical wavefronts.

## Kernels Used

| Kernel | Equation contribution |
|--------|-----------------------|
| `Diffusion` | `+∫ ∇E·∇v dΩ` (strong: `−∇²E`) |
| `ADMatReaction(k0sq)` | `−k₀²∫ E·v dΩ` (strong: `−k₀²E`) |
| `BodyForce(source_fn)` | `−∫ s·v dΩ` (strong: RHS = `+s`) |
| `EMRobinBC` | Absorbing boundary on all four walls |

## Connection to Radio Astronomy

This simulation models the core measurement of a two-element radio interferometer:

- The two "stars" are real astrophysical radio sources (e.g., Cygnus A and Cassiopeia A).
- The "observation plane" at x = 6 represents the baseline vector between two antennas.
- The fringe pattern `I(y)` is the cross-correlation of the signals from the two antennas.
- Fourier-transforming the fringe pattern recovers the source positions and separation.

The Very Large Array (VLA), ALMA, and the Event Horizon Telescope all operate on this principle, using many baseline pairs to sample the visibility function `V(b)` densely enough to reconstruct a high-resolution image via the inverse Fourier transform.

## Reference

Thompson, A.R., Moran, J.M., Swenson, G.W., *Interferometry and Synthesis in Radio Astronomy*, 3rd Ed. (Wiley-VCH, 2017), Chapters 1–2.
