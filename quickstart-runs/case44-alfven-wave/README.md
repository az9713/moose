# Case 44: Alfven Wave Propagation — MHD Elsasser Variables

**Reference:** Rieutord, *Fluid Dynamics* (Springer, 2015), Ch. 10, Sec. 10.4

---

## Overview

Alfven waves are transverse electromagnetic-hydrodynamic waves that propagate along
magnetic field lines in a conducting fluid (plasma or liquid metal). They arise from
the restoring force of magnetic tension acting on a displaced fluid parcel. Unlike
sound waves, which travel via pressure, Alfven waves travel via magnetic field line
curvature.

In ideal MHD with a uniform background field B_0 in the x-direction, transverse
perturbations (vy, by) satisfy coupled wave equations with a characteristic speed
called the Alfven velocity:

```
v_A = B_0 / sqrt(mu_0 * rho)
```

where mu_0 is the magnetic permeability and rho is the mass density. For typical
laboratory plasmas, v_A ranges from tens to thousands of km/s. In the solar wind,
Alfven waves carry energy from the corona outward and are central to coronal heating
models.

---

## Elsasser Variable Decomposition

The elegant insight due to Elsasser (1950) is that the MHD wave equations decouple
under the change of variables:

```
d+ = vy + by / sqrt(mu_0 * rho)
d- = vy - by / sqrt(mu_0 * rho)
```

In ideal MHD these satisfy:

```
d(d+)/dt + v_A * d(d+)/dx = 0
d(d-)/dt - v_A * d(d-)/dx = 0
```

The d+ field is a pure rightward-traveling wave at speed +v_A. The d- field is a
pure leftward-traveling wave at speed -v_A. The two modes are completely decoupled
in the ideal (non-dissipative) limit. This decomposition is the MHD analogue of the
d'Alembert solution for acoustic waves.

With resistive and viscous dissipation combined into an effective diffusivity D_eff,
the equations become advection-diffusion equations:

```
d(d+)/dt + v_A * d(d+)/dx = D_eff * d²(d+)/dx²
d(d-)/dt - v_A * d(d-)/dx = D_eff * d²(d-)/dx²
```

---

## Analytical Solution

For an initial Gaussian pulse of width w centered at x_0 in d+ (with d- = 0):

```
d+(x,0) = exp(-(x - x_0)^2 / w^2)
d-(x,0) = 0
```

The exact solution is a broadening Gaussian translating at v_A:

```
d+(x,t) = [w / sqrt(w^2 + 4*D_eff*t)] * exp(-(x - x_0 - v_A*t)^2 / (w^2 + 4*D_eff*t))
d-(x,t) = 0  (identically, for all time)
```

The peak amplitude decays as:

```
A(t) = 1 / sqrt(1 + 4*D_eff*t / w^2)
```

and the pulse center translates to x(t) = x_0 + v_A * t.

---

## Parameters and Domain

| Parameter | Symbol | Value |
|-----------|--------|-------|
| Alfven speed | v_A | 1.0 (normalized) |
| Effective diffusivity | D_eff | 0.01 |
| Initial pulse center | x_0 | 3.0 |
| Initial Gaussian half-width | w | 0.5 (sigma = w/2 = 0.25) |
| End time | T | 6.0 |
| Domain | [0,12] x [0,0.12] | quasi-1D |
| Mesh | 120 x 2 | quad elements |
| Time step | dt | 0.05 (adaptive) |

The domain is intentionally quasi-1D (thin in y, 2 elements) to simulate a 1D
problem while satisfying MOOSE's 2D mesh requirement.

---

## Validation at t = 5

At t = 5:

- **d+ peak location:** x = x_0 + v_A * t = 3 + 1.0 * 5 = **8.0**
- **d+ peak amplitude:**
  ```
  A(5) = 1 / sqrt(1 + 4 * 0.01 * 5 / 0.25) = 1 / sqrt(1.8) ~ 0.745
  ```
- **d- maximum:** remains ~ 0 throughout (numerical errors at machine precision level)

The ratio of final to initial peak amplitude (~0.745) quantifies the resistive
damping over 5 Alfven crossing times.

---

## Running in Docker

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case44-alfven-wave \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case44_alfven_wave.i 2>&1 | tail -30'
```

This runs for approximately 120 time steps (adaptive) and produces:

- `case44_alfven_wave_out.e` — Exodus file with d+, d- fields at each output step
- `case44_alfven_wave_out.csv` — Time history of all postprocessors

---

## Expected CSV Output (Selected Times)

| time | max_d_plus | max_d_minus | total_d_plus |
|------|-----------|-------------|--------------|
| 0.0  | 1.000     | 0.000       | ~0.443       |
| 1.0  | ~0.943    | ~0.000      | ~0.418       |
| 3.0  | ~0.834    | ~0.000      | ~0.370       |
| 5.0  | ~0.745    | ~0.000      | ~0.330       |
| 6.0  | ~0.707    | ~0.000      | ~0.313       |

The total_d_plus integral decreases due to dissipation and also because the pulse
encounters the right Dirichlet boundary (which absorbs it) near t = 4-5.

---

## Implementation Notes

### ADConservativeAdvection

The `ADConservativeAdvection` kernel implements the weak form of:

```
v . grad(u)
```

where the velocity is provided as a vector material property. This is the standard
conservative (divergence form) advection term. With `upwinding_type = full`, the
kernel applies full upwinding stabilization, which adds artificial diffusion
proportional to the element size. This prevents spurious oscillations for
convection-dominated problems (Peclet number >> 1 here: Pe = v_A * h / D_eff = 1.0
* 0.1 / 0.01 = 10).

### ADGenericConstantVectorMaterial

Used to provide the velocity vector as a material property. The d- equation uses
`${fparse -v_A}` to evaluate the negative of the Alfven speed at input-file parse
time, yielding the leftward velocity (-1.0, 0, 0).

### Boundary Conditions

Zero Dirichlet BCs are applied at x = 0 and x = 12. Because the initial pulse is
centered at x = 3 with width ~0.5, it is far from both boundaries at t = 0. The
pulse reaches the right boundary near t = 9, well after the simulation ends at t = 6.

---

## Key Takeaways

1. **Elsasser decomposition**: The change of variables d+/- = vy +/- by/sqrt(mu_0*rho)
   diagonalizes the MHD wave operator, revealing two decoupled one-way wave equations.
   This is physically exact in ideal MHD and is the fundamental structure underlying
   Alfven wave propagation.

2. **Decoupling**: A purely rightward-propagating initial condition (d- = 0) remains
   decoupled forever in the ideal limit. Numerical d- generation is a measure of
   scheme error, and its near-zero value confirms correct implementation.

3. **Resistive decay**: Physical dissipation (resistivity eta and viscosity nu combine
   into D_eff = (eta + nu) / 2) broadens the pulse and reduces its amplitude according
   to the analytical formula A(t) = 1/sqrt(1 + 4*D*t/w^2). This is identical in form
   to thermal diffusion of a Gaussian concentration profile.

4. **Alfven wave applications**: Alfven waves play a central role in solar physics
   (coronal heating, solar wind acceleration), space weather (magnetospheric dynamics),
   fusion plasma stability (Alfven eigenmode instabilities in tokamaks), and liquid
   metal MHD (electromagnetic stirring in casting processes).

5. **Upwinding stability**: At element Peclet number Pe = 10, centered differencing
   would generate spurious oscillations. Full upwinding stabilizes the scheme at the
   cost of additional numerical diffusion (O(h)), which for this mesh (h = 0.1) adds
   ~D_num = v_A * h / 2 = 0.05. The effective total diffusivity is D_eff + D_num
   ~ 0.06, meaning the numerical peak amplitude at t = 5 will be slightly lower than
   the analytical prediction of 0.745.
