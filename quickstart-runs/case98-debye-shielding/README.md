# Case 98: Debye Shielding — Linearized Poisson-Boltzmann Equation

## Overview

In a plasma or electrolyte, free mobile charges rearrange themselves in response to any external electric field. Positive ions drift away from positive potentials and accumulate near negative potentials, and vice versa for electrons or negative ions. This redistribution partially cancels the original field — an effect called Debye shielding (or Debye-Huckel screening). A test charge placed in a plasma does not produce the long-range 1/r Coulomb potential of free space; instead, the potential is "screened" and decays as exp(-r / lambda_D) / r, where lambda_D is the Debye length. For distances much larger than lambda_D the plasma looks electrically neutral and the test charge is invisible.

The mathematical description comes from combining the Poisson equation for the potential with the Boltzmann distribution for the ion densities. For small potentials (q Phi / k_B T << 1), the Boltzmann exponential linearizes and the result is the screened Poisson equation (also called the Yukawa equation or the linearized Poisson-Boltzmann equation): nabla^2 Phi - Phi / lambda_D^2 = -rho_free / epsilon. MIT 6.641 Lecture 7 derives this in the context of charge relaxation and electrolyte physics. The same equation appears in semiconductor physics as the Thomas-Fermi screening model, in colloid science as the Debye-Huckel component of DLVO theory, and in nuclear physics as the Yukawa meson-exchange potential.

This case is a steady-state 2D problem on a square domain. A narrow Gaussian source at the origin mimics a point charge, and the screening term causes the potential to decay exponentially away from the source rather than following the logarithmic 2D Coulomb potential. The combination of a `Diffusion` kernel (nabla^2 Phi) and a `CoefReaction` kernel (Phi / lambda_D^2) implements the screened Poisson operator, with a `BodyForce` kernel carrying the source. `PointValue` postprocessors at one, two, and five Debye lengths from the origin verify the exponential decay profile.

---

## The Physics

**Governing equation**

The linearized Poisson-Boltzmann (screened Poisson) equation:

```
nabla^2 Phi - Phi / lambda_D^2 = -f(x, y)
```

where lambda_D is the Debye length and f(x, y) is the source (localized charge distribution). Rearranging:

```
nabla^2 Phi - (1 / lambda_D^2) * Phi + f(x, y) = 0
```

**Analytic solution (3D)**

For a true point charge in 3D, the solution is the Yukawa potential:

```
Phi(r) = (Q / 4 pi epsilon r) * exp(-r / lambda_D)
```

**Solution in 2D**

In 2D (relevant to this simulation), the screened potential takes the form:

```
Phi(r) ~ K_0(r / lambda_D)
```

where K_0 is the modified Bessel function of the second kind of order zero. For r >> lambda_D this approaches:

```
Phi(r) ~ sqrt(pi / (2 r/lambda_D)) * exp(-r / lambda_D)
```

**Source term**

A narrow Gaussian centered at the origin replaces the singular point charge:

```
f(x, y) = 500 * exp(-(x^2 + y^2) / (2 * 0.01))
```

Width sigma_g = 0.1 m << lambda_D = 0.2 m, so the source appears point-like to the far field.

**Boundary conditions**

All four outer boundaries: Phi = 0. The domain extends 10 Debye lengths from the source (the boundary at |x| = 2 or |y| = 2 is at distance 10 lambda_D = 2 from the origin), where the screened potential has decayed to exp(-10) ≈ 5e-5. Setting Phi = 0 there introduces negligible error.

**Domain and discretization**

- 2D square: x in [-2, 2], y in [-2, 2]
- 40 x 40 bilinear quadrilateral elements (dx = dy = 0.1 m = 0.5 lambda_D)

**Parameters**

| Symbol | Value | Meaning |
|--------|-------|---------|
| lambda_D | 0.2 m | Debye length |
| 1 / lambda_D^2 | 25.0 m^-2 | Screening coefficient |
| Domain half-width | 2.0 m | = 10 lambda_D |
| Source amplitude | 500 | Calibrated so Phi_max ≈ 1 |
| Source width | 0.1 m | = 0.5 lambda_D (sub-Debye scale) |

---

## Input File Walkthrough

**Global parameter**

```
inv_lambda_D_sq = 25.0   # 1 / (0.2^2) = 25
```

This is the coefficient of the linear screening term. The Debye length itself is not directly used in the kernels — only its square inverse appears in the PDE.

**[Mesh]**

A centered 40 x 40 grid spanning [-2, 2] x [-2, 2]. The origin (source location) is at the center of the domain, at equal distance (10 lambda_D) from all four boundaries.

**[Variables]**

`phi` is the electric potential, initialized to zero. This is a steady-state problem but the initial condition is needed for the NEWTON solver starting guess.

**[Kernels]**

Three kernels implement the screened Poisson equation:

```
[laplacian]  type = Diffusion      — nabla^2 Phi
[screening]  type = CoefReaction   — (1/lambda_D^2) * Phi  (coefficient = 25)
[source]     type = BodyForce      — f(x, y)  (Gaussian source)
```

The weak form of the governing equation is:

```
integral( grad(phi) . grad(psi) + 25 * phi * psi - f * psi ) dV = 0
```

Note that `CoefReaction` adds a positive coefficient times phi to the residual. Because the screening term enters the PDE as -Phi/lambda_D^2 on the left-hand side (which moves to +Phi/lambda_D^2 when moved to the residual), the positive coefficient in `CoefReaction` is correct.

**[Functions]**

```
[source_function]  ParsedFunction  '500.0 * exp(-((x*x + y*y)) / (2 * 0.01))'
```

The amplitude 500 is chosen empirically so that the peak potential at the origin is approximately 1 V. The width parameter 0.01 is sigma_g^2 with sigma_g = 0.1 m.

**[BCs]**

```
[outer_zero]  DirichletBC  boundary = 'left right top bottom'  value = 0
```

A single BC block applies Phi = 0 to all four outer boundaries simultaneously using a space-separated boundary name list.

**[Postprocessors]**

| Name | Point | Expected value |
|------|-------|---------------|
| max_phi | (origin) | ~1.0 (set by source amplitude) |
| phi_at_1lambda | (0.2, 0) | K_0(1) * C ≈ 0.42 C |
| phi_at_2lambda | (0.4, 0) | K_0(2) * C ≈ 0.11 C |
| phi_at_5lambda | (1.0, 0) | K_0(5) * C ≈ 0.004 C |
| phi_integral | domain | Total electrostatic energy proxy |

The decay from phi_at_1lambda to phi_at_2lambda to phi_at_5lambda should show clearly faster-than-algebraic (exponential) falloff.

**[Executioner]**

Steady-state NEWTON with BoomerAMG multigrid. Tight tolerances (nl_rel_tol = 1e-10) are appropriate because the smooth Gaussian source and the exponential solution make the linear system well-conditioned.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case98-debye-shielding \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case98_debye_shielding.i 2>&1 | tail -30'
```

The steady-state solve converges in a handful of NEWTON iterations and completes in under one second. Exodus output is written to `case98_debye_shielding_out.e`; scalar postprocessors go to `case98_debye_shielding_out.csv`.

---

## Expected Results

The potential field should display a sharply peaked maximum at the origin surrounded by an exponentially decaying halo. The decay is noticeably faster than the 1/r (logarithmic in 2D) Coulomb profile: at a distance of one lambda_D the potential is already strongly suppressed compared to a bare Coulomb source, and by five lambda_D it is negligible.

Indicative values from the CSV (exact numbers depend on mesh resolution and source normalization):

| Distance from origin | r / lambda_D | Phi (approx.) |
|---------------------|--------------|---------------|
| 0.0 m (max_phi) | 0 | ~1.0 |
| 0.2 m | 1 | ~0.20-0.35 |
| 0.4 m | 2 | ~0.04-0.08 |
| 1.0 m | 5 | ~0.001-0.005 |

The ratio phi_at_2lambda / phi_at_1lambda is approximately K_0(2) / K_0(1) ≈ 0.114 / 0.421 ≈ 0.27, which can be verified directly from the CSV output. The domain integral phi_integral provides a check on the total "charge" screened by the plasma.

---

## Key Takeaways

- The linearized Poisson-Boltzmann (screened Poisson) equation nabla^2 Phi - Phi/lambda_D^2 = -f governs Debye shielding and appears identically in plasma physics, semiconductor physics (Thomas-Fermi), colloid science (DLVO), and nuclear physics (Yukawa).
- The Debye length lambda_D = sqrt(epsilon k_B T / (2 n_0 q^2)) is the shielding distance; beyond this scale the plasma appears electrically neutral.
- `Diffusion` + `CoefReaction` is the standard MOOSE pattern for the Helmholtz-type operator nabla^2 u - k^2 u.
- `BodyForce` with a `ParsedFunction` source provides a smooth, spatially localized excitation without the mesh singularity problems of a true point source.
- A single `DirichletBC` block can apply to multiple named boundaries simultaneously by listing them in a space-separated string: `boundary = 'left right top bottom'`.
- Placing the source at the domain center and extending the domain to 10 Debye lengths in all directions ensures that the zero far-field boundary condition introduces error of order exp(-10) ≈ 5e-5 — completely negligible relative to numerical discretization error.
