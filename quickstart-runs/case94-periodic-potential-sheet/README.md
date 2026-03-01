# Case 94: Spatially Periodic Potential Sheet — Laplace Equation with Sinusoidal Boundary Data

## Overview

A flat sheet of surface charge with a sinusoidal distribution Phi(x, y=0) = V0 sin(ax) is placed at y = 0. Above and below the sheet the medium is charge-free, so the electric potential satisfies Laplace's equation. Separation of variables yields the exact analytic solution immediately: the potential must oscillate in x with the same wavenumber a as the boundary, and must decay exponentially in y to satisfy the far-field condition. The result is Phi(x, y) = V0 sin(ax) exp(-ay), which is one of the most important closed-form results in electrostatics.

The central lesson from MIT 6.641 Lecture 10 is that spatial periodicity fixes the rate of decay. A charge sheet with short wavelength (large a) has fields that die out quickly in the transverse direction, while a charge sheet with long wavelength (small a) has fields that reach far into the surrounding space. This is why the electric field of a single infinite plane decays in the limit of vanishing wavenumber — or, rephrased, why a uniform surface charge produces a uniform E-field that does not decay.

This case introduces two new MOOSE ideas. First, the analytic solution is registered as a `ParsedFunction` and then projected onto an `AuxVariable` via `FunctionAux`, enabling pointwise comparison in the Exodus output. Second, `ElementL2Error` is used as a postprocessor to compute the global L2 norm of the numerical error, which provides a single scalar measure of how well the mesh resolves the exponential decay near y = 0.

---

## The Physics

**Governing equation**

Laplace's equation in the charge-free region:

```
nabla^2 Phi = d^2 Phi/dx^2 + d^2 Phi/dy^2 = 0
```

**Analytic solution (separation of variables)**

```
Phi(x, y) = V0 sin(a x) exp(-a y)
```

The wavenumber a = 2 pi corresponds to wavelength lambda = 1 m. The decay length (one e-folding) is 1/a = 1/(2 pi) ≈ 0.159 m.

**Boundary conditions**

| Boundary | Condition | Physical meaning |
|----------|-----------|-----------------|
| Bottom (y = 0) | Phi = V0 sin(ax) | Sinusoidal potential sheet |
| Top (y = L_y = 2) | Phi = 0 | Far-field: Phi/V0 = exp(-4 pi) ≈ 3.5e-6 |
| Left/right (x = 0, 1) | Natural (zero flux) | Consistent with periodicity: sin(a x) = 0 at endpoints |

**Domain and discretization**

- 2D rectangle: x in [0, 1] (one full wavelength), y in [0, 2] (12.6 decay lengths)
- 40 x 40 bilinear quadrilateral elements

**Parameters**

| Symbol | Value | Meaning |
|--------|-------|---------|
| V0 | 1.0 V | Peak sheet potential |
| a | 2 pi ≈ 6.283 rad/m | Spatial wavenumber |
| lambda | 1.0 m | Wavelength (= 2 pi / a) |
| 1/a | 0.159 m | Exponential decay length |

---

## Input File Walkthrough

**Global parameters**

```
V0   = 1.0
a    = 6.2831853   # 2pi
L_y  = 2.0
```

These are HIT-level constants substituted with `${...}` syntax throughout the file.

**[Mesh]**

A single `GeneratedMeshGenerator` produces a 40 x 40 grid. The x-extent is exactly one wavelength (xmax = 1.0 = 2 pi / a) and the y-extent is two meters — far enough that the analytic solution is negligible at the top.

**[Variables]**

`phi` is a first-order Lagrange scalar — the electric potential in volts.

**[Kernels]**

The `Diffusion` kernel provides the weak form of Laplace's equation:

```
integral( grad(phi) . grad(psi) dV ) = 0
```

No source term is needed because the domain is charge-free.

**[BCs]**

`FunctionDirichletBC` applies `bottom_potential = V0 sin(a x)` at y = 0. A standard `DirichletBC` enforces Phi = 0 at the top. The left and right boundaries receive the natural (zero-flux) Neumann condition implicitly — no explicit entry is required in MOOSE.

**[Functions]**

Two `ParsedFunction` objects are defined:

- `bottom_potential`: evaluates `V0 * sin(a * x)` — used by the bottom BC.
- `analytic_solution`: evaluates `V0 * sin(a * x) * exp(-a * y)` — used for validation.

**[AuxVariables] / [AuxKernels]**

`phi_exact` stores the analytic solution projected via `FunctionAux`. This is written to the Exodus file alongside `phi`, making a direct visual comparison straightforward in ParaView.

**[Postprocessors]**

| Name | Type | What it measures |
|------|------|-----------------|
| l2_error | ElementL2Error | Global L2 norm of Phi - Phi_exact |
| max_phi | ElementExtremeValue | Peak potential in domain (should be ~1 V) |
| phi_at_quarter | PointValue at (0.25, 0) | Should equal V0 sin(pi/2) = 1.0 |
| phi_at_one_decay | PointValue at (0.25, 0.159) | Should equal exp(-1) ≈ 0.368 |
| phi_at_half | PointValue at (0.25, 0.5) | Should equal exp(-pi) ≈ 0.043 |

**[Executioner]**

Steady-state NEWTON solve with BoomerAMG algebraic multigrid preconditioning. Tight tolerances (nl_rel_tol = 1e-10) ensure the L2 error is dominated by spatial discretization error rather than solver residual.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case94-periodic-potential-sheet \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case94_periodic_potential_sheet.i 2>&1 | tail -30'
```

The simulation runs in well under one second. Exodus output is written to `case94_periodic_potential_sheet_out.e` and scalar postprocessors to `case94_periodic_potential_sheet_out.csv`.

---

## Expected Results

The potential field should display the classic "evanescent wave" pattern: full-amplitude sinusoidal oscillations at y = 0, collapsing exponentially toward zero within the first ~0.5 m. The decay is so rapid that by y = 1 the field is essentially zero on the color scale.

Key quantitative checks from the CSV output:

| Postprocessor | Expected value |
|---------------|----------------|
| phi_at_quarter (y = 0) | 1.000 |
| phi_at_one_decay (y = 0.159) | 0.368 |
| phi_at_half (y = 0.5) | 0.043 |
| l2_error | < 1e-3 (40x40 mesh) |

In ParaView, loading both `phi` and `phi_exact` and plotting their difference confirms that the error is largest in the region near y = 0 where the gradient is steepest — the same place where adaptive mesh refinement would be most beneficial.

---

## Key Takeaways

- Separation of variables applied to Laplace's equation yields the exponential decay law Phi ~ exp(-a |y|) for any spatially periodic boundary condition with wavenumber a.
- The decay length 1/a equals the wavelength divided by 2 pi — short-wavelength charge patterns produce fields confined close to the surface.
- `FunctionDirichletBC` with a `ParsedFunction` is the standard MOOSE pattern for non-uniform, mathematically defined boundary data.
- Registering the analytic solution as a `ParsedFunction` and projecting it through `FunctionAux` enables clean, automated validation alongside every simulation run.
- `ElementL2Error` condenses the pointwise error over the entire domain into a single number, making mesh-convergence studies straightforward.
- The natural (zero-flux) Neumann boundary condition in MOOSE requires no input block — omitting a boundary from `[BCs]` automatically applies zero normal flux.
