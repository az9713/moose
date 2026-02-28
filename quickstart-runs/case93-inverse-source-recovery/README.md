# Case 93: Inverse Source Recovery — PDE-Constrained Estimation

## Overview

This case demonstrates PDE-constrained optimization to recover an unknown source distribution from sparse field measurements. Given observations of a scalar field u at 9 interior sensor locations, the optimizer recovers the amplitude of the source that generated those observations.

**Forward problem** (screened Poisson equation):

```
−∇²u + u = q(x,y)    on [0,1]²
u = 0                  on ∂Ω
```

**Inverse problem**: Find the source amplitude p₁ such that `q(x,y) = p₁ sin(πx) sin(πy)` best fits 9 measurements of u.

**True answer**: p₁ = 10. Starting guess: p₁ = 1.

## Physics: Screened Poisson Equation

The operator −∇² + αI (with α = 1) is the screened (or modified) Helmholtz / Yukawa operator. It governs:

- Steady-state diffusion with first-order decay
- Pressure correction in incompressible flow solvers
- Potential fields with exponential screening (Debye–Hückel)
- Source inversion problems in geophysics (ERT, gravity)

Unlike the pure Laplacian, the screening term αI makes the operator coercive (positive-definite) and its inverse (the Green's function) decays exponentially with characteristic length 1/√α rather than algebraically. This makes the forward problem well-conditioned and the inverse problem better-posed.

## Adjoint Method

The gradient dJ/dp₁ is computed analytically using the adjoint method — a single additional PDE solve rather than finite differences.

**Objective functional:**

```
J(p₁) = (1/2) Σᵢ [u(xᵢ, yᵢ; p₁) − u_obs,ᵢ]²
```

**Adjoint equation** (same operator, point sources at sensors):

```
−∇²λ + λ = Σᵢ eᵢ δ(x − xᵢ)    on [0,1]²
λ = 0                             on ∂Ω
```

where `eᵢ = u(xᵢ) − u_obs,ᵢ` is the measurement misfit.

**Gradient formula:**

```
dJ/dp₁ = ∫ λ(x,y) × sin(πx) sin(πy) dΩ
```

## Analytical Solution

For `q = p₁ sin(πx) sin(πy)` with zero Dirichlet BCs, the exact solution is:

```
u(x,y) = p₁ / (2π² + 1) × sin(πx) sin(πy)
```

With `2π² + 1 ≈ 20.739` and `p₁ = 10`:

```
u_exact(x,y) = 0.48218 × sin(πx) × sin(πy)
```

Measurement values at the 9 sensor locations:

| x    | y    | u_obs  |
|------|------|--------|
| 0.25 | 0.25 | 0.24109 |
| 0.50 | 0.25 | 0.34099 |
| 0.75 | 0.25 | 0.24109 |
| 0.25 | 0.50 | 0.34099 |
| 0.50 | 0.50 | 0.48218 |
| 0.75 | 0.50 | 0.34099 |
| 0.25 | 0.75 | 0.24109 |
| 0.50 | 0.75 | 0.34099 |
| 0.75 | 0.75 | 0.24109 |

## Three-File Architecture

```
case93_inverse_source_recovery.i   (MAIN — TAO optimizer)
         |
         |-- FORWARD execute_on --> case93_forward.i
         |         Solves: −∇²u + u = p₁ sin(πx)sin(πy)
         |         Returns: misfit eᵢ, objective J
         |
         |-- ADJOINT execute_on --> case93_adjoint.i
                   Solves: −∇²λ + λ = Σᵢ eᵢ δ(x−xᵢ)
                   Returns: gradient dJ/dp₁ = ∫ λ sin(πx)sin(πy) dΩ
```

## Files

| File | Purpose |
|------|---------|
| `case93_inverse_source_recovery.i` | Main optimizer (TAO L-BFGS, manages transfers) |
| `case93_forward.i` | Forward PDE solve (screened Poisson) |
| `case93_adjoint.i` | Adjoint sensitivity solve (misfit sources) |

## Mesh and Solver

| Parameter | Value |
|-----------|-------|
| Domain | [0,1] × [0,1] |
| Mesh | 30 × 30 elements |
| DOFs | ~841 free nodes per app |
| Solver | Newton + direct LU |
| Optimizer | TAO taoblmvm (L-BFGS) |
| Starting guess | p₁ = 1.0 |
| True value | p₁ = 10.0 |

## Expected Convergence

For this linear PDE with a single parameter, the objective J(p₁) is a quadratic function:

```
J(p₁) = (1/2) × K × (p₁ − 10)²
```

where K = Σᵢ [C sin(πxᵢ) sin(πyᵢ)]² with C = 1/(2π²+1). L-BFGS finds the exact minimum of a quadratic in one step. Expected behavior:

- **Iteration 0**: p₁ = 1.0, J is large, gradient is large and negative
- **Iteration 1**: p₁ ≈ 10.0, J ≈ 0, gradient ≈ 0 (converged)

## Running

```bash
# Run the main optimizer (drives forward and adjoint sub-apps automatically)
combined-opt -i case93_inverse_source_recovery.i

# With verbose TAO output to see iteration history:
# (verbose = true is already set in the main file)
combined-opt -i case93_inverse_source_recovery.i 2>&1 | grep -E "TAO|iter|Obj"
```

## Output Files

**Main app:**
- `case93_inverse_source_recovery_out.csv` — optimization history (p₁, J, |∇J| per iteration)

**Forward sub-app (created in sub-app working directory):**
- `case93_forward_out.e` — forward solution field u(x,y) at final p₁
- `case93_forward_out.csv` — misfit and objective at each step

**Adjoint sub-app:**
- `case93_adjoint_out.e` — adjoint field λ(x,y) at final iteration
- `case93_adjoint_out.csv` — gradient values per iteration

## Visualization

Open `case93_forward_out.e` in ParaView:
1. The converged field u(x,y) should show the sin(πx)sin(πy) mode shape.
2. Peak value at (0.5, 0.5) should be ≈ 0.48218.

Open `case93_adjoint_out.e` in ParaView:
1. The adjoint λ(x,y) shows the sensitivity pattern — peaked near the 9 sensors.
2. Each sensor appears as a rounded "bump" (screened Green's function, width ~ 1/√α = 1).

## Kernel and Reporter Summary

### Forward app kernels

| Kernel | Contribution (strong form) |
|--------|---------------------------|
| `Diffusion` | −∇²u |
| `CoefReaction(1.0)` | +u |
| `BodyForce(source_fn)` | RHS = p₁ sin(πx) sin(πy) |

### Adjoint app kernels

| Kernel/DiracKernel | Contribution (strong form) |
|--------------------|---------------------------|
| `Diffusion` | −∇²λ |
| `CoefReaction(1.0)` | +λ |
| `ReporterPointSource` | RHS = Σᵢ eᵢ δ(x−xᵢ) |

### Key reporters

| Reporter | Role |
|----------|------|
| `OptimizationReporter` (main) | Manages p₁, measurement data, J, ∇J |
| `src_values` (forward) | Receives p₁ from optimizer |
| `measure_data / OptimizationData` (forward) | Evaluates u at sensors, computes misfit |
| `misfit_data` (adjoint) | Receives misfit eᵢ and sensor coords |
| `gradient / ElementOptimizationSourceFunctionInnerProduct` (adjoint) | Computes dJ/dp₁ |

## Sign Convention Summary

The `CoefReaction` kernel (not `ADMatReaction` or `MatReaction`) is used here because:

- `CoefReaction` adds `+coefficient × u × test` → strong form `+αu` (positive screening)
- `ADMatReaction` adds `−rate × u × test` → strong form `−rate × u` (negative — wrong sign)

For the screened Poisson equation `−∇²u + αu = f`, `CoefReaction` with `coefficient = α` is correct.

## Connection to Case 47

This case follows the exact same three-file architecture as Case 47 (Heat Source Inversion), adapted from the heat equation to the screened Poisson equation. Key differences:

| Feature | Case 47 (Heat) | Case 93 (Screened Poisson) |
|---------|----------------|---------------------------|
| Forward PDE | −k∇²T = q | −∇²u + u = q |
| Operator | Laplacian only | Screened Helmholtz |
| Source param | Constant q | Modal: p₁ sin(πx)sin(πy) |
| Mesh | 10×10 | 30×30 |
| Sensors | 4 points | 9 points (3×3 grid) |
| BCs | Mixed (Dirichlet + Neumann) | Pure Dirichlet |

## References

- Vogel, C.R., *Computational Methods for Inverse Problems*, SIAM (2002), Ch. 1–3
- Gunzburger, M.D., *Perspectives in Flow Control and Optimization*, SIAM (2003), Ch. 1–2
- Biegler et al. (eds.), *Large-Scale PDE-Constrained Optimization*, Springer (2003), Ch. 1
- MOOSE Optimization Module documentation: https://mooseframework.inl.gov/modules/optimization/
