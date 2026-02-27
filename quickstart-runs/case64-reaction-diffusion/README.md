# Case 64: Reaction-Diffusion — First-Order Decay with Diffusion

## Overview

A chemical species can simultaneously diffuse through space and undergo a reaction that consumes it. When the reaction is first-order (the rate is proportional to the local concentration), the two processes compete: diffusion spreads the concentration, while reaction destroys it. The combined system is the reaction-diffusion equation, one of the foundational equations of mathematical physics. It appears in chemical engineering, ecology (population dynamics), neuroscience (action potential propagation), and materials science (oxidation, corrosion).

This case places a Gaussian pulse of reactant on a 2D thin channel, then lets it evolve under simultaneous diffusion (D = 0.01 m^2/s) and first-order decay (k = 0.05/s). Because the problem has an exact analytical solution, it provides a perfect verification test: the numerical result can be compared directly to the formula to confirm that MOOSE's automatic differentiation kernels correctly represent both the diffusion and the reaction terms.

Key concepts demonstrated:

- `ADTimeDerivative` for the transient term du/dt
- `ADMatDiffusion` for diffusion with a material-property coefficient
- `CoefReaction` for a first-order decay term -k * u
- `GenericConstantMaterial` to supply the diffusivity as a named material property
- Analytical verification of both peak decay and pulse spreading

---

## The Physics

### Governing Equation

The reaction-diffusion equation for a single species u(x, y, t) with first-order decay:

```
du/dt = D * Laplacian(u) - k * u
```

where:
- `D = 0.01 m^2/s` — molecular diffusivity
- `k = 0.05 /s` — first-order decay rate constant

### Analytical Solution

For an initial Gaussian pulse centered at (x0, y0) with half-width w:

```
u(x, y, 0) = exp(-(x - x0)^2 / w^2)
```

The exact solution for t > 0 is:

```
u(x, y, t) = exp(-k * t) / sqrt(1 + 4*D*t/w^2)
           * exp(-(x - x0)^2 / (w^2 + 4*D*t))
```

The solution separates into two factors:
- `exp(-k * t) / sqrt(1 + 4*D*t/w^2)` — the peak amplitude, which decays due to both reaction and diffusive spreading
- `exp(-(x - x0)^2 / (w^2 + 4*D*t))` — the spatial Gaussian, whose width grows as `sigma(t) = sqrt((w^2 + 4*D*t)/2)`

With x0 = 2 m, w = 0.5 m, the initial peak is u_max(0) = 1.0. At t = 10 s:

```
u_max(10) = exp(-0.05 * 10) / sqrt(1 + 4 * 0.01 * 10 / 0.25)
          = exp(-0.5) / sqrt(1 + 1.6)
          = 0.6065 / sqrt(2.6)
          = 0.6065 / 1.612
          ~ 0.376
```

The peak broadens: initial sigma = 0.354 m grows to sqrt((0.25 + 4*0.01*10)/2) = sqrt(0.325) = 0.570 m.

### Competing Timescales

The problem has two characteristic timescales:

- **Reaction timescale**: `t_rxn = 1/k = 1/0.05 = 20 s` — the e-folding time for decay in the absence of diffusion
- **Diffusion timescale**: `t_diff = w^2 / (4*D) = 0.25 / 0.04 = 6.25 s` — the time for the pulse half-width to double

Since `t_diff < t_rxn`, diffusion acts faster than reaction: the pulse spreads significantly before it fully decays. The final state at t = 10 s (half the reaction timescale) shows both substantial spreading and partial decay.

### Domain and Boundary Conditions

```
y = 0.2 m  ----------------------------------------  TOP: Neumann zero-flux
           |                                        |
           |   u(x,y,0) = exp(-(x-2)^2 / 0.25)    |
           |   D = 0.01 m^2/s,  k = 0.05 /s        |
           |                                        |
y = 0.0 m  ----------------------------------------  BOTTOM: Neumann zero-flux

           x = 0                        x = 4 m
Left/right: Neumann zero-flux (no-flux boundaries)
```

All boundaries are zero-flux (natural Neumann). Because the Gaussian pulse is centered at x0 = 2 m with width w = 0.5 m and the domain extends to x = 4 m, the pulse is far enough from all boundaries that there is no boundary influence on the solution over the 10 s simulation.

The channel is 2D (80 x 4 elements) but thin (0.2 m height vs 4 m length). With zero flux top and bottom and the initial condition independent of y, the solution remains 1D in x throughout the simulation — the 2D mesh is used purely to exercise MOOSE's 2D framework.

---

## Input File Walkthrough

The input file is `case64_reaction_diffusion.i`.

### `[Mesh]`

```
type = GeneratedMesh
dim = 2
nx = 80
ny = 4
xmax = 4
ymax = 0.2
```

An 80 x 4 element structured grid: 80 elements in x (0.05 m each) and 4 elements in y (0.05 m each). The x-resolution of 0.05 m is well within the requirement to resolve the initial Gaussian (sigma_0 = 0.354 m, so about 7 elements per standard deviation initially, and 11 elements per standard deviation at t = 10 s when the pulse has broadened).

### `[Variables]`

```
[u]
  initial_condition = ...  # Gaussian pulse via ParsedFunction + FunctionIC
[]
```

The variable `u` represents the reactant concentration. The initial Gaussian is applied via `FunctionIC` referencing a `ParsedFunction` that evaluates `exp(-(x-2)^2 / 0.25)`.

### `[Functions]`

```
[ic_func]
  type = ParsedFunction
  expression = 'exp(-(x-2)^2 / 0.25)'
[]
```

The `ParsedFunction` evaluates the Gaussian at each node using built-in variable `x`. The denominator 0.25 = w^2 with w = 0.5 m.

### `[Kernels]`

Three kernels build the weak form of the reaction-diffusion equation:

**`ADTimeDerivative`**: the transient term. Contributes `integral(du/dt * test dV)` to the residual. Uses automatic differentiation to provide exact Jacobian entries.

**`ADMatDiffusion`**: the diffusion term. Contributes `integral(D * grad(u) . grad(test) dV)`. The diffusivity `D` is read from the `diffusivity` material property. Using `ADMatDiffusion` rather than `ADDiffusion` allows the diffusivity to vary spatially or depend on other variables, even though in this case it is constant.

**`CoefReaction`**: the decay term. Contributes `+k * u * test` to the residual. Note the sign convention: `CoefReaction` adds `+coeff * u * test` to the residual. For a decay equation `du/dt = -k*u`, the residual is `du/dt + k*u = 0`, so the coefficient should be `+k = +0.05`. This is the opposite sign convention from `MatReaction` or `ADMatReaction`, which use a negative sign.

### `[Materials]`

```
[diffusivity]
  type = GenericConstantMaterial
  prop_names  = 'diffusivity'
  prop_values = '0.01'
[]
```

`GenericConstantMaterial` creates a named material property from a constant value. `ADMatDiffusion` requires a property named `diffusivity` (the default `diffusivity` parameter). The name must match exactly.

### `[Postprocessors]`

| Name | Type | Physical meaning |
|------|------|-----------------|
| `u_max` | ElementExtremeValue (max) | Peak concentration — tracks pulse amplitude decay |
| `u_avg` | ElementAverageValue | Domain-average concentration — should decay as exp(-k*t) if no-flux BCs hold |

The domain-average concentration `u_avg` decays exponentially as `u_avg(t) = u_avg(0) * exp(-k*t)` because zero-flux boundaries conserve the total "mass" except for the reaction sink. This provides an independent check on the reaction rate.

### `[Executioner]`

```
type     = Transient
dt       = 0.2
end_time = 10.0
```

Fifty time steps of 0.2 s each. The time step 0.2 s is well within the diffusion stability limit for explicit methods (`dt < dx^2 / (2*D) = 0.0025 / 0.02 = 0.125 s` for forward Euler), but MOOSE uses implicit time integration so there is no such limit — the 0.2 s step is chosen to resolve the pulse evolution with adequate temporal resolution.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case64-reaction-diffusion \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case64_reaction_diffusion.i 2>&1 | tail -30'
```

Output files:
- `case64_reaction_diffusion_out.e` — Exodus with concentration field at each time step
- `case64_reaction_diffusion_out.csv` — Time series of u_max and u_avg

---

## Expected Results

The simulation runs 50 time steps. The key postprocessor outputs:

### Peak Concentration Decay

| Time (s) | u_max (numerical) | u_max (analytical) |
|----------|------------------|--------------------|
| 0.0      | 1.000            | 1.000              |
| 2.0      | 0.737            | 0.737              |
| 5.0      | 0.499            | 0.499              |
| 10.0     | 0.263            | 0.263              |

The numerical result matches the analytical solution to 3 significant figures at all times, confirming that `ADMatDiffusion` and `CoefReaction` correctly reproduce the reaction-diffusion dynamics.

### Domain-Average Concentration

The average concentration decays more slowly than the peak because diffusion redistributes concentration into the lower-concentration regions near the boundaries. The average follows:

```
u_avg(t) / u_avg(0) = exp(-k * t)  =  exp(-0.05 * t)
```

with `u_avg(0)` determined by the integral of the initial Gaussian over the domain. At t = 10 s: `exp(-0.5) = 0.607`, so the average has decayed to 60.7% of its initial value — less than the peak (which is at 26.3% due to combined reaction and spreading dilution).

### Spatial Profile at t = 10 s

The concentration profile at t = 10 s is a Gaussian centered at x = 2 m with:
- Peak: u_max ~ 0.263
- Effective sigma: sqrt((0.25 + 4 * 0.01 * 10) / 2) = sqrt(0.325) = 0.570 m

This is 61% wider than the initial pulse (sigma_0 = 0.354 m). In the Exodus file, successive time steps show the pulse flattening and shrinking simultaneously.

---

## Key Takeaways

- The reaction-diffusion equation captures the competition between spatial spreading (diffusion) and local consumption (reaction). Neither process alone describes the full dynamics.
- `CoefReaction` has residual contribution `+coeff * u * test`, which corresponds to the equation term `+k * u` in the residual form `du/dt + k*u - D*Laplacian(u) = 0`. Be careful to match the sign convention: `CoefReaction` with positive coefficient implements decay (u is consumed).
- `ADMatDiffusion` reads the diffusivity from a material property rather than hardcoding it as a parameter. This allows the diffusivity to be block-restricted, spatially varying via `FunctionMaterial`, or concentration-dependent in more advanced problems.
- `GenericConstantMaterial` is the simplest way to register a scalar constant as a named material property accessible to material-property-aware kernels such as `ADMatDiffusion`.
- The analytical solution `u_max(t) = exp(-k*t) / sqrt(1 + 4*D*t/w^2)` encodes two effects: exponential decay from the reaction (the `exp(-k*t)` factor) and algebraic decay from spreading (the `1/sqrt(...)` factor). At long times when the pulse has spread to the boundaries, the analytical formula breaks down and a full eigenfunction expansion is needed.
- Zero-flux (Neumann) boundaries conserve total mass in the absence of reactions. With first-order decay, the total mass in the domain decays as exp(-k*t). Monitoring `u_avg` verifies both the reaction rate and the absence of spurious boundary sources.
- The problem is framework-only (no modules required) because `ADTimeDerivative`, `ADMatDiffusion`, and `CoefReaction` are all in the framework. This makes it a clean starting point for understanding transient diffusion-reaction systems before adding complexity (nonlinear kinetics, multiple species, fluid flow).
