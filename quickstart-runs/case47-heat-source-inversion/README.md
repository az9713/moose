# Case 47: Heat Source Inversion — PDE-Constrained Optimization

## Overview

In many engineering and scientific contexts, you cannot directly measure the quantity you care about most. You can place thermocouples in a furnace but cannot measure the internal heat source driving those temperatures. You can image brain activity with MRI but cannot directly observe the underlying neural currents. This class of problems — recovering an unknown input from observed outputs of a PDE — is called a PDE-constrained inverse problem.

This case demonstrates MOOSE's optimization module by recovering an unknown volumetric heat source q from four sparse temperature measurements. The forward problem is 2D steady heat conduction; the true source is q = 1000 W/m^3. Starting from an initial guess of q = 500, a TAO L-BFGS optimizer iterates by solving the forward PDE, computing a misfit between simulated and measured temperatures, then solving an adjoint PDE to compute the gradient of the misfit with respect to q. The optimizer converges to q = 1000 in a single iteration.

The three-file architecture (main + forward + adjoint) is the standard structure for adjoint-based PDE-constrained optimization in MOOSE. Understanding this structure unlocks a powerful class of problems: parameter identification, experimental design, and topology optimization all follow the same pattern.

## The Physics

**Forward problem**: Steady heat conduction on a 2x2 meter square domain:

```
-k * laplacian(T) = q    on [0,2] x [0,2]
T = 200 K                on bottom boundary (y = 0)
T = 100 K                on top boundary    (y = 2)
dT/dn = 0                on left and right boundaries (insulated)
k = 5 W/(m·K)
```

The boundary conditions model a wall heated from below and cooled from above, with insulated sides. The unknown source q represents additional internal heat generation (e.g., electrical resistive heating or chemical reaction).

**Measurement data**: Four temperature sensors at interior points:

| Sensor | x    | y    | T_measured |
|--------|------|------|------------|
| 1      | 0.20 | 0.20 | 226 K      |
| 2      | 0.80 | 0.60 | 254 K      |
| 3      | 0.20 | 1.40 | 214 K      |
| 4      | 0.80 | 1.80 | 146 K      |

**Objective function**: Sum of squared misfits:

```
J(q) = sum_i (T_sim(x_i, q) - T_meas_i)^2
```

**Adjoint problem**: To compute the gradient dJ/dq without finite differences, an adjoint equation is solved. For this forward PDE, the adjoint is:

```
-k * laplacian(lambda) = sum_i (T_sim(x_i) - T_meas_i) * delta(x - x_i)
lambda = 0    on all boundaries
```

where the point sources at measurement locations are weighted by the misfit. The gradient is then:

```
dJ/dq = integral(lambda) dV
```

This follows from the adjoint method: the gradient computation requires one additional linear solve (the adjoint), not one solve per parameter as finite differences would require.

## Input File Walkthrough

### Main Optimizer: `case47_main.i`

**`[Optimization]`** — Empty action block that activates the optimization module, analogous to `[StochasticTools]` in the UQ cases.

**`[OptimizationReporter]`** — Defines the optimization problem:

```
type = GeneralOptimization
objective_name = objective_value
parameter_names = 'parameter_results'
num_values = '1'
initial_condition = '500'
lower_bounds = '0.1'
upper_bounds = '10000'
```

`GeneralOptimization` is the most flexible type: it accepts gradient information from the adjoint app rather than computing it internally. The single parameter is the scalar q (one value). Bounds prevent the optimizer from exploring physically meaningless negative or astronomically large sources.

**`[Reporters]`** — The main app's `OptimizationData` reporter stores the measurement point locations and values, and accumulates the misfit after each forward solve:

```
measurement_points = '0.2 0.2 0
                      0.8 0.6 0
                      0.2 1.4 0
                      0.8 1.8 0'
measurement_values = '226 254 214 146'
```

**`[Executioner]`** — The TAO L-BFGS optimizer:

```
type = Optimize
tao_solver = taoblmvm
petsc_options_iname = '-tao_gatol -tao_max_it'
petsc_options_value = '1e-4       50'
```

TAO (Toolkit for Advanced Optimization) is PETSc's optimization library. `taoblmvm` is the limited-memory BFGS algorithm, which approximates the Hessian using gradient history. `-tao_gatol` is the gradient norm tolerance for convergence.

**`[MultiApps]`** — Two sub-apps execute at different phases of the optimization loop:

```
[forward]
  type = FullSolveMultiApp
  input_files = case47_forward.i
  execute_on = FORWARD
[]
[adjoint]
  type = FullSolveMultiApp
  input_files = case47_adjoint.i
  execute_on = ADJOINT
[]
```

The `FORWARD` and `ADJOINT` execute_on tokens are special to the Optimize executioner: it alternates between forward and adjoint solves on each iteration.

**`[Transfers]`** — Six transfers move data between the three apps:

- `toForward`: sends measurement coordinates, measurement values, and current q to the forward app
- `toAdjoint`: sends measurement coordinates, misfit values (T_sim - T_meas), and current q to the adjoint app
- `fromForward`: retrieves misfit values and the objective value back to the main app
- `fromAdjoint`: retrieves the gradient dJ/dq from the adjoint app

### Forward App: `case47_forward.i`

A standard heat conduction solve. The key addition is `ParsedOptimizationFunction`, which receives the optimizer's current q value through a reporter and exposes it as a MOOSE function:

```
[volumetric_heat_func]
  type = ParsedOptimizationFunction
  expression = q
  param_symbol_names = 'q'
  param_vector_name = 'params/q'
[]
```

The `BodyForce` kernel uses this function as its source. The `OptimizationData` reporter computes and stores the objective value and per-sensor misfit values after the solve completes.

### Adjoint App: `case47_adjoint.i`

The adjoint problem has the same diffusion operator as the forward problem but different source terms and boundary conditions. The boundary conditions are all homogeneous (zero Dirichlet everywhere), reflecting the fact that the adjoint of the forward operator with fixed Dirichlet BCs has zero adjoint BCs.

The point sources at measurement locations are applied via `ReporterPointSource`, a `DiracKernel` that reads point coordinates and weights from reporters:

```
[pt]
  type = ReporterPointSource
  variable = adjoint_T
  x_coord_name = misfit/measurement_xcoord
  y_coord_name = misfit/measurement_ycoord
  z_coord_name = misfit/measurement_zcoord
  value_name = misfit/misfit_values
[]
```

After the adjoint solve, `ElementOptimizationSourceFunctionInnerProduct` computes the gradient by integrating the adjoint solution against the derivative of the forward source with respect to q:

```
[gradient_vpp]
  type = ElementOptimizationSourceFunctionInnerProduct
  variable = adjoint_T
  function = volumetric_heat_func
[]
```

For a spatially constant q, the derivative of the source with respect to q is 1 everywhere, so the inner product reduces to the integral of the adjoint solution over the domain.

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" -w /work/case47-heat-source-inversion --entrypoint /bin/bash idaholab/moose:latest -c '/opt/moose/bin/combined-opt -i case47_main.i 2>&1 | tail -40'
```

The run is fast: one forward solve and one adjoint solve per optimizer iteration. Output files include `case47_main_out.csv` (objective history) and `case47_main_out_OptimizationReporter_0001.csv` (final parameter value).

## Expected Results

The optimizer converges in a single L-BFGS iteration:

| Iteration | Objective J(q)   | q (recovered) |
|-----------|------------------|---------------|
| 0 (initial) | ~0 (no solve yet) | 500.0       |
| 1           | 1.56e-23          | 1000.0      |

The recovered value q = 999.999... is essentially exact. The objective drops from a finite value at the initial guess to machine-zero after one step.

One-step convergence is not a coincidence. Because the forward PDE is linear in q, the simulated temperatures are linear functions of q: T_sim(x_i, q) = T_BC(x_i) + q * G(x_i), where G is the Green's function evaluated at the sensor locations. The objective J(q) is therefore exactly quadratic in the scalar q. L-BFGS initialized with a good Hessian approximation will minimize a quadratic exactly in one step. The adjoint gradient provides the exact gradient of this quadratic, so the step size is computed precisely.

## Key Takeaways

**The adjoint method computes gradients cheaply.** For a problem with N parameters, finite differences require N+1 forward solves to estimate the gradient. The adjoint method requires exactly 1 adjoint solve regardless of N. For problems with thousands of spatial degrees of freedom as parameters (image reconstruction, topology optimization), the adjoint is the only tractable approach.

**Three-file architecture separates concerns cleanly.** The main app knows about the optimization algorithm; the forward app knows about the physics; the adjoint app knows about the sensitivity calculation. None of the three needs to know the details of the others' implementations. This modularity makes it straightforward to swap in more complex physics.

**Linear PDEs give quadratic objectives, which optimize trivially.** The instant convergence here is a feature of the linear forward operator, not a general property. Nonlinear PDEs (combustion, fluid flow, phase field) give non-quadratic objectives that require many optimizer iterations and careful regularization.

**Ill-posedness is the hidden challenge.** This case has a single scalar unknown q and four measurements — an over-determined problem. In practice, recovering a spatially varying q(x,y) from sparse measurements is severely ill-posed: many different source distributions can fit the same four sensor readings. Tikhonov regularization (adding a penalty on the gradient of q) is the standard remedy.

**TAO connects MOOSE to PETSc's full optimization toolkit.** Beyond L-BFGS (`taoblmvm`), TAO supports gradient-free methods (`taonm`), bound-constrained methods, and nonlinear least-squares formulations. Switching algorithms requires changing only the `tao_solver` parameter.
