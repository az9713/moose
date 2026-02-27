# Case 45: Monte Carlo UQ — Uncertainty in Thermal Conductivity

## Overview

Physical simulations always carry uncertainty. Material properties measured in the lab have tolerances, manufacturing processes introduce variability, and field conditions are never perfectly known. Uncertainty Quantification (UQ) is the discipline of propagating those input uncertainties through a simulation to characterize the resulting uncertainty in outputs of interest.

This case introduces MOOSE's `stochastic_tools` module by running a Monte Carlo (MC) sampling study on a 2D steady heat conduction problem. The thermal conductivity k is modeled as a random variable drawn from a uniform distribution, and 30 independent simulations are run — one per sample — to produce a distribution of the domain-average and peak temperatures. No modification to the physics solver is required; the stochastic framework wraps the deterministic sub-application automatically.

The result reveals a key physical insight: because the steady-state temperature is proportional to q/k (heat source divided by conductivity), the output distribution is the reciprocal of a uniform distribution — right-skewed, with a long tail at high temperatures. This non-Gaussian behavior is invisible if you only run the simulation once with the nominal k value.

## The Physics

The sub-application solves steady heat conduction with a volumetric source on the unit square:

```
-k * laplacian(T) = q    on [0,1] x [0,1]
T = 0                    on all four boundaries (homogeneous Dirichlet)
```

Parameters:
- q = 100 W/m^3 (fixed volumetric heat source)
- k ~ Uniform(0.5, 3.0) W/(m·K) (uncertain thermal conductivity)

For the 1D analog, the analytical solution is T(x) = q/(2k) * x*(1-x), giving a peak temperature of q/(8k) at the center. In 2D the structure is similar: T scales as q/k, so halving k doubles all temperatures. This inverse dependence is what creates the right-skewed output distribution.

Outputs of interest collected from each sample:
- `avg_T`: domain-average nodal temperature (AverageNodalVariableValue)
- `max_T`: peak nodal temperature (NodalExtremeValue)

## Input File Walkthrough

### Main Orchestrator: `case45_monte_carlo_uq.i`

**`[StochasticTools]`** — An empty block that activates the stochastic tools module. Without this, none of the UQ objects below are registered. It signals to MOOSE that this application is a pure orchestrator with no PDE of its own.

**`[Distributions]`** — Declares the probability distribution for the uncertain parameter.

```
[k_dist]
  type = Uniform
  lower_bound = 0.5
  upper_bound = 3.0
[]
```

`Uniform` is the simplest non-trivial distribution. It assigns equal probability to all values in [0.5, 3.0]. Other available types include `Normal`, `Lognormal`, `Beta`, and `Weibull`.

**`[Samplers]`** — Defines how to draw samples from the distribution.

```
[mc]
  type = MonteCarlo
  num_rows = 30
  distributions = 'k_dist'
  seed = 2024
[]
```

`MonteCarlo` draws 30 independent random samples. Setting `seed` makes the run reproducible. Each "row" of the sampler is one vector of parameter values — here just a single scalar k per row.

**`[MultiApps]`** — Launches the sub-application for each sample.

```
[sub]
  type = SamplerFullSolveMultiApp
  sampler = mc
  input_files = 'case45_sub.i'
  mode = batch-restore
[]
```

`SamplerFullSolveMultiApp` is the workhorse of the MC workflow. It creates as many sub-app instances as there are sampler rows and runs each to completion. The `batch-restore` mode is memory-efficient: it runs samples sequentially, restoring the sub-app state between runs rather than holding all 30 solutions in memory simultaneously.

**`[Transfers]`** — Wires parameters into sub-apps and retrieves results.

```
[param]
  type = SamplerParameterTransfer
  to_multi_app = sub
  sampler = mc
  parameters = 'Materials/thermal/prop_values'
[]
[results]
  type = SamplerReporterTransfer
  from_multi_app = sub
  sampler = mc
  stochastic_reporter = storage
  from_reporter = 'avg_T/value max_T/value'
[]
```

`SamplerParameterTransfer` injects the sampled k value directly into the sub-app's input parameter `Materials/thermal/prop_values` before each solve. `SamplerReporterTransfer` collects the postprocessor values back into the main app's storage reporter.

**`[Reporters]`** — A `StochasticReporter` acts as a container that accumulates all 30 sets of results into vectors for output.

**`[Outputs]`** — `csv = true` writes one CSV file per timestep containing the full results table.

### Sub-Application: `case45_sub.i`

A standard steady-state FE solve. The `[Controls]` block containing `SamplerReceiver` is the only UQ-specific addition — it listens for parameter updates from the main app and applies them before the solve begins. Everything else is a standard MOOSE heat conduction problem.

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" -w /work/case45-monte-carlo-uq --entrypoint /bin/bash idaholab/moose:latest -c '/opt/moose/bin/combined-opt -i case45_monte_carlo_uq.i 2>&1 | tail -40'
```

The run completes in a few seconds. MOOSE reports 30 sub-app solves, each converging in a single Newton iteration (the problem is linear in T). Output is written to `case45_monte_carlo_uq_out_storage_0001.csv.0`.

## Expected Results

All 30 samples converge successfully. The output CSV contains columns `results:avg_T:value` and `results:max_T:value`:

| Quantity | Minimum | Maximum | Mean |
|----------|---------|---------|------|
| avg_T    | 1.086   | 6.265   | 2.272 |
| max_T    | 2.525   | 14.561  | 5.280 |

The ratio max_T / avg_T is nearly constant at 2.32 across all samples — this makes sense because both quantities scale identically as 1/k. The factor 2.32 reflects the geometry of the 2D domain (the peak is always 2.32 times the average for this boundary condition and source distribution).

The distribution of avg_T is right-skewed because avg_T ~ 1/k and k ~ Uniform(0.5, 3.0). The reciprocal of a uniform distribution has a long upper tail: the sample with k = 0.529 (near the lower bound) gives avg_T = 6.26, while the sample with k = 2.886 gives avg_T = 1.09. If you had only run the nominal case k = 1.75 (midpoint), you would have obtained avg_T ≈ 1.79 — which underestimates the mean of the distribution by 21% because the mean of 1/k is not 1/mean(k).

## Key Takeaways

**Monte Carlo is embarrassingly parallel.** Each sample is an independent solve; the framework handles the launch, data transfer, and collection automatically. In production, `mode = batch-restore` can be replaced with `mode = batch-reset` to run samples truly in parallel using MPI.

**The `[StochasticTools]` action transforms a deterministic app into a UQ driver.** The sub-application `case45_sub.i` is completely unmodified physics — the only addition is the `SamplerReceiver` control object that receives parameter updates.

**Output distributions reveal what point estimates hide.** The mean of the output distribution (2.27) differs significantly from the output at the mean input (1.79) because the input-output map is nonlinear (inverse). This Jensen's inequality effect is a fundamental reason to run UQ rather than simply perturbing the nominal case.

**Right-skewed distributions have heavier tails than normal.** The maximum observed avg_T (6.26) is almost three times the mean. Engineering safety margins based on a normal-distribution assumption would substantially underestimate the probability of high-temperature excursions.

**30 samples gives rough statistics, not tight confidence intervals.** For a serious analysis, 1000+ samples are typical for one-dimensional uncertain parameters. The `num_rows` parameter is the only change needed to refine the study.
