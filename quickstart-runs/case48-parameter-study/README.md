# Case 48: Latin Hypercube Parameter Study — Multi-Parameter UQ

## Overview

Case 45 showed the full manual assembly of a Monte Carlo UQ workflow: distributions, samplers, MultiApps, transfers, and reporters, each declared explicitly. For rapid exploratory analysis — where you want to run a parameter sweep quickly without constructing the full plumbing by hand — MOOSE provides the `[ParameterStudy]` action, which generates all of those objects automatically from a compact high-level specification.

This case runs a Latin Hypercube Sampling (LHS) study with three simultaneously uncertain parameters: thermal conductivity k, left boundary temperature T_left, and right boundary temperature T_right. The entire main input file is 15 lines. The same 50-sample study implemented manually (as in Case 45) would require approximately 70 lines covering distributions, samplers, MultiApps, transfers, and reporters.

Latin Hypercube Sampling is a stratified sampling strategy that guarantees better coverage of the parameter space than pure random Monte Carlo. The parameter space is divided into N equal-probability intervals along each axis, and exactly one sample is drawn from each interval — ensuring that no region of the input space is accidentally missed. For the same number of samples, LHS typically produces lower variance statistics than Monte Carlo.

## The Physics

The sub-application solves 2D steady heat conduction with a volumetric source on the unit square:

```
-k * laplacian(T) = 100    on [0,1] x [0,1]
T = T_left                 on left boundary  (x = 0)
T = T_right                on right boundary (x = 1)
dT/dn = 0                  on top and bottom (insulated)
```

Unlike Case 45 (which had T = 0 on all boundaries), this problem has non-zero and asymmetric Dirichlet conditions on left and right, with insulated top and bottom. The solution interpolates between T_left and T_right while the source term q = 100 W/m^3 adds a parabolic offset that depends on k.

The approximate analytical solution for the 1D version (ignoring the y-variation, which vanishes because top and bottom are insulated) is:

```
T(x) = T_left + (T_right - T_left)*x + (q/(2k)) * x*(1-x)
```

This shows that:
- The domain average T scales as (T_left + T_right)/2 + q/(12k)
- T_right has the dominant influence on the average because the right BC appears directly in the average
- k matters through the source term: lower k amplifies the parabolic source contribution

Uncertain parameters:
- k ~ Uniform(0.5, 5.0) W/(m·K)
- T_left ~ Uniform(0.0, 100.0) K
- T_right ~ Uniform(200.0, 400.0) K

Outputs of interest: domain-average temperature avg_T and peak temperature max_T.

## Input File Walkthrough

### Main File: `case48_parameter_study.i`

The entire main file is a single `[ParameterStudy]` block:

```
[ParameterStudy]
  input = case48_sub.i
  parameters = 'Materials/thermal/prop_values BCs/left/value BCs/right/value'
  quantities_of_interest = 'avg_T/value max_T/value'

  sampling_type = lhs
  num_samples = 50
  distributions = 'uniform uniform uniform'
  uniform_lower_bound = '0.5 0.0 200'
  uniform_upper_bound = '5.0 100.0 400'
  seed = 2024

  output_type = csv
[]
```

The `parameters` list maps directly to sub-app HIT paths, exactly as in the manual `SamplerParameterTransfer`. The action reads this list and creates the appropriate distributions, sampler, MultiApp, transfers, and reporters internally during input parsing — before the simulation even begins.

`sampling_type = lhs` selects Latin Hypercube Sampling. Other options include `mc` (pure Monte Carlo) and `sobol` (quasi-random Sobol sequences). The `distributions = 'uniform uniform uniform'` tells the action which distribution type to use for each parameter, and the `uniform_lower_bound` and `uniform_upper_bound` vectors provide the bounds.

`quantities_of_interest = 'avg_T/value max_T/value'` specifies which reporter values to collect from each sub-app run. The action automatically creates the `SamplerReporterTransfer` and `StochasticReporter` to accumulate these.

`output_type = csv` writes the results to a CSV file. The action handles output setup as well.

This compactness is the primary value of the `[ParameterStudy]` action: for standard parameter studies, there is no reason to write the underlying plumbing by hand.

### Sub-Application: `case48_sub.i`

The sub-app is a standard 2D heat conduction solve. Compared to Case 45's sub-app, the differences are:

1. The left and right boundaries now have their own `DirichletBC` objects (with separate names `left` and `right`) so the `ParameterStudy` can target `BCs/left/value` and `BCs/right/value` independently.
2. Top and bottom have no BC block — MOOSE's natural boundary condition (zero Neumann) applies by default, giving the insulated condition.
3. The `[Controls]` block with `SamplerReceiver` is still present — the `ParameterStudy` action still uses the same underlying transfer mechanism.

The postprocessors `avg_T` (AverageNodalVariableValue) and `max_T` (NodalExtremeValue) match the `quantities_of_interest` declared in the main file.

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" -w /work/case48-parameter-study --entrypoint /bin/bash idaholab/moose:latest -c '/opt/moose/bin/combined-opt -i case48_parameter_study.i 2>&1 | tail -40'
```

Output is written to `case48_parameter_study_csv_study_results_0001.csv`. Each row corresponds to one LHS sample and contains the three input parameter values and the two output quantities.

## Expected Results

All 50 samples converge. The output CSV contains columns for `BCs_left_value`, `BCs_right_value`, `Materials_thermal_prop_values`, `avg_T:value`, and `max_T:value`.

| Quantity | Minimum | Maximum | Mean |
|----------|---------|---------|------|
| avg_T    | 105.1   | 248.3   | 179.1 |
| max_T    | 203.0   | 399.7   | 300.1 |

The dominant observation is that max_T equals T_right in every single sample — the maximum temperature always occurs at the right boundary. This makes physical sense: T_right ~ Uniform(200, 400) sets the right boundary condition, and that boundary is always hotter than the interior (which is cooled toward T_left ~ Uniform(0, 100) on the left and has a modest source term). The peak of the domain is the right boundary point.

Examining the correlation structure:
- avg_T shows strong positive correlation with T_right (correlation ≈ 0.95), because T_right directly enters the domain-average through the boundary condition
- avg_T shows moderate positive correlation with T_left (correlation ≈ 0.5), because T_left also enters the average
- avg_T shows a secondary negative correlation with k (higher k reduces the source-driven temperature rise)
- The k effect is secondary because the source term q = 100 produces a source contribution of order q/(12k) ≈ 2–20 K — significant but smaller than the spread in T_left (0–100 K) and T_right (200–400 K)

The LHS strategy guarantees that the 50 samples cover the three-dimensional parameter space without gaps. You can verify this by examining the input columns: each dimension is divided into 50 equal intervals, and exactly one sample falls in each interval.

## Key Takeaways

**`[ParameterStudy]` condenses 70 lines to 15.** For standard parameter studies with uniform or normal distributions, the action replaces all manual plumbing. Use it for rapid exploration; switch to manual assembly only when you need non-standard transfers, custom samplers, or surrogate training.

**Latin Hypercube guarantees space coverage that Monte Carlo does not.** With 50 pure Monte Carlo samples in 3D, you might accidentally cluster samples in one corner of the parameter space. LHS stratification prevents this: each parameter's full range is sampled uniformly across all 50 draws.

**The dominant parameter is identifiable from the output CSV.** A scatter plot of avg_T vs. T_right shows a near-linear relationship with slope ≈ 0.5 (consistent with the 1D formula: the right BC contributes (T_right)/2 to the average). The k effect is visible only in the residuals after removing the T_right trend.

**Multiple uncertain parameters expose interactions.** In this case the three parameters affect the output approximately additively (no strong cross-terms), but in nonlinear problems the interaction terms can be the dominant source of output variability. LHS samples are sufficient to detect first-order effects; detecting interactions requires either more samples or a carefully designed factorial experiment.

**The same sub-app works with no modifications for any sampling strategy.** The `SamplerReceiver` control is agnostic to whether the upstream sampler is Monte Carlo, LHS, Sobol, or Quadrature. Switching between `[ParameterStudy]` with `sampling_type = lhs` and `sampling_type = mc` requires changing a single word in the main file.
