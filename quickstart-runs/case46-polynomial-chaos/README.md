# Case 46: Polynomial Chaos Expansion — Surrogate Modeling

## Overview

Monte Carlo UQ (Case 45) works well when each simulation is cheap, but realistic engineering models can take hours or days to run. If you need thousands of samples to characterize uncertainty, direct MC becomes infeasible. Surrogate modeling offers a way out: run a small, carefully chosen set of full simulations, fit a cheap algebraic approximation to the results, then evaluate the surrogate for as many samples as needed at negligible cost.

This case demonstrates Polynomial Chaos Expansion (PCE), one of the most mathematically principled surrogate methods. PCE expresses the model output as a sum of orthogonal polynomials in the uncertain inputs. For uniform distributions, the basis functions are Legendre polynomials. The key insight is that PCE coefficients can be computed exactly (not approximately) from a deterministic quadrature rule — no randomness involved in the training phase at all.

The problem is 1D diffusion-reaction with two uncertain parameters: diffusivity D and absorption cross-section sigma, both uniform on [2.5, 7.5]. A degree-5 PCE requires 36 quadrature points in 2D. After training, the surrogate evaluates 100 new Monte Carlo samples in microseconds rather than running 100 additional PDE solves.

## The Physics

The sub-application solves a 1D diffusion-reaction boundary value problem:

```
-D * u'' + sigma * u = 1    on [0, 10]
u(0) = u(10) = 0
```

This models neutron diffusion in a slab (or equivalently, heat conduction with linear absorption): D is the diffusion coefficient, sigma is the macroscopic absorption cross-section, and the unit source represents a spatially uniform neutron source or heat generation. The homogeneous Dirichlet conditions model a perfectly absorbing boundary.

The analytical solution is:

```
u(x) = (1/sigma) * (1 - cosh(sqrt(sigma/D) * (x-5)) / cosh(5 * sqrt(sigma/D)))
```

The domain-average of this solution is the quantity of interest. It decreases with increasing sigma (stronger absorption removes neutrons faster) and increases with increasing D (higher diffusivity spreads neutrons further from boundaries). The interplay of these two effects makes the average a non-trivially nonlinear function of (D, sigma).

Uncertain parameters:
- D ~ Uniform(2.5, 7.5) (diffusion coefficient)
- sigma ~ Uniform(2.5, 7.5) (absorption cross-section)

## Input File Walkthrough

### Main Orchestrator: `case46_polynomial_chaos.i`

**`[Distributions]`** — Declares two independent uniform distributions, one for D and one for sigma. The PCE basis functions are automatically matched to the distribution type (Legendre polynomials for uniform distributions, Hermite polynomials for normal distributions).

**`[Samplers]`** — Two samplers serve different roles:

```
[quadrature]
  type = Quadrature
  distributions = 'D_dist S_dist'
  execute_on = INITIAL
  order = 5
[]
[mc_eval]
  type = MonteCarlo
  num_rows = 100
  distributions = 'D_dist S_dist'
  execute_on = timestep_end
  seed = 2024
[]
```

`Quadrature` generates deterministic Gauss-Legendre quadrature points for training — 6 points per dimension at order 5, giving 6x6 = 36 total. These points are chosen to minimize the polynomial fitting error, not drawn randomly. `MonteCarlo` generates 100 random samples that will be evaluated by the trained surrogate. Note that `execute_on = INITIAL` for the quadrature sampler and `execute_on = timestep_end` for the MC sampler — they execute at different phases.

**`[MultiApps]` and `[Transfers]`** — Only the quadrature sampler drives actual PDE solves. The MC sampler is never connected to a MultiApp because it evaluates the surrogate, not the full model. The parameter transfer injects both D and sigma into the sub-app's material objects.

**`[Trainers]`** — The PCE training step:

```
[pc_trainer]
  type = PolynomialChaosTrainer
  execute_on = timestep_end
  order = 5
  distributions = 'D_dist S_dist'
  sampler = quadrature
  response = storage/data:avg:value
[]
```

`PolynomialChaosTrainer` collects the 36 quadrature responses and computes the PCE coefficients by numerical integration. The `order = 5` degree-5 expansion has (5+1)^2 = 36 terms in 2D, exactly matching the quadrature rule — this is a complete expansion at the Gauss-Legendre quadrature points.

**`[Surrogates]`** — The trained surrogate is stored as a `PolynomialChaos` object:

```
[pc_surrogate]
  type = PolynomialChaos
  trainer = pc_trainer
[]
```

**`[Reporters]`** — The `EvaluateSurrogate` reporter evaluates the trained PCE at the 100 MC sample points:

```
[eval]
  type = EvaluateSurrogate
  model = pc_surrogate
  sampler = mc_eval
  parallel_type = ROOT
  execute_on = final
[]
```

This runs entirely in memory — no sub-app launches, no PDE solves. The `execute_on = final` ensures the surrogate is fully trained before evaluation begins.

### Sub-Application: `case46_sub.i`

The sub-app implements the diffusion-reaction equation. A subtlety arises from MOOSE's kernel sign conventions:

`MatReaction` computes the residual contribution as `-rate * u * test`. To get the physically correct `+sigma * u * test` (the absorption term), the reaction rate must be negative. Rather than hard-coding a negative value (which cannot be sampled), a `DerivativeParsedMaterial` computes `neg_sigma = -sigma` from the material property `sigma`. This way the sampler sets the positive `sigma`, and the kernel uses the negated value:

```
[sigma_negate]
  type = DerivativeParsedMaterial
  property_name = neg_sigma
  expression = '-sigma'
  material_property_names = 'sigma'
  disable_fpoptimizer = true
  enable_jit = false
[]
```

The `disable_fpoptimizer = true` and `enable_jit = false` flags are required for Docker compatibility (the JIT compiler is unavailable in the container; see the Docker guide).

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" -w /work/case46-polynomial-chaos --entrypoint /bin/bash idaholab/moose:latest -c '/opt/moose/bin/combined-opt -i case46_polynomial_chaos.i 2>&1 | tail -40'
```

MOOSE runs 36 quadrature sub-app solves during the training phase, then evaluates the surrogate on 100 MC points. Output is written to `case46_polynomial_chaos_out_eval_0002.csv` (surrogate evaluations) and `case46_polynomial_chaos_out_storage_0002.csv.0` (training data).

## Expected Results

The surrogate evaluation on 100 Monte Carlo samples produces:

| Quantity | Value |
|----------|-------|
| Minimum avg(u) | 0.1091 |
| Maximum avg(u) | 0.3080 |
| Mean avg(u)    | 0.1772 |
| Std dev        | 0.0550 |

The output ranges from about 0.11 to 0.31. The asymmetry between minimum and maximum (0.11 to 0.31 around a mean of 0.18) reflects the nonlinear dependence on D and sigma: the domain average is more sensitive to sigma at low sigma values than at high sigma values.

The computational savings are substantial: 36 training solves replace what would otherwise require 1000+ MC solves to achieve comparable accuracy in estimating the mean and standard deviation. Moreover, once trained, the surrogate can be evaluated for any number of additional samples essentially for free.

To verify surrogate accuracy, you could compare the surrogate predictions against direct solves at a handful of the MC sample points. For a degree-5 PCE on two uniform parameters, the surrogate typically matches the full model to better than 0.1% for smooth response functions like this one.

## Key Takeaways

**Polynomial chaos is a spectral method for uncertainty.** Just as Fourier series represent periodic functions in terms of sine and cosine basis functions, PCE represents uncertain model outputs in terms of orthogonal polynomial basis functions matched to the input distributions. The coefficients carry all statistical information.

**Quadrature training is deterministic, not random.** The 36 training points are the Gauss-Legendre quadrature nodes for the 2D uniform distribution — the same nodes used in numerical integration. There is no sampling variability in the training step.

**The surrogate captures the full response surface, not just statistics.** Unlike simple Monte Carlo which gives you samples from the output distribution, the PCE surrogate gives you an explicit formula for the output as a function of inputs. You can extract mean, variance, sensitivity indices (Sobol indices), and conditional distributions all analytically from the coefficients.

**Two samplers, one MultiApp.** Only the quadrature sampler drives PDE solves; the MC evaluation sampler queries the trained surrogate in memory. This clean separation of training and evaluation is the fundamental architecture of surrogate-based UQ.

**The `[StochasticTools]` + `[Trainers]` + `[Surrogates]` pattern is composable.** The same pattern extends to Gaussian process surrogates, neural network surrogates, and active learning workflows. Once you understand this architecture, adding more sophisticated surrogate types requires only changing the `type` parameter.
