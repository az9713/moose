# ============================================================
# Case 46: Polynomial Chaos Expansion — Surrogate Modeling
# Build a polynomial chaos surrogate for a 1D diffusion-reaction
# problem with two uncertain parameters:
#   D     ~ Uniform(2.5, 7.5)   (diffusivity)
#   sigma ~ Uniform(2.5, 7.5)   (absorption cross-section)
#
# Sub app: -D*u'' + sigma*u = 1,  u(0)=u(L)=0,  L=10
# Analytical: u(x) = (1/sigma)*(1 - cosh(sqrt(sigma/D)*(x-5))/cosh(5*sqrt(sigma/D)))
#
# Training: Quadrature sampler (order 5) → 36 deterministic solves
# Evaluation: 100 MC samples evaluated via surrogate (zero PDE solves)
#
# This case shows how to train a cheap surrogate that replaces
# expensive full PDE solves for uncertainty propagation.
# ============================================================

[StochasticTools]
[]

[Distributions]
  [D_dist]
    type = Uniform
    lower_bound = 2.5
    upper_bound = 7.5
  []
  [S_dist]
    type = Uniform
    lower_bound = 2.5
    upper_bound = 7.5
  []
[]

[Samplers]
  # Deterministic quadrature points for training the surrogate
  [quadrature]
    type = Quadrature
    distributions = 'D_dist S_dist'
    execute_on = INITIAL
    order = 5
  []
  # Random MC samples to evaluate the trained surrogate
  [mc_eval]
    type = MonteCarlo
    num_rows = 100
    distributions = 'D_dist S_dist'
    execute_on = timestep_end
    seed = 2024
  []
[]

[MultiApps]
  [sub]
    type = SamplerFullSolveMultiApp
    input_files = 'case46_sub.i'
    sampler = quadrature
    mode = batch-restore
  []
[]

[Transfers]
  [param]
    type = SamplerParameterTransfer
    to_multi_app = sub
    sampler = quadrature
    parameters = 'Materials/diffusivity/prop_values Materials/sigma_input/prop_values'
  []
  [data]
    type = SamplerReporterTransfer
    from_multi_app = sub
    sampler = quadrature
    stochastic_reporter = storage
    from_reporter = 'avg/value'
  []
[]

[Reporters]
  [storage]
    type = StochasticReporter
  []
  # Evaluate the trained surrogate on MC samples (no PDE solves!)
  [eval]
    type = EvaluateSurrogate
    model = pc_surrogate
    sampler = mc_eval
    parallel_type = ROOT
    execute_on = final
  []
[]

[Surrogates]
  [pc_surrogate]
    type = PolynomialChaos
    trainer = pc_trainer
  []
[]

[Trainers]
  [pc_trainer]
    type = PolynomialChaosTrainer
    execute_on = timestep_end
    order = 5
    distributions = 'D_dist S_dist'
    sampler = quadrature
    response = storage/data:avg:value
  []
[]

[Outputs]
  [out]
    type = CSV
    execute_on = FINAL
  []
[]
