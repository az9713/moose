# ============================================================
# Case 45: Monte Carlo UQ — Uncertainty in Thermal Conductivity
# Propagate uncertainty in thermal conductivity through a 2D
# steady heat conduction problem using Monte Carlo sampling.
#
# Main app: orchestrates 30 MC samples via SamplerFullSolveMultiApp.
# Sub app (case45_sub.i): solves -k*laplacian(T) = 100 on [0,1]^2
# with T=0 on all boundaries. Solution T ~ q/(2k), so lower k
# gives higher temperatures.
#
# k ~ Uniform(0.5, 3.0)
#
# This case introduces the stochastic_tools module and the
# MultiApp sampling architecture.
# ============================================================

[StochasticTools]
[]

[Distributions]
  [k_dist]
    type = Uniform
    lower_bound = 0.5
    upper_bound = 3.0
  []
[]

[Samplers]
  [mc]
    type = MonteCarlo
    num_rows = 30
    distributions = 'k_dist'
    seed = 2024
  []
[]

[MultiApps]
  [sub]
    type = SamplerFullSolveMultiApp
    sampler = mc
    input_files = 'case45_sub.i'
    mode = batch-restore
  []
[]

[Transfers]
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
[]

[Reporters]
  [storage]
    type = StochasticReporter
  []
[]

[Outputs]
  csv = true
[]
