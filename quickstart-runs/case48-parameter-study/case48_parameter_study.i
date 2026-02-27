# ============================================================
# Case 48: Latin Hypercube Parameter Study — Multi-Parameter UQ
# Uses the high-level ParameterStudy action to run a Latin
# Hypercube sampling study with just a few lines of input.
#
# Sub app: 2D steady heat conduction on [0,1]x[0,1]
# Uncertain parameters:
#   k       ~ Uniform(0.5, 5.0)   thermal conductivity
#   T_left  ~ Uniform(0.0, 100.0) left boundary temperature
#   T_right ~ Uniform(200, 400)   right boundary temperature
#
# 50 LHS samples, collect avg_T and max_T.
#
# This case shows the simplest possible UQ workflow in MOOSE:
# the ParameterStudy action creates all samplers, multiapps,
# transfers, and reporters automatically.
# ============================================================

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
