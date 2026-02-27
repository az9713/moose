# ============================================================
# Case 47: Heat Source Inversion — PDE-Constrained Optimization
# Recover an unknown volumetric heat source q from sparse
# temperature measurements using adjoint-based optimization.
#
# Forward problem: -k*laplacian(T) = q  on [0,2]x[0,2]
#   T = 200 (bottom), T = 100 (top), insulated sides
#   k = 5 W/(m*K)
#
# Measurements: 4 sensors at interior points
# True source: q ~ 1000 W/m^3  (to be recovered)
# Initial guess: q = 500
#
# Architecture:
#   main.i    — TAO optimizer (L-BFGS), drives forward + adjoint
#   forward.i — solves heat equation with current q, computes misfit
#   adjoint.i — solves adjoint equation, computes gradient dJ/dq
#
# This case introduces the optimization module and adjoint methods.
# ============================================================

[Optimization]
[]

[OptimizationReporter]
  type = GeneralOptimization
  objective_name = objective_value
  parameter_names = 'parameter_results'
  num_values = '1'
  initial_condition = '500'
  lower_bounds = '0.1'
  upper_bounds = '10000'
[]

[Reporters]
  [main]
    type = OptimizationData
    measurement_points = '0.2 0.2 0
                          0.8 0.6 0
                          0.2 1.4 0
                          0.8 1.8 0'
    measurement_values = '226 254 214 146'
  []
[]

[Executioner]
  type = Optimize
  tao_solver = taoblmvm
  petsc_options_iname = '-tao_gatol -tao_max_it'
  petsc_options_value = '1e-4       50'
  verbose = true
[]

[MultiApps]
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
[]

[Transfers]
  # Send measurement data and current parameter to forward app
  [toForward]
    type = MultiAppReporterTransfer
    to_multi_app = forward
    from_reporters = 'main/measurement_xcoord
                      main/measurement_ycoord
                      main/measurement_zcoord
                      main/measurement_time
                      main/measurement_values
                      OptimizationReporter/parameter_results'
    to_reporters = 'measure_data/measurement_xcoord
                    measure_data/measurement_ycoord
                    measure_data/measurement_zcoord
                    measure_data/measurement_time
                    measure_data/measurement_values
                    params/q'
  []
  # Send misfit data and current parameter to adjoint app
  [toAdjoint]
    type = MultiAppReporterTransfer
    to_multi_app = adjoint
    from_reporters = 'main/measurement_xcoord
                      main/measurement_ycoord
                      main/measurement_zcoord
                      main/measurement_time
                      main/misfit_values
                      OptimizationReporter/parameter_results'
    to_reporters = 'misfit/measurement_xcoord
                    misfit/measurement_ycoord
                    misfit/measurement_zcoord
                    misfit/measurement_time
                    misfit/misfit_values
                    params/q'
  []
  # Collect objective value and misfit from forward app
  [fromForward]
    type = MultiAppReporterTransfer
    from_multi_app = forward
    from_reporters = 'measure_data/misfit_values measure_data/objective_value'
    to_reporters = 'main/misfit_values OptimizationReporter/objective_value'
  []
  # Collect gradient from adjoint app
  [fromAdjoint]
    type = MultiAppReporterTransfer
    from_multi_app = adjoint
    from_reporters = 'gradient_vpp/inner_product'
    to_reporters = 'OptimizationReporter/grad_parameter_results'
  []
[]

[Outputs]
  csv = true
[]
