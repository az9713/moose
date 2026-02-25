# ============================================================
# Case 13: Comprehensive Postprocessor-Driven Analysis
# Transient heat equation with multiple diagnostic quantities
# rho*cp*dT/dt = div(k*grad T) + Q,  T=0 on walls
# k=2, rho*cp=1, Q=5
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 30
  ny   = 30
[]

[Variables]
  [T]
  []
[]

[Kernels]
  [time_deriv]
    type     = TimeDerivative
    variable = T
  []
  [heat_conduction]
    type        = MatDiffusion
    variable    = T
    diffusivity = k
  []
  [heat_source]
    type     = BodyForce
    variable = T
    value    = 5.0   # Q = 5 W/m^3
  []
[]

[BCs]
  [zero_temp_walls]
    type     = DirichletBC
    variable = T
    boundary = 'left right top bottom'
    value    = 0
  []
[]

[Materials]
  [thermal_properties]
    type        = GenericConstantMaterial
    prop_names  = 'k'
    prop_values = '2.0'   # k = 2 W/(m K)
  []
[]

[Postprocessors]
  # Maximum temperature over the entire domain.
  # Useful for checking against safety limits.
  [max_temp]
    type       = ElementExtremeValue
    variable   = T
    value_type = max
  []

  # Volume-averaged temperature.
  # For a unit square: avg_T = (1/1) * int(T) dV
  [avg_temp]
    type     = ElementAverageValue
    variable = T
  []

  # Total thermal energy stored = int( rho*cp*T ) dV.
  # With rho*cp=1: this equals the integral of T.
  # ElementIntegralVariablePostprocessor computes int( T ) dV.
  [total_energy]
    type     = ElementIntegralVariablePostprocessor
    variable = T
  []

  # The L2 norm of the temperature field: sqrt( int(T^2) dV ).
  # Measures the "strength" of the field beyond just its average.
  [T_L2_norm]
    type     = ElementL2Norm
    variable = T
  []

  # Time step size (for reference in the CSV).
  [dt_size]
    type = TimestepSize
  []
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'

  [TimeStepper]
    type           = IterationAdaptiveDT
    dt             = 0.005
    optimal_iterations = 4
    growth_factor  = 2.0
    cutback_factor = 0.5
  []

  start_time = 0.0
  end_time   = 1.0
  nl_rel_tol = 1e-8
[]

[Outputs]
  # Write Exodus for spatial visualization.
  exodus = true

  # Write all postprocessors to a CSV file.
  # Columns: time, max_temp, avg_temp, total_energy, T_L2_norm, dt_size
  csv = true
[]
