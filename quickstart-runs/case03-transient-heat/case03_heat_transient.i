# ============================================================
# Case 3: Transient Heat Equation
# rho*cp * dT/dt = div(k*grad T) + Q
# rho=cp=k=1, Q=1, T=0 on all walls, T_0=0
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 20
  ny   = 20
[]

[Variables]
  # T is temperature (units: K in a real problem; dimensionless here)
  [T]
  []
[]

[Kernels]
  # TimeDerivative contributes  rho*cp * dT/dt  to the residual.
  # Because rho=cp=1 we do not need a coefficient.
  [time_deriv]
    type     = TimeDerivative
    variable = T
  []

  # MatDiffusion reads the conductivity from a material property named 'k'.
  # This is more general than bare Diffusion: it allows k to vary in space.
  [heat_conduction]
    type        = MatDiffusion
    variable    = T
    diffusivity = k    # material property name (defined below)
  []

  # BodyForce applies a volumetric source term:  int( phi_i * Q ) dV
  [heat_source]
    type     = BodyForce
    variable = T
    value    = 1.0   # Q = 1 W/m^3 (uniform, constant)
  []
[]

[BCs]
  # Zero temperature on all four walls
  [walls]
    type     = DirichletBC
    variable = T
    boundary = 'left right top bottom'
    value    = 0
  []
[]

[Materials]
  # GenericConstantMaterial declares one or more scalar material properties
  # with constant values.  'prop_names' and 'prop_values' must be the same length.
  [thermal_props]
    type        = GenericConstantMaterial
    prop_names  = 'k'     # thermal conductivity
    prop_values = '1.0'   # W/(m K)
  []
[]

[Postprocessors]
  # ElementAverageValue computes the spatial average of a variable:
  #   (1/|Omega|) * int( T ) dV
  # Printed to screen and written to the CSV file each timestep.
  [avg_temperature]
    type     = ElementAverageValue
    variable = T
  []

  # ElementExtremeValue tracks the maximum of T over all quadrature points.
  [max_temperature]
    type          = ElementExtremeValue
    variable      = T
    value_type    = max   # 'max' or 'min'
  []
[]

[Executioner]
  type = Transient   # time-marching solve

  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'

  # ConstantDT uses a fixed timestep size.
  [TimeStepper]
    type = ConstantDT
    dt   = 0.01    # 10 ms per step
  []

  # March from t=0 to t=0.5.
  start_time = 0.0
  end_time   = 0.5

  # Tolerances for Newton convergence.
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true   # spatial fields at each timestep
  csv    = true   # postprocessor history (avg_temp, max_temp vs. time)
[]
