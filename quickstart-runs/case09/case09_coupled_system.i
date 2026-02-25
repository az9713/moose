# ============================================================
# Case 9: Coupled Two-Variable Reaction-Diffusion
# du/dt = Du*div(grad u) + v
# dv/dt = Dv*div(grad v) - u
# Du=1, Dv=0.5, zero initial conditions, u=1 on left face
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 30
  ny   = 30
[]

[Variables]
  [u]
  []
  [v]
  []
[]

[Kernels]
  # --- Equation for u ---
  [u_time]
    type     = TimeDerivative
    variable = u
  []
  [u_diff]
    type     = MatDiffusion
    variable = u
    diffusivity = Du   # material property
  []
  # CoupledForce adds  int( phi_i * v ) dV  to the equation for u.
  # 'v' here is the MOOSE variable v defined in [Variables].
  [u_source]
    type     = CoupledForce
    variable = u
    v        = v       # the coupling variable
  []

  # --- Equation for v ---
  [v_time]
    type     = TimeDerivative
    variable = v
  []
  [v_diff]
    type     = MatDiffusion
    variable = v
    diffusivity = Dv
  []
  # CoupledForce with a negative coefficient adds  -int( phi_i * u ) dV.
  [v_sink]
    type        = CoupledForce
    variable    = v
    v           = u    # couples u into the v equation
    coef        = -1.0 # negative: u acts as a sink for v
  []
[]

[BCs]
  # u = 1 on the left boundary (drives the system).
  [u_left]
    type     = DirichletBC
    variable = u
    boundary = left
    value    = 1.0
  []
  # u = 0 on the right boundary.
  [u_right]
    type     = DirichletBC
    variable = u
    boundary = right
    value    = 0.0
  []
  # v = 0 on all walls.
  [v_walls]
    type     = DirichletBC
    variable = v
    boundary = 'left right top bottom'
    value    = 0.0
  []
[]

[Materials]
  [diffusivities]
    type        = GenericConstantMaterial
    prop_names  = 'Du  Dv'
    prop_values = '1.0 0.5'
  []
[]

[Postprocessors]
  [avg_u]
    type     = ElementAverageValue
    variable = u
  []
  [avg_v]
    type     = ElementAverageValue
    variable = v
  []
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'

  dt       = 0.01
  end_time = 2.0

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 15

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    growth_factor = 2.0
    cutback_factor = 0.5
    optimal_iterations = 8
  []
[]

[Outputs]
  exodus = true
  csv    = true
[]
