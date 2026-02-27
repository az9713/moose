# ============================================================
# Case 73: Level Set Bubble Advection
# A smooth circular bubble defined by a level set function
# is advected by a solid-body rotation velocity field.
#
# Level set equation:
#   dphi/dt + v . nabla(phi) = 0
#
# Velocity field (solid-body rotation about (0.5, 0.5)):
#   v_x = -(y - 0.5)
#   v_y =  (x - 0.5)
#
# The bubble starts at (0.5, 0.75) with radius 0.15.
# After a half-rotation (t = pi), it reaches (0.5, 0.25).
# SUPG stabilization prevents oscillations in the pure
# advection equation.
#
# Tracked quantities:
#   total_phi: integral of phi (should be conserved)
#   max_phi: peak value (should stay near 1.0)
#   min_phi: trough value (should stay near 0.0)
#
# Domain: [0, 1]^2, 40x40 Q1 elements
# Time: t in [0, pi], dt = 0.1 (~31 steps)
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 40
  ny = 40
  xmin = 0.0
  xmax = 1.0
  ymin = 0.0
  ymax = 1.0
[]

[Variables]
  [phi]
  []
[]

[AuxVariables]
  [velocity]
    family = LAGRANGE_VEC
  []
[]

[ICs]
  # Smoothed bubble: tanh profile, radius 0.15 at (0.5, 0.75)
  [phi_ic]
    type = FunctionIC
    variable = phi
    function = phi_init
  []
  # Solid-body rotation velocity
  [vel_ic]
    type = VectorFunctionIC
    variable = velocity
    function = vel_func
  []
[]

[Functions]
  [phi_init]
    type = ParsedFunction
    expression = '0.5*(1.0 + tanh((0.15 - sqrt((x-0.5)^2+(y-0.75)^2))/0.02))'
  []
  [vel_func]
    type = ParsedVectorFunction
    expression_x = '-(y - 0.5)'
    expression_y = '(x - 0.5)'
  []
[]

[Kernels]
  # Standard Galerkin terms
  [time]
    type = TimeDerivative
    variable = phi
  []
  [advection]
    type = LevelSetAdvection
    variable = phi
    velocity = velocity
  []
  # SUPG stabilization terms
  [advection_supg]
    type = LevelSetAdvectionSUPG
    variable = phi
    velocity = velocity
  []
  [time_supg]
    type = LevelSetTimeDerivativeSUPG
    variable = phi
    velocity = velocity
  []
[]

[Postprocessors]
  # Total mass (integral of phi) — should be conserved
  [total_phi]
    type = ElementIntegralVariablePostprocessor
    variable = phi
    execute_on = 'initial timestep_end'
  []
  # Peak value — should stay near 1.0
  [max_phi]
    type = ElementExtremeValue
    variable = phi
    value_type = max
    execute_on = 'initial timestep_end'
  []
  # Minimum value — should stay near 0.0
  [min_phi]
    type = ElementExtremeValue
    variable = phi
    value_type = min
    execute_on = 'initial timestep_end'
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  dt = 0.1
  end_time = 3.14159
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-8
[]

[Outputs]
  exodus = true
  csv = true
[]
