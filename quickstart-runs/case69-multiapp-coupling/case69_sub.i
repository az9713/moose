# ============================================================
# Case 69 Sub-App: Heat Source with Temperature Feedback
# Computes heat source field q that responds to temperature T
# transferred from the main app.
#
# Governing equation:
#   dq/dt = nabla^2 q - q + S(x,y) - 0.5*T
#
# S(x,y) = 10 * exp(-((x-0.5)^2 + (y-0.5)^2) / 0.05)
#   Gaussian source centered at (0.5, 0.5)
#
# The -q term (first-order decay) and -0.5*T (feedback)
# ensure the source reaches a bounded steady state.
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 20
  ny = 20
[]

[Variables]
  [q]
    initial_condition = 0
  []
[]

[AuxVariables]
  [T_from_main]
  []
[]

[Kernels]
  [time]
    type = TimeDerivative
    variable = q
  []
  [diffusion]
    type = Diffusion
    variable = q
  []
  # First-order decay: -q (Reaction kernel adds +q to residual)
  [decay]
    type = Reaction
    variable = q
  []
  # Gaussian heat generation
  [source]
    type = BodyForce
    variable = q
    function = source_func
  []
  # Negative feedback from temperature
  # CoupledForce adds -coef*v to residual => coef=-0.5 gives +0.5*T in residual
  # This means the RHS gets -0.5*T (higher T reduces q)
  [feedback]
    type = CoupledForce
    variable = q
    v = T_from_main
    coef = -0.5
  []
[]

[Functions]
  [source_func]
    type = ParsedFunction
    expression = '10*exp(-((x-0.5)^2+(y-0.5)^2)/0.05)'
  []
[]

[BCs]
  [all]
    type = DirichletBC
    variable = q
    boundary = 'left right top bottom'
    value = 0
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  dt = 0.2
  end_time = 2.0
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
[]
