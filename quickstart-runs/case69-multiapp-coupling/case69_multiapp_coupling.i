# ============================================================
# Case 69: MultiApp Coupled Diffusion — Bidirectional Picard
# Demonstrates MOOSE's MultiApp system for multi-physics
# coupling via Picard (fixed-point) iteration.
#
# Main app solves heat conduction with a source from sub-app:
#   dT/dt = nabla^2 T + q(x,y)
#
# Sub-app solves a reaction-diffusion for heat source q:
#   dq/dt = nabla^2 q - q + S(x,y) - 0.5*T
#
# The negative feedback (-0.5*T) ensures the coupled system
# reaches a stable equilibrium: as T rises, q drops.
#
# Transfer mechanism:
#   Main -> Sub: T field (drives feedback term)
#   Sub -> Main: q field (drives heat source)
#
# Domain: [0, 1] x [0, 1], 20x20 elements
# Time: t in [0, 2], dt = 0.2
# Picard tolerance: 1e-6 (up to 10 iterations per step)
#
# Expected: T peaks at center (~0.8) where Gaussian source
# is strongest; q settles to feedback-limited steady value.
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 20
  ny = 20
[]

[Variables]
  [T]
    initial_condition = 0
  []
[]

[AuxVariables]
  [source_from_sub]
  []
[]

[Kernels]
  [time]
    type = TimeDerivative
    variable = T
  []
  [diffusion]
    type = Diffusion
    variable = T
  []
  # Heat source transferred from sub-app
  [source]
    type = CoupledForce
    variable = T
    v = source_from_sub
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = T
    boundary = left
    value = 0
  []
  [right]
    type = DirichletBC
    variable = T
    boundary = right
    value = 0
  []
  [top]
    type = DirichletBC
    variable = T
    boundary = top
    value = 0
  []
  [bottom]
    type = DirichletBC
    variable = T
    boundary = bottom
    value = 0
  []
[]

[MultiApps]
  [sub]
    type = TransientMultiApp
    positions = '0 0 0'
    input_files = case69_sub.i
    execute_on = 'timestep_end'
  []
[]

[Transfers]
  # q from sub-app becomes the heat source in main
  [source_from_sub]
    type = MultiAppNearestNodeTransfer
    from_multi_app = sub
    source_variable = q
    variable = source_from_sub
  []
  # T from main drives the feedback in sub-app
  [T_to_sub]
    type = MultiAppNearestNodeTransfer
    to_multi_app = sub
    source_variable = T
    variable = T_from_main
  []
[]

[Postprocessors]
  [T_max]
    type = ElementExtremeValue
    variable = T
    execute_on = 'initial timestep_end'
  []
  [T_avg]
    type = ElementAverageValue
    variable = T
    execute_on = 'initial timestep_end'
  []
  [source_avg]
    type = ElementAverageValue
    variable = source_from_sub
    execute_on = 'initial timestep_end'
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
  nl_rel_tol = 1e-8
  fixed_point_max_its = 10
  fixed_point_rel_tol = 1e-6
  fixed_point_abs_tol = 1e-8
[]

[Outputs]
  exodus = true
  csv = true
[]
