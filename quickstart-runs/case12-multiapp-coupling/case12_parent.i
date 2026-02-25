# ============================================================
# Case 12 Parent Application
# Solves thermal problem, then hands T to the sub-application.
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 20
  ny   = 20
[]

[Variables]
  [T]
  []
[]

[Kernels]
  [diffusion]
    type     = Diffusion
    variable = T
  []
  [source]
    type     = BodyForce
    variable = T
    value    = 1.0
  []
[]

[BCs]
  [walls]
    type     = DirichletBC
    variable = T
    boundary = 'left right top bottom'
    value    = 0
  []
[]

[Executioner]
  type       = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]

[MultiApps]
  # FullSolveMultiApp runs the sub-application to convergence
  # once at the specified execute_on point.
  [thermal_sub]
    type        = FullSolveMultiApp
    input_files = case12_sub.i    # path to the sub-app input file
    execute_on  = timestep_end    # run after the parent solves
    positions   = '0 0 0'        # sub-app origin in parent coordinates
  []
[]

[Transfers]
  # MultiAppCopyTransfer copies nodal values from one app to another.
  # The source variable (T in parent) must match the target mesh topology.
  [send_temperature]
    type          = MultiAppCopyTransfer
    to_multi_app  = thermal_sub    # direction: parent -> sub
    source_variable = T            # variable in the parent
    variable        = T_from_parent  # AuxVariable in the sub
  []
[]

[Outputs]
  exodus = true
[]
