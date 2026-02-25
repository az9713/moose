# ============================================================
# Case 12 Sub-Application
# Receives temperature T from parent as an AuxVariable,
# then solves its own diffusion problem for field 'phi'.
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 20
  ny   = 20
[]

[Variables]
  [phi]
  []
[]

[AuxVariables]
  # 'T_from_parent' is populated by the MultiAppCopyTransfer
  # in the parent input file.  It is read-only from the sub's perspective.
  [T_from_parent]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  # Diffusion equation for phi.
  [diffusion]
    type     = Diffusion
    variable = phi
  []

  # CoupledForce uses T_from_parent as a source term for phi.
  # This creates one-way coupling: T influences phi but not vice versa.
  [coupling]
    type = CoupledForce
    variable = phi
    v        = T_from_parent
    coef     = 0.1
  []
[]

[BCs]
  [walls]
    type     = DirichletBC
    variable = phi
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

[Outputs]
  exodus = true
[]
