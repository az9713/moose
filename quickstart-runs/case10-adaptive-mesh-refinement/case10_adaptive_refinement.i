# ============================================================
# Case 10: Adaptive Mesh Refinement (AMR)
# Laplace equation with a boundary discontinuity that
# creates a steep gradient -> drives refinement near the corner
# ============================================================

[Mesh]
  # Start coarse; AMR will add resolution where needed.
  [gen]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 8
    ny   = 8
  []
[]

[Variables]
  [u]
  []
[]

[Kernels]
  [diffusion]
    type     = Diffusion
    variable = u
  []
[]

[BCs]
  # u=0 on left and bottom; u=1 on right.
  # The corner at (1,0) is a singular point and drives refinement.
  [left]
    type     = DirichletBC
    variable = u
    boundary = left
    value    = 0
  []
  [bottom]
    type     = DirichletBC
    variable = u
    boundary = bottom
    value    = 0
  []
  [right]
    type     = DirichletBC
    variable = u
    boundary = right
    value    = 1
  []
  # top: zero-flux Neumann (natural, no entry needed)
[]

[Adaptivity]
  # 'marker' is the name of the Marker to use for steady-state adaptation.
  marker = err_marker

  # Perform up to 4 cycles of refinement before the solve.
  initial_steps  = 4
  initial_marker = err_marker

  # After solving, refine once more.
  steps = 1

  # Maximum number of refinement levels relative to the base mesh.
  max_h_level = 5

  [Indicators]
    # GradientJumpIndicator estimates error from the jump in grad(u)
    # across element faces. Large jumps indicate large error.
    [jump_indicator]
      type     = GradientJumpIndicator
      variable = u
    []
  []

  [Markers]
    # ErrorFractionMarker refines the top 'refine' fraction of elements
    # by error and coarsens the bottom 'coarsen' fraction.
    [err_marker]
      type      = ErrorFractionMarker
      indicator = jump_indicator
      refine    = 0.5    # refine the worst 50% of elements
      coarsen   = 0.05   # coarsen the best 5%
    []
  []
[]

[Executioner]
  type       = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]

[Outputs]
  exodus = true   # output includes refined mesh geometry
[]
