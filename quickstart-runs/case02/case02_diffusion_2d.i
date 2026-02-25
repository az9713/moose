# ============================================================
# Case 2: Steady-State 2-D Diffusion
# Solves -div(grad u) = 0 on the unit square
# Exact solution: u(x,y) = x
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2    # now 2-D
  nx   = 20   # elements in x
  ny   = 20   # elements in y
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
  # Named boundaries on a 2-D GeneratedMesh:
  #   left   -> x = xmin
  #   right  -> x = xmax
  #   bottom -> y = ymin
  #   top    -> y = ymax
[]

[Variables]
  [u]
  []
[]

[Kernels]
  [diffusion]
    type     = Diffusion   # same kernel, works in any dimension
    variable = u
  []
[]

[BCs]
  [left]
    type     = DirichletBC
    variable = u
    boundary = left
    value    = 0
  []
  [right]
    type     = DirichletBC
    variable = u
    boundary = right
    value    = 1
  []
  # No BCs on 'top' and 'bottom' => zero-flux Neumann (the natural condition).
[]

[Executioner]
  type       = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]

[Outputs]
  # Exodus output — open with ParaView for 2-D visualization.
  exodus = true
[]
