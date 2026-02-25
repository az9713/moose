# ============================================================
# Case 4: Poisson with Manufactured Solution (MMS)
# -div(grad u) = f,  u = 0 on boundary
# u_exact = sin(pi*x)*sin(pi*y)
# f = 2*pi^2*sin(pi*x)*sin(pi*y)
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 20
  ny   = 20
[]

[Variables]
  [u]
  []
[]

[Functions]
  # ParsedFunction evaluates a mathematical expression at (x, y, z, t).
  # Built-in constants: pi.
  [exact_fn]
    type       = ParsedFunction
    expression = 'sin(pi*x)*sin(pi*y)'
  []

  [forcing_fn]
    type       = ParsedFunction
    expression = '2*pi^2*sin(pi*x)*sin(pi*y)'
  []
[]

[Kernels]
  [diffusion]
    type     = Diffusion
    variable = u
  []

  # BodyForce applies  int( phi_i * f ) dV.
  # The 'function' parameter evaluates a Function object at each quadrature point.
  [source]
    type     = BodyForce
    variable = u
    function = forcing_fn
  []
[]

[BCs]
  # FunctionDirichletBC evaluates a Function object at each boundary node.
  # Since sin(pi*x)*sin(pi*y) = 0 on all four sides of [0,1]^2,
  # this is equivalent to a zero Dirichlet BC here.
  [all_walls]
    type     = FunctionDirichletBC
    variable = u
    boundary = 'left right top bottom'
    function = exact_fn
  []
[]

[Postprocessors]
  # ElementL2Error computes || u_h - u_exact ||_{L2}
  # where u_h is the FE solution and u_exact is a Function.
  [L2_error]
    type     = ElementL2Error
    variable = u
    function = exact_fn
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
  csv    = true   # captures L2_error value
[]
