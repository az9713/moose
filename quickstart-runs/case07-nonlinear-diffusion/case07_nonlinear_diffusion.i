# ============================================================
# Case 7: Nonlinear Diffusion
# -div(k(T)*grad T) = Q,  T=0 on all walls
# k(T) = 1 + T,   Q = 10
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
  # ADMatDiffusion is the AD (automatic-differentiation) version of MatDiffusion.
  # It computes the exact Jacobian via dual-number arithmetic, which gives
  # quadratic Newton convergence even for nonlinear k(T).
  [diffusion]
    type        = ADMatDiffusion
    variable    = T
    diffusivity = k   # references the AD material property declared below
  []

  [source]
    type     = BodyForce
    variable = T
    value    = 10.0
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

[Materials]
  # ADPiecewiseLinearInterpolationMaterial defines k as a piecewise-linear
  # function of the coupled variable T.  For k(T) = 1 + T we just need two
  # points that span the expected range.  AD gives automatic derivatives.
  [conductivity]
    type        = ADPiecewiseLinearInterpolationMaterial
    property    = k
    variable    = T
    xy_data     = '0 1
                   10 11'   # k(T) = 1 + T, linear between T=0 and T=10
  []
[]

[Postprocessors]
  [max_T]
    type       = ElementExtremeValue
    variable   = T
    value_type = max
  []
[]

[Executioner]
  type       = Steady
  solve_type = 'NEWTON'   # use full Newton (not PJFNK) since AD provides exact Jacobian

  # Tighter tolerances to demonstrate quadratic convergence.
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12

  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]

[Outputs]
  exodus = true
  csv    = true
[]
