# ============================================================
# Case 5: Spatially Varying Conductivity
# -div((1+x)*grad u) = 0,  u=0 left, u=1 right
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
[]

[Functions]
  # A function used to define the conductivity field k(x) = 1 + x.
  [k_fn]
    type       = ParsedFunction
    expression = '1 + x'
  []
[]

[Materials]
  # GenericFunctionMaterial maps Function objects to material properties.
  # At each quadrature point it evaluates k_fn(x,y,z,t) and stores the
  # result in the material property 'k'.
  [conductivity]
    type        = GenericFunctionMaterial
    prop_names  = 'k'
    prop_values = 'k_fn'   # references the [Functions] block entry above
  []
[]

[Kernels]
  # MatDiffusion uses the material property 'k' at each quadrature point.
  # This is the term:  int( k(x) * grad(phi_i) . grad(u) ) dV
  [diffusion]
    type        = MatDiffusion
    variable    = u
    diffusivity = k
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
[]

[Postprocessors]
  # Track the average solution value for a quick sanity check.
  [avg_u]
    type     = ElementAverageValue
    variable = u
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
  csv    = true
[]
