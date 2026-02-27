# ============================================================
# Case 48 Sub App: 2D Steady Heat Conduction
# Solves -k*laplacian(T) = 100 on [0,1]x[0,1]
# T = T_left on left, T = T_right on right, insulated top/bottom.
# k, T_left, T_right are overridden by the ParameterStudy action.
# The source term makes the solution depend on k (not just BCs).
# ============================================================

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 20
    ny = 20
  []
[]

[Variables]
  [T]
  []
[]

[Kernels]
  [diffusion]
    type = MatDiffusion
    variable = T
    diffusivity = k
  []
  [source]
    type = BodyForce
    variable = T
    value = 100.0  # volumetric heat source makes solution depend on k
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = T
    boundary = left
    value = 50.0  # default; overridden by sampler
  []
  [right]
    type = DirichletBC
    variable = T
    boundary = right
    value = 300.0  # default; overridden by sampler
  []
[]

[Materials]
  [thermal]
    type = GenericConstantMaterial
    prop_names = k
    prop_values = 2.0  # default; overridden by sampler
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Postprocessors]
  [avg_T]
    type = AverageNodalVariableValue
    variable = T
  []
  [max_T]
    type = NodalExtremeValue
    variable = T
    value_type = max
  []
[]

[Outputs]
[]
