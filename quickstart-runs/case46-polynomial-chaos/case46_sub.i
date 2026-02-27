# ============================================================
# Case 46 Sub App: 1D Diffusion-Reaction
# Solves -D*u'' + sigma*u = 1 on [0,10]
# u(0) = u(10) = 0
# D and sigma are overridden by the sampler.
#
# MatReaction residual = -rate * u * test
# To get +sigma*u*test, we create neg_sigma = -sigma via
# DerivativeParsedMaterial, so -(-sigma)*u*test = +sigma*u*test.
# ============================================================

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 100
    xmax = 10
  []
[]

[Variables]
  [u]
  []
[]

[Kernels]
  [diffusion]
    type = MatDiffusion
    variable = u
    diffusivity = D
  []
  [absorption]
    type = MatReaction
    variable = u
    reaction_rate = neg_sigma  # -(-sigma)*u*test = +sigma*u*test
  []
  [source]
    type = BodyForce
    variable = u
    value = 1.0
  []
[]

[Materials]
  [diffusivity]
    type = GenericConstantMaterial
    prop_names = D
    prop_values = 5.0  # default; overridden by sampler
  []
  [sigma_input]
    type = GenericConstantMaterial
    prop_names = sigma
    prop_values = 5.0  # positive sigma; overridden by sampler
  []
  [sigma_negate]
    type = DerivativeParsedMaterial
    property_name = neg_sigma
    expression = '-sigma'
    material_property_names = 'sigma'
    disable_fpoptimizer = true
    enable_jit = false
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 0
  []
  [right]
    type = DirichletBC
    variable = u
    boundary = right
    value = 0
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Controls]
  [stochastic]
    type = SamplerReceiver
  []
[]

[Postprocessors]
  [avg]
    type = AverageNodalVariableValue
    variable = u
  []
[]

[Outputs]
[]
