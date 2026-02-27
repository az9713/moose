# ============================================================
# Case 45 Sub App: 2D Steady Heat Conduction with Source
# Solves -k*laplacian(T) = 100 on [0,1]^2
# T=0 on all boundaries (homogeneous Dirichlet).
# k is overridden by the sampler via SamplerReceiver.
#
# With a volumetric source, T ~ q/(2k), so the solution
# depends inversely on k — lower k means higher temperatures.
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
    value = 100.0  # volumetric heat source W/m^3
  []
[]

[BCs]
  [all]
    type = DirichletBC
    variable = T
    boundary = 'left right top bottom'
    value = 0
  []
[]

[Materials]
  [thermal]
    type = GenericConstantMaterial
    prop_names = k
    prop_values = 1.0  # default; overridden by sampler
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
